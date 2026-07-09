# Construct matched microbiome-metabolome cohort and generate QC outputs

source("R/00_utils_io.R")

ensure_project_dirs()

message_header("Load processed objects")

metadata <- readRDS("data/processed/metadata_clean.rds")
microbiome <- readRDS("data/processed/microbiome_raw.rds")
metabolomics <- readRDS("data/processed/metabolomics_raw.rds")

micro_mat <- microbiome$abundance
metab_mat <- metabolomics$abundance

metadata_ids <- unique(metadata$sample_id_clean)
metadata_ids <- metadata_ids[!is.na(metadata_ids) & metadata_ids != ""]

micro_ids <- rownames(micro_mat)
metab_ids <- rownames(metab_mat)

message("Metadata unique sample IDs: ", length(metadata_ids))
message("Microbiome sample IDs: ", length(micro_ids))
message("Metabolomics sample IDs: ", length(metab_ids))

message_header("Match samples across metadata, microbiome, and metabolomics")

matched_ids <- Reduce(intersect, list(metadata_ids, micro_ids, metab_ids))

if (length(matched_ids) == 0) {
  stop(
    "No matched samples found across metadata, microbiome, and metabolomics.\n",
    "Check sample ID suffixes and metadata External ID column."
  )
}

matched_metadata_raw <- metadata %>%
  filter(sample_id_clean %in% matched_ids) %>%
  arrange(sample_id_clean) %>%
  group_by(sample_id_clean) %>%
  slice(1) %>%
  ungroup()

matched_metadata <- matched_metadata_raw %>%
  filter(
    !is.na(sample_id_clean),
    sample_id_clean != "",
    !is.na(participant_id_clean),
    participant_id_clean != "",
    !is.na(diagnosis_clean),
    diagnosis_clean %in% c("CD", "UC", "non_IBD")
  )

matched_metadata <- matched_metadata %>%
  mutate(
    week_sort = if_else(is.na(week_num_clean), Inf, week_num_clean),
    baseline_priority = if_else(!is.na(week_num_clean) & week_num_clean == 0, 0, 1)
  ) %>%
  arrange(participant_id_clean, baseline_priority, week_sort, sample_id_clean)

primary_independent_ids <- matched_metadata %>%
  group_by(participant_id_clean) %>%
  slice(1) %>%
  ungroup() %>%
  pull(sample_id_clean)

matched_metadata <- matched_metadata %>%
  mutate(
    is_primary_independent = sample_id_clean %in% primary_independent_ids
  ) %>%
  select(-week_sort, -baseline_priority)

matched_sample_ids <- matched_metadata$sample_id_clean

micro_matched <- micro_mat[matched_sample_ids, , drop = FALSE]
metab_matched <- metab_mat[matched_sample_ids, , drop = FALSE]

saveRDS(matched_metadata, "data/processed/matched_metadata.rds")
saveRDS(micro_matched, "data/processed/microbiome_matched_raw.rds")
saveRDS(metab_matched, "data/processed/metabolomics_matched_raw.rds")

message("Matched samples after metadata filter: ", nrow(matched_metadata))
message("Matched participants: ", n_distinct(matched_metadata$participant_id_clean))
message("Primary independent samples: ", sum(matched_metadata$is_primary_independent))

message_header("Build sample matching flow table")

count_participants_for_ids <- function(ids) {
  tmp <- metadata %>%
    filter(sample_id_clean %in% ids)
  
  count_unique_nonmissing(tmp$participant_id_clean)
}

micro_metab_ids <- intersect(micro_ids, metab_ids)
matched_all_raw_ids <- matched_ids

sample_matching_flow <- tibble(
  step = c(
    "metadata_total",
    "microbiome_available",
    "metabolomics_available",
    "matched_microbiome_metabolome",
    "matched_metadata_microbiome_metabolome",
    "after_missing_metadata_filter",
    "primary_independent_subset"
  ),
  n_samples = c(
    length(metadata_ids),
    length(micro_ids),
    length(metab_ids),
    length(micro_metab_ids),
    length(matched_all_raw_ids),
    nrow(matched_metadata),
    sum(matched_metadata$is_primary_independent)
  ),
  n_participants = c(
    count_unique_nonmissing(metadata$participant_id_clean),
    count_participants_for_ids(intersect(metadata_ids, micro_ids)),
    count_participants_for_ids(intersect(metadata_ids, metab_ids)),
    count_participants_for_ids(intersect(metadata_ids, micro_metab_ids)),
    count_participants_for_ids(matched_all_raw_ids),
    n_distinct(matched_metadata$participant_id_clean),
    n_distinct(matched_metadata$participant_id_clean[matched_metadata$is_primary_independent])
  )
)

write_tsv_safe(sample_matching_flow, "results/tables/sample_matching_flow.tsv")
print(sample_matching_flow)

message_header("Build group distribution table")

group_distribution <- bind_rows(
  matched_metadata %>%
    count(diagnosis_clean, name = "n_samples") %>%
    mutate(
      analysis_set = "longitudinal_matched_set",
      n_participants = map_int(diagnosis_clean, function(g) {
        matched_metadata %>%
          filter(diagnosis_clean == g) %>%
          pull(participant_id_clean) %>%
          n_distinct()
      })
    ),
  matched_metadata %>%
    filter(is_primary_independent) %>%
    count(diagnosis_clean, name = "n_samples") %>%
    mutate(
      analysis_set = "primary_independent_set",
      n_participants = n_samples
    )
) %>%
  select(analysis_set, diagnosis_clean, n_samples, n_participants) %>%
  arrange(analysis_set, diagnosis_clean)

write_tsv_safe(group_distribution, "results/tables/group_distribution.tsv")
print(group_distribution)

message_header("Create QC figures")

flow_plot_data <- sample_matching_flow %>%
  mutate(step = factor(step, levels = step))

p1 <- ggplot(flow_plot_data, aes(x = step, y = n_samples)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Sample retention during multi-omics cohort construction",
    x = NULL,
    y = "Number of samples"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = "results/figures/sample_matching_flow.png",
  plot = p1,
  width = 8,
  height = 5,
  dpi = 300
)

p2 <- ggplot(
  group_distribution,
  aes(x = diagnosis_clean, y = n_samples, fill = diagnosis_clean)
) +
  geom_col() +
  facet_wrap(~ analysis_set) +
  labs(
    title = "Diagnosis group distribution after sample matching",
    x = "Diagnosis group",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = "results/figures/group_distribution_barplot.png",
  plot = p2,
  width = 8,
  height = 4.5,
  dpi = 300
)

message_header("Write analysis decisions")

column_map <- attr(metadata, "column_map")

n_metadata_no_microbiome <- length(setdiff(metadata_ids, micro_ids))
n_metadata_no_metabolomics <- length(setdiff(metadata_ids, metab_ids))
n_micro_metab_not_metadata <- length(setdiff(micro_metab_ids, metadata_ids))
n_missing_metadata_filter <- nrow(matched_metadata_raw) - nrow(matched_metadata)
n_repeated_participants <- matched_metadata %>%
  count(participant_id_clean) %>%
  filter(n > 1) %>%
  nrow()

decision_lines <- c(
  "# Analysis decisions",
  "",
  "## Cohort construction",
  "",
  paste0("- Metadata sample identifier source column: `", column_map$sample_id_col, "`."),
  paste0("- Participant identifier source column: `", column_map$participant_id_col, "`."),
  paste0("- Diagnosis source column: `", column_map$diagnosis_col, "`."),
  paste0("- Week number source column: `", ifelse(is.na(column_map$week_col), "not available", column_map$week_col), "`."),
  "- Omics sample identifiers were normalised by removing technical suffixes such as `_taxonomic_profile`, `_taxonomic`, and `_metabolomics`.",
  "- Matched samples were defined as sample IDs present in metadata, microbiome taxonomic profiles, and stool metabolomics profiles.",
  "- Retained diagnosis groups: CD, UC, and non_IBD.",
  "",
  "## Repeated samples",
  "",
  "- Repeated samples were retained in `data/processed/matched_metadata.rds` for longitudinal sensitivity analysis.",
  "- A participant-level independent subset was created using the column `is_primary_independent`.",
  "- For participant-level exploratory integration, one stool sample per participant was retained.",
  "- Selection rule for the independent subset: prefer baseline week when available; otherwise use the earliest available week; ties are broken by sample ID.",
  "",
  "## Exclusion summary",
  "",
  paste0("- Metadata samples without microbiome profile: ", n_metadata_no_microbiome, "."),
  paste0("- Metadata samples without metabolomics profile: ", n_metadata_no_metabolomics, "."),
  paste0("- Microbiome-metabolomics samples not found in metadata: ", n_micro_metab_not_metadata, "."),
  paste0("- Matched samples removed due to missing participant ID or diagnosis: ", n_missing_metadata_filter, "."),
  "",
  "## Analysis sets",
  "",
  paste0("- Longitudinal matched set: `data/processed/matched_metadata.rds`, all rows where `is_primary_independent` can be TRUE or FALSE; n samples = ", nrow(matched_metadata), "."),
  paste0("- Primary independent set: rows where `is_primary_independent == TRUE`; n samples = ", sum(matched_metadata$is_primary_independent), "."),
  paste0("- Number of participants with repeated matched samples: ", n_repeated_participants, ".")
)

writeLines(decision_lines, "docs/analysis_decisions.md")

message_header("QC completed")

message("Main outputs written:")
message("- data/processed/metadata_clean.rds")
message("- data/processed/microbiome_raw.rds")
message("- data/processed/metabolomics_raw.rds")
message("- data/processed/matched_metadata.rds")
message("- results/tables/raw_table_dimensions.tsv")
message("- results/tables/sample_id_overlap_before_matching.tsv")
message("- results/tables/sample_matching_flow.tsv")
message("- results/tables/group_distribution.tsv")
message("- results/figures/sample_matching_flow.png")
message("- results/figures/group_distribution_barplot.png")
message("- docs/analysis_decisions.md")