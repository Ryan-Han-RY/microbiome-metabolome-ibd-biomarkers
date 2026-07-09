# Microbiome taxonomic preprocessing: filtering, pseudocount, CLR transformation

source("R/00_utils_io.R")

ensure_project_dirs()

message_header("Load matched microbiome data")

metadata <- readRDS("data/processed/matched_metadata.rds")
microbiome_raw <- readRDS("data/processed/microbiome_raw.rds")
micro_mat <- readRDS("data/processed/microbiome_matched_raw.rds")

feature_info <- microbiome_raw$feature_info

metadata <- metadata %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD")
  )

standardise_relative_abundance <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  
  row_sums <- rowSums(mat, na.rm = TRUE)
  median_sum <- median(row_sums[row_sums > 0], na.rm = TRUE)
  
  if (is.finite(median_sum) && median_sum > 2) {
    mat <- mat / 100
  }
  
  mat
}

close_rows <- function(mat) {
  rs <- rowSums(mat, na.rm = TRUE)
  mat <- mat[rs > 0, , drop = FALSE]
  rs <- rowSums(mat, na.rm = TRUE)
  sweep(mat, 1, rs, "/")
}

clr_transform <- function(mat, pseudocount = NULL) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  positive_values <- mat[mat > 0 & is.finite(mat)]
  
  if (length(positive_values) == 0) {
    stop("No positive values found for CLR transformation.")
  }
  
  if (is.null(pseudocount)) {
    pseudocount <- min(positive_values, na.rm = TRUE) / 2
  }
  
  mat_pc <- mat + pseudocount
  mat_pc <- close_rows(mat_pc)
  
  log_mat <- log(mat_pc)
  gm <- rowMeans(log_mat)
  clr <- sweep(log_mat, 1, gm, "-")
  
  attr(clr, "pseudocount") <- pseudocount
  clr
}

remove_unclassified_taxa <- function(feature_ids) {
  x <- str_to_lower(feature_ids)
  
  keep <- !str_detect(
    x,
    "unclassified|unknown|uncharacterized|unassigned|noname|metagenome|g__$|s__$|\\|g__$|\\|s__$"
  )
  
  keep
}

prepare_tax_rank <- function(rank_name, output_prefix) {
  message_header(paste("Prepare", rank_name, "taxonomic table"))
  
  rank_features <- feature_info %>%
    filter(tax_rank == rank_name) %>%
    pull(feature_id)
  
  rank_features <- intersect(rank_features, colnames(micro_mat))
  
  if (length(rank_features) < 5) {
    stop("Too few features detected for rank: ", rank_name)
  }
  
  mat_rank <- micro_mat[, rank_features, drop = FALSE]
  mat_rank <- standardise_relative_abundance(mat_rank)
  
  keep_classified <- remove_unclassified_taxa(colnames(mat_rank))
  mat_classified <- mat_rank[, keep_classified, drop = FALSE]
  
  prevalence <- colMeans(mat_classified > 0, na.rm = TRUE)
  mean_abundance <- colMeans(mat_classified, na.rm = TRUE)
  
  keep_prevalence <- prevalence >= 0.10
  keep_abundance <- mean_abundance >= 0.0001
  keep_final <- keep_prevalence & keep_abundance
  
  mat_filtered <- mat_classified[, keep_final, drop = FALSE]
  mat_filtered <- close_rows(mat_filtered)
  
  mat_clr <- clr_transform(mat_filtered)
  
  filtering_summary <- tibble(
    taxonomic_rank = rank_name,
    step = c(
      "raw_rank_features",
      "after_removing_unclassified_taxa",
      "after_prevalence_filter",
      "after_mean_abundance_filter",
      "final_filtered_features"
    ),
    n_features = c(
      ncol(mat_rank),
      ncol(mat_classified),
      sum(keep_classified & prevalence >= 0.10),
      sum(keep_classified & mean_abundance >= 0.0001),
      ncol(mat_filtered)
    ),
    n_samples = c(
      nrow(mat_rank),
      nrow(mat_classified),
      nrow(mat_classified),
      nrow(mat_classified),
      nrow(mat_filtered)
    )
  )
  
  saveRDS(mat_filtered, paste0("data/processed/taxa_", output_prefix, "_filtered.rds"))
  saveRDS(mat_clr, paste0("data/processed/taxa_", output_prefix, "_clr.rds"))
  
  list(
    filtered = mat_filtered,
    clr = mat_clr,
    summary = filtering_summary
  )
}

genus_result <- prepare_tax_rank("genus", "genus")
species_result <- prepare_tax_rank("species", "species")

taxa_filtering_summary <- bind_rows(
  genus_result$summary,
  species_result$summary
)

write_tsv_safe(
  taxa_filtering_summary,
  "results/tables/taxa_filtering_summary.tsv"
)

message_header("Create top taxa relative abundance figure")

genus_rel <- genus_result$filtered

top_taxa <- sort(colMeans(genus_rel, na.rm = TRUE), decreasing = TRUE)
top_taxa <- names(top_taxa)[seq_len(min(15, length(top_taxa)))]

plot_df <- genus_rel[, top_taxa, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id_clean") %>%
  pivot_longer(
    cols = -sample_id_clean,
    names_to = "taxon",
    values_to = "relative_abundance"
  ) %>%
  left_join(
    metadata %>% select(sample_id_clean, diagnosis_clean),
    by = "sample_id_clean"
  ) %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD"),
    taxon_label = str_replace(taxon, "^.*\\|", ""),
    taxon_label = str_replace(taxon_label, "^g__", "")
  ) %>%
  group_by(diagnosis_plot, taxon_label) %>%
  summarise(
    mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(
  plot_df,
  aes(x = reorder(taxon_label, mean_relative_abundance), y = mean_relative_abundance)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ diagnosis_plot) +
  labs(
    title = "Top genus-level taxa by diagnosis group",
    x = "Taxon",
    y = "Mean relative abundance"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/top_taxa_relative_abundance_barplot.png",
  p,
  width = 9,
  height = 6,
  dpi = 300
)

message_header("Microbiome preprocessing completed")

print(taxa_filtering_summary)

message("Outputs written:")
message("- data/processed/taxa_genus_filtered.rds")
message("- data/processed/taxa_genus_clr.rds")
message("- data/processed/taxa_species_filtered.rds")
message("- data/processed/taxa_species_clr.rds")
message("- results/tables/taxa_filtering_summary.tsv")
message("- results/figures/top_taxa_relative_abundance_barplot.png")