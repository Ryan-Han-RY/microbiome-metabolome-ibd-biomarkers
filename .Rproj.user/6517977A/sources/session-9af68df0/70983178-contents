# Differential abundance analysis for genus-level taxa using MaasLin2

source("R/00_utils_io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  stop("Package `Maaslin2` is not installed. Run: BiocManager::install('Maaslin2')")
}

ensure_project_dirs()
ensure_dir("results/models/maaslin2_taxa")

message_header("Load genus-level microbiome data and metadata")

metadata_all <- readRDS("data/processed/matched_metadata.rds")
taxa_genus_rel <- readRDS("data/processed/taxa_genus_filtered.rds")

taxa_genus_rel <- as.data.frame(taxa_genus_rel)
taxa_genus_rel[] <- lapply(taxa_genus_rel, as.numeric)

metadata_primary <- metadata_all %>%
  filter(is_primary_independent) %>%
  mutate(
    diagnosis_clean = factor(diagnosis_clean, levels = c("non_IBD", "CD", "UC")),
    diagnosis_plot = recode(as.character(diagnosis_clean), "non_IBD" = "non-IBD")
  )

common_ids <- intersect(metadata_primary$sample_id_clean, rownames(taxa_genus_rel))

metadata_primary <- metadata_primary %>%
  filter(sample_id_clean %in% common_ids) %>%
  arrange(sample_id_clean)

taxa_genus_rel <- taxa_genus_rel[metadata_primary$sample_id_clean, , drop = FALSE]

message("Samples used for MaasLin2 taxa analysis: ", nrow(metadata_primary))
message("Genus-level taxa used: ", ncol(taxa_genus_rel))

message_header("Prepare MaasLin2 feature names")

make_taxa_label <- function(x) {
  y <- x
  y <- str_replace(y, "^.*\\|g__", "")
  y <- str_replace(y, "^g__", "")
  y <- str_replace_all(y, "[^A-Za-z0-9]+", "_")
  y <- str_replace_all(y, "^_+|_+$", "")
  y[y == "" | is.na(y)] <- "unknown_taxon"
  make.unique(y, sep = "_")
}

taxa_original_names <- colnames(taxa_genus_rel)
taxa_maaslin_names <- make_taxa_label(taxa_original_names)

taxa_feature_map <- tibble(
  original_feature_id = taxa_original_names,
  maaslin_feature_id = taxa_maaslin_names
)

colnames(taxa_genus_rel) <- taxa_maaslin_names

write_tsv_safe(
  taxa_feature_map,
  "results/tables/maaslin2_taxa_feature_name_map.tsv"
)

message_header("Prepare metadata for MaasLin2")

clean_covariate <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

metadata_model_base <- metadata_primary %>%
  transmute(
    sample_id_clean = sample_id_clean,
    diagnosis_clean = as.character(diagnosis_clean),
    participant_id_clean = participant_id_clean,
    age_model = age_clean,
    sex_model = clean_covariate(sex_clean),
    antibiotic_model = clean_covariate(antibiotic_use_clean)
  )

select_covariates <- function(meta) {
  covariates <- c()
  
  if ("age_model" %in% names(meta)) {
    age_ok <- mean(!is.na(meta$age_model)) >= 0.80 &&
      dplyr::n_distinct(meta$age_model, na.rm = TRUE) > 5
    if (age_ok) covariates <- c(covariates, "age_model")
  }
  
  if ("sex_model" %in% names(meta)) {
    sex_ok <- mean(!is.na(meta$sex_model)) >= 0.80 &&
      dplyr::n_distinct(meta$sex_model, na.rm = TRUE) >= 2
    if (sex_ok) covariates <- c(covariates, "sex_model")
  }
  
  if ("antibiotic_model" %in% names(meta)) {
    antibiotic_ok <- mean(!is.na(meta$antibiotic_model)) >= 0.80 &&
      dplyr::n_distinct(meta$antibiotic_model, na.rm = TRUE) >= 2
    if (antibiotic_ok) covariates <- c(covariates, "antibiotic_model")
  }
  
  covariates
}

prepare_comparison_metadata <- function(comparison_name) {
  meta <- metadata_model_base
  
  if (comparison_name == "IBD_vs_non_IBD") {
    meta <- meta %>%
      mutate(
        model_group = if_else(diagnosis_clean == "non_IBD", "non_IBD", "IBD")
      ) %>%
      filter(model_group %in% c("non_IBD", "IBD")) %>%
      mutate(model_group = factor(model_group, levels = c("non_IBD", "IBD")))
  }
  
  if (comparison_name == "CD_vs_non_IBD") {
    meta <- meta %>%
      filter(diagnosis_clean %in% c("non_IBD", "CD")) %>%
      mutate(model_group = factor(diagnosis_clean, levels = c("non_IBD", "CD")))
  }
  
  if (comparison_name == "UC_vs_non_IBD") {
    meta <- meta %>%
      filter(diagnosis_clean %in% c("non_IBD", "UC")) %>%
      mutate(model_group = factor(diagnosis_clean, levels = c("non_IBD", "UC")))
  }
  
  if (comparison_name == "CD_vs_UC") {
    meta <- meta %>%
      filter(diagnosis_clean %in% c("UC", "CD")) %>%
      mutate(model_group = factor(diagnosis_clean, levels = c("UC", "CD")))
  }
  
  covariates <- select_covariates(meta)
  
  meta <- meta %>%
    select(sample_id_clean, model_group, all_of(covariates)) %>%
    filter(!is.na(model_group))
  
  for (covar in covariates) {
    meta <- meta %>% filter(!is.na(.data[[covar]]))
  }
  
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$sample_id_clean
  meta$sample_id_clean <- NULL
  
  list(
    metadata = meta,
    fixed_effects = c("model_group", covariates),
    covariates_used = covariates
  )
}

extract_maaslin_results <- function(output_dir, comparison_name) {
  result_path <- file.path(output_dir, "all_results.tsv")
  
  if (!file.exists(result_path)) {
    stop("MaasLin2 did not create all_results.tsv in: ", output_dir)
  }
  
  result <- readr::read_tsv(result_path, show_col_types = FALSE)
  
  # MaasLin2 usually returns columns:
  # metadata, feature, value, coef, stderr, N, N.not.zero, pval, qval
  result <- result %>%
    as_tibble() %>%
    filter(metadata == "model_group") %>%
    mutate(
      comparison = comparison_name,
      feature_id = feature,
      feature_label = feature,
      beta = coef,
      standard_error = stderr,
      p_value = pval,
      q_value = qval,
      n_samples = N,
      method = "MaasLin2_LM_CLR"
    ) %>%
    select(
      feature_id,
      comparison,
      metadata,
      value,
      beta,
      standard_error,
      p_value,
      q_value,
      n_samples,
      method,
      feature_label,
      everything()
    ) %>%
    arrange(q_value, p_value)
  
  result
}

run_maaslin2_comparison <- function(comparison_name) {
  message_header(paste("Run MaasLin2:", comparison_name))
  
  cmp <- prepare_comparison_metadata(comparison_name)
  
  meta <- cmp$metadata
  sample_ids <- rownames(meta)
  
  data <- taxa_genus_rel[sample_ids, , drop = FALSE]
  data <- as.data.frame(data)
  
  output_dir <- file.path("results/models/maaslin2_taxa", comparison_name)
  
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  
  message("Samples: ", nrow(data))
  message("Features: ", ncol(data))
  message("Fixed effects: ", paste(cmp$fixed_effects, collapse = ", "))
  message("Covariates used: ", ifelse(length(cmp$covariates_used) == 0, "none", paste(cmp$covariates_used, collapse = ", ")))
  
  set.seed(20260707)
  
  Maaslin2::Maaslin2(
    input_data = data,
    input_metadata = meta,
    output = output_dir,
    fixed_effects = cmp$fixed_effects,
    random_effects = NULL,
    reference = NULL,
    normalization = "CLR",
    transform = "NONE",
    analysis_method = "LM",
    min_abundance = 0,
    min_prevalence = 0,
    min_variance = 0,
    correction = "BH",
    max_significance = 0.10,
    standardize = TRUE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = max(1, parallel::detectCores() - 1)
  )
  
  extract_maaslin_results(output_dir, comparison_name)
}

comparisons <- c(
  "IBD_vs_non_IBD",
  "CD_vs_non_IBD",
  "UC_vs_non_IBD",
  "CD_vs_UC"
)

taxa_results <- map_dfr(comparisons, run_maaslin2_comparison) %>%
  left_join(
    taxa_feature_map,
    by = c("feature_id" = "maaslin_feature_id")
  ) %>%
  mutate(
    original_feature_id = if_else(is.na(original_feature_id), feature_id, original_feature_id),
    direction = case_when(
      beta > 0 ~ "higher_in_case_group",
      beta < 0 ~ "lower_in_case_group",
      TRUE ~ "no_direction"
    )
  ) %>%
  arrange(comparison, q_value, p_value)

taxa_significant <- taxa_results %>%
  filter(!is.na(q_value), q_value < 0.10) %>%
  arrange(q_value, p_value)

write_tsv_safe(
  taxa_results,
  "results/tables/maaslin2_taxa_all_results.tsv"
)

write_tsv_safe(
  taxa_significant,
  "results/tables/maaslin2_taxa_significant_results.tsv"
)

message_header("Create MaasLin2 taxa volcano plot")

volcano_df <- taxa_results %>%
  mutate(
    neg_log10_q = -log10(pmax(q_value, 1e-300)),
    significant = if_else(!is.na(q_value) & q_value < 0.10, "FDR < 0.10", "Not significant")
  )

p_volcano <- ggplot(volcano_df, aes(x = beta, y = neg_log10_q, shape = significant)) +
  geom_point(alpha = 0.75, size = 1.8) +
  facet_wrap(~ comparison, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Genus-level taxa differential abundance using MaasLin2",
    x = "MaasLin2 coefficient after CLR normalization",
    y = "-log10(FDR)",
    shape = "Result"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/taxa_effect_size_volcano.png",
  p_volcano,
  width = 10,
  height = 7,
  dpi = 300
)

message_header("Create MaasLin2 taxa dotplot")

dot_df <- taxa_results %>%
  group_by(comparison) %>%
  arrange(q_value, p_value) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(
    feature_label = str_trunc(feature_label, width = 45),
    neg_log10_p = -log10(pmax(p_value, 1e-300))
  )

p_dot <- ggplot(dot_df, aes(x = beta, y = reorder(feature_label, beta))) +
  geom_point(aes(size = neg_log10_p, shape = q_value < 0.10), alpha = 0.85) +
  facet_wrap(~ comparison, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Top disease-associated genus-level taxa from MaasLin2",
    x = "MaasLin2 coefficient after CLR normalization",
    y = "Taxon",
    size = "-log10(p)",
    shape = "FDR < 0.10"
  ) +
  theme_bw(base_size = 11)

ggsave(
  "results/figures/taxa_differential_abundance_dotplot.png",
  p_dot,
  width = 11,
  height = 8,
  dpi = 300
)

message_header("Update analysis decisions")

decision_addition <- c(
  "",
  "## MaasLin2 differential abundance modelling",
  "",
  "- Genus-level microbiome differential abundance was performed using MaasLin2.",
  "- Input features were filtered genus-level relative abundance tables from the primary independent set.",
  "- MaasLin2 was run with LM analysis, CLR normalization, no additional transform, and BH FDR correction.",
  "- The primary independent set was used to avoid repeated-measure leakage.",
  "- Comparisons tested separately: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.",
  "- Candidate covariates were age, sex, and antibiotic use when sufficiently available.",
  "- Significance threshold was FDR < 0.10."
)

cat(paste(decision_addition, collapse = "\n"), file = "docs/analysis_decisions.md", append = TRUE)

message_header("MaasLin2 taxa analysis completed")

message("Total MaasLin2 taxa results: ", nrow(taxa_results))
message("FDR < 0.10 MaasLin2 taxa results: ", nrow(taxa_significant))

print(taxa_significant %>% count(comparison))

message("Outputs written:")
message("- results/tables/maaslin2_taxa_all_results.tsv")
message("- results/tables/maaslin2_taxa_significant_results.tsv")
message("- results/tables/maaslin2_taxa_feature_name_map.tsv")
message("- results/models/maaslin2_taxa/")
message("- results/figures/taxa_effect_size_volcano.png")
message("- results/figures/taxa_differential_abundance_dotplot.png")

