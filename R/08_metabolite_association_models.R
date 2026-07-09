# Metabolite association analysis using MaAsLin2
#
# Input:
#   data/processed/matched_metadata.rds
#   data/processed/metabolites_filtered.rds
#   data/processed/metabolites_log_scaled.rds
#
# Output:
#   results/tables/metabolite_association_all_results.tsv
#   results/tables/metabolite_association_significant_results.tsv
#   results/tables/maaslin2_metabolite_feature_name_map.tsv
#   results/models/maaslin2_metabolites/
#   results/figures/metabolite_volcano.png
#   results/figures/top_metabolites_heatmap.png

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
  stop(
    "Package `Maaslin2` is not installed.\n",
    "Run this first:\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install('Maaslin2')"
  )
}

ensure_project_dirs()
ensure_dir("results/models/maaslin2_metabolites")
ensure_dir("results/tables")
ensure_dir("results/figures")

message_header("Load metabolomics data and metadata")

metadata_all <- readRDS("data/processed/matched_metadata.rds")
metabolites_filtered <- readRDS("data/processed/metabolites_filtered.rds")
metabolites_log_scaled <- readRDS("data/processed/metabolites_log_scaled.rds")

metabolites_filtered <- as.data.frame(metabolites_filtered)
metabolites_filtered[] <- lapply(metabolites_filtered, as.numeric)

metabolites_log_scaled <- as.data.frame(metabolites_log_scaled)
metabolites_log_scaled[] <- lapply(metabolites_log_scaled, as.numeric)

metadata_primary <- metadata_all %>%
  filter(is_primary_independent) %>%
  mutate(
    diagnosis_clean = as.character(diagnosis_clean),
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD")
  )

common_ids <- Reduce(
  intersect,
  list(
    metadata_primary$sample_id_clean,
    rownames(metabolites_filtered),
    rownames(metabolites_log_scaled)
  )
)

metadata_primary <- metadata_primary %>%
  filter(sample_id_clean %in% common_ids) %>%
  arrange(sample_id_clean)

metabolites_filtered <- metabolites_filtered[metadata_primary$sample_id_clean, , drop = FALSE]
metabolites_log_scaled <- metabolites_log_scaled[metadata_primary$sample_id_clean, , drop = FALSE]

message("Primary independent samples available in metabolomics table: ", nrow(metadata_primary))
message("Metabolite features available: ", ncol(metabolites_filtered))

if (nrow(metadata_primary) < 20) {
  stop("Too few primary independent samples available for MaAsLin2 metabolite analysis.")
}

if (ncol(metabolites_filtered) < 10) {
  stop("Too few metabolite features available for MaAsLin2 metabolite analysis.")
}

message_header("Prepare Maaslin2 metabolite feature names")

make_safe_feature_name <- function(x) {
  y <- as.character(x)
  y <- str_replace_all(y, "[^A-Za-z0-9]+", "_")
  y <- str_replace_all(y, "^_+|_+$", "")
  y[y == "" | is.na(y)] <- "metabolite_feature"
  make.unique(y, sep = "_")
}

metabolite_original_names <- colnames(metabolites_filtered)
metabolite_maaslin_names <- make_safe_feature_name(metabolite_original_names)

metabolite_feature_map <- tibble(
  original_metabolite_id = metabolite_original_names,
  maaslin_metabolite_id = metabolite_maaslin_names
)

colnames(metabolites_filtered) <- metabolite_maaslin_names

write_tsv_safe(
  metabolite_feature_map,
  "results/tables/maaslin2_metabolite_feature_name_map.tsv"
)

message_header("Prepare metadata for Maaslin2")

clean_covariate <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

metadata_model_base <- metadata_primary %>%
  transmute(
    sample_id_clean = sample_id_clean,
    diagnosis_clean = diagnosis_clean,
    age_model = age_clean,
    sex_model = clean_covariate(sex_clean),
    antibiotic_model = clean_covariate(antibiotic_use_clean)
  )

select_covariates <- function(meta) {
  covariates <- c()
  
  if ("age_model" %in% names(meta)) {
    age_ok <- mean(!is.na(meta$age_model)) >= 0.80 &&
      dplyr::n_distinct(meta$age_model, na.rm = TRUE) > 5
    if (age_ok) {
      covariates <- c(covariates, "age_model")
    }
  }
  
  if ("sex_model" %in% names(meta)) {
    sex_ok <- mean(!is.na(meta$sex_model)) >= 0.80 &&
      dplyr::n_distinct(meta$sex_model, na.rm = TRUE) >= 2
    if (sex_ok) {
      covariates <- c(covariates, "sex_model")
    }
  }
  
  if ("antibiotic_model" %in% names(meta)) {
    antibiotic_ok <- mean(!is.na(meta$antibiotic_model)) >= 0.80 &&
      dplyr::n_distinct(meta$antibiotic_model, na.rm = TRUE) >= 2
    if (antibiotic_ok) {
      covariates <- c(covariates, "antibiotic_model")
    }
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
  
  if (!"model_group" %in% names(meta)) {
    stop("Unknown comparison name: ", comparison_name)
  }
  
  if (dplyr::n_distinct(meta$model_group) < 2) {
    stop("Comparison has fewer than two groups: ", comparison_name)
  }
  
  covariates <- select_covariates(meta)
  
  meta <- meta %>%
    select(sample_id_clean, model_group, all_of(covariates)) %>%
    filter(!is.na(model_group))
  
  for (covariate_name in covariates) {
    meta <- meta %>%
      filter(!is.na(.data[[covariate_name]]))
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
    stop("Maaslin2 did not create all_results.tsv in: ", output_dir)
  }
  
  raw_result <- readr::read_tsv(result_path, show_col_types = FALSE)
  
  required_cols <- c(
    "feature", "metadata", "value", "coef", "stderr", "N", "pval", "qval"
  )
  
  missing_cols <- setdiff(required_cols, names(raw_result))
  
  if (length(missing_cols) > 0) {
    stop(
      "Maaslin2 output is missing expected columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  raw_result %>%
    as_tibble() %>%
    filter(metadata == "model_group") %>%
    mutate(
      comparison = comparison_name,
      metabolite_id = feature,
      beta = coef,
      standard_error = stderr,
      p_value = pval,
      q_value = qval,
      n_samples = N,
      method = "Maaslin2_LM_LOG",
      direction = case_when(
        beta > 0 ~ "higher_in_case_group",
        beta < 0 ~ "lower_in_case_group",
        TRUE ~ "no_direction"
      )
    ) %>%
    select(
      metabolite_id,
      comparison,
      metadata,
      value,
      beta,
      standard_error,
      p_value,
      q_value,
      n_samples,
      method,
      direction,
      everything()
    ) %>%
    arrange(q_value, p_value)
}

run_maaslin2_comparison <- function(comparison_name) {
  message_header(paste("Run Maaslin2 metabolite model:", comparison_name))
  
  comparison_setup <- prepare_comparison_metadata(comparison_name)
  
  model_metadata <- comparison_setup$metadata
  sample_ids <- rownames(model_metadata)
  
  model_data <- metabolites_filtered[sample_ids, , drop = FALSE]
  model_data <- as.data.frame(model_data)
  
  output_dir <- file.path("results/models/maaslin2_metabolites", comparison_name)
  
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  
  message("Samples: ", nrow(model_data))
  message("Features: ", ncol(model_data))
  message("Fixed effects: ", paste(comparison_setup$fixed_effects, collapse = ", "))
  message(
    "Covariates used: ",
    ifelse(
      length(comparison_setup$covariates_used) == 0,
      "none",
      paste(comparison_setup$covariates_used, collapse = ", ")
    )
  )
  
  set.seed(20260707)
  
  Maaslin2::Maaslin2(
    input_data = model_data,
    input_metadata = model_metadata,
    output = output_dir,
    fixed_effects = comparison_setup$fixed_effects,
    random_effects = NULL,
    reference = NULL,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    min_abundance = 0,
    min_prevalence = 0,
    min_variance = 0,
    correction = "BH",
    max_significance = 0.10,
    standardize = TRUE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = 1
  )
  
  extract_maaslin_results(output_dir, comparison_name)
}

comparisons <- c(
  "IBD_vs_non_IBD",
  "CD_vs_non_IBD",
  "UC_vs_non_IBD",
  "CD_vs_UC"
)

message_header("Run all Maaslin2 metabolite comparisons")

metabolite_results <- purrr::map_dfr(
  comparisons,
  run_maaslin2_comparison
) %>%
  left_join(
    metabolite_feature_map,
    by = c("metabolite_id" = "maaslin_metabolite_id")
  ) %>%
  mutate(
    original_metabolite_id = if_else(
      is.na(original_metabolite_id),
      metabolite_id,
      original_metabolite_id
    )
  ) %>%
  arrange(comparison, q_value, p_value)

metabolite_significant <- metabolite_results %>%
  filter(!is.na(q_value), q_value < 0.10) %>%
  arrange(q_value, p_value)

write_tsv_safe(
  metabolite_results,
  "results/tables/metabolite_association_all_results.tsv"
)

write_tsv_safe(
  metabolite_significant,
  "results/tables/metabolite_association_significant_results.tsv"
)

message_header("Create Maaslin2 metabolite volcano plot")

volcano_df <- metabolite_results %>%
  mutate(
    neg_log10_q = -log10(pmax(q_value, 1e-300)),
    significant = if_else(
      !is.na(q_value) & q_value < 0.10,
      "FDR < 0.10",
      "Not significant"
    )
  )

p_volcano <- ggplot(
  volcano_df,
  aes(x = beta, y = neg_log10_q, shape = significant)
) +
  geom_point(alpha = 0.65, size = 1.2) +
  facet_wrap(~ comparison, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Metabolite association models using MaAsLin2",
    x = "MaAsLin2 coefficient after log transform",
    y = "-log10(FDR)",
    shape = "Result"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/metabolite_volcano.png",
  p_volcano,
  width = 10,
  height = 7,
  dpi = 300
)

message_header("Create top metabolites heatmap")

top_metabolites_original <- metabolite_results %>%
  arrange(q_value, p_value) %>%
  filter(!is.na(p_value)) %>%
  distinct(original_metabolite_id, .keep_all = TRUE) %>%
  slice_head(n = 30) %>%
  pull(original_metabolite_id)

top_metabolites_original <- intersect(
  top_metabolites_original,
  colnames(metabolites_log_scaled)
)

if (length(top_metabolites_original) >= 2) {
  heatmap_meta <- metadata_primary %>%
    arrange(diagnosis_clean, sample_id_clean) %>%
    mutate(sample_order = row_number())
  
  heatmap_mat <- metabolites_log_scaled[
    heatmap_meta$sample_id_clean,
    top_metabolites_original,
    drop = FALSE
  ]
  
  heatmap_df <- as.data.frame(heatmap_mat) %>%
    rownames_to_column("sample_id_clean") %>%
    pivot_longer(
      cols = -sample_id_clean,
      names_to = "metabolite_id",
      values_to = "scaled_abundance"
    ) %>%
    left_join(
      heatmap_meta %>%
        select(sample_id_clean, diagnosis_plot, sample_order),
      by = "sample_id_clean"
    ) %>%
    mutate(
      metabolite_id = factor(metabolite_id, levels = rev(top_metabolites_original))
    )
  
  p_heatmap <- ggplot(
    heatmap_df,
    aes(x = sample_order, y = metabolite_id, fill = scaled_abundance)
  ) +
    geom_tile() +
    facet_grid(. ~ diagnosis_plot, scales = "free_x", space = "free_x") +
    labs(
      title = "Top disease-associated metabolite features",
      x = "Samples ordered by diagnosis",
      y = "Metabolite feature",
      fill = "Scaled abundance"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ggsave(
    "results/figures/top_metabolites_heatmap.png",
    p_heatmap,
    width = 11,
    height = 7,
    dpi = 300
  )
} else {
  message("Too few metabolites available for heatmap.")
}

message_header("Update analysis decisions")

decision_addition <- c(
  "",
  "## MaAsLin2 metabolite association modelling",
  "",
  "- Metabolite association models were performed using MaAsLin2.",
  "- Input metabolite features were half-minimum imputed positive abundance values from the primary independent set.",
  "- MaAsLin2 was run with LM analysis, no normalization, LOG transform, and BH FDR correction.",
  "- The primary independent set was used to avoid repeated-measure leakage.",
  "- Comparisons tested separately: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.",
  "- Candidate covariates were age, sex, and antibiotic use when sufficiently available.",
  "- Significance threshold was FDR < 0.10.",
  "- Significant metabolite features were treated as candidate disease-associated features for downstream integration, not as clinically validated biomarkers."
)

cat(
  paste(decision_addition, collapse = "\n"),
  file = "docs/analysis_decisions.md",
  append = TRUE
)

message_header("Maaslin2 metabolite analysis completed")

message("Total Maaslin2 metabolite results: ", nrow(metabolite_results))
message("FDR < 0.10 Maaslin2 metabolite results: ", nrow(metabolite_significant))

print(metabolite_significant %>% count(comparison))

message("Outputs written:")
message("- results/tables/metabolite_association_all_results.tsv")
message("- results/tables/metabolite_association_significant_results.tsv")
message("- results/tables/maaslin2_metabolite_feature_name_map.tsv")
message("- results/models/maaslin2_metabolites/")
message("- results/figures/metabolite_volcano.png")
message("- results/figures/top_metabolites_heatmap.png")