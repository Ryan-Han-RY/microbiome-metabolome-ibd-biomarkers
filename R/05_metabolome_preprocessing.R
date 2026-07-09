# Metabolomics preprocessing: missingness filtering, half-minimum imputation, log2 transform, z-score scaling

source("R/00_utils_io.R")

ensure_project_dirs()

message_header("Load matched metabolomics data")

metadata <- readRDS("data/processed/matched_metadata.rds")
metab_mat <- readRDS("data/processed/metabolomics_matched_raw.rds")

metadata <- metadata %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD")
  )

metab_mat <- as.matrix(metab_mat)
storage.mode(metab_mat) <- "numeric"

metab_mat[!is.finite(metab_mat)] <- NA
metab_mat[metab_mat <= 0] <- NA

message("Input metabolomics samples: ", nrow(metab_mat))
message("Input metabolomics features: ", ncol(metab_mat))

message_header("Assess metabolite missingness")

missingness <- colMeans(is.na(metab_mat))

missingness_table <- tibble(
  metabolite_id = names(missingness),
  missingness = as.numeric(missingness)
) %>%
  arrange(desc(missingness))

keep_missingness <- missingness <= 0.30
metab_filtered_missing <- metab_mat[, keep_missingness, drop = FALSE]

message("Features retained after missingness filter: ", ncol(metab_filtered_missing))

half_min_impute <- function(mat) {
  mat <- as.matrix(mat)
  
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    positive <- x[!is.na(x) & x > 0]
    
    if (length(positive) == 0) {
      mat[, j] <- NA_real_
    } else {
      replacement <- min(positive, na.rm = TRUE) / 2
      x[is.na(x)] <- replacement
      mat[, j] <- x
    }
  }
  
  mat
}

metab_imputed <- half_min_impute(metab_filtered_missing)

keep_non_missing <- colSums(is.na(metab_imputed)) == 0
metab_imputed <- metab_imputed[, keep_non_missing, drop = FALSE]

metab_log <- log2(metab_imputed)

feature_sd <- apply(metab_log, 2, sd, na.rm = TRUE)
keep_variance <- is.finite(feature_sd) & feature_sd > 1e-8

metab_log <- metab_log[, keep_variance, drop = FALSE]
metab_scaled <- scale(metab_log)

saveRDS(metab_imputed, "data/processed/metabolites_filtered.rds")
saveRDS(metab_scaled, "data/processed/metabolites_log_scaled.rds")

metabolite_filtering_summary <- tibble(
  step = c(
    "raw_metabolite_features",
    "after_missingness_filter",
    "after_complete_imputation",
    "after_near_zero_variance_filter",
    "final_log_scaled_features"
  ),
  n_features = c(
    ncol(metab_mat),
    ncol(metab_filtered_missing),
    ncol(metab_imputed),
    sum(keep_variance),
    ncol(metab_scaled)
  ),
  n_samples = c(
    nrow(metab_mat),
    nrow(metab_filtered_missing),
    nrow(metab_imputed),
    nrow(metab_log),
    nrow(metab_scaled)
  )
)

write_tsv_safe(
  metabolite_filtering_summary,
  "results/tables/metabolite_filtering_summary.tsv"
)

write_tsv_safe(
  missingness_table,
  "results/tables/metabolite_missingness_summary.tsv"
)

message_header("Create metabolite missingness figure")

p_missing <- ggplot(missingness_table, aes(x = missingness)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = 0.30, linetype = "dashed") +
  labs(
    title = "Metabolite missingness distribution",
    x = "Missingness proportion",
    y = "Number of metabolites"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/metabolite_missingness_distribution.png",
  p_missing,
  width = 7,
  height = 5,
  dpi = 300
)

message_header("Create metabolomics PCA QC figure")

primary_ids <- metadata %>%
  filter(is_primary_independent) %>%
  pull(sample_id_clean)

primary_ids <- intersect(primary_ids, rownames(metab_scaled))

pca_mat <- metab_scaled[primary_ids, , drop = FALSE]
pca_meta <- metadata %>%
  filter(sample_id_clean %in% primary_ids) %>%
  arrange(match(sample_id_clean, primary_ids))

pca <- prcomp(pca_mat, center = FALSE, scale. = FALSE)

pca_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- tibble(
  sample_id_clean = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(
    pca_meta %>% select(sample_id_clean, diagnosis_clean),
    by = "sample_id_clean"
  ) %>%
  mutate(diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, shape = diagnosis_plot)) +
  geom_point(size = 2.2, alpha = 0.85) +
  labs(
    title = "Metabolomics PCA QC by diagnosis",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)"),
    shape = "Diagnosis"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/metabolomics_pca_qc.png",
  p_pca,
  width = 7,
  height = 5.5,
  dpi = 300
)

message_header("Metabolomics preprocessing completed")

print(metabolite_filtering_summary)

message("Outputs written:")
message("- data/processed/metabolites_filtered.rds")
message("- data/processed/metabolites_log_scaled.rds")
message("- results/tables/metabolite_filtering_summary.tsv")
message("- results/tables/metabolite_missingness_summary.tsv")
message("- results/figures/metabolite_missingness_distribution.png")
message("- results/figures/metabolomics_pca_qc.png")