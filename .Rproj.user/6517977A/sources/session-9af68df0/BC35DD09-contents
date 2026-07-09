# Exploratory multi-omics integration using mixOmics multiblock sPLS-DA
#
# Inputs:
#   data/processed/integration_taxa_matrix.rds
#   data/processed/integration_metabolite_matrix.rds
#   data/processed/integration_metadata.rds
#   results/tables/taxa_metabolite_spearman_all_pairs.tsv
#
# Outputs:
#   results/tables/mixomics_selected_microbiome_features.tsv
#   results/tables/mixomics_selected_metabolite_features.tsv
#   results/tables/mixomics_model_performance.tsv
#   results/figures/mixomics_individuals_plot.png
#   results/figures/mixomics_microbiome_loading_plot.png
#   results/figures/mixomics_metabolite_loading_plot.png
#   results/figures/mixomics_circos_plot.png
#   results/models/mixomics_splsda_fit.rds

source("R/00_utils_io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(readr)
})

if (!requireNamespace("mixOmics", quietly = TRUE)) {
  stop(
    "Package `mixOmics` is not installed.\n",
    "Run this first:\n",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
    "BiocManager::install('mixOmics')"
  )
}

ensure_project_dirs()
ensure_dir("results/tables")
ensure_dir("results/figures")
ensure_dir("results/models")

message_header("Load selected integration features")

taxa_mat <- readRDS("data/processed/integration_taxa_matrix.rds")
metab_mat <- readRDS("data/processed/integration_metabolite_matrix.rds")
metadata <- readRDS("data/processed/integration_metadata.rds")

taxa_mat <- as.matrix(taxa_mat)
metab_mat <- as.matrix(metab_mat)

storage.mode(taxa_mat) <- "numeric"
storage.mode(metab_mat) <- "numeric"

common_ids <- Reduce(
  intersect,
  list(
    metadata$sample_id_clean,
    rownames(taxa_mat),
    rownames(metab_mat)
  )
)

metadata <- metadata %>%
  filter(sample_id_clean %in% common_ids) %>%
  arrange(sample_id_clean)

taxa_mat <- taxa_mat[metadata$sample_id_clean, , drop = FALSE]
metab_mat <- metab_mat[metadata$sample_id_clean, , drop = FALSE]

Y <- factor(metadata$diagnosis_clean, levels = c("non_IBD", "CD", "UC"))

message("Samples used for mixOmics: ", nrow(metadata))
message("Microbiome block features: ", ncol(taxa_mat))
message("Metabolome block features: ", ncol(metab_mat))
message("Diagnosis groups:")
print(table(Y))

if (nrow(metadata) < 30) {
  stop("Too few samples for sPLS-DA.")
}

if (ncol(taxa_mat) < 2 || ncol(metab_mat) < 2) {
  stop("Too few selected features for mixOmics integration.")
}

message_header("Run multiblock sPLS-DA")

X <- list(
  microbiome = taxa_mat,
  metabolome = metab_mat
)

design <- matrix(
  0.1,
  ncol = length(X),
  nrow = length(X),
  dimnames = list(names(X), names(X))
)
diag(design) <- 0

ncomp <- 2

keepX <- list(
  microbiome = rep(min(10, ncol(taxa_mat)), ncomp),
  metabolome = rep(min(25, ncol(metab_mat)), ncomp)
)

set.seed(20260707)

fit <- mixOmics::block.splsda(
  X = X,
  Y = Y,
  ncomp = ncomp,
  keepX = keepX,
  design = design,
  scale = TRUE,
  near.zero.var = TRUE
)

saveRDS(fit, "results/models/mixomics_splsda_fit.rds")

message_header("Extract selected features")

extract_loadings <- function(fit, block_name, feature_type) {
  loading_mat <- fit$loadings[[block_name]]
  
  if (is.null(loading_mat)) {
    return(tibble())
  }
  
  loading_df <- as.data.frame(loading_mat) %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(
      cols = -feature_id,
      names_to = "component",
      values_to = "loading"
    ) %>%
    mutate(
      feature_type = feature_type,
      abs_loading = abs(loading),
      selected_by_spls = abs_loading > 0
    ) %>%
    filter(selected_by_spls) %>%
    arrange(component, desc(abs_loading))
  
  loading_df
}

selected_microbiome <- extract_loadings(
  fit = fit,
  block_name = "microbiome",
  feature_type = "microbiome_taxon"
) %>%
  mutate(
    feature_label = feature_id %>%
      str_replace("^.*\\|g__", "") %>%
      str_replace("^g__", "")
  )

selected_metabolites <- extract_loadings(
  fit = fit,
  block_name = "metabolome",
  feature_type = "metabolite"
) %>%
  mutate(
    feature_label = str_trunc(feature_id, 60)
  )

write_tsv_safe(
  selected_microbiome,
  "results/tables/mixomics_selected_microbiome_features.tsv"
)

write_tsv_safe(
  selected_metabolites,
  "results/tables/mixomics_selected_metabolite_features.tsv"
)

message("sPLS-selected microbiome features: ", nrow(selected_microbiome))
message("sPLS-selected metabolite features: ", nrow(selected_metabolites))

message_header("Create individuals plot")

scores_micro <- as.data.frame(fit$variates$microbiome[, 1:2, drop = FALSE]) %>%
  rownames_to_column("sample_id_clean")

scores_metab <- as.data.frame(fit$variates$metabolome[, 1:2, drop = FALSE]) %>%
  rownames_to_column("sample_id_clean")

names(scores_micro)[2:3] <- c("micro_comp1", "micro_comp2")
names(scores_metab)[2:3] <- c("metab_comp1", "metab_comp2")

scores_df <- metadata %>%
  select(sample_id_clean, diagnosis_clean) %>%
  left_join(scores_micro, by = "sample_id_clean") %>%
  left_join(scores_metab, by = "sample_id_clean") %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD"),
    comp1 = rowMeans(cbind(scale(micro_comp1), scale(metab_comp1)), na.rm = TRUE),
    comp2 = rowMeans(cbind(scale(micro_comp2), scale(metab_comp2)), na.rm = TRUE)
  )

p_individuals <- ggplot(
  scores_df,
  aes(x = comp1, y = comp2, shape = diagnosis_plot)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    title = "mixOmics multiblock sPLS-DA sample scores",
    x = "Integrated component 1",
    y = "Integrated component 2",
    shape = "Diagnosis"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/mixomics_individuals_plot.png",
  p_individuals,
  width = 7,
  height = 5.5,
  dpi = 300
)

message_header("Estimate exploratory apparent classification performance")

centroids <- scores_df %>%
  group_by(diagnosis_clean) %>%
  summarise(
    centroid_comp1 = mean(comp1, na.rm = TRUE),
    centroid_comp2 = mean(comp2, na.rm = TRUE),
    .groups = "drop"
  )

predict_nearest_centroid <- function(x1, x2, centroids) {
  d <- (centroids$centroid_comp1 - x1)^2 + (centroids$centroid_comp2 - x2)^2
  centroids$diagnosis_clean[which.min(d)]
}

performance_df <- scores_df %>%
  rowwise() %>%
  mutate(
    predicted_group = predict_nearest_centroid(comp1, comp2, centroids)
  ) %>%
  ungroup()

apparent_accuracy <- mean(performance_df$predicted_group == performance_df$diagnosis_clean)

model_performance <- tibble(
  model = "mixOmics_block_sPLS_DA",
  performance_type = "apparent_nearest_centroid_accuracy_on_integrated_scores",
  n_samples = nrow(scores_df),
  n_components = ncomp,
  n_microbiome_features_input = ncol(taxa_mat),
  n_metabolite_features_input = ncol(metab_mat),
  keepX_microbiome_per_component = paste(keepX$microbiome, collapse = ";"),
  keepX_metabolome_per_component = paste(keepX$metabolome, collapse = ";"),
  apparent_accuracy = apparent_accuracy,
  interpretation_note = "Exploratory training-set metric only; not clinical validation."
)

write_tsv_safe(
  model_performance,
  "results/tables/mixomics_model_performance.tsv"
)

message_header("Create loading plots")

plot_loading <- function(df, output_path, title_text) {
  if (nrow(df) == 0) {
    message("No selected features for: ", title_text)
    return(invisible(NULL))
  }
  
  plot_df <- df %>%
    group_by(component) %>%
    arrange(desc(abs_loading), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup() %>%
    mutate(
      feature_label = str_trunc(feature_label, 45)
    )
  
  p <- ggplot(
    plot_df,
    aes(x = loading, y = reorder(feature_label, loading))
  ) +
    geom_col() +
    facet_wrap(~ component, scales = "free_y") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = title_text,
      x = "sPLS loading",
      y = "Selected feature"
    ) +
    theme_bw(base_size = 11)
  
  ggsave(
    output_path,
    p,
    width = 9,
    height = 7,
    dpi = 300
  )
}

plot_loading(
  selected_microbiome,
  "results/figures/mixomics_microbiome_loading_plot.png",
  "mixOmics selected microbiome features"
)

plot_loading(
  selected_metabolites,
  "results/figures/mixomics_metabolite_loading_plot.png",
  "mixOmics selected metabolite features"
)

message_header("Create exploratory circos-like cross-omics plot")

spearman_all <- readr::read_tsv(
  "results/tables/taxa_metabolite_spearman_all_pairs.tsv",
  show_col_types = FALSE
)

selected_micro_ids <- selected_microbiome %>%
  distinct(feature_id) %>%
  pull(feature_id)

selected_metab_ids <- selected_metabolites %>%
  distinct(feature_id) %>%
  pull(feature_id)

circos_pairs <- spearman_all %>%
  filter(
    taxon %in% selected_micro_ids,
    metabolite %in% selected_metab_ids
  ) %>%
  arrange(correlation_fdr, desc(abs_rho)) %>%
  slice_head(n = 80)

if (nrow(circos_pairs) > 0) {
  taxa_nodes <- circos_pairs %>%
    distinct(node = taxon, label = taxon_label) %>%
    arrange(label) %>%
    mutate(
      node_type = "microbiome",
      angle = seq(pi * 0.65, pi * 1.35, length.out = n())
    )
  
  metabolite_nodes <- circos_pairs %>%
    distinct(node = metabolite, label = metabolite_label) %>%
    arrange(label) %>%
    mutate(
      node_type = "metabolite",
      angle = seq(-pi * 0.35, pi * 0.35, length.out = n())
    )
  
  nodes <- bind_rows(taxa_nodes, metabolite_nodes) %>%
    mutate(
      x = cos(angle),
      y = sin(angle)
    )
  
  edges <- circos_pairs %>%
    left_join(
      nodes %>% select(taxon = node, x_taxon = x, y_taxon = y),
      by = "taxon"
    ) %>%
    left_join(
      nodes %>% select(metabolite = node, x_metab = x, y_metab = y),
      by = "metabolite"
    )
  
  p_circos <- ggplot() +
    geom_segment(
      data = edges,
      aes(
        x = x_taxon,
        y = y_taxon,
        xend = x_metab,
        yend = y_metab,
        linewidth = abs_rho,
        linetype = correlation_direction
      ),
      alpha = 0.45
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, shape = node_type),
      size = 2.6
    ) +
    geom_text(
      data = nodes,
      aes(x = x * 1.12, y = y * 1.12, label = str_trunc(label, 28)),
      size = 2.5
    ) +
    coord_equal() +
    labs(
      title = "Cross-omics links among sPLS-selected features",
      x = NULL,
      y = NULL,
      linewidth = "|rho|",
      linetype = "Direction",
      shape = "Feature type"
    ) +
    theme_void(base_size = 11)
  
  ggsave(
    "results/figures/mixomics_circos_plot.png",
    p_circos,
    width = 9,
    height = 9,
    dpi = 300
  )
} else {
  message("No overlapping sPLS-selected taxa-metabolite correlation pairs available for circos plot.")
  
  empty_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No overlapping sPLS-selected taxa-metabolite pairs",
      size = 5
    ) +
    theme_void()
  
  ggsave(
    "results/figures/mixomics_circos_plot.png",
    empty_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
}

message_header("Update analysis decisions")

decision_addition <- c(
  "",
  "## mixOmics multiblock sPLS-DA integration",
  "",
  "- Multiblock sPLS-DA was used as an exploratory supervised feature-selection approach.",
  "- The outcome was diagnosis group, and the two blocks were microbiome genus-level CLR features and log-scaled metabolite features.",
  "- Only the primary independent sample set was used to avoid repeated-measure leakage.",
  "- The design matrix used a low block connection weight of 0.1 to balance cross-omics correlation and diagnosis discrimination.",
  "- The sPLS model was used for exploratory feature selection, not as a validated clinical classifier.",
  "- Apparent classification accuracy was reported only as a training-set descriptive metric."
)

cat(
  paste(decision_addition, collapse = "\n"),
  file = "docs/analysis_decisions.md",
  append = TRUE
)

message_header("mixOmics integration completed")

message("Outputs written:")
message("- results/tables/mixomics_selected_microbiome_features.tsv")
message("- results/tables/mixomics_selected_metabolite_features.tsv")
message("- results/tables/mixomics_model_performance.tsv")
message("- results/figures/mixomics_individuals_plot.png")
message("- results/figures/mixomics_microbiome_loading_plot.png")
message("- results/figures/mixomics_metabolite_loading_plot.png")
message("- results/figures/mixomics_circos_plot.png")
message("- results/models/mixomics_splsda_fit.rds")