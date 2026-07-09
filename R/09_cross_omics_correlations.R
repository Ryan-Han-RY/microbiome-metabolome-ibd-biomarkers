# Cross-omics feature selection and taxa-metabolite correlation analysis


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

ensure_project_dirs()
ensure_dir("results/tables")
ensure_dir("results/figures")
ensure_dir("data/processed")

message_header("Load matched metadata and processed omics matrices")

metadata <- readRDS("data/processed/matched_metadata.rds")
taxa_clr <- readRDS("data/processed/taxa_genus_clr.rds")
taxa_rel <- readRDS("data/processed/taxa_genus_filtered.rds")
metab_scaled <- readRDS("data/processed/metabolites_log_scaled.rds")

taxa_clr <- as.matrix(taxa_clr)
taxa_rel <- as.matrix(taxa_rel)
metab_scaled <- as.matrix(metab_scaled)

storage.mode(taxa_clr) <- "numeric"
storage.mode(taxa_rel) <- "numeric"
storage.mode(metab_scaled) <- "numeric"

primary_meta <- metadata %>%
  filter(is_primary_independent) %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD")
  )

common_ids <- Reduce(
  intersect,
  list(
    primary_meta$sample_id_clean,
    rownames(taxa_clr),
    rownames(taxa_rel),
    rownames(metab_scaled)
  )
)

primary_meta <- primary_meta %>%
  filter(sample_id_clean %in% common_ids) %>%
  arrange(sample_id_clean)

taxa_clr <- taxa_clr[primary_meta$sample_id_clean, , drop = FALSE]
taxa_rel <- taxa_rel[primary_meta$sample_id_clean, , drop = FALSE]
metab_scaled <- metab_scaled[primary_meta$sample_id_clean, , drop = FALSE]

message("Primary independent samples used: ", nrow(primary_meta))
message("Available genus-level taxa: ", ncol(taxa_clr))
message("Available metabolite features: ", ncol(metab_scaled))

message_header("Load disease association results")

taxa_all <- readr::read_tsv(
  "results/tables/maaslin2_taxa_all_results.tsv",
  show_col_types = FALSE
)

taxa_sig <- readr::read_tsv(
  "results/tables/maaslin2_taxa_significant_results.tsv",
  show_col_types = FALSE
)

metab_all <- readr::read_tsv(
  "results/tables/metabolite_association_all_results.tsv",
  show_col_types = FALSE
)

metab_sig <- readr::read_tsv(
  "results/tables/metabolite_association_significant_results.tsv",
  show_col_types = FALSE
)

message("Taxa association results: ", nrow(taxa_all))
message("FDR-significant taxa results: ", nrow(taxa_sig))
message("Metabolite association results: ", nrow(metab_all))
message("FDR-significant metabolite results: ", nrow(metab_sig))

message_header("Select integration features")

clean_taxon_label <- function(x) {
  x %>%
    str_replace("^.*\\|g__", "") %>%
    str_replace("^g__", "") %>%
    str_replace_all("_", " ")
}

clean_metabolite_label <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_trunc(50)
}

taxa_ranked <- taxa_all %>%
  mutate(
    original_feature_id = if_else(
      is.na(original_feature_id),
      feature_id,
      original_feature_id
    ),
    abs_beta = abs(beta),
    q_value_safe = if_else(is.na(q_value), 1, q_value),
    p_value_safe = if_else(is.na(p_value), 1, p_value)
  ) %>%
  group_by(original_feature_id) %>%
  arrange(q_value_safe, p_value_safe, desc(abs_beta), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    selection_source = if_else(q_value_safe < 0.10, "FDR_lt_0.10", "top_by_effect_size"),
    taxon_label = clean_taxon_label(original_feature_id)
  ) %>%
  arrange(q_value_safe, p_value_safe, desc(abs_beta))

# Your MaAsLin2 taxa FDR result is currently empty, so this intentionally keeps top taxa by effect size.
if (any(taxa_ranked$q_value_safe < 0.10, na.rm = TRUE)) {
  selected_taxa <- taxa_ranked %>%
    filter(q_value_safe < 0.10)
} else {
  selected_taxa <- taxa_ranked %>%
    arrange(desc(abs_beta), p_value_safe) %>%
    slice_head(n = 50)
}

selected_taxa <- selected_taxa %>%
  filter(original_feature_id %in% colnames(taxa_clr)) %>%
  distinct(original_feature_id, .keep_all = TRUE)

metab_ranked <- metab_all %>%
  mutate(
    original_metabolite_id = if_else(
      is.na(original_metabolite_id),
      metabolite_id,
      original_metabolite_id
    ),
    abs_beta = abs(beta),
    q_value_safe = if_else(is.na(q_value), 1, q_value),
    p_value_safe = if_else(is.na(p_value), 1, p_value)
  ) %>%
  group_by(original_metabolite_id) %>%
  arrange(q_value_safe, p_value_safe, desc(abs_beta), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    selection_source = if_else(q_value_safe < 0.10, "FDR_lt_0.10", "top_by_effect_size"),
    metabolite_label = clean_metabolite_label(original_metabolite_id)
  ) %>%
  arrange(q_value_safe, p_value_safe, desc(abs_beta))

if (nrow(metab_sig) > 0) {
  selected_metabolites <- metab_ranked %>%
    filter(q_value_safe < 0.10) %>%
    slice_head(n = 100)
} else {
  selected_metabolites <- metab_ranked %>%
    arrange(desc(abs_beta), p_value_safe) %>%
    slice_head(n = 100)
}

selected_metabolites <- selected_metabolites %>%
  filter(original_metabolite_id %in% colnames(metab_scaled)) %>%
  distinct(original_metabolite_id, .keep_all = TRUE)

message("Selected taxa for integration: ", nrow(selected_taxa))
message("Selected metabolites for integration: ", nrow(selected_metabolites))

taxa_integration_matrix <- taxa_clr[, selected_taxa$original_feature_id, drop = FALSE]
metab_integration_matrix <- metab_scaled[, selected_metabolites$original_metabolite_id, drop = FALSE]

saveRDS(taxa_integration_matrix, "data/processed/integration_taxa_matrix.rds")
saveRDS(metab_integration_matrix, "data/processed/integration_metabolite_matrix.rds")
saveRDS(primary_meta, "data/processed/integration_metadata.rds")

write_tsv_safe(
  selected_taxa,
  "results/tables/integration_selected_taxa.tsv"
)

write_tsv_safe(
  selected_metabolites,
  "results/tables/integration_selected_metabolites.tsv"
)

integration_feature_selection_summary <- tibble(
  feature_type = c("microbiome_taxa", "metabolites"),
  n_total_features = c(ncol(taxa_clr), ncol(metab_scaled)),
  n_after_filtering = c(nrow(selected_taxa), nrow(selected_metabolites)),
  selection_criterion = c(
    if_else(any(taxa_ranked$q_value_safe < 0.10, na.rm = TRUE),
            "MaAsLin2 FDR < 0.10",
            "Top taxa by absolute MaAsLin2 effect size because no taxa survived FDR < 0.10"),
    if_else(nrow(metab_sig) > 0,
            "MaAsLin2 FDR < 0.10, capped at top 100 by q-value/effect size",
            "Top metabolites by absolute MaAsLin2 effect size because no metabolites survived FDR < 0.10")
  ),
  used_for_correlation = c(TRUE, TRUE),
  used_for_spls = c(TRUE, TRUE)
)

write_tsv_safe(
  integration_feature_selection_summary,
  "results/tables/integration_feature_selection_summary.tsv"
)

message_header("Spearman taxa-metabolite correlation")

run_spearman_pair <- function(taxon_id, metabolite_id, taxa_mat, metab_mat) {
  x <- taxa_mat[, taxon_id]
  y <- metab_mat[, metabolite_id]
  
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 10 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(
      taxon = taxon_id,
      metabolite = metabolite_id,
      spearman_rho = NA_real_,
      p_value = NA_real_,
      n_samples = sum(ok)
    ))
  }
  
  test <- suppressWarnings(
    cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  )
  
  tibble(
    taxon = taxon_id,
    metabolite = metabolite_id,
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    n_samples = sum(ok)
  )
}

pair_grid <- tidyr::expand_grid(
  taxon = colnames(taxa_integration_matrix),
  metabolite = colnames(metab_integration_matrix)
)

spearman_results <- purrr::pmap_dfr(
  pair_grid,
  function(taxon, metabolite) {
    run_spearman_pair(
      taxon_id = taxon,
      metabolite_id = metabolite,
      taxa_mat = taxa_integration_matrix,
      metab_mat = metab_integration_matrix
    )
  }
) %>%
  mutate(
    correlation_fdr = p.adjust(p_value, method = "BH"),
    abs_rho = abs(spearman_rho),
    correlation_direction = case_when(
      spearman_rho > 0 ~ "positive",
      spearman_rho < 0 ~ "negative",
      TRUE ~ "none"
    ),
    taxon_label = clean_taxon_label(taxon),
    metabolite_label = clean_metabolite_label(metabolite)
  ) %>%
  arrange(correlation_fdr, desc(abs_rho), p_value)

spearman_significant <- spearman_results %>%
  filter(
    !is.na(correlation_fdr),
    abs_rho >= 0.30,
    correlation_fdr < 0.10
  ) %>%
  arrange(correlation_fdr, desc(abs_rho))

write_tsv_safe(
  spearman_results,
  "results/tables/taxa_metabolite_spearman_all_pairs.tsv"
)

write_tsv_safe(
  spearman_significant,
  "results/tables/taxa_metabolite_spearman_significant_pairs.tsv"
)

message("All taxa-metabolite pairs tested: ", nrow(spearman_results))
message("Significant taxa-metabolite pairs: ", nrow(spearman_significant))

message_header("Create correlation heatmap")

heatmap_pairs <- if (nrow(spearman_significant) > 0) {
  spearman_significant %>%
    slice_head(n = 100)
} else {
  spearman_results %>%
    arrange(desc(abs_rho), p_value) %>%
    slice_head(n = 100)
}

heatmap_taxa <- heatmap_pairs %>%
  distinct(taxon, taxon_label) %>%
  pull(taxon)

heatmap_metabolites <- heatmap_pairs %>%
  distinct(metabolite, metabolite_label) %>%
  pull(metabolite)

heatmap_df <- spearman_results %>%
  filter(taxon %in% heatmap_taxa, metabolite %in% heatmap_metabolites) %>%
  mutate(
    taxon_label = factor(taxon_label, levels = unique(clean_taxon_label(heatmap_taxa))),
    metabolite_label = factor(metabolite_label, levels = rev(unique(clean_metabolite_label(heatmap_metabolites))))
  )

p_heatmap <- ggplot(
  heatmap_df,
  aes(x = taxon_label, y = metabolite_label, fill = spearman_rho)
) +
  geom_tile() +
  labs(
    title = "Taxa-metabolite Spearman correlation",
    x = "Genus-level taxon",
    y = "Metabolite feature",
    fill = "Spearman rho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "results/figures/taxa_metabolite_correlation_heatmap.png",
  p_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)

message_header("Create correlation network plot")

network_pairs <- if (nrow(spearman_significant) > 0) {
  spearman_significant %>%
    arrange(correlation_fdr, desc(abs_rho)) %>%
    slice_head(n = 60)
} else {
  spearman_results %>%
    arrange(desc(abs_rho), p_value) %>%
    slice_head(n = 60)
}

taxa_nodes <- network_pairs %>%
  distinct(node = taxon, label = taxon_label) %>%
  arrange(label) %>%
  mutate(
    node_type = "taxon",
    x = 0,
    y = row_number()
  )

metabolite_nodes <- network_pairs %>%
  distinct(node = metabolite, label = metabolite_label) %>%
  arrange(label) %>%
  mutate(
    node_type = "metabolite",
    x = 1,
    y = row_number()
  )

nodes <- bind_rows(taxa_nodes, metabolite_nodes)

edges <- network_pairs %>%
  left_join(
    taxa_nodes %>% select(taxon = node, taxon_y = y),
    by = "taxon"
  ) %>%
  left_join(
    metabolite_nodes %>% select(metabolite = node, metabolite_y = y),
    by = "metabolite"
  )

p_network <- ggplot() +
  geom_segment(
    data = edges,
    aes(
      x = 0,
      xend = 1,
      y = taxon_y,
      yend = metabolite_y,
      linewidth = abs_rho,
      linetype = correlation_direction
    ),
    alpha = 0.55
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y, shape = node_type),
    size = 2.5
  ) +
  geom_text(
    data = nodes %>% filter(node_type == "taxon"),
    aes(x = x - 0.03, y = y, label = label),
    hjust = 1,
    size = 3
  ) +
  geom_text(
    data = nodes %>% filter(node_type == "metabolite"),
    aes(x = x + 0.03, y = y, label = label),
    hjust = 0,
    size = 2.5
  ) +
  scale_x_continuous(limits = c(-0.45, 1.45), breaks = c(0, 1), labels = c("Taxa", "Metabolites")) +
  labs(
    title = "Taxa-metabolite correlation network",
    x = NULL,
    y = NULL,
    linewidth = "|rho|",
    linetype = "Direction",
    shape = "Feature type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  "results/figures/taxa_metabolite_correlation_network.png",
  p_network,
  width = 12,
  height = 8,
  dpi = 300
)

message_header("Covariate-adjusted correlation by residualisation")

clean_covariate <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

adjust_meta <- primary_meta %>%
  transmute(
    sample_id_clean,
    age_model = age_clean,
    sex_model = clean_covariate(sex_clean),
    antibiotic_model = clean_covariate(antibiotic_use_clean)
  )

select_covariates <- function(meta) {
  covariates <- c()
  
  if (mean(!is.na(meta$age_model)) >= 0.80 &&
      n_distinct(meta$age_model, na.rm = TRUE) > 5) {
    covariates <- c(covariates, "age_model")
  }
  
  if (mean(!is.na(meta$sex_model)) >= 0.80 &&
      n_distinct(meta$sex_model, na.rm = TRUE) >= 2) {
    covariates <- c(covariates, "sex_model")
  }
  
  if (mean(!is.na(meta$antibiotic_model)) >= 0.80 &&
      n_distinct(meta$antibiotic_model, na.rm = TRUE) >= 2) {
    covariates <- c(covariates, "antibiotic_model")
  }
  
  covariates
}

covariates_used <- select_covariates(adjust_meta)

message(
  "Covariates used for residualisation: ",
  ifelse(length(covariates_used) == 0, "none", paste(covariates_used, collapse = ", "))
)

if (length(covariates_used) > 0) {
  complete_adjust_ids <- adjust_meta %>%
    filter(if_all(all_of(covariates_used), ~ !is.na(.x))) %>%
    pull(sample_id_clean)
} else {
  complete_adjust_ids <- adjust_meta$sample_id_clean
}

adjust_meta <- adjust_meta %>%
  filter(sample_id_clean %in% complete_adjust_ids) %>%
  arrange(sample_id_clean)

taxa_adjust_mat <- taxa_integration_matrix[adjust_meta$sample_id_clean, , drop = FALSE]
metab_adjust_mat <- metab_integration_matrix[adjust_meta$sample_id_clean, , drop = FALSE]

residualize_matrix <- function(mat, meta, covariates) {
  mat <- as.matrix(mat)
  
  if (length(covariates) == 0) {
    return(scale(mat))
  }
  
  meta_df <- meta %>%
    select(all_of(covariates)) %>%
    mutate(across(where(is.character), factor))
  
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  formula_text <- paste("y ~", paste(covariates, collapse = " + "))
  
  for (j in seq_len(ncol(mat))) {
    df <- meta_df %>%
      mutate(y = mat[, j])
    
    fit <- tryCatch(
      lm(as.formula(formula_text), data = df),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      out[, j] <- scale(mat[, j])[, 1]
    } else {
      out[, j] <- residuals(fit)
    }
  }
  
  out
}

taxa_residuals <- residualize_matrix(taxa_adjust_mat, adjust_meta, covariates_used)
metab_residuals <- residualize_matrix(metab_adjust_mat, adjust_meta, covariates_used)

adjusted_pair_grid <- tidyr::expand_grid(
  taxon = colnames(taxa_residuals),
  metabolite = colnames(metab_residuals)
)

adjusted_results <- purrr::pmap_dfr(
  adjusted_pair_grid,
  function(taxon, metabolite) {
    run_spearman_pair(
      taxon_id = taxon,
      metabolite_id = metabolite,
      taxa_mat = taxa_residuals,
      metab_mat = metab_residuals
    )
  }
) %>%
  mutate(
    adjusted_correlation_fdr = p.adjust(p_value, method = "BH"),
    abs_adjusted_rho = abs(spearman_rho),
    adjusted_direction = case_when(
      spearman_rho > 0 ~ "positive",
      spearman_rho < 0 ~ "negative",
      TRUE ~ "none"
    ),
    covariates_used = ifelse(length(covariates_used) == 0, "none", paste(covariates_used, collapse = ";")),
    taxon_label = clean_taxon_label(taxon),
    metabolite_label = clean_metabolite_label(metabolite)
  ) %>%
  arrange(adjusted_correlation_fdr, desc(abs_adjusted_rho), p_value)

write_tsv_safe(
  adjusted_results,
  "results/tables/covariate_adjusted_correlation_pairs.tsv"
)

adjusted_heatmap_pairs <- adjusted_results %>%
  arrange(adjusted_correlation_fdr, desc(abs_adjusted_rho)) %>%
  slice_head(n = 100)

adjusted_heatmap_taxa <- adjusted_heatmap_pairs %>%
  distinct(taxon) %>%
  pull(taxon)

adjusted_heatmap_metabolites <- adjusted_heatmap_pairs %>%
  distinct(metabolite) %>%
  pull(metabolite)

adjusted_heatmap_df <- adjusted_results %>%
  filter(
    taxon %in% adjusted_heatmap_taxa,
    metabolite %in% adjusted_heatmap_metabolites
  ) %>%
  mutate(
    taxon_label = factor(taxon_label, levels = unique(clean_taxon_label(adjusted_heatmap_taxa))),
    metabolite_label = factor(metabolite_label, levels = rev(unique(clean_metabolite_label(adjusted_heatmap_metabolites))))
  )

p_adjusted_heatmap <- ggplot(
  adjusted_heatmap_df,
  aes(x = taxon_label, y = metabolite_label, fill = spearman_rho)
) +
  geom_tile() +
  labs(
    title = "Covariate-adjusted taxa-metabolite correlation",
    x = "Genus-level taxon",
    y = "Metabolite feature",
    fill = "Residual Spearman rho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "results/figures/covariate_adjusted_correlation_heatmap.png",
  p_adjusted_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)

message_header("Update analysis decisions")

decision_addition <- c(
  "",
  "## Cross-omics correlation and integration feature selection",
  "",
  "- Cross-omics correlation was restricted to selected disease-associated or top-ranked features rather than all taxa-metabolite combinations.",
  "- Because no genus-level taxa survived MaAsLin2 FDR < 0.10, taxa were selected by top absolute MaAsLin2 effect size and nominal association strength.",
  "- Metabolite features were selected from MaAsLin2 FDR-significant disease-associated features, capped at the top 100 to keep integration interpretable.",
  "- Spearman correlation was used between genus-level CLR abundance and log-scaled metabolite abundance.",
  "- Correlation p-values were corrected using Benjamini-Hochberg FDR.",
  "- Strong correlation pairs were defined as absolute Spearman rho >= 0.30 and FDR < 0.10.",
  "- Covariate-adjusted correlation was performed by residualising taxa and metabolite features against available covariates before correlation."
)

cat(
  paste(decision_addition, collapse = "\n"),
  file = "docs/analysis_decisions.md",
  append = TRUE
)

message_header("Cross-omics correlation completed")

message("Outputs written:")
message("- results/tables/integration_feature_selection_summary.tsv")
message("- results/tables/integration_selected_taxa.tsv")
message("- results/tables/integration_selected_metabolites.tsv")
message("- results/tables/taxa_metabolite_spearman_all_pairs.tsv")
message("- results/tables/taxa_metabolite_spearman_significant_pairs.tsv")
message("- results/tables/covariate_adjusted_correlation_pairs.tsv")
message("- results/figures/taxa_metabolite_correlation_heatmap.png")
message("- results/figures/taxa_metabolite_correlation_network.png")
message("- results/figures/covariate_adjusted_correlation_heatmap.png")
