# Final multi-evidence candidate biomarker ranking
#
# Inputs:
#   results/tables/maaslin2_taxa_all_results.tsv
#   results/tables/metabolite_association_all_results.tsv
#   results/tables/taxa_metabolite_spearman_all_pairs.tsv
#   results/tables/taxa_metabolite_spearman_significant_pairs.tsv
#   results/tables/covariate_adjusted_correlation_pairs.tsv
#   results/tables/mixomics_selected_microbiome_features.tsv
#   results/tables/mixomics_selected_metabolite_features.tsv
#   data/processed/taxa_genus_filtered.rds
#   data/processed/metabolites_filtered.rds
#
# Outputs:
#   results/tables/final_candidate_biomarker_ranking.tsv
#   results/figures/top_candidate_biomarker_dotplot.png
#   results/figures/top_candidate_biomarker_heatmap.png

source("R/00_utils_io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(readr)
})

ensure_project_dirs()
ensure_dir("results/tables")
ensure_dir("results/figures")

message_header("Load evidence tables")

taxa_all <- readr::read_tsv(
  "results/tables/maaslin2_taxa_all_results.tsv",
  show_col_types = FALSE
)

metab_all <- readr::read_tsv(
  "results/tables/metabolite_association_all_results.tsv",
  show_col_types = FALSE
)

spearman_all <- readr::read_tsv(
  "results/tables/taxa_metabolite_spearman_all_pairs.tsv",
  show_col_types = FALSE
)

spearman_sig <- readr::read_tsv(
  "results/tables/taxa_metabolite_spearman_significant_pairs.tsv",
  show_col_types = FALSE
)

adjusted_pairs <- readr::read_tsv(
  "results/tables/covariate_adjusted_correlation_pairs.tsv",
  show_col_types = FALSE
)

mixomics_taxa <- readr::read_tsv(
  "results/tables/mixomics_selected_microbiome_features.tsv",
  show_col_types = FALSE
)

mixomics_metabolites <- readr::read_tsv(
  "results/tables/mixomics_selected_metabolite_features.tsv",
  show_col_types = FALSE
)

taxa_rel <- readRDS("data/processed/taxa_genus_filtered.rds")
metab_filtered <- readRDS("data/processed/metabolites_filtered.rds")

taxa_rel <- as.matrix(taxa_rel)
metab_filtered <- as.matrix(metab_filtered)

storage.mode(taxa_rel) <- "numeric"
storage.mode(metab_filtered) <- "numeric"

message("Taxa association rows: ", nrow(taxa_all))
message("Metabolite association rows: ", nrow(metab_all))
message("Spearman all pairs: ", nrow(spearman_all))
message("Spearman significant pairs: ", nrow(spearman_sig))

message_header("Prepare feature-level evidence")

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
    str_trunc(60)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(abs(x))
}

pick_direction <- function(beta) {
  beta <- beta[is.finite(beta)]
  if (length(beta) == 0) return("unknown")
  b <- beta[which.max(abs(beta))]
  if (b > 0) "higher_in_case_group" else if (b < 0) "lower_in_case_group" else "no_direction"
}

taxa_evidence <- taxa_all %>%
  mutate(
    original_feature_id = if_else(
      is.na(original_feature_id),
      feature_id,
      original_feature_id
    )
  ) %>%
  group_by(original_feature_id) %>%
  summarise(
    taxon_fdr = safe_min(q_value),
    taxon_min_p = safe_min(p_value),
    taxon_effect_size = safe_max_abs(beta),
    taxon_effect_direction = pick_direction(beta),
    taxon_best_comparison = comparison[which.min(if_else(is.na(q_value), 1, q_value))][1],
    .groups = "drop"
  ) %>%
  rename(taxon = original_feature_id) %>%
  mutate(
    taxon_label = clean_taxon_label(taxon)
  )

metab_evidence <- metab_all %>%
  mutate(
    original_metabolite_id = if_else(
      is.na(original_metabolite_id),
      metabolite_id,
      original_metabolite_id
    )
  ) %>%
  group_by(original_metabolite_id) %>%
  summarise(
    metabolite_fdr = safe_min(q_value),
    metabolite_min_p = safe_min(p_value),
    metabolite_effect_size = safe_max_abs(beta),
    metabolite_effect_direction = pick_direction(beta),
    metabolite_best_comparison = comparison[which.min(if_else(is.na(q_value), 1, q_value))][1],
    .groups = "drop"
  ) %>%
  rename(metabolite = original_metabolite_id) %>%
  mutate(
    metabolite_label = clean_metabolite_label(metabolite)
  )

taxa_prevalence <- tibble(
  taxon = colnames(taxa_rel),
  prevalence = colMeans(taxa_rel > 0, na.rm = TRUE)
)

metabolite_missingness <- tibble(
  metabolite = colnames(metab_filtered),
  missingness = colMeans(is.na(metab_filtered) | !is.finite(metab_filtered))
)

mixomics_taxa_score <- mixomics_taxa %>%
  group_by(feature_id) %>%
  summarise(
    taxon_spls_component = paste(unique(component), collapse = ";"),
    taxon_spls_loading = max(abs_loading, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(taxon = feature_id)

mixomics_metab_score <- mixomics_metabolites %>%
  group_by(feature_id) %>%
  summarise(
    metabolite_spls_component = paste(unique(component), collapse = ";"),
    metabolite_spls_loading = max(abs_loading, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(metabolite = feature_id)

adjusted_best <- adjusted_pairs %>%
  group_by(taxon, metabolite) %>%
  arrange(adjusted_correlation_fdr, desc(abs_adjusted_rho), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    taxon,
    metabolite,
    adjusted_spearman_rho = spearman_rho,
    adjusted_correlation_fdr = adjusted_correlation_fdr
  )

message_header("Build candidate pair table")

candidate_pairs <- if (nrow(spearman_sig) > 0) {
  spearman_sig
} else {
  spearman_all %>%
    arrange(correlation_fdr, desc(abs_rho), p_value) %>%
    slice_head(n = 300)
}

ranking_base <- candidate_pairs %>%
  select(
    taxon,
    metabolite,
    spearman_rho,
    correlation_fdr,
    abs_rho,
    correlation_direction
  ) %>%
  left_join(taxa_evidence, by = "taxon") %>%
  left_join(metab_evidence, by = "metabolite") %>%
  left_join(taxa_prevalence, by = "taxon") %>%
  left_join(metabolite_missingness, by = "metabolite") %>%
  left_join(mixomics_taxa_score, by = "taxon") %>%
  left_join(mixomics_metab_score, by = "metabolite") %>%
  left_join(adjusted_best, by = c("taxon", "metabolite")) %>%
  mutate(
    taxon_spls_loading = if_else(is.na(taxon_spls_loading), 0, taxon_spls_loading),
    metabolite_spls_loading = if_else(is.na(metabolite_spls_loading), 0, metabolite_spls_loading),
    spls_loading = pmax(taxon_spls_loading, metabolite_spls_loading, na.rm = TRUE),
    spls_component = paste(
      if_else(is.na(taxon_spls_component), "taxon_not_selected", taxon_spls_component),
      if_else(is.na(metabolite_spls_component), "metabolite_not_selected", metabolite_spls_component),
      sep = "|"
    ),
    prevalence = if_else(is.na(prevalence), 0, prevalence),
    missingness = if_else(is.na(missingness), 0, missingness),
    taxon_fdr = if_else(is.na(taxon_fdr), 1, taxon_fdr),
    metabolite_fdr = if_else(is.na(metabolite_fdr), 1, metabolite_fdr),
    correlation_fdr = if_else(is.na(correlation_fdr), 1, correlation_fdr),
    taxon_effect_size = if_else(is.na(taxon_effect_size), 0, taxon_effect_size),
    metabolite_effect_size = if_else(is.na(metabolite_effect_size), 0, metabolite_effect_size),
    abs_rho = if_else(is.na(abs_rho), 0, abs_rho)
  )

scale01 <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  
  if (!is.finite(min_x) || !is.finite(max_x) || max_x == min_x) {
    return(rep(0.5, length(x)))
  }
  
  (x - min_x) / (max_x - min_x)
}

ranking_scored <- ranking_base %>%
  mutate(
    taxon_effect_score = scale01(taxon_effect_size),
    taxon_fdr_score = scale01(-log10(pmax(taxon_fdr, 1e-300))),
    metabolite_effect_score = scale01(metabolite_effect_size),
    metabolite_fdr_score = scale01(-log10(pmax(metabolite_fdr, 1e-300))),
    correlation_score = scale01(abs_rho),
    correlation_fdr_score = scale01(-log10(pmax(correlation_fdr, 1e-300))),
    spls_score = scale01(spls_loading),
    prevalence_score = prevalence,
    missingness_score = 1 - missingness,
    adjusted_correlation_score = scale01(abs(adjusted_spearman_rho)),
    overall_score =
      0.12 * taxon_effect_score +
      0.08 * taxon_fdr_score +
      0.15 * metabolite_effect_score +
      0.15 * metabolite_fdr_score +
      0.15 * correlation_score +
      0.10 * correlation_fdr_score +
      0.10 * spls_score +
      0.05 * prevalence_score +
      0.05 * missingness_score +
      0.05 * adjusted_correlation_score,
    biological_note = case_when(
      correlation_fdr < 0.10 & abs_rho >= 0.30 & spls_loading > 0 &
        metabolite_fdr < 0.10 ~
        "High-priority exploratory pair supported by metabolite disease association, cross-omics correlation, and sPLS selection; requires external validation.",
      correlation_fdr < 0.10 & abs_rho >= 0.30 & metabolite_fdr < 0.10 ~
        "Disease-associated metabolite paired with correlated microbial feature; candidate for downstream biological interpretation.",
      spls_loading > 0 ~
        "Selected by exploratory sPLS integration but requires stronger univariate/correlation support.",
      TRUE ~
        "Exploratory candidate retained by multi-evidence ranking; not clinically validated."
    )
  ) %>%
  arrange(desc(overall_score), correlation_fdr, desc(abs_rho)) %>%
  mutate(rank = row_number()) %>%
  transmute(
    rank,
    taxon,
    metabolite,
    taxon_label,
    metabolite_label,
    taxon_effect_direction,
    metabolite_effect_direction,
    taxon_fdr,
    metabolite_fdr,
    spearman_rho,
    correlation_fdr,
    adjusted_spearman_rho,
    adjusted_correlation_fdr,
    spls_component,
    spls_loading,
    prevalence,
    missingness,
    biological_note,
    overall_score,
    taxon_effect_size,
    metabolite_effect_size,
    taxon_best_comparison,
    metabolite_best_comparison
  )

write_tsv_safe(
  ranking_scored,
  "results/tables/final_candidate_biomarker_ranking.tsv"
)

message("Final candidate pairs ranked: ", nrow(ranking_scored))

message_header("Create top candidate biomarker dotplot")

top_dot <- ranking_scored %>%
  slice_head(n = 25) %>%
  mutate(
    pair_label = paste0(
      str_trunc(taxon_label, 25),
      " | ",
      str_trunc(metabolite_label, 35)
    ),
    pair_label = factor(pair_label, levels = rev(pair_label))
  )

p_dot <- ggplot(
  top_dot,
  aes(x = overall_score, y = pair_label)
) +
  geom_point(aes(size = abs(spearman_rho), shape = metabolite_fdr < 0.10), alpha = 0.85) +
  labs(
    title = "Top candidate microbiome-metabolome biomarker pairs",
    x = "Overall multi-evidence score",
    y = "Taxon | metabolite",
    size = "|Spearman rho|",
    shape = "Metabolite FDR < 0.10"
  ) +
  theme_bw(base_size = 11)

ggsave(
  "results/figures/top_candidate_biomarker_dotplot.png",
  p_dot,
  width = 10,
  height = 8,
  dpi = 300
)

message_header("Create top candidate biomarker heatmap")

top_heatmap <- ranking_scored %>%
  slice_head(n = 50)

heatmap_df <- top_heatmap %>%
  mutate(
    taxon_label = factor(taxon_label, levels = unique(taxon_label)),
    metabolite_label = factor(metabolite_label, levels = rev(unique(metabolite_label)))
  )

p_heatmap <- ggplot(
  heatmap_df,
  aes(x = taxon_label, y = metabolite_label, fill = spearman_rho)
) +
  geom_tile() +
  labs(
    title = "Top candidate biomarker pair correlation structure",
    x = "Genus-level taxon",
    y = "Metabolite feature",
    fill = "Spearman rho"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "results/figures/top_candidate_biomarker_heatmap.png",
  p_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)

message_header("Update analysis decisions")

decision_addition <- c(
  "",
  "## Final candidate biomarker ranking",
  "",
  "- Candidate biomarker ranking was performed at the taxa-metabolite pair level.",
  "- The ranking combined taxa MaAsLin2 effect size/FDR, metabolite MaAsLin2 effect size/FDR, taxa-metabolite correlation strength/FDR, sPLS loading, prevalence, missingness robustness, and covariate-adjusted correlation evidence.",
  "- Because genus-level taxa did not survive MaAsLin2 FDR correction, taxa evidence was interpreted as exploratory and weighted together with cross-omics and metabolite evidence.",
  "- The final ranking represents candidate biomarker discovery, not clinical validation."
)

cat(
  paste(decision_addition, collapse = "\n"),
  file = "docs/analysis_decisions.md",
  append = TRUE
)

message_header("Final candidate biomarker ranking completed")

message("Outputs written:")
message("- results/tables/final_candidate_biomarker_ranking.tsv")
message("- results/figures/top_candidate_biomarker_dotplot.png")
message("- results/figures/top_candidate_biomarker_heatmap.png")