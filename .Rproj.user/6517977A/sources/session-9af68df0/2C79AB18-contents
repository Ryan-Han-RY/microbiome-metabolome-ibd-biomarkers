# Diversity and ordination analysis for microbiome and metabolomics

source("R/00_utils_io.R")

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package `vegan` is not installed. Run: install.packages('vegan')")
}

ensure_project_dirs()

message_header("Load preprocessed data")

metadata <- readRDS("data/processed/matched_metadata.rds")
taxa_genus_rel <- readRDS("data/processed/taxa_genus_filtered.rds")
taxa_genus_clr <- readRDS("data/processed/taxa_genus_clr.rds")
metab_scaled <- readRDS("data/processed/metabolites_log_scaled.rds")

metadata <- metadata %>%
  mutate(
    diagnosis_plot = recode(diagnosis_clean, "non_IBD" = "non-IBD")
  )

primary_meta <- metadata %>%
  filter(is_primary_independent) %>%
  arrange(sample_id_clean)

primary_ids <- primary_meta$sample_id_clean

common_ids <- Reduce(
  intersect,
  list(primary_ids, rownames(taxa_genus_rel), rownames(taxa_genus_clr), rownames(metab_scaled))
)

if (length(common_ids) < 10) {
  stop("Too few common primary independent samples for ordination.")
}

primary_meta <- primary_meta %>%
  filter(sample_id_clean %in% common_ids) %>%
  arrange(match(sample_id_clean, common_ids))

taxa_rel <- taxa_genus_rel[primary_meta$sample_id_clean, , drop = FALSE]
taxa_clr <- taxa_genus_clr[primary_meta$sample_id_clean, , drop = FALSE]
metab <- metab_scaled[primary_meta$sample_id_clean, , drop = FALSE]

message("Primary independent samples used: ", nrow(primary_meta))

message_header("Alpha diversity")

alpha_df <- tibble(
  sample_id_clean = rownames(taxa_rel),
  shannon = vegan::diversity(taxa_rel, index = "shannon"),
  observed_richness = rowSums(taxa_rel > 0)
) %>%
  left_join(
    primary_meta %>% select(sample_id_clean, diagnosis_clean, diagnosis_plot),
    by = "sample_id_clean"
  )

write_tsv_safe(alpha_df, "results/tables/alpha_diversity_values.tsv")

run_kruskal <- function(df, value_col) {
  kw_result <- kruskal.test(df[[value_col]] ~ df$diagnosis_clean)
  
  tibble(
    metric = value_col,
    method = "Kruskal-Wallis",
    statistic = unname(kw_result$statistic),
    p_value = kw_result$p.value
  )
}

alpha_results <- bind_rows(
  run_kruskal(alpha_df, "shannon"),
  run_kruskal(alpha_df, "observed_richness")
) %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))

write_tsv_safe(alpha_results, "results/tables/alpha_diversity_results.tsv")

p_shannon <- ggplot(alpha_df, aes(x = diagnosis_plot, y = shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.7) +
  labs(
    title = "Shannon diversity by diagnosis",
    x = "Diagnosis",
    y = "Shannon diversity"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/alpha_diversity_shannon_boxplot.png",
  p_shannon,
  width = 6.5,
  height = 5,
  dpi = 300
)

p_observed <- ggplot(alpha_df, aes(x = diagnosis_plot, y = observed_richness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.7) +
  labs(
    title = "Observed genus richness by diagnosis",
    x = "Diagnosis",
    y = "Observed richness"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/alpha_diversity_observed_boxplot.png",
  p_observed,
  width = 6.5,
  height = 5,
  dpi = 300
)

message_header("Bray-Curtis PCoA")

bray_dist <- vegan::vegdist(taxa_rel, method = "bray")

bray_pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)

bray_var <- round(100 * bray_pcoa$eig[1:2] / sum(bray_pcoa$eig[bray_pcoa$eig > 0]), 1)

bray_df <- tibble(
  sample_id_clean = rownames(taxa_rel),
  PCoA1 = bray_pcoa$points[, 1],
  PCoA2 = bray_pcoa$points[, 2]
) %>%
  left_join(
    primary_meta %>% select(sample_id_clean, diagnosis_clean, diagnosis_plot),
    by = "sample_id_clean"
  )

p_bray <- ggplot(bray_df, aes(x = PCoA1, y = PCoA2, shape = diagnosis_plot)) +
  geom_point(size = 2.3, alpha = 0.85) +
  labs(
    title = "Bray-Curtis PCoA by diagnosis",
    x = paste0("PCoA1 (", bray_var[1], "%)"),
    y = paste0("PCoA2 (", bray_var[2], "%)"),
    shape = "Diagnosis"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/bray_pcoa_by_diagnosis.png",
  p_bray,
  width = 7,
  height = 5.5,
  dpi = 300
)

message_header("Aitchison PCA")

aitchison_pca <- prcomp(taxa_clr, center = TRUE, scale. = FALSE)
aitchison_var <- round(100 * (aitchison_pca$sdev^2 / sum(aitchison_pca$sdev^2)), 1)

aitchison_df <- tibble(
  sample_id_clean = rownames(aitchison_pca$x),
  PC1 = aitchison_pca$x[, 1],
  PC2 = aitchison_pca$x[, 2]
) %>%
  left_join(
    primary_meta %>% select(sample_id_clean, diagnosis_clean, diagnosis_plot),
    by = "sample_id_clean"
  )

p_aitchison <- ggplot(aitchison_df, aes(x = PC1, y = PC2, shape = diagnosis_plot)) +
  geom_point(size = 2.3, alpha = 0.85) +
  labs(
    title = "Aitchison PCA by diagnosis",
    x = paste0("PC1 (", aitchison_var[1], "%)"),
    y = paste0("PC2 (", aitchison_var[2], "%)"),
    shape = "Diagnosis"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/aitchison_pca_by_diagnosis.png",
  p_aitchison,
  width = 7,
  height = 5.5,
  dpi = 300
)

message_header("PERMANOVA and betadisper")

permanova_bray <- vegan::adonis2(
  bray_dist ~ diagnosis_clean,
  data = primary_meta,
  permutations = 999
)

aitchison_dist <- dist(taxa_clr)

permanova_aitchison <- vegan::adonis2(
  aitchison_dist ~ diagnosis_clean,
  data = primary_meta,
  permutations = 999
)

permanova_results <- bind_rows(
  as.data.frame(permanova_bray) %>%
    rownames_to_column("term") %>%
    mutate(distance = "Bray-Curtis"),
  as.data.frame(permanova_aitchison) %>%
    rownames_to_column("term") %>%
    mutate(distance = "Aitchison")
) %>%
  as_tibble()

write_tsv_safe(permanova_results, "results/tables/permanova_results.tsv")

bray_betadisper <- vegan::betadisper(bray_dist, primary_meta$diagnosis_clean)
bray_betadisper_perm <- vegan::permutest(bray_betadisper, permutations = 999)

aitchison_betadisper <- vegan::betadisper(aitchison_dist, primary_meta$diagnosis_clean)
aitchison_betadisper_perm <- vegan::permutest(aitchison_betadisper, permutations = 999)

extract_betadisper <- function(obj, distance_name) {
  tab <- as.data.frame(obj$tab)
  tab %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(distance = distance_name)
}

betadisper_results <- bind_rows(
  extract_betadisper(bray_betadisper_perm, "Bray-Curtis"),
  extract_betadisper(aitchison_betadisper_perm, "Aitchison")
)

write_tsv_safe(betadisper_results, "results/tables/betadisper_results.tsv")

message_header("Metabolomics PCA by diagnosis")

metab_pca <- prcomp(metab, center = FALSE, scale. = FALSE)
metab_var <- round(100 * (metab_pca$sdev^2 / sum(metab_pca$sdev^2)), 1)

metab_pca_df <- tibble(
  sample_id_clean = rownames(metab_pca$x),
  PC1 = metab_pca$x[, 1],
  PC2 = metab_pca$x[, 2]
) %>%
  left_join(
    primary_meta %>% select(sample_id_clean, diagnosis_clean, diagnosis_plot),
    by = "sample_id_clean"
  )

p_metab <- ggplot(metab_pca_df, aes(x = PC1, y = PC2, shape = diagnosis_plot)) +
  geom_point(size = 2.3, alpha = 0.85) +
  labs(
    title = "Metabolomics PCA by diagnosis",
    x = paste0("PC1 (", metab_var[1], "%)"),
    y = paste0("PC2 (", metab_var[2], "%)"),
    shape = "Diagnosis"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "results/figures/metabolomics_pca_by_diagnosis.png",
  p_metab,
  width = 7,
  height = 5.5,
  dpi = 300
)

message_header("Diversity and ordination completed")

message("Outputs written:")
message("- results/figures/alpha_diversity_shannon_boxplot.png")
message("- results/figures/alpha_diversity_observed_boxplot.png")
message("- results/figures/bray_pcoa_by_diagnosis.png")
message("- results/figures/aitchison_pca_by_diagnosis.png")
message("- results/figures/metabolomics_pca_by_diagnosis.png")
message("- results/tables/alpha_diversity_results.tsv")
message("- results/tables/permanova_results.tsv")
message("- results/tables/betadisper_results.tsv")