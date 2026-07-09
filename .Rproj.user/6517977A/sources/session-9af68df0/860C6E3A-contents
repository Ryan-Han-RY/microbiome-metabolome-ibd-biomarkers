# Annotate metabolite feature IDs using available HMP2/IBDMDB annotation files
#
# Main goal:
#   Map feature IDs such as C18n_QIxxxx to available metabolite names when possible.
#
# Outputs:
#   data/metadata/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx
#   results/tables/c18n_metabolite_annotation_raw.tsv
#   results/tables/final_candidate_biomarker_ranking_annotated.tsv
#   results/tables/metabolite_annotation_coverage_summary.tsv

source("R/00_utils_io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
})

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

ensure_project_dirs()
ensure_dir("data/metadata")
ensure_dir("results/tables")

message_header("Download official HMP2 C18n metabolite annotation file")

annotation_url <- "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MBX/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"

annotation_file <- "data/metadata/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"

if (!file.exists(annotation_file)) {
  download.file(
    annotation_url,
    destfile = annotation_file,
    mode = "wb"
  )
}

message("Annotation file exists: ", file.exists(annotation_file))

message_header("Read annotation file")

sheets <- readxl::excel_sheets(annotation_file)
message("Sheets:")
print(sheets)

ann_raw <- readxl::read_excel(annotation_file, sheet = sheets[1])

ann <- ann_raw %>%
  as.data.frame() %>%
  as_tibble()

message("Annotation table columns:")
print(names(ann))

write_tsv_safe(
  ann %>% mutate(across(everything(), as.character)),
  "results/tables/c18n_metabolite_annotation_raw.tsv"
)

message_header("Auto-detect ID and name columns")

clean_colnames <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    tolower()
}

ann_clean <- ann %>%
  rename_with(clean_colnames) %>%
  mutate(across(everything(), as.character))

message("Cleaned annotation columns:")
print(names(ann_clean))

id_candidates <- names(ann_clean)[
  str_detect(names(ann_clean), "qi|feature|metabolite.*id|compound.*id|id")
]

name_candidates <- names(ann_clean)[
  str_detect(names(ann_clean), "name|metabolite|compound|annotation|identity")
]

message("Possible ID columns:")
print(id_candidates)

message("Possible name columns:")
print(name_candidates)

# The auto-detection below is intentionally conservative.
# If the wrong columns are selected, manually set id_col and name_col after checking the printed column names.
id_col <- id_candidates[1]
name_col <- setdiff(name_candidates, id_col)[1]

if (is.na(id_col) || is.na(name_col)) {
  stop(
    "Could not automatically detect ID/name columns. Check the printed column names and manually set id_col/name_col."
  )
}

message("Selected ID column: ", id_col)
message("Selected name column: ", name_col)

make_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[ _.-]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    toupper()
}

annotation_map <- ann_clean %>%
  transmute(
    annotation_id_raw = .data[[id_col]],
    metabolite_name_candidate = .data[[name_col]]
  ) %>%
  mutate(
    annotation_key = make_key(annotation_id_raw),
    metabolite_name_candidate = str_squish(metabolite_name_candidate),
    metabolite_name_candidate = na_if(metabolite_name_candidate, "")
  ) %>%
  filter(!is.na(annotation_key), !is.na(metabolite_name_candidate)) %>%
  distinct(annotation_key, .keep_all = TRUE)

message("Annotation rows available after cleaning: ", nrow(annotation_map))

message_header("Annotate final biomarker ranking")

ranking <- readr::read_tsv(
  "results/tables/final_candidate_biomarker_ranking.tsv",
  show_col_types = FALSE
)

ranking_annotated <- ranking %>%
  mutate(
    metabolite_key_full = make_key(metabolite),
    metabolite_key_without_platform = str_replace(metabolite_key_full, "^(C18N|C8P|HILP)_", ""),
    metabolite_platform = str_extract(metabolite_key_full, "^(C18N|C8P|HILP)")
  ) %>%
  left_join(
    annotation_map %>%
      rename(
        metabolite_key_full = annotation_key,
        c18n_annotation_by_full_id = metabolite_name_candidate
      ),
    by = "metabolite_key_full"
  ) %>%
  left_join(
    annotation_map %>%
      rename(
        metabolite_key_without_platform = annotation_key,
        c18n_annotation_by_qi_id = metabolite_name_candidate
      ),
    by = "metabolite_key_without_platform"
  ) %>%
  mutate(
    metabolite_annotation = coalesce(
      c18n_annotation_by_full_id,
      c18n_annotation_by_qi_id
    ),
    metabolite_display_name = case_when(
      !is.na(metabolite_annotation) ~ metabolite_annotation,
      TRUE ~ metabolite_label
    ),
    annotation_status = case_when(
      !is.na(metabolite_annotation) ~ "annotated_from_HMP2_C18n_table",
      metabolite_platform == "C18N" ~ "C18n_feature_without_matched_annotation",
      TRUE ~ "unannotated_feature_no_matching_platform_annotation"
    )
  ) %>%
  select(
    rank,
    taxon,
    taxon_label,
    metabolite,
    metabolite_label,
    metabolite_platform,
    metabolite_display_name,
    metabolite_annotation,
    annotation_status,
    everything()
  )

write_tsv_safe(
  ranking_annotated,
  "results/tables/final_candidate_biomarker_ranking_annotated.tsv"
)

coverage <- ranking_annotated %>%
  count(metabolite_platform, annotation_status, name = "n_features") %>%
  arrange(metabolite_platform, annotation_status)

write_tsv_safe(
  coverage,
  "results/tables/metabolite_annotation_coverage_summary.tsv"
)

message_header("Annotation coverage summary")
print(coverage)

message_header("Top annotated ranking preview")

ranking_annotated %>%
  select(
    rank,
    taxon_label,
    metabolite,
    metabolite_display_name,
    annotation_status,
    metabolite_fdr,
    spearman_rho,
    correlation_fdr,
    overall_score
  ) %>%
  head(20) %>%
  print(n = 20, width = Inf)

message_header("Metabolite annotation completed")

message("Outputs written:")
message("- results/tables/c18n_metabolite_annotation_raw.tsv")
message("- results/tables/final_candidate_biomarker_ranking_annotated.tsv")
message("- results/tables/metabolite_annotation_coverage_summary.tsv")