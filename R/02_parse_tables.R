# Parse metadata, microbiome taxonomic profiles, and metabolomics BIOM table

source("R/00_utils_io.R")

ensure_project_dirs()

message_header("Parse HMP2 metadata")

metadata_path <- "data/raw/hmp2_metadata_2018-08-20.csv"
microbiome_path <- "data/raw/taxonomic_profiles.tsv.gz"
metabolomics_path <- "data/raw/HMP2_metabolomics_w_metadata.biom.gz"
metabolite_anno_path <- "data/raw/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"

if (!file.exists(metadata_path)) stop("Missing file: ", metadata_path)
if (!file.exists(microbiome_path)) stop("Missing file: ", microbiome_path)
if (!file.exists(metabolomics_path)) stop("Missing file: ", metabolomics_path)

metadata_raw <- readr::read_csv(
  metadata_path,
  show_col_types = FALSE,
  progress = FALSE
)

metadata_clean <- metadata_raw %>%
  clean_colnames()

sample_id_col <- choose_col(
  metadata_clean,
  candidates = c(
    "external_id", "external id", "external.id",
    "sample_id", "sample id", "sampleid",
    "site_sub_coll", "site_sub_coll_id"
  ),
  required = TRUE,
  purpose = "metadata sample ID"
)

participant_id_col <- choose_col(
  metadata_clean,
  candidates = c(
    "participant_id", "participant id", "participant.id",
    "subject_id", "subject id", "subjectid",
    "host_subject_id"
  ),
  required = TRUE,
  purpose = "participant ID"
)

diagnosis_col <- choose_col(
  metadata_clean,
  candidates = c(
    "diagnosis", "diagnosis_study", "disease", "disease_state",
    "ibd_diagnosis", "study_diagnosis"
  ),
  required = TRUE,
  purpose = "diagnosis"
)

week_col <- choose_col(
  metadata_clean,
  candidates = c(
    "week_num", "week_number", "week", "visit_week",
    "week_num_clean", "interval_week"
  ),
  required = FALSE,
  purpose = "week number"
)

age_col <- choose_col(
  metadata_clean,
  candidates = c("age", "age_at_collection", "age_years"),
  required = FALSE,
  purpose = "age"
)

sex_col <- choose_col(
  metadata_clean,
  candidates = c("sex", "gender"),
  required = FALSE,
  purpose = "sex"
)

antibiotic_col <- choose_col(
  metadata_clean,
  candidates = c("antibiotics", "antibiotic", "antibiotic_use", "antibiotics_use"),
  required = FALSE,
  purpose = "antibiotic use"
)

metadata_clean <- metadata_clean %>%
  mutate(
    sample_id_clean = normalize_sample_id(.data[[sample_id_col]]),
    participant_id_clean = as.character(.data[[participant_id_col]]),
    diagnosis_clean = normalize_diagnosis(.data[[diagnosis_col]]),
    week_num_clean = if (!is.na(week_col)) {
      readr::parse_number(as.character(.data[[week_col]]))
    } else {
      NA_real_
    },
    age_clean = if (!is.na(age_col)) {
      suppressWarnings(as.numeric(.data[[age_col]]))
    } else {
      NA_real_
    },
    sex_clean = if (!is.na(sex_col)) {
      as.character(.data[[sex_col]])
    } else {
      NA_character_
    },
    antibiotic_use_clean = if (!is.na(antibiotic_col)) {
      as.character(.data[[antibiotic_col]])
    } else {
      NA_character_
    }
  )

attr(metadata_clean, "column_map") <- list(
  sample_id_col = sample_id_col,
  participant_id_col = participant_id_col,
  diagnosis_col = diagnosis_col,
  week_col = week_col,
  age_col = age_col,
  sex_col = sex_col,
  antibiotic_col = antibiotic_col
)

saveRDS(metadata_clean, "data/processed/metadata_clean.rds")

metadata_column_map <- tibble(
  standard_variable = c(
    "sample_id_clean",
    "participant_id_clean",
    "diagnosis_clean",
    "week_num_clean",
    "age_clean",
    "sex_clean",
    "antibiotic_use_clean"
  ),
  source_column = c(
    sample_id_col,
    participant_id_col,
    diagnosis_col,
    week_col,
    age_col,
    sex_col,
    antibiotic_col
  )
)

write_tsv_safe(metadata_column_map, "results/tables/metadata_column_map.tsv")

message("Metadata rows: ", nrow(metadata_clean))
message("Metadata columns: ", ncol(metadata_clean))
message("Sample ID column used: ", sample_id_col)
message("Participant ID column used: ", participant_id_col)
message("Diagnosis column used: ", diagnosis_col)

message_header("Parse microbiome taxonomic profile")

microbiome_raw <- read_taxonomic_profile(microbiome_path)
saveRDS(microbiome_raw, "data/processed/microbiome_raw.rds")

message("Microbiome samples: ", nrow(microbiome_raw$abundance))
message("Microbiome features: ", ncol(microbiome_raw$abundance))

message_header("Parse metabolomics BIOM table")

metabolomics_raw <- read_biom_matrix(metabolomics_path)
saveRDS(metabolomics_raw, "data/processed/metabolomics_raw.rds")

message("Metabolomics samples: ", nrow(metabolomics_raw$abundance))
message("Metabolomics features: ", ncol(metabolomics_raw$abundance))

message_header("Parse metabolite annotation if available")

if (file.exists(metabolite_anno_path)) {
  metabolite_annotation <- readxl::read_excel(metabolite_anno_path)
  metabolite_annotation <- metabolite_annotation %>%
    as.data.frame() %>%
    clean_colnames()
  
  saveRDS(metabolite_annotation, "data/metadata/metabolite_annotation.rds")
  message("Metabolite annotation rows: ", nrow(metabolite_annotation))
} else {
  message("Metabolite annotation file not found. Skipping.")
}

message_header("Write raw table dimensions")

raw_table_dimensions <- tibble(
  object = c("metadata_clean", "microbiome_raw", "metabolomics_raw"),
  n_samples_or_rows = c(
    nrow(metadata_clean),
    nrow(microbiome_raw$abundance),
    nrow(metabolomics_raw$abundance)
  ),
  n_features_or_columns = c(
    ncol(metadata_clean),
    ncol(microbiome_raw$abundance),
    ncol(metabolomics_raw$abundance)
  ),
  source_file = c(
    metadata_path,
    microbiome_path,
    metabolomics_path
  )
)

write_tsv_safe(raw_table_dimensions, "results/tables/raw_table_dimensions.tsv")

message_header("Write sample ID overlap before matching")

metadata_ids <- unique(metadata_clean$sample_id_clean)
metadata_ids <- metadata_ids[!is.na(metadata_ids) & metadata_ids != ""]

microbiome_ids <- rownames(microbiome_raw$abundance)
metabolomics_ids <- rownames(metabolomics_raw$abundance)

sample_id_overlap_before_matching <- bind_rows(
  tibble(
    comparison = "metadata_vs_microbiome",
    n_left = length(metadata_ids),
    n_right = length(microbiome_ids),
    n_overlap = length(intersect(metadata_ids, microbiome_ids))
  ),
  tibble(
    comparison = "metadata_vs_metabolomics",
    n_left = length(metadata_ids),
    n_right = length(metabolomics_ids),
    n_overlap = length(intersect(metadata_ids, metabolomics_ids))
  ),
  tibble(
    comparison = "microbiome_vs_metabolomics",
    n_left = length(microbiome_ids),
    n_right = length(metabolomics_ids),
    n_overlap = length(intersect(microbiome_ids, metabolomics_ids))
  ),
  tibble(
    comparison = "metadata_vs_microbiome_vs_metabolomics",
    n_left = length(metadata_ids),
    n_right = NA_integer_,
    n_overlap = length(Reduce(intersect, list(metadata_ids, microbiome_ids, metabolomics_ids)))
  )
)

write_tsv_safe(
  sample_id_overlap_before_matching,
  "results/tables/sample_id_overlap_before_matching.tsv"
)

print(sample_id_overlap_before_matching)

message_header("Parsing completed")