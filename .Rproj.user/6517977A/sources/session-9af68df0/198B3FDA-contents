# Download and document public IBDMDB/HMP2 input files

source("R/00_utils_io.R")

message_header("Create project folders")
ensure_project_dirs()

message_header("Write project-level helper files")

gitignore_lines <- c(
  ".Rhistory",
  ".RData",
  ".Ruserdata",
  ".DS_Store",
  "",
  "# Do not upload raw public data files to GitHub",
  "data/raw/*",
  "!data/raw/.gitkeep",
  "",
  "# Large intermediate objects",
  "data/processed/*.rds",
  "",
  "# Rendered reports can be regenerated",
  "*.html",
  "*.pdf"
)

writeLines(gitignore_lines, ".gitignore")

if (!file.exists("README.md")) {
  writeLines(c(
    "# Microbiome-Metabolome Integration Reveals Candidate Biomarker Signatures Associated with IBD Disease State and Activity",
    "",
    "This repository contains a reproducible exploratory multi-omics analysis using public IBDMDB/HMP2 data.",
    "",
    "The aim is candidate biomarker discovery, not clinical validation.",
    "",
    "Raw data are not included in this repository due to file size.",
    "Scripts reproduce the download, preprocessing, sample matching, quality control, and downstream analysis steps from public IBDMDB/HMP2 resources.",
    "",
    "## Main workflow",
    "",
    "1. Download public metadata, metagenomic taxonomic profiles, and stool metabolomics profiles.",
    "2. Parse input tables and standardise sample identifiers.",
    "3. Construct a matched microbiome-metabolome cohort.",
    "4. Generate quality-control tables and figures.",
    "5. Perform downstream microbiome-metabolome integration and candidate biomarker ranking.",
    "",
    "## Data policy",
    "",
    "Raw files are stored locally under `data/raw/` and are excluded from GitHub.",
    "Processed intermediate `.rds` objects are also excluded because they can be regenerated."
  ), "README.md")
}

if (!file.exists("analysis_report.qmd")) {
  writeLines(c(
    "---",
    "title: \"Microbiome-Metabolome Integration in IBD\"",
    "format: html",
    "editor: visual",
    "---",
    "",
    "## Aim",
    "",
    "This project investigates whether paired gut microbial taxa and stool metabolite profiles can identify candidate biomarker signatures associated with IBD diagnosis and disease activity using matched metagenomic and metabolomic samples from the IBDMDB/HMP2 cohort.",
    "",
    "## Reproducibility",
    "",
    "All analysis code is stored in the `R/` folder. Raw data are downloaded from public IBDMDB/HMP2 resources and are not included in this repository."
  ), "analysis_report.qmd")
}

if (!file.exists("docs/data_dictionary.md")) {
  writeLines(c(
    "# Data dictionary",
    "",
    "This file records the main analysis-ready variables created during preprocessing.",
    "",
    "| variable | meaning | source |",
    "|---|---|---|",
    "| sample_id_clean | Standardised sample identifier used for matching | HMP2 metadata / omics sample names |",
    "| participant_id_clean | Participant-level identifier | HMP2 metadata |",
    "| diagnosis_clean | Standardised diagnosis group: CD, UC, non_IBD | HMP2 metadata |",
    "| week_num_clean | Numeric study week when available | HMP2 metadata |",
    "| is_primary_independent | TRUE for the one-sample-per-participant subset | Derived |"
  ), "docs/data_dictionary.md")
}

if (!file.exists("docs/analysis_decisions.md")) {
  writeLines(c(
    "# Analysis decisions",
    "",
    "This file records practical analysis decisions made during cohort construction and preprocessing.",
    "",
    "It is updated by the sample matching script after the matched cohort is created."
  ), "docs/analysis_decisions.md")
}

message_header("Download input files")

sources <- tibble::tribble(
  ~file_name, ~data_type, ~source, ~url, ~destination, ~used_in_analysis, ~notes,
  "hmp2_metadata_2018-08-20.csv",
  "metadata",
  "IBDMDB/HMP2",
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv",
  "data/raw/hmp2_metadata_2018-08-20.csv",
  "yes",
  "Main HMP2 sample metadata.",
  
  "taxonomic_profiles.tsv.gz",
  "metagenomics_taxa",
  "IBDMDB/HMP2 MGX products",
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz",
  "data/raw/taxonomic_profiles.tsv.gz",
  "yes",
  "Merged metagenomic taxonomic profile table.",
  
  "HMP2_metabolomics_w_metadata.biom.gz",
  "metabolomics",
  "IBDMDB/HMP2 MBX products",
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MBX/HMP2_metabolomics_w_metadata.biom.gz",
  "data/raw/HMP2_metabolomics_w_metadata.biom.gz",
  "yes",
  "Merged stool metabolomics BIOM table with metadata.",
  
  "C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx",
  "metabolite_annotation",
  "IBDMDB/HMP2 MBX products",
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MBX/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx",
  "data/raw/C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx",
  "optional",
  "Metabolite annotation table, used if needed for interpretation."
)

walk2(sources$url, sources$destination, download_file_checked)

manifest <- pmap_dfr(
  sources,
  function(file_name, data_type, source, url, destination, used_in_analysis, notes) {
    file_manifest_row(
      file_name = file_name,
      data_type = data_type,
      source = source,
      url = url,
      destination = destination,
      used_in_analysis = used_in_analysis,
      notes = notes
    )
  }
)

write_tsv_safe(manifest, "data_manifest.tsv")

message_header("Downloaded files")
print(manifest %>% select(file_name, data_type, file_size_mb, md5, used_in_analysis))

message_header("Download script completed")