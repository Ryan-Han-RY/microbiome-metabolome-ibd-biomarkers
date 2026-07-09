# Utility functions for IBDMDB/HMP2 microbiome-metabolome analysis

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tibble)
  library(fs)
})

message_header <- function(text) {
  cat("\n", strrep("=", 72), "\n", text, "\n", strrep("=", 72), "\n", sep = "")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

ensure_project_dirs <- function() {
  dirs <- c(
    "R",
    "data",
    "data/raw",
    "data/processed",
    "data/metadata",
    "results",
    "results/figures",
    "results/tables",
    "results/models",
    "docs"
  )
  
  walk(dirs, ensure_dir)
  
  gitkeep_dirs <- c("data/raw", "data/processed", "data/metadata",
                    "results/figures", "results/tables", "results/models")
  walk(file.path(gitkeep_dirs, ".gitkeep"), ~ if (!file.exists(.x)) file.create(.x))
  
  invisible(TRUE)
}

clean_names_vec <- function(x) {
  x <- gsub("\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  make.unique(x, sep = "_")
}

clean_colnames <- function(df) {
  names(df) <- clean_names_vec(names(df))
  df
}

normalize_sample_id <- function(x) {
  x %>%
    as.character() %>%
    str_replace("^#", "") %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace("_taxonomic_profile_3$", "") %>%
    str_replace("_taxonomic_profile$", "") %>%
    str_replace("_taxonomic$", "") %>%
    str_replace("_metabolomics$", "") %>%
    str_replace("\\.metabolomics$", "")
}

normalize_diagnosis <- function(x) {
  y <- as.character(x)
  y <- str_trim(str_to_lower(y))
  
  case_when(
    str_detect(y, "^cd$|crohn") ~ "CD",
    str_detect(y, "^uc$|ulcerative") ~ "UC",
    str_detect(y, "non.*ibd|nonibd|control|healthy") ~ "non_IBD",
    TRUE ~ NA_character_
  )
}

choose_col <- function(df, candidates, required = TRUE, purpose = "column") {
  candidates_clean <- clean_names_vec(candidates)
  hit <- candidates_clean[candidates_clean %in% names(df)]
  
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  
  if (required) {
    stop(
      "Could not find ", purpose, ". Tried: ",
      paste(candidates_clean, collapse = ", "),
      "\nAvailable columns include:\n",
      paste(head(names(df), 80), collapse = ", ")
    )
  }
  
  NA_character_
}

write_tsv_safe <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(as_tibble(x), path, na = "NA")
}

download_file_checked <- function(url, dest) {
  ensure_dir(dirname(dest))
  
  if (file.exists(dest) && file.size(dest) > 0) {
    message("File already exists: ", dest)
    return(invisible(TRUE))
  }
  
  message("Downloading: ", basename(dest))
  
  ok <- tryCatch(
    {
      download.file(
        url = url,
        destfile = dest,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )
      TRUE
    },
    error = function(e) {
      message("Download failed: ", conditionMessage(e))
      FALSE
    }
  )
  
  if (!ok || !file.exists(dest) || file.size(dest) == 0) {
    stop(
      "Download failed for: ", basename(dest), "\n",
      "Open the IBDMDB results page manually and download this file into data/raw/."
    )
  }
  
  invisible(TRUE)
}

file_manifest_row <- function(file_name, data_type, source, url,
                              destination, used_in_analysis = "yes",
                              notes = "") {
  checksum <- if (file.exists(destination)) {
    unname(tools::md5sum(destination))
  } else {
    NA_character_
  }
  
  size_mb <- if (file.exists(destination)) {
    round(file.size(destination) / 1024^2, 3)
  } else {
    NA_real_
  }
  
  tibble(
    file_name = file_name,
    data_type = data_type,
    source = source,
    download_date = as.character(Sys.Date()),
    original_format = tools::file_ext(file_name),
    destination = destination,
    file_size_mb = size_mb,
    md5 = checksum,
    used_in_analysis = used_in_analysis,
    notes = notes,
    url = url
  )
}

read_taxonomic_profile <- function(path) {
  message("Reading microbiome taxonomic profile: ", path)
  
  con <- if (str_detect(path, "\\.gz$")) gzfile(path, "rt") else file(path, "rt")
  first_line <- readLines(con, n = 1, warn = FALSE)
  close(con)
  
  if (length(first_line) == 0) {
    stop("Empty file: ", path)
  }
  
  if (startsWith(first_line, "#")) {
    col_names <- str_split(str_replace(first_line, "^#", ""), "\t", simplify = TRUE)
    col_names <- as.character(col_names)
    df <- readr::read_tsv(
      path,
      skip = 1,
      col_names = col_names,
      show_col_types = FALSE,
      progress = FALSE
    )
  } else {
    df <- readr::read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE
    )
  }
  
  df <- as.data.frame(df, check.names = FALSE)
  names(df)[1] <- "feature_id"
  
  candidate_cols <- names(df)[-1]
  
  sample_cols <- candidate_cols[
    str_detect(candidate_cols, "^[A-Z]{3}[A-Z0-9]+")
  ]
  
  if (length(sample_cols) < 10) {
    numeric_score <- map_dbl(candidate_cols, function(z) {
      v <- suppressWarnings(as.numeric(df[[z]]))
      mean(!is.na(v))
    })
    
    sample_cols <- candidate_cols[numeric_score > 0.8]
  }
  
  if (length(sample_cols) < 10) {
    stop(
      "Could not detect enough microbiome sample columns. ",
      "Check the taxonomic profile table format."
    )
  }
  
  abundance_df <- df[, sample_cols, drop = FALSE]
  abundance_df[] <- lapply(abundance_df, function(v) suppressWarnings(as.numeric(v)))
  
  mat_feature_sample <- as.matrix(abundance_df)
  rownames(mat_feature_sample) <- df$feature_id
  
  mat_sample_feature <- t(mat_feature_sample)
  rownames(mat_sample_feature) <- normalize_sample_id(rownames(mat_sample_feature))
  
  if (any(duplicated(rownames(mat_sample_feature)))) {
    mat_sample_feature <- rowsum(
      mat_sample_feature,
      group = rownames(mat_sample_feature),
      reorder = FALSE
    )
  }
  
  feature_info <- tibble(feature_id = colnames(mat_sample_feature)) %>%
    mutate(
      tax_rank = case_when(
        str_detect(feature_id, "\\|t__|^t__") ~ "strain",
        str_detect(feature_id, "\\|s__|^s__") ~ "species",
        str_detect(feature_id, "\\|g__|^g__") ~ "genus",
        str_detect(feature_id, "\\|f__|^f__") ~ "family",
        str_detect(feature_id, "\\|o__|^o__") ~ "order",
        str_detect(feature_id, "\\|c__|^c__") ~ "class",
        str_detect(feature_id, "\\|p__|^p__") ~ "phylum",
        str_detect(feature_id, "\\|k__|^k__") ~ "kingdom",
        TRUE ~ "unclassified"
      ),
      feature_label = str_replace(feature_id, "^.*\\|", "")
    )
  
  list(
    abundance = mat_sample_feature,
    feature_info = feature_info,
    source_file = path
  )
}

read_biom_matrix <- function(path) {
  message("Reading metabolomics BIOM file: ", path)
  
  if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("Package `biomformat` is not installed. Run: BiocManager::install('biomformat')")
  }
  
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package `R.utils` is not installed. Run: install.packages('R.utils')")
  }
  
  read_path <- path
  
  if (str_detect(path, "\\.gz$")) {
    read_path <- tempfile(fileext = ".biom")
    R.utils::gunzip(
      filename = path,
      destname = read_path,
      overwrite = TRUE,
      remove = FALSE
    )
  }
  
  biom_obj <- biomformat::read_biom(read_path)
  mat <- as.matrix(biomformat::biom_data(biom_obj))
  
  score_sample_like <- function(ids) {
    ids2 <- normalize_sample_id(ids)
    mean(str_detect(ids2, "^[A-Z]{3}[A-Z0-9]+"), na.rm = TRUE)
  }
  
  row_score <- score_sample_like(rownames(mat))
  col_score <- score_sample_like(colnames(mat))
  
  if (col_score >= row_score) {
    mat_sample_feature <- t(mat)
  } else {
    mat_sample_feature <- mat
  }
  
  rownames(mat_sample_feature) <- normalize_sample_id(rownames(mat_sample_feature))
  colnames(mat_sample_feature) <- make.unique(as.character(colnames(mat_sample_feature)))
  
  if (any(duplicated(rownames(mat_sample_feature)))) {
    mat_sample_feature <- rowsum(
      mat_sample_feature,
      group = rownames(mat_sample_feature),
      reorder = FALSE
    )
  }
  
  storage.mode(mat_sample_feature) <- "numeric"
  
  feature_info <- tibble(
    feature_id = colnames(mat_sample_feature)
  )
  
  list(
    abundance = mat_sample_feature,
    feature_info = feature_info,
    source_file = path
  )
}

count_unique_nonmissing <- function(x) {
  n_distinct(x[!is.na(x) & x != ""])
}