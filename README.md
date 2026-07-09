# Microbiome-Metabolome Integration Reveals Candidate Biomarker Signatures Associated with IBD Disease State and Activity

This repository contains a reproducible exploratory multi-omics analysis using public IBDMDB/HMP2 data.

The aim is candidate biomarker discovery, not clinical validation.

Raw data are not included in this repository due to file size.
Scripts reproduce the download, preprocessing, sample matching, quality control, and downstream analysis steps from public IBDMDB/HMP2 resources.

## Main workflow

1. Download public metadata, metagenomic taxonomic profiles, and stool metabolomics profiles.
2. Parse input tables and standardise sample identifiers.
3. Construct a matched microbiome-metabolome cohort.
4. Generate quality-control tables and figures.
5. Perform downstream microbiome-metabolome integration and candidate biomarker ranking.

## Data policy

Raw files are stored locally under `data/raw/` and are excluded from GitHub.
Processed intermediate `.rds` objects are also excluded because they can be regenerated.
