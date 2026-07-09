# Analysis decisions

## Cohort construction

- Metadata sample identifier source column: `external_id`.
- Participant identifier source column: `participant_id`.
- Diagnosis source column: `diagnosis`.
- Week number source column: `week_num`.
- Omics sample identifiers were normalised by removing technical suffixes such as `_taxonomic_profile`, `_taxonomic`, and `_metabolomics`.
- Matched samples were defined as sample IDs present in metadata, microbiome taxonomic profiles, and stool metabolomics profiles.
- Retained diagnosis groups: CD, UC, and non_IBD.

## Repeated samples

- Repeated samples were retained in `data/processed/matched_metadata.rds` for longitudinal sensitivity analysis.
- A participant-level independent subset was created using the column `is_primary_independent`.
- For participant-level exploratory integration, one stool sample per participant was retained.
- Selection rule for the independent subset: prefer baseline week when available; otherwise use the earliest available week; ties are broken by sample ID.

## Exclusion summary

- Metadata samples without microbiome profile: 1254.
- Metadata samples without metabolomics profile: 2346.
- Microbiome-metabolomics samples not found in metadata: 0.
- Matched samples removed due to missing participant ID or diagnosis: 0.

## Analysis sets

- Longitudinal matched set: `data/processed/matched_metadata.rds`, all rows where `is_primary_independent` can be TRUE or FALSE; n samples = 388.
- Primary independent set: rows where `is_primary_independent == TRUE`; n samples = 105.
- Number of participants with repeated matched samples: 99.

## Differential abundance modelling

- Genus-level microbiome differential abundance was performed on CLR-transformed taxa abundance.
- The primary independent set was used to avoid repeated-measure leakage.
- Comparisons tested: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.
- Linear models were used as the main robust workflow for this stage.
- Candidate covariates were age, sex, and antibiotic use when sufficiently available.
- FDR correction was applied within each comparison using the Benjamini-Hochberg method.
## Metabolite association modelling

- Metabolite association models were performed on log2-transformed, z-score scaled metabolomics features.
- The primary independent set was used to avoid repeated-measure leakage.
- Comparisons tested: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.
- Candidate covariates were age, sex, and antibiotic use when sufficiently available.
- FDR correction was applied within each comparison using the Benjamini-Hochberg method.
## MaasLin2 differential abundance modelling

- Genus-level microbiome differential abundance was performed using MaasLin2.
- Input features were filtered genus-level relative abundance tables from the primary independent set.
- MaasLin2 was run with LM analysis, CLR normalization, no additional transform, and BH FDR correction.
- The primary independent set was used to avoid repeated-measure leakage.
- Comparisons tested separately: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.
- Candidate covariates were age, sex, and antibiotic use when sufficiently available.
- Significance threshold was FDR < 0.10.
## MaAsLin2 metabolite association modelling

- Metabolite association models were performed using MaAsLin2.
- Input metabolite features were half-minimum imputed positive abundance values from the primary independent set.
- MaAsLin2 was run with LM analysis, no normalization, LOG transform, and BH FDR correction.
- The primary independent set was used to avoid repeated-measure leakage.
- Comparisons tested separately: IBD vs non_IBD, CD vs non_IBD, UC vs non_IBD, and CD vs UC.
- Candidate covariates were age, sex, and antibiotic use when sufficiently available.
- Significance threshold was FDR < 0.10.
- Significant metabolite features were treated as candidate disease-associated features for downstream integration, not as clinically validated biomarkers.
## Cross-omics correlation and integration feature selection

- Cross-omics correlation was restricted to selected disease-associated or top-ranked features rather than all taxa-metabolite combinations.
- Because no genus-level taxa survived MaAsLin2 FDR < 0.10, taxa were selected by top absolute MaAsLin2 effect size and nominal association strength.
- Metabolite features were selected from MaAsLin2 FDR-significant disease-associated features, capped at the top 100 to keep integration interpretable.
- Spearman correlation was used between genus-level CLR abundance and log-scaled metabolite abundance.
- Correlation p-values were corrected using Benjamini-Hochberg FDR.
- Strong correlation pairs were defined as absolute Spearman rho >= 0.30 and FDR < 0.10.
- Covariate-adjusted correlation was performed by residualising taxa and metabolite features against available covariates before correlation.
## mixOmics multiblock sPLS-DA integration

- Multiblock sPLS-DA was used as an exploratory supervised feature-selection approach.
- The outcome was diagnosis group, and the two blocks were microbiome genus-level CLR features and log-scaled metabolite features.
- Only the primary independent sample set was used to avoid repeated-measure leakage.
- The design matrix used a low block connection weight of 0.1 to balance cross-omics correlation and diagnosis discrimination.
- The sPLS model was used for exploratory feature selection, not as a validated clinical classifier.
- Apparent classification accuracy was reported only as a training-set descriptive metric.
## Final candidate biomarker ranking

- Candidate biomarker ranking was performed at the taxa-metabolite pair level.
- The ranking combined taxa MaAsLin2 effect size/FDR, metabolite MaAsLin2 effect size/FDR, taxa-metabolite correlation strength/FDR, sPLS loading, prevalence, missingness robustness, and covariate-adjusted correlation evidence.
- Because genus-level taxa did not survive MaAsLin2 FDR correction, taxa evidence was interpreted as exploratory and weighted together with cross-omics and metabolite evidence.
- The final ranking represents candidate biomarker discovery, not clinical validation.