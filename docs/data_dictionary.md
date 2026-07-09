# Data dictionary

This file records the main analysis-ready variables created during preprocessing.

| variable | meaning | source |
|---|---|---|
| sample_id_clean | Standardised sample identifier used for matching | HMP2 metadata / omics sample names |
| participant_id_clean | Participant-level identifier | HMP2 metadata |
| diagnosis_clean | Standardised diagnosis group: CD, UC, non_IBD | HMP2 metadata |
| week_num_clean | Numeric study week when available | HMP2 metadata |
| is_primary_independent | TRUE for the one-sample-per-participant subset | Derived |
