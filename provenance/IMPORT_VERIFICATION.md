> Historical record; current completeness is documented in the root REPRODUCTION_STATUS.md. Pending statements below describe the earlier snapshot.

# Import verification — 2026-11-25

The user's browser upload flattened the staging directory. Commit f333115a5e776b12fbf835f13d641201d58f348d restores the exact paths in STAGING_MANIFEST.csv for 77 files (54 R scripts and 23 scientific outputs/metadata/environment records). The existing checkpoint was compared after newline normalization and replaced by the archived uploaded copy. All 77 relocated Git blob IDs were verified unchanged; other repository files were preserved, including the original staging ZIP. All 23 unchanged scientific/output/metadata files match the manifest byte lengths. This verifies relocation integrity, not independent SHA-256 recomputation against the original project archive.

## Final synthesis preflight
Script: scripts/63_FINAL_FIGURE8_D1_D8.R. It uses ggplot2 and patchwork and expects to run from the repository root (or PARPI_PROJECT_ROOT).

The following inputs are present at the paths required by the script:
- results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv
- results/GSE163854/D8_analysis/D8_primary_Hallmark_GSEA.csv
- results/GSE165914/D7_analysis/D7A_M468_Hallmark_GSEA.csv
- results/GSE165914/D7_analysis/D7B_SUM1315_Hallmark_GSEA.csv
- results/GSE202434/D5_integration/D5_H3_P3_Hallmark_master_comparison.csv

Each input has exactly 50 distinct Hallmark identifiers; the five identifier sets are identical. This is an input-level preflight only. No R scripts, DE, GSEA or integration were executed in this check.

## Remaining
Full raw-data acquisition and upstream execution-order audit, missing or ambiguous historical scripts, historical environment records, final publication supplementary tables and artwork, data redistribution and code license review. Keep private until the author approves release.
