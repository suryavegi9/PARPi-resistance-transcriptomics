# Source archive audit — 2026-11-25

Inspected the supplied PARPi_resistance_project.zip archive without executing any biological analysis.

## Inventory

- 966 files; 190,735,231 bytes compressed; 275,971,335 bytes uncompressed.
- 54 top-level R script files covering 53 distinct script numbers.
- Top-level numbered scripts 07, 12, 16, 18, 23, 29, 40, 45, 50 and 52 are absent. Numbering gaps do not alone establish missing analysis steps; verify dependencies.
- The D1 root script 01 differs from scripts/01: the root version contains appended exploratory commands. Preserve both variants for review.
- The two top-level copies of script 25 are byte-identical.
- Original scripts use workstation-specific absolute paths. A narrow automated path adaptation is not equivalent to a tested reproduction.
- R history, temporary R data, raw data, duplicate export trees, nested archives and manuscript drafts are not suitable for indiscriminate repository import.
- The top-level R script credential-pattern scan found no matches; this is not a comprehensive secret audit.
- Only D2 and D3 sessionInfo files in the canonical project results tree are nonempty; the final synthesis sessionInfo file is empty.

## Final synthesis

The archived Figure8_FINAL_checkpoint.csv contains accession counts 0/4/24/12/10, model-unit counts 0/5/35/5/5, treatment Spearman rho -0.255702280912365 and 24 concordant directions. SHA-256 of original file: 366d59d018ea8530a69eb4490ab039151d93e7aaded7b74df3f8645cf489f2ac.

The five frozen source CSVs for script 63 are present:
- results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv
- results/GSE163854/D8_analysis/D8_primary_Hallmark_GSEA.csv
- results/GSE165914/D7_analysis/D7A_M468_Hallmark_GSEA.csv
- results/GSE165914/D7_analysis/D7B_SUM1315_Hallmark_GSEA.csv
- results/GSE202434/D5_integration/D5_H3_P3_Hallmark_master_comparison.csv

The six-accession master is a historical dependency, not the final scientific endpoint. The archive contains the original Figure 8 artwork; final renumbered composite artwork and publication-build S1–S3 files are not established by this archive.

No DE, normalization, GSEA or synthesis script was executed. Repository remains private; release is not authorized by this audit.
