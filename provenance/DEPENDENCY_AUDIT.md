> Historical record; current completeness is documented in the root RELEASE_AUDIT.md. Pending statements below describe the earlier snapshot.

# Static reproducibility dependency audit — 2026-11-25

Scope: inspected the original project ZIP, archived scripts and input/output paths; no R expression model, enrichment or integration execution. This is not an end-to-end reproducibility certification.

## Final synthesis inputs
- Script 62 (final accession integration) reads the frozen D1–D6 master, searches D7 for a 50-row merged table with pathway, M468_NES, M468_FDR, SUM1315_NES and SUM1315_FDR, and reads a D8 primary-vote or primary Hallmark table. The previously omitted archived D7_accession_Hallmark_integration.csv has been restored in results/GSE165914/D7_analysis/. The Git blob matches the original archived bytes exactly (17b81ef9610d9fb1e973bbf641e1d500f22d709b; SHA-256 e31d26e93fc2f6e1d5149386aaf0b6c56740a0d53da311ec472943ac184ce40c). The file contains 50 distinct pathways, 33 concordant directions and 6 concordant significant pathways.
- Script 63 (final synthesis figure) reads the frozen D1–D6 master, separate D7A/D7B primary Hallmark tables, D8 primary Hallmark table and D5 direct exposure comparison; all five are present (see IMPORT_VERIFICATION.md). It requires ggplot2 and patchwork. The extra D7 integration table has a different basename and does not conflict with script 63's exact filename search.
- Both scripts write outputs. Test them in a separate working copy, not over the frozen publication outputs. These are static input checks, not successful execution tests.

## Upstream workflow map (script numbers)
- D1: 01–06. Preserve the shorter 01 variant separately. The root 01 includes exploratory additions; verify intended execution boundary.
- D2: 08–11.
- D3: 13–15.
- D4: 17, missing primary regional model, 19–22. Scripts 17, 19 and 21 require interactive file.choose(). Use the verified GSE342642_Data.txt.gz input. Script 19 consumes primary_DEG/GSE342642_primary_DEG_results.csv; script 20 consumes primary_DEG/GSE342642_moderated_t_ranking.csv. Script 21 also reads the primary DE and primary Hallmark outputs. Script 21 documents the corresponding two-pass blocked voom and duplicateCorrelation model for the low-depth-region exclusion, not the missing original primary script.
- Historical D1–D4: 24–28. Historical only, not the final conclusion.
- D5: 30–39.
- D6: 41–44.
- Historical D1–D6: 46–49 consume saved historical integration materials. No original writer of FINAL_D1_D6_Hallmark_master.csv was identified in the numbered scripts; do not infer which missing script number generated it.
- D7: 51, 53–55.
- D8: 56–60.
- D9 eligibility investigation: 61, excluded from primary synthesis.
- Final: 62–63, separate final endpoints based on frozen pathway outputs.

## Source gaps
Missing numbered scripts: 07, 12, 16, 18, 23, 29, 40, 45, 50, 52. A gap in numbering is not itself evidence of a missing required analysis. No original D4 primary script was found in root scripts, scripts/ copies, export script folders or the two archived R histories (11 lines and empty). Prior conversation retrieval did not recover the exact code either. Original primary results, moderated-t ranking, tumor-aware filtering and both voom plots are present in the original ZIP. The primary analysis summary records 15,623 before filtering, 14,277 after filtering, within-tumor correlation 0.39477637734058, 1,984 FDR-significant genes, 588 up and 453 down large-effect genes. Input accessibility and identity verified separately (D4_PUBLIC_INPUT_VERIFICATION.md). Do not present a newly constructed model script as historical code.

## Environment and raw-data access
The archive contains raw input files for the eight included accessions; these are not all included in the GitHub staging import. Upstream scripts require the corresponding data_raw paths and, in some cases, intermediate results. The original project archive must be retained. Nonempty canonical sessionInfo records are available for D2/D3 only; other sessionInfo-named records in the archive are empty, apart from duplicate exports and the export status document. Dependencies observed in library()/require() calls include DESeq2, edgeR, limma, clusterProfiler, msigdbr, org.Hs.eg.db, enrichplot, ggplot2, ggrepel, patchwork, dplyr, tidyr, readxl, data.table, R.utils and BiocManager; this is not an exhaustive or pinned dependency list. Rscript is not installed in the audit runtime.

## Release gate
Final integration input preflight is substantially complete. Raw-to-results reproduction remains unverified, and two important historical source gaps remain. Recover original source code where possible, or with explicit author approval prepare and independently test clearly labeled reconstruction scripts against saved outputs. Archive the historical outputs and disclose limitations. Do not run all scripts indiscriminately, overwrite frozen outputs, or make the repository public solely on the basis of this audit.
