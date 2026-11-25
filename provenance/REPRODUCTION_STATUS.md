# Reproduction status — 25 November 2026

This record supersedes the pending-test status in the earlier static DEPENDENCY_AUDIT.md and IMPORT_VERIFICATION.md. Those records remain historical audit evidence.

## Verified stages

| Stage | Evidence | Scope |
|---|---|---|
| D4 public source identity | D4_PUBLIC_INPUT_VERIFICATION.md | Author-provided fresh GEO download matched the archived input byte-for-byte |
| D4 primary model | D4_reconstruction/20261125_132322/VERIFICATION.md | Author-run reconstruction; independently compared exported tables against historical outputs, 14,277 retained genes |
| Historical D1–D6 master | D1_D6_reconstruction/20261125_133657/VERIFICATION.md | Author-run reconstruction from 11 pathway tables; 50 pathways, 49 columns, independent exported-table comparison |
| Final D1–D8 synthesis | Author execution report and release checklist | Isolated executions of scripts 62 and 63, 46/46 column comparisons; complete final evidence inventory and transfer remain pending |

Both missing historical scripts are addressed by clearly labeled reconstructions, not recovered original source. Script 62 required 13 indexing syntax repairs. Original archived script and frozen inputs were preserved. The reported script 62 comparison covers 11 overlapping accession-table columns (including pathway ID), not the 12 mistakenly stated in the console report. The complete historical script 62 master is unavailable, so full historical master equivalence is not claimed. Script 63 comparison covers 35 columns across five frozen CSVs. No DE or GSEA was performed during the final integration test.

## Expected frozen final checkpoints

- Eight accessions: strong recurrent 0; directional 4; context-dependent 24; dataset-specific 12; inconsistent 10.
- Twelve model/contrast units: 0, 5, 35, 5, 5, respectively. Units are not independent studies.
- D7: 33/50 concordant directions; 6 concordant-significant pathways.
- D8: 15 primary significant Hallmarks.
- Direct versus observational treatment-state comparison: Spearman rho approximately -0.2557; 24/50 matching directions.

## Reproduction instructions and limitations

The validated reconstruction scripts are under scripts/. They require the original input files and frozen references at their documented relative paths. They produce timestamped validation outputs rather than overwrite historical results. Do not execute all numbered scripts as a single pipeline.

For final synthesis, use isolated copies and the exact validated wrapper, with the recorded input hashes and syntax-repair log. The precise portable wrapper/corrected-script paths and remaining validation bundle files must be confirmed in the final file inventory before publishing runnable instructions.

Complete upstream raw-data acquisition, all DE and GSEA reruns, historical package and gene-set version recovery, and clean-clone execution have not been established. Do not represent this repository as fully reproducible from raw sequencing data. The author retains original archives. Before public release, inspect repository contents and history for private paths, credentials, unnecessary raw data, and third-party rights; document source terms and exclusions. No DOI has been minted and no public release has been made.
