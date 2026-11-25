# PARPi resistance transcriptomics

Cross-dataset reproducibility of transcriptional programs in acquired PARP-inhibitor resistance across independent cancer models.

**Status:** Private publication staging repository. The manuscript is in preparation. This repository is not yet a complete, independently runnable reproduction package.

## Study design

Eight GEO accessions were analyzed separately using expression-data-appropriate workflows. Ten Hallmark hypotheses were fixed in the discovery dataset before external validation. The primary synthesis uses one conservative vote per accession; a secondary analysis uses 12 model/contrast units. Treatment-state comparisons and the D8 serial-selection trajectory are complementary analyses, not additional independent accession votes.

| ID | GEO accession | System |
|---|---|---|
| D1 | [GSE235980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235980) | UWB1.289, BRCA1-deficient and restored backgrounds |
| D2 | [GSE117765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117765) | PEO1 |
| D3 | [GSE153867](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153867) | A2780 |
| D4 | [GSE342642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342642) | HGSOC PDX; public input identity verified |
| D5 | [GSE202434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202434) | HCC1806 and PSN1, talazoparib maintenance/withdrawal |
| D6 | [GSE191231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191231) | Recurrent HGSOC PDX after olaparib selection |
| D7 | [GSE165914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165914) | MDA-MB-468 and SUM1315 |
| D8 | [GSE163854](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163854) | PH039 PDX, serial niraparib selection |

GSE287274 was excluded because the required sample-level metadata and expression data could not be recovered reproducibly at the time of eligibility assessment.

## Planned repository layout

- `scripts/`: original analysis scripts and reviewed portable copies
- `metadata/`: exact sample mappings, contrasts and source records
- `results/`: frozen derived analysis outputs
- `figures/`: final main Figures 1–5 and supplementary Figures S1–S4
- `supplementary/`: Tables S1–S3
- `environment/`: archived session information and dependency records
- `provenance/`: source inventory, validation records and release checklist

See [release checklist](provenance/RELEASE_CHECKLIST.md) for the status of these materials. A listed directory does not imply its scientific contents have already been uploaded.

## Reproduction

The final synthesis scripts and their frozen pathway inputs are present. See [dependency audit](provenance/DEPENDENCY_AUDIT.md) and [D4 public input verification](provenance/D4_PUBLIC_INPUT_VERIFICATION.md). Full upstream execution instructions require historical source recovery and testing. Do not treat the historical six-accession integration as the final result. No analysis has been rerun during repository initialization.

## Citation and licensing

Manuscript in preparation; no article DOI or archived code DOI has been assigned. Code licensing and release metadata remain subject to author approval. Original data remain subject to their source terms.

Author: Vegi Surya Vara Prasad.
