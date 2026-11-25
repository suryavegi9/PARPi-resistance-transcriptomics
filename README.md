# PARPi resistance transcriptomics

Cross-dataset reproducibility of transcriptional programs in acquired PARP-inhibitor resistance across independent cancer models.

**Status:** Private publication staging repository. The manuscript is in preparation. The central synthesis has passed numerical execution verification from frozen Hallmark inputs. Full upstream raw-data-to-results execution has not been established; release preparation is ongoing.

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

The following targeted reconstruction and execution checks have been completed:
 
- D4 primary differential expression: author-executed [reconstructed R script](scripts/18_GSE342642_primary_RECONSTRUCTED.R) reproduced filtering and all six gene-level statistics for 14,277 genes within numerical tolerance. [Verification record](provenance/D4_reconstruction/20261125_132322/VERIFICATION.md).
- Historical D1–D6 Hallmark master: author-executed [reconstruction](scripts/FINAL_D1_D6_master_RECONSTRUCTED.R) reproduced 50 pathways and 49 columns from 11 archived Hallmark tables. [Verification record](provenance/D1_D6_reconstruction/20261125_133657/VERIFICATION.md).
- Final synthesis: the author's isolated execution of script 62 (13 documented syntax repairs) and script 63 passed all 46 column comparisons against frozen reference tables, including five script 63 tables and 11 overlapping script 62 accession fields. Original input fingerprints remained unchanged. The complete historical script 62 output was not available for comparison. The release checklist records the execution pass; exact final evidence bundle transfer and inventory remain pending.
 
Final classification of the 50 Hallmark pathways:
 
| Classification | Eight accessions | Twelve model/contrast units |
|---|---:|---:|
| Strong recurrent | 0 | 0 |
| Directional | 4 | 5 |
| Context-dependent | 24 | 35 |
| Dataset/unit-specific | 12 | 5 |
| Inconsistent | 10 | 5 |
 
The 12 units are not independent studies. These counts describe the specified recurrence criteria, not the absence of shared molecular resistance mechanisms. The historical six-accession integration is not the final result.
 
See the [reproduction status](provenance/REPRODUCTION_STATUS.md), [historical dependency audit](provenance/DEPENDENCY_AUDIT.md), and [D4 public input verification](provenance/D4_PUBLIC_INPUT_VERIFICATION.md). Full upstream execution order, acquisition instructions, environment records and portable validation files require further review. Do not run every numbered script indiscriminately or overwrite frozen outputs.

## Citation and licensing

Manuscript in preparation; no article DOI or archived code DOI has been assigned. Cite the software using [CITATION.cff](CITATION.cff). Original code is licensed under [MIT](LICENSE); eligible original derived tables and figures are CC BY 4.0, subject to rights held by the author. See [LICENSING.md](LICENSING.md) for scope and exclusions. Third-party materials remain subject to their original terms.

Author: Vegi Surya Vara Prasad.
