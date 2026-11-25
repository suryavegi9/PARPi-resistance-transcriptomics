# Analysis scripts
 
The staged copies of 53 numbered R scripts and one alternative D1 script were imported with their original and staged fingerprints documented. Historical scripts are retained; portability changes and later syntax repairs must be distinguished from changes to scientific logic.
 
Two missing historical sources have been addressed by explicitly labeled and numerically validated reconstructions:
 
- [18_GSE342642_primary_RECONSTRUCTED.R](18_GSE342642_primary_RECONSTRUCTED.R): D4 primary regional model; validated on all 14,277 retained genes.
- [FINAL_D1_D6_master_RECONSTRUCTED.R](FINAL_D1_D6_master_RECONSTRUCTED.R): historical master; 50 pathways, 49 columns, 11 original pathway inputs.
 
The author executed the final synthesis scripts in isolated copies: script 62 required 13 syntax repairs, while script 63 executed successfully. The comparison passed 46/46 columns. Preserve the original 62 separately from its repaired version. The complete historical 62 output was unavailable; the comparison uses 11 overlapping fields of the frozen accession table. Exact final validation bundle archival and path inventory remain outstanding.
 
See [reproduction status](../provenance/REPRODUCTION_STATUS.md) and [release checklist](../provenance/RELEASE_CHECKLIST.md). These checks do not certify full upstream raw-to-results reproduction. Historical environment records are incomplete. Do not run all scripts indiscriminately, or write over frozen outputs.
