# Release candidate audit â€” 26 November 2026

**NOT READY â€” the author must approve historical-email exposure or choose a separate clean-history release repository.** Technical checks below distinguish execution evidence from historical limitations. The repository remains private.

| Item | Status | Evidence / limitation |
|---|---|---|
| Original research project | PASS | All 1,018 baseline paths, sizes and SHA-256 values unchanged; ORIGINAL_PROJECT_INTEGRITY.json. |
| Missing numbered scripts | PASS | Eleven candidates recovered and fingerprinted; script 64 and exact verification wrapper recovered separately. Candidates 18/45 are not represented as proven historical executed sources. |
| D4 reconstruction | PASS | Exact executed reconstruction and evidence; 15,623/14,277 genes, 1,984 FDR-significant, 588 up/453 down; six statistics and correlation 0.39477637734058. No release DE rerun. |
| Historical D1-D6 | PASS | Eleven frozen Hallmark inputs; all 50 pathways/49 columns verified; complete nonempty evidence. |
| Final author synthesis | PASS | Successful 20261125_165053 wrapper/report/fingerprints/repair log/tables/artwork/session; 46/46 comparisons. Failed run excluded. |
| Script-62 identity and erratum | PASS | Exact path-bearing repaired source retained privately; public portable source fingerprinted. Eleven overlap fields plus 35 script-63 columns; original report preserved. |
| Portable synthesis | PASS | Clean committed clone: 46/46; all five script-63 tables byte-identical; unchanged inputs and clean worktree. Both overwrite guards tested. See provenance/portable_synthesis/clean_clone_20261126/TEST_RECORD.json. |
| Final Figures 1-5 and S1-S4 | PASS | Nine exact PNG objects from latest available author-review manuscript; source SHA-256/member map and captions. All visually inspected. |
| Native final artwork assembly source | PENDING, disclosed limitation | Not recovered after local source/archive search. Final PNG assets are present; extraction source provided. Historical plotting scripts are not claimed to recreate exact final compositions. |
| Graphical abstract | PASS | Exact author-generated PDF/PNG/TIFF, R generator and nonempty session; isolated repeat succeeded and PNG pixels match the original exactly; prototype excluded. |
| S1-S3 | PASS | Exact publication CSVs/workbook; every workbook data cell matches CSV; isolated table assembly agrees within 1e-13. S3 derived hierarchy explicitly documented. |
| Locked checkpoints | PASS | Accession 0/4/24/12/10; model 0/5/35/5/5; treatment 0/9/9/2/30; D7/D8, transitions, treatment concordance, hypotheses and glycolysis checked. |
| Raw-data acquisition | PASS | All eleven official downloads across eight accessions matched original archived SHA-256 values. Download-only script; raw matrices excluded. |
| Environment records | PASS within scope | Six recovered nonempty sessions indexed; observed versions recorded. Historical database/environment completeness not claimed. |
| R syntax | PASS within scope | 69/70 parse; sole failure is deliberately retained original script 62. Portable repaired 62 passes. |
| Licensing | PASS within scope | MIT and eligible-output CC BY 4.0 scopes preserved; Hallmark attribution; raw matrices/external libraries excluded. |
| Citation metadata | PASS | CFF 1.2 schema validation; candidate version 1.0.0; no invented DOI/date. |
| Current-tree and history privacy scan | PASS within scope | No credential, personal-path or personal-email matches in scanned file contents/history blobs; binary metadata inspected. Commit email is a separate pending gate. See provenance/SECURITY_REVIEW.json. |
| Historical commit email | PENDING â€” PUBLIC RELEASE BLOCKER | Personal email in 45 original commits. New no-reply settings do not remove it. Explicit author choice required. |
| Redundant archive | PASS | Original archive retained privately; its members preserved in release tree, including original staging note. Removed from current tracked tree, still present in historical Git objects. |
| Clean clone, full manifest and link check | PASS | Scientific checkpoint and relative Markdown link checks pass. Manifest covers every release file except itself; canonical Git bytes used to prevent checkout conversion mismatches. Final remote checkout is verified after push. |
| Visibility, DOI, history rewrite | PASS | No public action, archival deposit, force push or history rewrite authorized/performed. |

Priority 1: author history privacy decision in provenance/PRIVACY_DECISION.md. Technical execution checks pass; private release branch is retained for review while the privacy gate remains open. The absent native artwork assembly source is disclosed as a reproduction limitation, not missing publication PNGs. Historical raw-to-results execution and full historical script-62 master equivalence are not claimed.

Required exact public source identities and inventories are in RELEASE_MANIFEST.csv, provenance/SOURCE_TRANSFER_MANIFEST.csv, provenance/RECOVERED_SOURCES.csv and provenance/SCIENTIFIC_OUTPUT_MANIFEST.csv. The release manifest covers every release file except itself; Git preserves the manifest. Private exact source paths and original input inventories are intentionally not distributed.


## Exact missing or privately retained material
- Public executable/figure/table deliverables are present. Exact original path-bearing repaired 62 is intentionally private: `private_archive/final_synthesis/62_SYNTAX_REPAIRED.R`; its SHA-256 and public counterpart are recorded in HISTORICAL_SOURCE_RECONCILIATION.md.
- Final editable/vector assembly source for manuscript Figures 1-5/S1-S4 was not located; no authoritative filename is known, so none is invented. Recovered PNGs and extraction map are authoritative available artwork.
- Full historical script-62 accession master comparator was not available. The regenerated master is present under `provenance/final_synthesis/20261125_165053/script62/results/FINAL_D1_D8_integration/accession_level/FINAL_D1_D8_accession_Hallmark_master.csv`; equivalence is limited to eleven overlapping fields.
- No complete historical lockfile/database-version record was recovered. Nonempty sessions are indexed; empty session files are not treated as evidence.

## Release disposition
Prepared v1.0.0 metadata and release notes are not a public release. Candidate changes are committed on `release/v1.0.0-preparation`; main is held unchanged while the author decision is pending. This preserves the original private history and leaves the choice of public repository identity open. No force push, tag, public visibility change or DOI publication is performed.
