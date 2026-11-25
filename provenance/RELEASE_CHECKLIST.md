# Publication release checklist

This is a staging checklist, not a certificate of reproducibility.

## Before public release

- [x] Inventory and preserve available D1–D8 scripts; missing historical scripts documented.
- [ ] Review portable copies; replace machine-specific paths without altering scientific logic.
- [ ] Verify exact sample mappings and input acquisition instructions.
- [x] Verify D4 public accessibility and archived input identity (D4_PUBLIC_INPUT_VERIFICATION.md).
- [x] Document D6 article/accession discrepancy; cite the actual GEO record separately.
- [x] Verify final synthesis input paths and restore missing D7 merged table (DEPENDENCY_AUDIT.md).
- [x] Transparently address both missing historical scripts: author-executed D4 primary DE and D1–D6 master reconstructions numerically validated; exact executed code and verification records archived.
- [x] Execute final synthesis scripts 62 (13 syntax repairs) and 63 in isolation; 46/46 columns passed, five script 63 CSVs byte-identical; verification record and portable corrected script archived.
- [ ] Transfer remaining exact final synthesis evidence bundle after privacy review.
- [x] Document remaining upstream execution limitations in REPRODUCTION_STATUS.md.
- [ ] Assess whether further upstream execution checks are needed and feasible before release; do not claim end-to-end reproduction.
- [ ] Archive frozen outputs and verified supplementary Tables S1–S3.
- [ ] Assemble and visually audit final Figures 1–5 and supplementary Figures S1–S4.
- [ ] Archive available R session information; explicitly document missing historical records.
- [ ] Audit final manuscript claims, citations and AI-use disclosure.
- [x] Author approved MIT code license and CC BY 4.0 for eligible original outputs; LICENSE and LICENSING.md added.
- [x] Add CITATION.cff with author and ORCID; version and DOI to be added at release.
- [ ] Complete per-file third-party redistribution rights audit, including gene-set resources.
- [ ] Check for credentials, personal filesystem paths and unintended large files.
- [ ] Author approves public release; create tagged release and archive DOI if desired.

The D4 primary model was rerun during reconstruction validation. No new DE or GSEA was executed during the final synthesis validation. Frozen results were not overwritten.
