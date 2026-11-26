# Portable execution evidence
`INPUT_MANIFEST.csv` fingerprints the two portable scripts, six frozen inputs and five frozen reference tables. It is separate from the original wrapper's MD5 manifest.

`clean_clone_20261126/TEST_RECORD.json` identifies the exact committed code executed from a local clean clone without hardlinks. All 46 columns passed; all five script-63 CSVs were byte-identical; sources and clone working tree remained unchanged. Both attempts to use an existing output directory or an in-repository output were rejected. The run executed no differential expression, GSEA or model fitting.

The historical script-63 PDF device emitted Greek kappa/rho conversion warnings. This limits validation PDF typography, not numeric equivalence; final publication PNGs are separate recovered assets. An initial manifest check identified line-ending differences in six unchanged documentation files. Their canonical Git bytes were restored; the final manifest is checked against a fresh committed checkout.

The runner and its scientific inputs are unchanged after the tested commit. Subsequent commits add documentation, integrity records and manifests. The full release is checked again after push; the final response records the remote commit and counts.
