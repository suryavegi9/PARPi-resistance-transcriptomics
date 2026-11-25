# Historical D1-D6 master reconstruction verification

Author execution: 2026-11-25 13:36:57 (local timestamp). This is reconstructed source, not recovered historical code.

The uploaded validation ZIP passed CRC testing. All eleven Hallmark input MD5s and the archived reference master MD5 were independently checked against the original project ZIP; all twelve matched. Reconstructed and archived master tables were independently read and aligned by unique pathway identifier. All 50 pathways and 49 columns matched, including all numerical values at parsed CSV precision and exact categorical and logical values. The author's R report records maximum pre-export numerical difference 4.88498130835069e-15, with a tolerance of 1e-12 for NES/FDR. All 49 columns pass. Historical class counts: strong 1, directional 12, context-dependent 15, dataset-specific significant 12, inconsistent/weak 10. These are historical six-accession results, not the final eight-accession classifications.

The uploaded R script and originally supplied candidate have identical lines after stripping leading/trailing whitespace. The exact author-uploaded script is archived at scripts/FINAL_D1_D6_master_RECONSTRUCTED.R (Git blob 00d3cae2d510276dfddb82fe6fb0267e412455f5). The report, comparison, class counts and input fingerprints were committed with original Git blob fingerprints. The reconstructed master and complete sessionInfo remain in the author's evidence ZIP, not uploaded in this commit. R 4.6.1 on Windows 11 x64 recorded for this reconstruction.

## SHA-256 evidence inventory

- Uploaded ZIP: af6fe4906c2456544464b3b52861cf1167f23555c0b2aa5ffc7b0ac4573ceb7c
- Executed R script: c80233697864a74475c91e0556e8fc4cbc18924c54101367c87c9d8aa1402acd
- classification_counts.csv: 514b8cf01a9efcef81c47158816ecb1898c49c0621c03c9f020c39471b5c902d
- column_comparison.csv: 49c3478f213cff0af4aa1fe093caf54b9e4fa2155026e35eb0548e4a441f2261
- D1_D6_reconstruction_report.txt: 0b2b7e034a99605644d90bce526388b7a8bde032a6c63525b2b7a5bed69abd9e
- input_fingerprints.csv: 7c3994b8080cc55135224c39e5907f915a805894ee40d2564b6a8a71e2039ee7
- reconstructed_Hallmark_master.csv: 4c8fa2dfab0b386696234369b1b834b2086a770510023395bb8adc9755babb34
- sessionInfo.txt: c02e8c85713e9dd9676caa04392822830e140aadb0ad87095f0aaa75d5ad414b

No DE or GSEA executed. Frozen original outputs unchanged according to author MD5 checks. Both historical source gaps have validated reconstruction scripts, but historical originals remain missing. Final scripts 62 and 63 still require isolated execution checks. Complete raw-to-results execution is not established. Repository remains private.
