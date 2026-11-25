# D4 primary-model reconstruction verification

Author execution: 2026-11-25 13:23:22 (local timestamp). This is a reconstruction, not recovered historical source code.

The uploaded validation ZIP passed CRC verification. Independently read both exported reconstructed CSVs and the original project ZIP. Joined tables by unique Gene_ID; 15,623 filter IDs and all 14,277 retained DE IDs matched. All four median CPM columns were identical at parsed exported precision; all retained flags matched. All six DE statistic columns (logFC, AveExpr, t, P.Value, adj.P.Val, B) matched at parsed exported precision; all status and FDR membership matched. FDR count 1,984; 588 up and 453 down large-effect genes. The author's report records matching consensus correlation 0.39477637734058 and original files unchanged by MD5. No GSEA was executed.

Execution environment reported in the supplied sessionInfo: R 4.6.1, Windows 11 x64, edgeR 4.10.7, limma 3.68.5, statmod 1.5.2 (loaded namespace). These versions refer to the reconstruction run, not necessarily the original historical analysis.

## Evidence inventory (SHA-256)

- Uploaded 20261125_132322.zip: 2385dbce7064493a1ac8dabab0748e9eb965794ee43c852e4c862c73c9905c47
- D4_reconstruction_report.txt: 553c123ae7fb4f8258deafd553ae8257f87d169e524208d85c029d917a289266
- numeric_comparison.csv: acc28020b05657f60e1030baa61fa6f5541e838e43a1122f949d2c0649a06928
- reconstructed_DE.csv: 3ab750d69f43190579a99b90fbc2d7f61bcee8d015d8feb505662b92fcaada09
- reconstructed_filter.csv: 945f38ac3bc0337adeb9e4bd8036e86b90a697f7a073658732d887735e316e35
- sessionInfo.txt: 4fa6af28db5767533363b24f6fa7cd7363ea57f63517e24a8bcdf71e71df2e21

The report and numeric comparison were committed with exact original Git blob IDs. The two large reconstructed CSVs and full sessionInfo are retained in the author's evidence ZIP, not yet uploaded to GitHub. The exact executed R script was not included in the ZIP. A previously supplied candidate reconstruction exists, but the executed copy needs confirmation before archiving it as the validated source.

D4 primary DE numerical reproduction passes. D4 downstream enrichment execution, historical D1–D6 master reconstruction and final integration execution remain separate tasks. Repository stays private. No locked results or manuscript claims changed.
