# Final synthesis execution verification

Author run: 2026-11-25 16:50:53 local time. The supplied ZIP includes an earlier failed run (16:42:29, syntax failure) and the successful run (16:50:53). Only the latter is used as successful evidence.

Independent checks against the original project ZIP: ZIP CRC passed, all 13 original input/reference MD5 fingerprints matched, all isolated input copies matched. Exactly 13 malformed split double-bracket indices in script 62 were repaired; the resulting source matches the uploaded repaired source line-for-line. Working directory redirected only in the isolated execution wrapper. The portable repaired script is archived at scripts/62_FINAL_D1_D8_accession_integration_SYNTAX_REPAIRED.R (Git blob fea21f290a9d8033cc70b07d80340972defe3a80); its only additional change from the author repair is portable working-directory setup.

All five regenerated script 63 CSVs are byte-for-byte identical to their archived originals. Script 62's 11 overlapping fields match the archived accession source table exactly after pathway alignment. All 46 column comparisons passed; maximum exported numerical difference zero. The original full historical script 62 output is unavailable. The report says 12 overlapping fields, a reporting typo; the actual mapping includes 11 and 46 total comparisons. The regenerated historical Figure 8 has the previously noted legend overlap and is not a replacement for the corrected publication Figure 5. Session record: R 4.6.1, Windows 11 x64, ggplot2 4.0.3, patchwork 1.3.2.

## Evidence SHA-256

- Original uploaded ZIP: d08e13b67744702c7807fdec012eb2174f8de3784f29e1d1f5fc1bbd55b1e7ad
- Exact executed wrapper: 5e25642ee558e91bb08c4ad2c7377028232cb0269f23e10292dea81f9dead009
- Repaired script 62 (original workstation path): 93639ab13374351ef76dadfa69f9973bbfdcb6ff402a09b9c57fc975f44a8b4b
- Portable repaired script 62: f417988855582dfddf3f5ecc859aa6a95b52d874d322a201a8ef806eab260a99
- final_comparison.csv: b92e155f56a5f57b2507ca868ca0bf5b2ce2444e168b1f76c95979b89116efd3
- FINAL_execution_report.txt: d164e3a83c0840375d81644dbf4bafb36c3171d5a843eeb62018282db5bf5d99
- input_fingerprints.csv: cdddf78ace1911013e9a52b5b7933016724e022260546984a05cc743955cd25f
- script62_syntax_repair_log.csv: a2db2580ba343e7cbf651262cef890eb7b0f11d1739e1d05ed64c2272ae0d028
- sessionInfo.txt: 8c120616433dd776feeb637324f066557c50df7740d6e52343514bc074eaa592

Exact comparison and report plus portable corrected script are committed. Exact wrapper, full evidence files and regenerated figure are preserved in a separate downloadable evidence bundle, pending repository transfer and review for personal paths. Original scripts and frozen results unchanged. Central synthesis numerical reproduction demonstrated, not complete raw-to-results DE/GSEA reproduction. No new biological analyses or public release.
