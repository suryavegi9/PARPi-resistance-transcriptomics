# D4 public input verification

Date: 2026-11-25.

Source: NCBI GEO GSE342642, publicly released 2026-10-07 according to the author-provided official GEO screenshots. BioProject PRJNA1508690; 12 regional specimens from four independent PDX tumors (two control, two resistant). No linked article is shown in the GEO record.

Public supplementary file: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE342nnn/GSE342642/suppl/GSE342642_Data.txt.gz

Compared the newly downloaded file supplied by the author, `GSE342642_Data (1).txt.gz`, with `PARPi_resistance_project/data_raw/GSE342642_Data.txt.gz` inside the original project ZIP. Both compressed files are exactly 4,907,012 bytes, SHA-256 `28610d9c5aadcdb40e11e88d2bff3546da0bb3712c2f78969a3178992512034a`. Both decompressed files are exactly 11,302,731 bytes, SHA-256 `77631b6688c6205fc144e4e5a75d28b623d429e9195319eb1f3543fae5c024fc`. Direct byte equality passed at both levels. The header contains FPKM, TPM and read-count columns for the 12 specimens.

This verifies the input file against the author's newly downloaded public copy; the assistant did not independently download the file from NCBI. No expression analysis was rerun. D4 public accessibility and input identity are resolved; historical primary-model script recovery and end-to-end execution verification remain separate open tasks. Earlier staging warnings are retained as historical records and superseded for public accessibility by this record.
