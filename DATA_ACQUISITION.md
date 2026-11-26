# Raw-data acquisition
The [machine-readable manifest](metadata/DATA_ACQUISITION_MANIFEST.csv) identifies all 11 required supplementary inputs for eight accessions: exact official filenames, HTTPS URLs, local destinations, original SHA-256 values and decompression rules. Official GEO FTP listings and SOFT metadata were retrieved on 25 November 2026. Raw matrices and full GEO contact metadata are not redistributed.

Run `python scripts/acquire_data.py --output ../PARPi_raw_inputs` from the repository. The target must be new and outside this repository. The script downloads only; it does not execute analysis. A failed checksum is a stop condition, not permission to update a frozen reference. Hashes refer to destination bytes; three files are decompressed before comparison. Retain gzip files for inputs marked `none`; the historical R readers handle these. D7 uses the author's counts, not the alternative NCBI generated matrix.

Sample mappings, contrasts, preprocessing and historical execution ranges are in [DATASET_CONTRAST_MANIFEST.csv](metadata/DATASET_CONTRAST_MANIFEST.csv). D6 reads sheet 1 of the author workbook and retains the expected-count/TPM distinction. D5 supplied normalized counts and TPM are not raw sequencing counts. No conversion to a common expression scale is warranted across datasets. Actual sample IDs are preserved in the linked frozen metadata files.

GSE287274 was excluded under the original eligibility assessment because necessary expression and sample information could not be recovered reproducibly. Its recovery script is historical evidence, not an extra included dataset. D8 serial-selection trajectory is exploratory and does not add accession votes. The D6 candidate article's GSE131231/GSE191231 discrepancy remains a citation limitation; the analyzed workbook and acquisition accession are GSE191231.

See [rights and attribution](THIRD_PARTY_RIGHTS.md). Download access does not transfer ownership of source data.
