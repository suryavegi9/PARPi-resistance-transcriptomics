# ============================================================
# GSE191231 - DATASET 6
# EXPRESSION-SPECIFIC AUDIT
#
# Goals:
# 1. Audit RSEM expected counts separately
# 2. Audit TPM separately
# 3. Check library sizes and TPM sums
# 4. Inspect sample correlations / PCA
# 5. Investigate duplicated Ensembl IDs and gene symbols
#
# NO DIFFERENTIAL EXPRESSION YET
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(readxl)


# ============================================================
# 1. LOAD DATA
# ============================================================

file_path <- paste0(
  "data_raw/",
  "GSE191231_Hyun_Analysis_tpm_edger.xlsx"
)

dat <- read_excel(
  file_path,
  sheet = 1
)

dat <- as.data.frame(
  dat,
  check.names = FALSE
)


# ============================================================
# 2. IDENTIFY EXPRESSION COLUMNS
# ============================================================

count_cols <- grep(
  "expected_count$",
  colnames(dat),
  value = TRUE
)

tpm_cols <- grep(
  "\\.TPM$",
  colnames(dat),
  value = TRUE
)


cat("\n========================================\n")
cat("EXPRESSION COLUMN QC\n")
cat("========================================\n\n")

cat(
  "Expected-count columns:",
  length(count_cols),
  "\n"
)

print(count_cols)

cat(
  "\nTPM columns:",
  length(tpm_cols),
  "\n"
)

print(tpm_cols)


stopifnot(
  length(count_cols) == 4,
  length(tpm_cols) == 4
)


count_prefix <- sub(
  "\\.expected_count$",
  "",
  count_cols
)

tpm_prefix <- sub(
  "\\.TPM$",
  "",
  tpm_cols
)


stopifnot(
  identical(
    count_prefix,
    tpm_prefix
  )
)


# ============================================================
# 3. BUILD MATRICES
# ============================================================

counts <- as.matrix(
  dat[
    ,
    count_cols
  ]
)

tpms <- as.matrix(
  dat[
    ,
    tpm_cols
  ]
)


storage.mode(counts) <- "numeric"
storage.mode(tpms) <- "numeric"


sample_names <- c(
  "Control_834",
  "Control_846",
  "OlapR_827",
  "OlapR_878"
)


colnames(counts) <- sample_names
colnames(tpms) <- sample_names


# ============================================================
# 4. EXPECTED-COUNT QC
# ============================================================

cat("\n========================================\n")
cat("RSEM EXPECTED-COUNT QC\n")
cat("========================================\n\n")


cat(
  "NA:",
  sum(is.na(counts)),
  "\n"
)

cat(
  "Negative:",
  sum(counts < 0, na.rm = TRUE),
  "\n"
)

cat(
  "Zero:",
  sum(counts == 0, na.rm = TRUE),
  "\n"
)


count_integer_fraction <- mean(
  abs(
    counts - round(counts)
  ) < 1e-8,
  na.rm = TRUE
)


cat(
  "Integer-like fraction:",
  round(
    count_integer_fraction,
    6
  ),
  "\n"
)


cat("\nExpected-count quantiles:\n")

print(
  quantile(
    counts,
    probs = c(
      0,
      0.01,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


count_library_sizes <- colSums(
  counts,
  na.rm = TRUE
)


cat("\nEstimated-count library sizes:\n")

print(
  count_library_sizes
)


cat(
  "\nLibrary-size max/min ratio:",
  round(
    max(count_library_sizes) /
      min(count_library_sizes),
    3
  ),
  "\n"
)


# ============================================================
# 5. TPM QC
# ============================================================

cat("\n========================================\n")
cat("TPM QC\n")
cat("========================================\n\n")


cat(
  "NA:",
  sum(is.na(tpms)),
  "\n"
)

cat(
  "Negative:",
  sum(tpms < 0, na.rm = TRUE),
  "\n"
)

cat(
  "Zero:",
  sum(tpms == 0, na.rm = TRUE),
  "\n"
)


cat("\nTPM quantiles:\n")

print(
  quantile(
    tpms,
    probs = c(
      0,
      0.01,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


tpm_sums <- colSums(
  tpms,
  na.rm = TRUE
)


cat("\nTPM sums:\n")

print(
  tpm_sums
)


# ============================================================
# 6. COUNT-vs-TPM AGREEMENT WITHIN EACH SAMPLE
# ============================================================

cat("\n========================================\n")
cat("EXPECTED COUNT vs TPM AGREEMENT\n")
cat("========================================\n\n")


for (i in seq_len(ncol(counts))) {
  
  rho <- cor(
    log2(counts[, i] + 1),
    log2(tpms[, i] + 1),
    method = "spearman"
  )
  
  cat(
    sample_names[i],
    ": rho =",
    round(rho, 4),
    "\n"
  )
}


# ============================================================
# 7. SAMPLE-SAMPLE CORRELATION
# ============================================================

log_tpm <- log2(
  tpms + 1
)


tpm_cor <- cor(
  log_tpm,
  method = "spearman"
)


cat("\n========================================\n")
cat("SAMPLE-SAMPLE SPEARMAN CORRELATION")
cat("\nlog2(TPM + 1)\n")
cat("========================================\n\n")

print(
  round(
    tpm_cor,
    4
  )
)


# ============================================================
# 8. PCA QC
#
# Use moderately expressed genes only for visualization.
# This is NOT DE filtering.
# ============================================================

pca_keep <- rowSums(
  tpms >= 1
) >= 2


cat("\n========================================\n")
cat("PCA QC\n")
cat("========================================\n\n")

cat(
  "Genes used for PCA:",
  sum(pca_keep),
  "/",
  nrow(tpms),
  "\n"
)


pca_input <- t(
  log_tpm[
    pca_keep,
    ,
    drop = FALSE
  ]
)


pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)


variance <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )


pca_table <- data.frame(
  
  sample = rownames(
    pca$x
  ),
  
  group = c(
    "Vehicle",
    "Vehicle",
    "Olaparib_selected",
    "Olaparib_selected"
  ),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  stringsAsFactors = FALSE
)


cat(
  "PC1 variance:",
  round(
    variance[1],
    2
  ),
  "%\n"
)

cat(
  "PC2 variance:",
  round(
    variance[2],
    2
  ),
  "%\n\n"
)


print(
  pca_table,
  row.names = FALSE
)


# ============================================================
# 9. EXACT ENSEMBL-ID DUPLICATES
# ============================================================

exact_dup_index <- duplicated(
  dat$gene
) |
  duplicated(
    dat$gene,
    fromLast = TRUE
  )


exact_dups <- dat[
  exact_dup_index,
  c(
    "gene",
    "id",
    count_cols,
    tpm_cols
  )
]


cat("\n========================================\n")
cat("EXACT ENSEMBL-ID DUPLICATES\n")
cat("========================================\n\n")

cat(
  "Rows involved:",
  nrow(exact_dups),
  "\n"
)

cat(
  "Unique duplicated IDs:",
  length(
    unique(
      exact_dups$gene
    )
  ),
  "\n\n"
)


if (nrow(exact_dups) > 0) {
  
  print(
    head(
      exact_dups,
      30
    ),
    row.names = FALSE
  )
}


# ============================================================
# 10. ENSEMBL IDs AFTER REMOVING VERSION NUMBER
# ============================================================

dat$ensembl_base <- sub(
  "\\..*$",
  "",
  dat$gene
)


base_dup_index <- duplicated(
  dat$ensembl_base
) |
  duplicated(
    dat$ensembl_base,
    fromLast = TRUE
  )


base_dups <- dat[
  base_dup_index,
  c(
    "gene",
    "ensembl_base",
    "id",
    count_cols,
    tpm_cols
  )
]


cat("\n========================================\n")
cat("ENSEMBL DUPLICATES AFTER VERSION REMOVAL\n")
cat("========================================\n\n")

cat(
  "Rows involved:",
  nrow(base_dups),
  "\n"
)

cat(
  "Unique duplicated base IDs:",
  length(
    unique(
      base_dups$ensembl_base
    )
  ),
  "\n"
)


# ============================================================
# 11. DUPLICATED GENE SYMBOLS
# ============================================================

symbol_dup_index <- duplicated(
  dat$id
) |
  duplicated(
    dat$id,
    fromLast = TRUE
  )


symbol_dups <- dat[
  symbol_dup_index,
  c(
    "gene",
    "ensembl_base",
    "id",
    count_cols,
    tpm_cols
  )
]


cat("\n========================================\n")
cat("DUPLICATED GENE SYMBOLS\n")
cat("========================================\n\n")

cat(
  "Rows involved:",
  nrow(symbol_dups),
  "\n"
)

cat(
  "Unique duplicated symbols:",
  length(
    unique(
      symbol_dups$id
    )
  ),
  "\n\n"
)


print(
  head(
    symbol_dups,
    30
  ),
  row.names = FALSE
)


# ============================================================
# 12. PROVIDED AUTHOR STATISTICS
#
# Audit only.
# We will NOT use these as our primary results.
# ============================================================

cat("\n========================================\n")
cat("SUBMITTER-PROVIDED STATISTICS AUDIT\n")
cat("========================================\n\n")


cat(
  "logFC range:",
  range(
    dat$logFC,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "PValue range:",
  range(
    dat$PValue,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR range:",
  range(
    dat$FDR,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Absolute logFC > 100:",
  sum(
    abs(dat$logFC) > 100,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 13. SAVE AUDIT TABLES
# ============================================================

outdir <- "results/GSE191231/audit"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  pca_table,
  file.path(
    outdir,
    "GSE191231_PCA_audit_coordinates.csv"
  ),
  row.names = FALSE
)


write.csv(
  exact_dups,
  file.path(
    outdir,
    "GSE191231_exact_Ensembl_duplicates.csv"
  ),
  row.names = FALSE
)


write.csv(
  base_dups,
  file.path(
    outdir,
    "GSE191231_base_Ensembl_duplicates.csv"
  ),
  row.names = FALSE
)


write.csv(
  symbol_dups,
  file.path(
    outdir,
    "GSE191231_symbol_duplicates.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    sample = sample_names,
    estimated_count_sum = count_library_sizes,
    TPM_sum = tpm_sums
  ),
  file.path(
    outdir,
    "GSE191231_library_QC.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("GSE191231 EXPRESSION AUDIT COMPLETE\n")
cat("========================================\n")