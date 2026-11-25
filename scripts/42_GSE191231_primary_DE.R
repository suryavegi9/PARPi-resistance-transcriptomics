# ============================================================
# GSE191231 - DATASET 6
# PRIMARY DIFFERENTIAL EXPRESSION + QUANTIFICATION SENSITIVITY
#
# PRIMARY:
#   TPM >= 1 in >= 2 samples
#   log2(TPM + 1)
#   limma-trend
#
# SENSITIVITY:
#   RSEM expected counts
#   TMM + voom + limma
#
# Biological contrast:
#   Olaparib-selected recurrent/resistant PDX
#   versus
#   Vehicle-control PDX
#
# IMPORTANT:
# - Submitter logFC/PValue/FDR are NOT used.
# - Exact duplicated Ensembl rows are deduplicated.
# - Gene-symbol duplication is resolved only later for GSEA.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


library(readxl)
library(limma)
library(edgeR)


# ============================================================
# 1. INPUT
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
# 2. EXPRESSION COLUMNS
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


stopifnot(
  length(count_cols) == 4,
  length(tpm_cols) == 4
)


sample_names <- c(
  "Control_834",
  "Control_846",
  "OlapR_827",
  "OlapR_878"
)


group <- factor(
  c(
    "Control",
    "Control",
    "OlapR",
    "OlapR"
  ),
  levels = c(
    "Control",
    "OlapR"
  )
)


# ============================================================
# 3. ENSEMBL IDENTIFIERS
# ============================================================

dat$ensembl_id <- dat$gene

dat$ensembl_base <- sub(
  "\\..*$",
  "",
  dat$gene
)


# ============================================================
# 4. VERIFY EXACT DUPLICATED ENSEMBL FEATURES
# ============================================================

dup_ids <- unique(
  dat$ensembl_id[
    duplicated(dat$ensembl_id) |
      duplicated(
        dat$ensembl_id,
        fromLast = TRUE
      )
  ]
)


cat("\n========================================\n")
cat("EXACT ENSEMBL DUPLICATE AUDIT\n")
cat("========================================\n\n")

cat(
  "Duplicated Ensembl IDs:",
  length(dup_ids),
  "\n"
)


# Expression columns used to check duplicate equality

expr_cols <- c(
  count_cols,
  tpm_cols
)


dup_check <- data.frame(
  ensembl_id = character(),
  n_rows = integer(),
  identical_expression = logical(),
  stringsAsFactors = FALSE
)


for (g in dup_ids) {
  
  z <- dat[
    dat$ensembl_id == g,
    ,
    drop = FALSE
  ]
  
  
  m <- as.matrix(
    z[
      ,
      expr_cols,
      drop = FALSE
    ]
  )
  
  
  # Compare every duplicated row against row 1
  
  identical_expression <- all(
    apply(
      m,
      2,
      function(x) {
        max(
          abs(
            x - x[1]
          ),
          na.rm = TRUE
        ) < 1e-8
      }
    )
  )
  
  
  dup_check <- rbind(
    dup_check,
    data.frame(
      ensembl_id = g,
      n_rows = nrow(z),
      identical_expression = identical_expression
    )
  )
}


cat(
  "All duplicated IDs have identical expression:",
  all(dup_check$identical_expression),
  "\n"
)


# STOP rather than silently deduplicate
# if the expression values differ.

stopifnot(
  all(
    dup_check$identical_expression
  )
)


# ============================================================
# 5. DEDUPLICATE EXACT ENSEMBL IDs
#
# Preference:
# 1. non-PAR annotation
# 2. otherwise first occurrence
#
# Expression is identical within every exact duplicate group,
# so this changes annotation redundancy, not expression.
# ============================================================

dat$annotation_priority <- ifelse(
  dat$id == "PAR",
  1,
  0
)


dat <- dat[
  order(
    dat$ensembl_id,
    dat$annotation_priority
  ),
]


dat <- dat[
  !duplicated(
    dat$ensembl_id
  ),
]


rownames(dat) <- NULL


cat("\nAfter exact-ID deduplication:\n")

cat(
  "Rows:",
  nrow(dat),
  "\n"
)

cat(
  "Unique Ensembl IDs:",
  length(
    unique(
      dat$ensembl_id
    )
  ),
  "\n"
)


stopifnot(
  nrow(dat) ==
    length(
      unique(
        dat$ensembl_id
      )
    )
)


# ============================================================
# 6. BUILD TPM MATRIX
# ============================================================

tpm <- as.matrix(
  dat[
    ,
    tpm_cols,
    drop = FALSE
  ]
)


storage.mode(tpm) <- "numeric"

rownames(tpm) <- dat$ensembl_id

colnames(tpm) <- sample_names


# ============================================================
# 7. PRIMARY FILTER
#
# TPM >= 1 in at least 2 of 4 samples
#
# Reason:
# n = 2 independent mice per biological condition.
# ============================================================

keep_primary <- rowSums(
  tpm >= 1
) >= 2


cat("\n========================================\n")
cat("PRIMARY TPM FILTER\n")
cat("========================================\n\n")

cat(
  "Before filtering:",
  nrow(tpm),
  "\n"
)

cat(
  "Retained:",
  sum(keep_primary),
  "\n"
)

cat(
  "Removed:",
  sum(!keep_primary),
  "\n"
)


tpm_f <- tpm[
  keep_primary,
  ,
  drop = FALSE
]


annotation_f <- dat[
  keep_primary,
  c(
    "ensembl_id",
    "ensembl_base",
    "id"
  )
]


# ============================================================
# 8. LOG TRANSFORMATION
# ============================================================

log_tpm <- log2(
  tpm_f + 1
)


# ============================================================
# 9. PRIMARY PCA
# ============================================================

pca <- prcomp(
  t(log_tpm),
  center = TRUE,
  scale. = FALSE
)


pca_var <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )


pca_table <- data.frame(
  
  sample = sample_names,
  
  group = group,
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("PRIMARY PCA\n")
cat("========================================\n\n")

cat(
  "PC1:",
  round(
    pca_var[1],
    2
  ),
  "%\n"
)

cat(
  "PC2:",
  round(
    pca_var[2],
    2
  ),
  "%\n\n"
)


print(
  pca_table,
  row.names = FALSE
)


# ============================================================
# 10. PRIMARY LIMMA-TREND MODEL
# ============================================================

design <- model.matrix(
  ~ 0 + group
)


colnames(design) <- c(
  "Control",
  "OlapR"
)


fit <- lmFit(
  log_tpm,
  design
)


contrast_matrix <- makeContrasts(
  OlapR_vs_Control =
    OlapR - Control,
  levels = design
)


fit2 <- contrasts.fit(
  fit,
  contrast_matrix
)


fit2 <- eBayes(
  fit2,
  trend = TRUE
)


primary <- topTable(
  fit2,
  coef = "OlapR_vs_Control",
  number = Inf,
  sort.by = "none"
)


primary$ensembl_id <-
  annotation_f$ensembl_id

primary$ensembl_base <-
  annotation_f$ensembl_base

primary$symbol <-
  annotation_f$id


primary <- primary[
  ,
  c(
    "ensembl_id",
    "ensembl_base",
    "symbol",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


# ============================================================
# 11. PRIMARY DEG CLASSIFICATION
# ============================================================

primary$DEG <- "NS"


primary$DEG[
  primary$adj.P.Val < 0.05 &
    primary$logFC >= 1
] <- "UP"


primary$DEG[
  primary$adj.P.Val < 0.05 &
    primary$logFC <= -1
] <- "DOWN"


cat("\n========================================\n")
cat("PRIMARY DIFFERENTIAL EXPRESSION\n")
cat("========================================\n\n")


cat(
  "BH FDR < 0.05:",
  sum(
    primary$adj.P.Val < 0.05
  ),
  "\n"
)


cat(
  "BH FDR < 0.05 and |log2FC| >= 1:",
  sum(
    primary$DEG != "NS"
  ),
  "\n"
)


cat(
  "UP:",
  sum(
    primary$DEG == "UP"
  ),
  "\n"
)


cat(
  "DOWN:",
  sum(
    primary$DEG == "DOWN"
  ),
  "\n"
)


cat("\nTop primary DE results:\n")

print(
  head(
    primary[
      order(
        primary$adj.P.Val
      ),
    ],
    20
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================
# 12. SENSITIVITY ANALYSIS:
# RSEM EXPECTED COUNTS + TMM + VOOM
#
# IMPORTANT:
# This is NOT the primary analysis.
#
# The workbook does not provide the RSEM effective-length
# matrix needed for the preferred tximport count+offset method.
#
# This therefore tests robustness to an alternative
# expression scale / statistical model.
# ============================================================
# ============================================================


counts <- as.matrix(
  dat[
    ,
    count_cols,
    drop = FALSE
  ]
)


storage.mode(counts) <- "numeric"

rownames(counts) <- dat$ensembl_id

colnames(counts) <- sample_names


# ------------------------------------------------------------
# 12.1 DGEList
# ------------------------------------------------------------

dge <- DGEList(
  counts = counts,
  group = group
)


# n = 2 per group

keep_count <- filterByExpr(
  dge,
  group = group,
  min.count = 10
)


dge <- dge[
  keep_count,
  ,
  keep.lib.sizes = FALSE
]


dge <- calcNormFactors(
  dge,
  method = "TMM"
)


cat("\n========================================\n")
cat("EXPECTED-COUNT SENSITIVITY FILTER\n")
cat("========================================\n\n")

cat(
  "Retained:",
  nrow(dge),
  "\n"
)


# ------------------------------------------------------------
# 12.2 VOOM
# ------------------------------------------------------------

v <- voom(
  dge,
  design,
  plot = FALSE
)


fit_count <- lmFit(
  v,
  design
)


fit_count <- contrasts.fit(
  fit_count,
  contrast_matrix
)


fit_count <- eBayes(
  fit_count
)


sens <- topTable(
  fit_count,
  coef = "OlapR_vs_Control",
  number = Inf,
  sort.by = "none"
)


sens$ensembl_id <- rownames(
  sens
)


sens <- sens[
  ,
  c(
    "ensembl_id",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


colnames(sens)[
  colnames(sens) != "ensembl_id"
] <- paste0(
  "sensitivity_",
  colnames(sens)[
    colnames(sens) != "ensembl_id"
  ]
)


# ============================================================
# 13. PRIMARY-vs-SENSITIVITY COMPARISON
# ============================================================

comparison <- merge(
  primary,
  sens,
  by = "ensembl_id"
)


logfc_rho <- cor(
  comparison$logFC,
  comparison$sensitivity_logFC,
  method = "spearman"
)


t_rho <- cor(
  comparison$t,
  comparison$sensitivity_t,
  method = "spearman"
)


direction_same <- mean(
  sign(comparison$logFC) ==
    sign(
      comparison$sensitivity_logFC
    )
)


primary_sig <- comparison$adj.P.Val < 0.05

sens_sig <-
  comparison$sensitivity_adj.P.Val < 0.05


cat("\n========================================\n")
cat("PRIMARY vs SENSITIVITY CONCORDANCE\n")
cat("========================================\n\n")


cat(
  "Genes compared:",
  nrow(comparison),
  "\n"
)


cat(
  "Spearman logFC rho:",
  round(
    logfc_rho,
    4
  ),
  "\n"
)


cat(
  "Spearman t-statistic rho:",
  round(
    t_rho,
    4
  ),
  "\n"
)


cat(
  "Same logFC direction:",
  round(
    100 * direction_same,
    2
  ),
  "%\n"
)


cat(
  "Primary FDR < 0.05:",
  sum(primary_sig),
  "\n"
)


cat(
  "Sensitivity FDR < 0.05:",
  sum(sens_sig),
  "\n"
)


cat(
  "Significant in both:",
  sum(
    primary_sig &
      sens_sig
  ),
  "\n"
)


# ============================================================
# 14. ROBUSTNESS OF PRIMARY LARGE-EFFECT DEGs
# ============================================================

primary_large <- comparison[
  comparison$adj.P.Val < 0.05 &
    abs(
      comparison$logFC
    ) >= 1,
]


cat("\n========================================\n")
cat("PRIMARY LARGE-EFFECT DEG ROBUSTNESS\n")
cat("========================================\n\n")


cat(
  "Primary large-effect DEGs:",
  nrow(primary_large),
  "\n"
)


if (
  nrow(primary_large) > 0
) {
  
  cat(
    "Same direction in sensitivity:",
    round(
      100 *
        mean(
          sign(
            primary_large$logFC
          ) ==
            sign(
              primary_large$sensitivity_logFC
            )
        ),
      2
    ),
    "%\n"
  )
  
  
  cat(
    "Sensitivity FDR < 0.05:",
    sum(
      primary_large$sensitivity_adj.P.Val <
        0.05
    ),
    "/",
    nrow(primary_large),
    "\n"
  )
}


# ============================================================
# 15. DUPLICATED SYMBOL INFORMATION FOR FUTURE GSEA
# ============================================================

symbol_counts <- table(
  primary$symbol
)


duplicate_symbols <- names(
  symbol_counts[
    symbol_counts > 1
  ]
)


cat("\n========================================\n")
cat("GSEA SYMBOL AUDIT\n")
cat("========================================\n\n")


cat(
  "Unique symbols in primary filtered set:",
  length(
    unique(
      primary$symbol
    )
  ),
  "\n"
)


cat(
  "Duplicated symbols requiring rank collapse:",
  length(
    duplicate_symbols
  ),
  "\n"
)


# ============================================================
# 16. SAVE RESULTS
# ============================================================

outdir <- "results/GSE191231"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  primary,
  file.path(
    outdir,
    "GSE191231_primary_limma_TPM_DE.csv"
  ),
  row.names = FALSE
)


write.csv(
  sens,
  file.path(
    outdir,
    "GSE191231_sensitivity_RSEM_voom_DE.csv"
  ),
  row.names = FALSE
)


write.csv(
  comparison,
  file.path(
    outdir,
    "GSE191231_primary_vs_sensitivity_DE.csv"
  ),
  row.names = FALSE
)


write.csv(
  pca_table,
  file.path(
    outdir,
    "GSE191231_primary_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


write.csv(
  dup_check,
  file.path(
    outdir,
    "GSE191231_exact_duplicate_validation.csv"
  ),
  row.names = FALSE
)


# Save filtered TPM matrix

write.csv(
  data.frame(
    ensembl_id =
      rownames(tpm_f),
    tpm_f,
    check.names = FALSE
  ),
  file.path(
    outdir,
    "GSE191231_filtered_TPM_matrix.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("GSE191231 PRIMARY DE COMPLETE\n")
cat("========================================\n")