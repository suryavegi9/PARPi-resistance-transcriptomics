# ============================================================
# GSE202434 - PSN1
# REPLICATE QC + MAINTAIN-2 SENSITIVITY
#
# Goal:
# Determine whether PSN1-maintain-2 is an influential outlier
# and whether primary DE conclusions depend on that sample.
#
# IMPORTANT:
# We do NOT automatically exclude maintain-2.
# Exclusion analysis is sensitivity analysis only.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)


# ============================================================
# 1. LOAD ORIGINAL MATRIX
# ============================================================

file_path <- paste0(
  "data_raw/",
  "GSE202434_Exp43_Tala_Normalized_Counts_ExpressedGenesONLY_ProteinCoding_ONLY.csv.gz"
)

dat <- read.csv(
  gzfile(file_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(dat)[1] <- "gene"

expr <- as.matrix(
  dat[, -1]
)

rownames(expr) <- dat$gene

storage.mode(expr) <- "numeric"


# ============================================================
# 2. SAME FILTER AS PRIMARY ANALYSIS
# ============================================================

expr_back <- 2^expr - 1

keep <- rowSums(
  expr_back >= 1
) >= 3

expr_f <- expr[
  keep,
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("FILTER QC\n")
cat("========================================\n\n")

cat(
  "Genes retained:",
  nrow(expr_f),
  "\n"
)


# ============================================================
# 3. SAMPLE METADATA
# ============================================================

sample_names <- colnames(expr_f)

group <- ifelse(
  
  grepl(
    "parental",
    sample_names,
    ignore.case = TRUE
  ),
  
  "Parental",
  
  ifelse(
    
    grepl(
      "maintain",
      sample_names,
      ignore.case = TRUE
    ),
    
    "TalaR_M",
    
    "Withdraw"
  )
)


metadata <- data.frame(
  sample = sample_names,
  group = group,
  stringsAsFactors = FALSE
)


# ============================================================
# 4. SAMPLE-SAMPLE SPEARMAN CORRELATION
# ============================================================

sample_cor <- cor(
  expr_f,
  method = "spearman"
)


cat("\n========================================\n")
cat("SAMPLE SPEARMAN CORRELATIONS\n")
cat("========================================\n\n")

print(
  round(
    sample_cor,
    4
  )
)


cat("\n========================================\n")
cat("WITHIN-GROUP CORRELATIONS\n")
cat("========================================\n\n")


for (g in unique(metadata$group)) {
  
  samples_g <- metadata$sample[
    metadata$group == g
  ]
  
  sub_cor <- sample_cor[
    samples_g,
    samples_g,
    drop = FALSE
  ]
  
  vals <- sub_cor[
    upper.tri(sub_cor)
  ]
  
  cat(
    g,
    ": mean rho =",
    round(mean(vals), 4),
    "; range =",
    paste(
      round(range(vals), 4),
      collapse = " to "
    ),
    "\n"
  )
}


# ============================================================
# 5. LIMMA SAMPLE QUALITY WEIGHTS
#
# arrayWeights estimates how reliable each sample is relative
# to the model.
#
# Weight < 1 = downweighted
# Weight > 1 = relatively higher consistency
#
# This is QC evidence, not an automatic exclusion rule.
# ============================================================

group_factor <- factor(
  metadata$group,
  levels = c(
    "Parental",
    "Withdraw",
    "TalaR_M"
  )
)

design <- model.matrix(
  ~ 0 + group_factor
)

colnames(design) <- c(
  "Parental",
  "Withdraw",
  "TalaR_M"
)


weights <- arrayWeights(
  expr_f,
  design = design
)


weight_table <- data.frame(
  sample = colnames(expr_f),
  group = metadata$group,
  array_weight = weights,
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("LIMMA SAMPLE QUALITY WEIGHTS\n")
cat("========================================\n\n")

print(
  weight_table[
    order(weight_table$array_weight),
  ],
  row.names = FALSE
)


# ============================================================
# 6. LEAVE-MAINTAIN-2-OUT PCA
# ============================================================

exclude_sample <- "PSN1-maintain-2"

keep_samples <- colnames(expr_f) !=
  exclude_sample


expr_noM2 <- expr_f[
  ,
  keep_samples,
  drop = FALSE
]


pca_noM2 <- prcomp(
  t(expr_noM2),
  center = TRUE,
  scale. = FALSE
)


variance_noM2 <- 100 *
  pca_noM2$sdev^2 /
  sum(pca_noM2$sdev^2)


metadata_noM2 <- metadata[
  metadata$sample != exclude_sample,
]


pca_noM2_data <- data.frame(
  
  sample = rownames(pca_noM2$x),
  
  group = metadata_noM2$group[
    match(
      rownames(pca_noM2$x),
      metadata_noM2$sample
    )
  ],
  
  PC1 = pca_noM2$x[, 1],
  
  PC2 = pca_noM2$x[, 2],
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("PCA WITHOUT PSN1-maintain-2\n")
cat("========================================\n\n")

cat(
  "PC1 variance:",
  round(variance_noM2[1], 2),
  "%\n"
)

cat(
  "PC2 variance:",
  round(variance_noM2[2], 2),
  "%\n\n"
)

print(
  pca_noM2_data,
  row.names = FALSE
)


# ============================================================
# 7. PRIMARY P3 RESULTS
# ============================================================

primary_P3 <- read.csv(
  "results/GSE202434/PSN1/PSN1_P3_M_vs_WD_all.csv",
  stringsAsFactors = FALSE
)


# ============================================================
# 8. P3 SENSITIVITY WITHOUT MAINTAIN-2
# ============================================================

group_noM2 <- factor(
  
  metadata_noM2$group,
  
  levels = c(
    "Parental",
    "Withdraw",
    "TalaR_M"
  )
)


design_noM2 <- model.matrix(
  ~ 0 + group_noM2
)

colnames(design_noM2) <- c(
  "Parental",
  "Withdraw",
  "TalaR_M"
)


fit_noM2 <- lmFit(
  expr_noM2,
  design_noM2
)


contrast_noM2 <- makeContrasts(
  
  P3_M_vs_WD =
    TalaR_M - Withdraw,
  
  levels = design_noM2
)


fit_noM2 <- contrasts.fit(
  fit_noM2,
  contrast_noM2
)


fit_noM2 <- eBayes(
  fit_noM2,
  trend = TRUE
)


P3_noM2 <- topTable(
  fit_noM2,
  coef = "P3_M_vs_WD",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)


P3_noM2$gene <- rownames(
  P3_noM2
)


# ============================================================
# 9. ALIGN PRIMARY + SENSITIVITY RESULTS
# ============================================================

comparison <- merge(
  
  primary_P3[
    ,
    c(
      "gene",
      "logFC",
      "t",
      "adj.P.Val"
    )
  ],
  
  P3_noM2[
    ,
    c(
      "gene",
      "logFC",
      "t",
      "adj.P.Val"
    )
  ],
  
  by = "gene",
  
  suffixes = c(
    "_primary",
    "_noM2"
  )
)


# ============================================================
# 10. ROBUSTNESS METRICS
# ============================================================

logFC_rho <- cor(
  
  comparison$logFC_primary,
  
  comparison$logFC_noM2,
  
  method = "spearman"
)


t_rho <- cor(
  
  comparison$t_primary,
  
  comparison$t_noM2,
  
  method = "spearman"
)


same_direction <- mean(
  
  sign(comparison$logFC_primary) ==
    sign(comparison$logFC_noM2)
)


primary_DEG <-
  comparison$adj.P.Val_primary < 0.05 &
  abs(comparison$logFC_primary) >= 1


noM2_DEG <-
  comparison$adj.P.Val_noM2 < 0.05 &
  abs(comparison$logFC_noM2) >= 1


cat("\n========================================\n")
cat("P3 SENSITIVITY: PRIMARY VS NO MAINTAIN-2\n")
cat("========================================\n\n")


cat(
  "Genes compared:",
  nrow(comparison),
  "\n"
)

cat(
  "Spearman rho - logFC:",
  round(logFC_rho, 4),
  "\n"
)

cat(
  "Spearman rho - moderated t:",
  round(t_rho, 4),
  "\n"
)

cat(
  "Same logFC direction:",
  round(
    100 * same_direction,
    2
  ),
  "%\n"
)

cat(
  "Primary large-effect DEGs:",
  sum(primary_DEG),
  "\n"
)

cat(
  "No-maintain-2 large-effect DEGs:",
  sum(noM2_DEG),
  "\n"
)


# ============================================================
# 11. PRIMARY DEGs: DIRECTION ROBUSTNESS
# ============================================================

primary_deg_table <- comparison[
  primary_DEG,
]


cat("\n========================================\n")
cat("PRIMARY P3 DEG ROBUSTNESS\n")
cat("========================================\n\n")


cat(
  "Primary DEGs:",
  nrow(primary_deg_table),
  "\n"
)

cat(
  "Primary DEGs same direction without maintain-2:",
  round(
    100 *
      mean(
        sign(primary_deg_table$logFC_primary) ==
          sign(primary_deg_table$logFC_noM2)
      ),
    2
  ),
  "%\n"
)


cat(
  "Primary DEGs still FDR<0.05 without maintain-2:",
  round(
    100 *
      mean(
        primary_deg_table$adj.P.Val_noM2 < 0.05
      ),
    2
  ),
  "%\n"
)


# ============================================================
# 12. SAVE RESULTS
# ============================================================

outdir <- "results/GSE202434/PSN1/QC"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  sample_cor,
  file.path(
    outdir,
    "PSN1_sample_Spearman_correlations.csv"
  )
)


write.csv(
  weight_table,
  file.path(
    outdir,
    "PSN1_arrayWeights.csv"
  ),
  row.names = FALSE
)


write.csv(
  pca_noM2_data,
  file.path(
    outdir,
    "PSN1_PCA_without_maintain2.csv"
  ),
  row.names = FALSE
)


write.csv(
  comparison,
  file.path(
    outdir,
    "PSN1_P3_primary_vs_without_maintain2.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("PSN1 REPLICATE QC COMPLETE\n")
cat("========================================\n")