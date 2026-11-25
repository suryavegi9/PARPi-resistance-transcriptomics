# ============================================================
# GSE117765 - PRIMARY VALIDATION ANALYSIS
#
# Dataset 2:
# PEO1 parental (n = 2)
# PEO1 independently derived olaparib-resistant clones (n = 4)
#
# PRE-SPECIFIED PIPELINE:
# TPM
# -> NA treated as 0
# -> duplicate gene symbols collapsed
# -> TPM >= 1 in >= 2 samples
# -> log2(TPM + 1)
# -> limma
# -> moderated t-statistic
# -> Hallmark GSEA
#
# IMPORTANT:
# Dataset 1 hypotheses were frozen BEFORE this analysis.
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

required_packages <- c(
  "limma",
  "ggplot2",
  "clusterProfiler",
  "msigdbr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste(
        "Package", pkg,
        "is not installed."
      )
    )
  }
}

library(limma)
library(ggplot2)
library(clusterProfiler)
library(msigdbr)

dir.create(
  "results/GSE117765",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/GSE117765",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. READ TPM MATRIX
# ============================================================

candidate_files <- c(
  "data_raw/GSE117765/GSE117765_matrix.txt.gz",
  "data_raw/GSE117765/GSE117765_matrix.txt"
)

existing_files <- candidate_files[file.exists(candidate_files)]

if (length(existing_files) == 0) {
  stop("GSE117765 matrix file not found.")
}

file_path <- existing_files[1]

if (grepl("\\.gz$", file_path)) {
  
  raw <- read.delim(
    gzfile(file_path),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  raw <- read.delim(
    file_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


cat("\nRaw dimensions:")
cat(nrow(raw), "genes x", ncol(raw) - 1, "samples\n")


# ============================================================
# 2. NA -> ZERO
# ============================================================
#
# Rationale:
# - TPM matrix contains no explicit zero values.
# - Missing values therefore most plausibly represent
#   non-detected / zero-expression entries after merging
#   per-sample expression tables.
#
# This assumption will later be tested by sensitivity analysis.
# ============================================================

expression <- raw

n_missing_before <- sum(
  is.na(expression[, -1])
)

expression[, -1][
  is.na(expression[, -1])
] <- 0

cat(
  "\nMissing TPM values converted to zero:",
  n_missing_before,
  "\n"
)


# ============================================================
# 3. COLLAPSE DUPLICATE GENE SYMBOLS
# ============================================================

cat(
  "Duplicated gene symbols before collapse:",
  sum(duplicated(expression$gene)),
  "\n"
)

expression_collapsed <- aggregate(
  . ~ gene,
  data = expression,
  FUN = sum
)

rownames(expression_collapsed) <-
  expression_collapsed$gene

tpm <- as.matrix(
  expression_collapsed[, -1]
)

storage.mode(tpm) <- "numeric"

cat(
  "Unique genes after collapse:",
  nrow(tpm),
  "\n"
)


# ============================================================
# 4. METADATA
# ============================================================

metadata <- data.frame(
  
  sample = colnames(tpm),
  
  phenotype = factor(
    c(
      "Sensitive",
      "Sensitive",
      "Resistant",
      "Resistant",
      "Resistant",
      "Resistant"
    ),
    levels = c(
      "Sensitive",
      "Resistant"
    )
  ),
  
  clone = c(
    "Parent_1",
    "Parent_2",
    "Clone10",
    "Clone18",
    "Clone4",
    "Clone5"
  ),
  
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample

stopifnot(
  identical(
    rownames(metadata),
    colnames(tpm)
  )
)

write.csv(
  metadata,
  "results/GSE117765/GSE117765_metadata.csv",
  row.names = FALSE
)


# ============================================================
# 5. PRE-SPECIFIED LOW-EXPRESSION FILTER
# ============================================================
#
# Keep genes with TPM >= 1
# in at least 2 samples.
#
# This rule is fixed before viewing validation results.
# ============================================================

keep <- rowSums(
  tpm >= 1
) >= 2

tpm_filtered <- tpm[
  keep,
  ,
  drop = FALSE
]

cat(
  "\nGenes before filtering:",
  nrow(tpm),
  "\n"
)

cat(
  "Genes after TPM filter:",
  nrow(tpm_filtered),
  "\n"
)


write.csv(
  tpm_filtered,
  "results/GSE117765/GSE117765_filtered_TPM.csv"
)


# ============================================================
# 6. LOG2 TRANSFORMATION
# ============================================================

log_expr <- log2(
  tpm_filtered + 1
)

write.csv(
  log_expr,
  "results/GSE117765/GSE117765_log2TPMplus1.csv"
)


# ============================================================
# 7. PCA - QC BEFORE INTERPRETATION
# ============================================================

pca <- prcomp(
  t(log_expr),
  center = TRUE,
  scale. = FALSE
)

variance <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pca_df <- data.frame(
  
  sample = rownames(pca$x),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  phenotype = metadata[
    rownames(pca$x),
    "phenotype"
  ],
  
  clone = metadata[
    rownames(pca$x),
    "clone"
  ]
)


p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = phenotype
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text(
    aes(label = clone),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  labs(
    title = "GSE117765 transcriptomic structure",
    subtitle = "PEO1 parental cells and independently derived olaparib-resistant clones",
    x = paste0(
      "PC1 (",
      round(variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance[2], 1),
      "%)"
    ),
    color = "Phenotype"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_pca)

ggsave(
  "figures/GSE117765/GSE117765_PCA.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 8. LIMMA DIFFERENTIAL EXPRESSION
# ============================================================

design <- model.matrix(
  ~ 0 + phenotype,
  data = metadata
)

colnames(design) <- c(
  "Sensitive",
  "Resistant"
)

fit <- lmFit(
  log_expr,
  design
)

contrast_matrix <- makeContrasts(
  Resistant_vs_Sensitive =
    Resistant - Sensitive,
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


limma_results <- topTable(
  fit2,
  coef = "Resistant_vs_Sensitive",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

limma_results$gene <- rownames(
  limma_results
)

limma_results <- limma_results[
  ,
  c(
    "gene",
    setdiff(
      colnames(limma_results),
      "gene"
    )
  )
]

write.csv(
  limma_results,
  "results/GSE117765/GSE117765_limma_all_results.csv",
  row.names = FALSE
)


cat("\n========================================\n")
cat("LIMMA ANALYSIS\n")
cat("========================================\n")

cat(
  "\nGenes tested:",
  nrow(limma_results),
  "\n"
)

cat(
  "BH-adjusted P < 0.05:",
  sum(
    limma_results$adj.P.Val < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "BH-adjusted P < 0.05 AND |logFC| >= 1:",
  sum(
    limma_results$adj.P.Val < 0.05 &
      abs(limma_results$logFC) >= 1,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 9. RANK GENES USING MODERATED t-STATISTIC
# ============================================================

ranked_genes <- fit2$t[, 1]

names(ranked_genes) <- rownames(
  fit2$t
)

ranked_genes <- ranked_genes[
  is.finite(ranked_genes)
]

ranked_genes <- sort(
  ranked_genes,
  decreasing = TRUE
)


cat(
  "\nGenes entering GSEA:",
  length(ranked_genes),
  "\n"
)


# ============================================================
# 10. OBTAIN MSigDB HALLMARK SETS
# ============================================================

hallmark_db <- tryCatch(
  
  msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ),
  
  error = function(e) {
    
    msigdbr(
      species = "Homo sapiens",
      category = "H"
    )
  }
)

hallmark_t2g <- unique(
  hallmark_db[
    ,
    c(
      "gs_name",
      "gene_symbol"
    )
  ]
)

colnames(hallmark_t2g) <- c(
  "term",
  "gene"
)

cat(
  "Hallmark gene sets:",
  length(
    unique(hallmark_t2g$term)
  ),
  "\n"
)


# ============================================================
# 11. HALLMARK GSEA
# ============================================================

set.seed(123)

gsea_hallmark <- GSEA(
  
  geneList = ranked_genes,
  
  TERM2GENE = hallmark_t2g,
  
  minGSSize = 10,
  
  maxGSSize = 500,
  
  pvalueCutoff = 1,
  
  pAdjustMethod = "BH",
  
  eps = 0,
  
  seed = TRUE,
  
  verbose = FALSE
)


gsea_results <- as.data.frame(
  gsea_hallmark
)

write.csv(
  gsea_results,
  "results/GSE117765/GSE117765_Hallmark_GSEA_all.csv",
  row.names = FALSE
)


gsea_sig <- subset(
  gsea_results,
  p.adjust < 0.05
)

write.csv(
  gsea_sig,
  "results/GSE117765/GSE117765_Hallmark_GSEA_significant.csv",
  row.names = FALSE
)


cat(
  "\nSignificant Hallmark pathways:",
  nrow(gsea_sig),
  "\n"
)


# ============================================================
# 12. TEST FROZEN DATASET-1 PREDICTIONS
# ============================================================

validation <- data.frame(
  
  Hallmark = c(
    
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    
    "HALLMARK_G2M_CHECKPOINT",
    
    "HALLMARK_E2F_TARGETS",
    
    "HALLMARK_MTORC1_SIGNALING",
    
    "HALLMARK_MYC_TARGETS_V1",
    
    "HALLMARK_GLYCOLYSIS",
    
    "HALLMARK_INFLAMMATORY_RESPONSE",
    
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    
    "HALLMARK_UV_RESPONSE_DN"
  ),
  
  expected_direction = c(
    
    "UP",
    "UP",
    
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN"
  ),
  
  stringsAsFactors = FALSE
)


validation$NES <-
  gsea_results$NES[
    match(
      validation$Hallmark,
      gsea_results$ID
    )
  ]


validation$p_adjust <-
  gsea_results$p.adjust[
    match(
      validation$Hallmark,
      gsea_results$ID
    )
  ]


validation$observed_direction <- ifelse(
  
  validation$NES > 0,
  
  "UP",
  
  ifelse(
    validation$NES < 0,
    "DOWN",
    NA
  )
)


validation$direction_concordant <-
  validation$observed_direction ==
  validation$expected_direction


validation$significant <-
  validation$p_adjust < 0.05


validation$significant_validation <-
  validation$direction_concordant &
  validation$significant


write.csv(
  validation,
  "results/GSE117765/GSE117765_frozen_hypothesis_validation.csv",
  row.names = FALSE
)


# ============================================================
# 13. PRINT VALIDATION TABLE
# ============================================================

cat("\n========================================\n")
cat("FROZEN HYPOTHESIS VALIDATION\n")
cat("========================================\n\n")

print(
  validation,
  row.names = FALSE
)


cat("\nDirectional concordance:")
cat(
  sum(
    validation$direction_concordant,
    na.rm = TRUE
  ),
  "/",
  nrow(validation),
  "\n"
)


cat(
  "Directionally concordant AND significant:",
  sum(
    validation$significant_validation,
    na.rm = TRUE
  ),
  "/",
  nrow(validation),
  "\n"
)


# ============================================================
# 14. SAVE SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/GSE117765/sessionInfo_GSE117765.txt"
)


# ============================================================
# DONE
# ============================================================

cat("\n========================================\n")
cat("GSE117765 PRIMARY ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nOutputs saved under:\n",
  "results/GSE117765/\n"
)

cat(
  "\nPCA saved under:\n",
  "figures/GSE117765/\n"
)