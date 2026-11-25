# ============================================================
# GSE153867 - DATASET 3 PRIMARY ANALYSIS
#
# PRIMARY COMPARISON:
# A2780 olaparib-resistant C-1...C-8
# versus
# A2780 parental O-1...O-8
#
# CS and R groups belong to the separate C13* C/EBP-beta
# knockdown experiment and are EXCLUDED.
#
# PRE-SPECIFIED PIPELINE:
# FPKM
# -> collapse duplicate symbols
# -> FPKM >= 1 in >= 4 of 16 samples
# -> log2(FPKM + 1)
# -> limma-trend
# -> moderated t statistic
# -> Hallmark GSEA
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(msigdbr)

dir.create(
  "results/GSE153867",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/GSE153867",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. WRITE PRE-SPECIFIED ANALYSIS PLAN
# ============================================================

analysis_plan <- c(
  
  "GSE153867 PRE-SPECIFIED PRIMARY ANALYSIS",
  
  "",
  
  "Primary biological comparison:",
  "A2780 olaparib-resistant C-1 through C-8",
  "versus",
  "A2780 parental O-1 through O-8",
  
  "",
  
  "Excluded samples:",
  "CS samples = C13* shCON",
  "R samples = C13* shCEBPB",
  
  "",
  
  "Primary pipeline:",
  "1. Collapse duplicate gene symbols by summing FPKM.",
  "2. Retain genes with FPKM >= 1 in at least 4/16 samples.",
  "3. Transform using log2(FPKM + 1).",
  "4. Differential expression using limma with trend=TRUE.",
  "5. Rank genes by moderated t statistic.",
  "6. Hallmark GSEA.",
  "7. Hallmark significance threshold: BH-adjusted P < 0.05.",
  
  "",
  
  "Frozen Dataset-1 predictions remain unchanged.",
  
  "",
  
  "Planned sensitivity analysis:",
  "Repeat using stricter FPKM >= 1 in at least 8/16 samples."
)

writeLines(
  analysis_plan,
  "results/GSE153867/GSE153867_prespecified_analysis_plan.txt"
)


# ============================================================
# 2. READ FPKM MATRIX
# ============================================================

file_path <- "data_raw/GSE153867_fpkm.txt"

if (!file.exists(file_path)) {
  stop("GSE153867_fpkm.txt not found.")
}

raw <- read.delim(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 3. SELECT ONLY PRIMARY RESISTANCE COMPARISON
# ============================================================

parental_samples <- paste0(
  "O-",
  1:8
)

resistant_samples <- paste0(
  "C-",
  1:8
)

selected_samples <- c(
  parental_samples,
  resistant_samples
)

stopifnot(
  all(
    selected_samples %in%
      colnames(raw)
  )
)


expression <- raw[
  ,
  c(
    "GeneID",
    selected_samples
  )
]


cat("\n========================================\n")
cat("PRIMARY COMPARISON\n")
cat("========================================\n")

cat("\nParental samples:\n")
print(parental_samples)

cat("\nResistant samples:\n")
print(resistant_samples)

cat(
  "\nSamples included:",
  length(selected_samples),
  "\n"
)

cat(
  "Samples excluded from full matrix:",
  ncol(raw) - 1 - length(selected_samples),
  "\n"
)


# ============================================================
# 4. COLLAPSE DUPLICATE GENE SYMBOLS
# ============================================================

cat(
  "\nDuplicated gene symbols before collapse:",
  sum(
    duplicated(
      expression$GeneID
    )
  ),
  "\n"
)


expression_collapsed <- aggregate(
  . ~ GeneID,
  data = expression,
  FUN = sum
)

rownames(expression_collapsed) <-
  expression_collapsed$GeneID

fpkm <- as.matrix(
  expression_collapsed[, -1]
)

storage.mode(fpkm) <- "numeric"


cat(
  "Unique genes after collapse:",
  nrow(fpkm),
  "\n"
)


# ============================================================
# 5. PRIMARY EXPRESSION FILTER
# ============================================================
#
# At least FPKM >= 1 in >= 25% of samples.
#
# 25% of 16 = 4 samples.
# ============================================================

keep <- rowSums(
  fpkm >= 1
) >= 4

fpkm_filtered <- fpkm[
  keep,
  ,
  drop = FALSE
]


cat(
  "\nGenes before filtering:",
  nrow(fpkm),
  "\n"
)

cat(
  "Genes after FPKM filter:",
  nrow(fpkm_filtered),
  "\n"
)


write.csv(
  fpkm_filtered,
  "results/GSE153867/GSE153867_primary_filtered_FPKM.csv"
)


# ============================================================
# 6. LOG2 TRANSFORMATION
# ============================================================

log_expr <- log2(
  fpkm_filtered + 1
)

write.csv(
  log_expr,
  "results/GSE153867/GSE153867_log2FPKMplus1.csv"
)


# ============================================================
# 7. METADATA
# ============================================================

metadata <- data.frame(
  
  sample = colnames(log_expr),
  
  phenotype = factor(
    ifelse(
      grepl(
        "^C-",
        colnames(log_expr)
      ),
      "Resistant",
      "Parental"
    ),
    levels = c(
      "Parental",
      "Resistant"
    )
  ),
  
  stringsAsFactors = FALSE
)

rownames(metadata) <-
  metadata$sample


stopifnot(
  identical(
    rownames(metadata),
    colnames(log_expr)
  )
)


write.csv(
  metadata,
  "results/GSE153867/GSE153867_primary_metadata.csv",
  row.names = FALSE
)


# ============================================================
# 8. PCA
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
  
  sample = rownames(
    pca$x
  ),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  phenotype = metadata[
    rownames(pca$x),
    "phenotype"
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
  geom_text_repel(
    aes(
      label = sample
    ),
    size = 3.3,
    show.legend = FALSE
  ) +
  labs(
    title = "GSE153867 transcriptomic structure",
    subtitle =
      "A2780 parental and olaparib-resistant samples",
    x = paste0(
      "PC1 (",
      round(
        variance[1],
        1
      ),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance[2],
        1
      ),
      "% variance)"
    ),
    color = "Phenotype"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )


print(
  p_pca
)


ggsave(
  "figures/GSE153867/GSE153867_primary_PCA.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 9. LIMMA DIFFERENTIAL EXPRESSION
# ============================================================

design <- model.matrix(
  ~ 0 + phenotype,
  data = metadata
)

colnames(design) <- c(
  "Parental",
  "Resistant"
)


fit <- lmFit(
  log_expr,
  design
)


contrast_matrix <- makeContrasts(
  
  Resistant_vs_Parental =
    Resistant - Parental,
  
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
  coef = "Resistant_vs_Parental",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)


limma_results$gene <-
  rownames(
    limma_results
  )


limma_results <-
  limma_results[
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
  "results/GSE153867/GSE153867_limma_all_results.csv",
  row.names = FALSE
)


cat("\n========================================\n")
cat("LIMMA RESULTS\n")
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
      abs(
        limma_results$logFC
      ) >= 1,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 10. MODERATED t RANKING
# ============================================================

ranked_genes <-
  fit2$t[, 1]

names(ranked_genes) <-
  rownames(
    fit2$t
  )

ranked_genes <-
  ranked_genes[
    is.finite(
      ranked_genes
    )
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
# 11. MSigDB HALLMARK SETS
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


colnames(
  hallmark_t2g
) <- c(
  "term",
  "gene"
)


cat(
  "Hallmark gene sets:",
  length(
    unique(
      hallmark_t2g$term
    )
  ),
  "\n"
)


# ============================================================
# 12. HALLMARK GSEA
# ============================================================

set.seed(123)


gsea <- GSEA(
  
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


gsea_results <-
  as.data.frame(
    gsea
  )


write.csv(
  gsea_results,
  "results/GSE153867/GSE153867_Hallmark_GSEA_all.csv",
  row.names = FALSE
)


gsea_sig <- subset(
  gsea_results,
  p.adjust < 0.05
)


write.csv(
  gsea_sig,
  "results/GSE153867/GSE153867_Hallmark_GSEA_significant.csv",
  row.names = FALSE
)


cat(
  "\nSignificant Hallmark pathways:",
  nrow(gsea_sig),
  "\n"
)


# ============================================================
# 13. TEST FROZEN DATASET-1 HYPOTHESIS
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


validation$observed_direction <-
  ifelse(
    validation$NES > 0,
    "UP",
    "DOWN"
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
  "results/GSE153867/GSE153867_frozen_hypothesis_validation.csv",
  row.names = FALSE
)


# ============================================================
# 14. PRINT FROZEN VALIDATION
# ============================================================

cat("\n========================================\n")
cat("FROZEN HYPOTHESIS VALIDATION\n")
cat("========================================\n\n")


print(
  validation,
  row.names = FALSE
)


cat(
  "\nDirectional concordance:",
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
# 15. PRINT ALL SIGNIFICANT HALLMARKS
# ============================================================

cat("\n========================================\n")
cat("SIGNIFICANT HALLMARK PATHWAYS\n")
cat("========================================\n\n")


if (nrow(gsea_sig) > 0) {
  
  print(
    gsea_sig[
      order(
        gsea_sig$p.adjust
      ),
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No Hallmark pathways reached BH-adjusted P < 0.05.\n"
  )
}


# ============================================================
# 16. SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  "results/GSE153867/sessionInfo_GSE153867.txt"
)


# ============================================================
# DONE
# ============================================================

cat("\n========================================\n")
cat("GSE153867 PRIMARY ANALYSIS COMPLETE\n")
cat("========================================\n")