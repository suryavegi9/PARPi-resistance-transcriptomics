# ============================================================
# GSE163854 - D8
# PRIMARY DIFFERENTIAL EXPRESSION + HALLMARK GSEA
# + INDEPENDENT edgeR-QL SENSITIVITY ANALYSIS
#
# FROZEN PRIMARY COMPARISON:
#
# Sensitive:
#   P2A1, P2A2, P2A3
#
# Resistant:
#   F451, F452, F453
#
# Frozen expression filter:
#   CPM >= 1 in >= 3 of 6 primary samples
#
# PRIMARY:
#   TMM -> voom -> limma
#   Hallmark GSEA ranked by moderated t
#
# SENSITIVITY:
#   edgeR quasi-likelihood
#   Hallmark GSEA ranked by sign(logFC) * sqrt(QLF)
#
# DEG reporting threshold:
#   BH FDR < 0.05
#
# Large-effect DEG threshold:
#   BH FDR < 0.05 AND |log2FC| >= 1
#
# NO post-hoc filtering or tuning.
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

qc_dir <-
  "results/GSE163854/QC"

analysis_dir <-
  "results/GSE163854/D8_analysis"

dir.create(
  analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("GSE163854 - D8 DE + HALLMARK GSEA\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

cran_packages <- c(
  "msigdbr"
)

bioc_packages <- c(
  "edgeR",
  "limma",
  "clusterProfiler"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    install.packages(pkg)
    
  }
}


if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager"
  )
}


for (pkg in bioc_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
}


library(edgeR)
library(limma)
library(clusterProfiler)
library(msigdbr)


cat(
  "edgeR version:",
  as.character(
    packageVersion(
      "edgeR"
    )
  ),
  "\n"
)

cat(
  "limma version:",
  as.character(
    packageVersion(
      "limma"
    )
  ),
  "\n"
)

cat(
  "clusterProfiler version:",
  as.character(
    packageVersion(
      "clusterProfiler"
    )
  ),
  "\n"
)

cat(
  "msigdbr version:",
  as.character(
    packageVersion(
      "msigdbr"
    )
  ),
  "\n\n"
)


# ============================================================
# 2. LOAD FROZEN QC OBJECT
# ============================================================

qc_file <- file.path(
  qc_dir,
  "D8_QC_objects.rds"
)


if (!file.exists(
  qc_file
)) {
  
  stop(
    "D8 QC object from Script 57 not found."
  )
}


qc <- readRDS(
  qc_file
)


# ============================================================
# 3. VERIFY FROZEN PRIMARY DESIGN
# ============================================================

expected_samples <- c(
  "P2A1",
  "P2A2",
  "P2A3",
  "F451",
  "F452",
  "F453"
)


if (!identical(
  colnames(
    qc$primary_dge
  ),
  expected_samples
)) {
  
  cat(
    "Observed primary sample order:\n"
  )
  
  print(
    colnames(
      qc$primary_dge
    )
  )
  
  
  stop(
    "Primary samples differ from frozen D8 design."
  )
}


if (
  nrow(
    qc$primary_dge
  ) != 13410
) {
  
  stop(
    paste0(
      "Expected 13,410 frozen genes, observed ",
      nrow(
        qc$primary_dge
      ),
      "."
    )
  )
}


cat(
  "Frozen primary samples verified.\n"
)

cat(
  "Frozen genes verified:",
  nrow(
    qc$primary_dge
  ),
  "\n"
)


# ============================================================
# 4. COPY FROZEN DGE OBJECT
# ============================================================

dge <- qc$primary_dge


# Recalculate TMM using the current edgeR function name.
# This should reproduce the frozen QC normalization factors.

dge <- edgeR::normLibSizes(
  dge,
  method = "TMM"
)


cat(
  "\nTMM normalization factors:\n"
)


print(
  dge$samples[
    ,
    c(
      "lib.size",
      "norm.factors"
    )
  ]
)


# ============================================================
# 5. DESIGN MATRIX
#
# Sensitive is reference conceptually.
#
# Contrast:
#   Resistant - Sensitive
#
# Therefore:
#   positive logFC = higher in resistant
#   negative logFC = lower in resistant
# ============================================================

group <- factor(
  c(
    "Sensitive",
    "Sensitive",
    "Sensitive",
    "Resistant",
    "Resistant",
    "Resistant"
  ),
  levels = c(
    "Sensitive",
    "Resistant"
  )
)


design <- model.matrix(
  ~ 0 + group
)


colnames(
  design
) <- c(
  "Sensitive",
  "Resistant"
)


rownames(
  design
) <- colnames(
  dge
)


cat(
  "\nDesign matrix:\n\n"
)


print(
  design
)


contrast_matrix <- limma::makeContrasts(
  Resistant_minus_Sensitive =
    Resistant - Sensitive,
  levels = design
)


cat(
  "\nContrast matrix:\n\n"
)


print(
  contrast_matrix
)


# ============================================================
# 6. PRIMARY VOOM
# ============================================================

pdf(
  file.path(
    analysis_dir,
    "D8_voom_mean_variance.pdf"
  ),
  width = 7,
  height = 6
)


v <- limma::voom(
  dge,
  design,
  plot = TRUE
)


dev.off()


# ============================================================
# 7. PRIMARY LIMMA MODEL
# ============================================================

fit <- limma::lmFit(
  v,
  design
)


fit <- limma::contrasts.fit(
  fit,
  contrast_matrix
)


fit <- limma::eBayes(
  fit
)


de_limma <- limma::topTable(
  fit,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)


de_limma$gene <- rownames(
  de_limma
)


de_limma <- de_limma[
  ,
  c(
    "gene",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


rownames(
  de_limma
) <- NULL


# ============================================================
# 8. DEG CLASSIFICATION
# ============================================================

de_limma$FDR_significant <-
  de_limma$adj.P.Val < 0.05


de_limma$large_DEG <-
  de_limma$adj.P.Val < 0.05 &
  abs(
    de_limma$logFC
  ) >= 1


de_limma$direction <- "NS"


de_limma$direction[
  de_limma$large_DEG &
    de_limma$logFC > 0
] <- "UP"


de_limma$direction[
  de_limma$large_DEG &
    de_limma$logFC < 0
] <- "DOWN"


write.csv(
  de_limma,
  file.path(
    analysis_dir,
    "D8_primary_limma_DE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. PRIMARY DE SUMMARY
# ============================================================

n_fdr <- sum(
  de_limma$FDR_significant
)


n_large <- sum(
  de_limma$large_DEG
)


n_up <- sum(
  de_limma$direction == "UP"
)


n_down <- sum(
  de_limma$direction == "DOWN"
)


cat("\n")
cat("============================================================\n")
cat("PRIMARY LIMMA DIFFERENTIAL EXPRESSION\n")
cat("============================================================\n\n")


cat(
  "Genes tested:",
  nrow(
    de_limma
  ),
  "\n"
)


cat(
  "BH FDR < 0.05:",
  n_fdr,
  "\n"
)


cat(
  "BH FDR < 0.05 and |logFC| >= 1:",
  n_large,
  "\n"
)


cat(
  "Large-effect UP:",
  n_up,
  "\n"
)


cat(
  "Large-effect DOWN:",
  n_down,
  "\n"
)


# ============================================================
# 10. TOP PRIMARY DE GENES
# ============================================================

top_de <- de_limma[
  order(
    de_limma$adj.P.Val,
    -abs(
      de_limma$logFC
    )
  ),
  ,
  drop = FALSE
]


cat(
  "\nTop 20 genes by FDR:\n\n"
)


print(
  head(
    top_de,
    20
  ),
  row.names = FALSE
)


# ============================================================
# 11. PRIMARY GSEA RANKING
#
# Moderated t-statistic.
# ============================================================

rank_limma <- de_limma$t


names(
  rank_limma
) <- de_limma$gene


rank_limma <- rank_limma[
  is.finite(
    rank_limma
  )
]


rank_limma <- sort(
  rank_limma,
  decreasing = TRUE
)


limma_tie_fraction <-
  1 -
  length(
    unique(
      rank_limma
    )
  ) /
  length(
    rank_limma
  )


cat(
  "\nPrimary ranking genes:",
  length(
    rank_limma
  ),
  "\n"
)


cat(
  "Primary exact-tie fraction:",
  round(
    limma_tie_fraction,
    6
  ),
  "\n"
)


# ============================================================
# 12. HALLMARK GENE SETS
# ============================================================

hallmark_msig <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)


required_msig_columns <- c(
  "gs_name",
  "gene_symbol"
)


if (!all(
  required_msig_columns %in%
  colnames(
    hallmark_msig
  )
)) {
  
  stop(
    "Expected msigdbr columns gs_name/gene_symbol are missing."
  )
}


term2gene <- unique(
  hallmark_msig[
    ,
    c(
      "gs_name",
      "gene_symbol"
    )
  ]
)


term2gene <- term2gene[
  !is.na(
    term2gene$gene_symbol
  ) &
    term2gene$gene_symbol != "",
  ,
  drop = FALSE
]


cat(
  "\nHallmark gene sets:",
  length(
    unique(
      term2gene$gs_name
    )
  ),
  "\n"
)


if (
  "db_version" %in%
  colnames(
    hallmark_msig
  )
) {
  
  cat(
    "MSigDB version:",
    paste(
      unique(
        hallmark_msig$db_version
      ),
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# 13. GSEA FUNCTION
# ============================================================

run_hallmark_gsea <- function(
    ranked_vector,
    label
) {
  
  
  set.seed(
    20261110
  )
  
  
  gsea_object <- clusterProfiler::GSEA(
    
    geneList =
      ranked_vector,
    
    TERM2GENE =
      term2gene,
    
    exponent =
      1,
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    eps =
      0,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    verbose =
      FALSE,
    
    by =
      "fgsea"
    
  )
  
  
  result <- as.data.frame(
    gsea_object@result
  )
  
  
  if (
    nrow(
      result
    ) == 0
  ) {
    
    stop(
      paste(
        "No Hallmark GSEA results returned for",
        label
      )
    )
  }
  
  
  result$analysis <-
    label
  
  
  return(
    result
  )
}


# ============================================================
# 14. PRIMARY HALLMARK GSEA
# ============================================================

gsea_limma <- run_hallmark_gsea(
  ranked_vector = rank_limma,
  label = "limma_voom"
)


write.csv(
  gsea_limma,
  file.path(
    analysis_dir,
    "D8_primary_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)


gsea_limma_sorted <- gsea_limma[
  order(
    gsea_limma$p.adjust,
    -abs(
      gsea_limma$NES
    )
  ),
  ,
  drop = FALSE
]


cat("\n")
cat("============================================================\n")
cat("PRIMARY HALLMARK GSEA\n")
cat("============================================================\n\n")


cat(
  "Hallmarks returned:",
  nrow(
    gsea_limma
  ),
  "\n"
)


cat(
  "Hallmarks FDR < 0.05:",
  sum(
    gsea_limma$p.adjust < 0.05
  ),
  "\n\n"
)


cat(
  "Significant Hallmarks:\n\n"
)


print(
  gsea_limma_sorted[
    gsea_limma_sorted$p.adjust < 0.05,
    c(
      "ID",
      "NES",
      "pvalue",
      "p.adjust",
      "setSize"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. edgeR QUASI-LIKELIHOOD SENSITIVITY
# ============================================================

dge_ql <- dge


dge_ql <- edgeR::estimateDisp(
  dge_ql,
  design
)


fit_ql <- edgeR::glmQLFit(
  dge_ql,
  design,
  robust = TRUE
)


contrast_vector <- c(
  -1,
  1
)


ql_test <- edgeR::glmQLFTest(
  fit_ql,
  contrast = contrast_vector
)


de_edgeR <- edgeR::topTags(
  ql_test,
  n = Inf,
  sort.by = "none"
)$table


de_edgeR$gene <- rownames(
  de_edgeR
)


de_edgeR <- de_edgeR[
  ,
  c(
    "gene",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR"
  )
]


rownames(
  de_edgeR
) <- NULL


de_edgeR$QL_signed_stat <-
  sign(
    de_edgeR$logFC
  ) *
  sqrt(
    pmax(
      de_edgeR$F,
      0
    )
  )


write.csv(
  de_edgeR,
  file.path(
    analysis_dir,
    "D8_edgeR_QL_sensitivity_DE.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. LIMMA VS edgeR GENE-LEVEL VALIDATION
# ============================================================

de_compare <- merge(
  
  de_limma[
    ,
    c(
      "gene",
      "logFC",
      "t",
      "adj.P.Val",
      "large_DEG"
    )
  ],
  
  de_edgeR[
    ,
    c(
      "gene",
      "logFC",
      "QL_signed_stat",
      "FDR"
    )
  ],
  
  by = "gene",
  
  suffixes = c(
    "_limma",
    "_edgeR"
  )
  
)


logfc_rho <- cor(
  de_compare$logFC_limma,
  de_compare$logFC_edgeR,
  method = "spearman"
)


nonzero_compare <- de_compare[
  de_compare$logFC_limma != 0 &
    de_compare$logFC_edgeR != 0,
  ,
  drop = FALSE
]


genomewide_sign_agreement <- mean(
  sign(
    nonzero_compare$logFC_limma
  ) ==
    sign(
      nonzero_compare$logFC_edgeR
    )
)


large_compare <- de_compare[
  de_compare$large_DEG,
  ,
  drop = FALSE
]


large_sign_agreement <- mean(
  sign(
    large_compare$logFC_limma
  ) ==
    sign(
      large_compare$logFC_edgeR
    )
)


large_logfc_rho <- cor(
  large_compare$logFC_limma,
  large_compare$logFC_edgeR,
  method = "spearman"
)


both_sig <- de_compare[
  de_compare$adj.P.Val < 0.05 &
    de_compare$FDR < 0.05,
  ,
  drop = FALSE
]


both_sig_sign_agreement <- mean(
  sign(
    both_sig$logFC_limma
  ) ==
    sign(
      both_sig$logFC_edgeR
    )
)


write.csv(
  de_compare,
  file.path(
    analysis_dir,
    "D8_limma_vs_edgeR_gene_comparison.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("edgeR-QL GENE-LEVEL SENSITIVITY\n")
cat("============================================================\n\n")


cat(
  "Matched genes:",
  nrow(
    de_compare
  ),
  "\n"
)


cat(
  "Genome-wide logFC Spearman rho:",
  round(
    logfc_rho,
    5
  ),
  "\n"
)


cat(
  "Genome-wide direction agreement:",
  round(
    100 *
      genomewide_sign_agreement,
    2
  ),
  "%\n"
)


cat(
  "Primary large DEGs:",
  nrow(
    large_compare
  ),
  "\n"
)


cat(
  "Large-DEG logFC Spearman rho:",
  round(
    large_logfc_rho,
    5
  ),
  "\n"
)


cat(
  "Large-DEG direction agreement:",
  round(
    100 *
      large_sign_agreement,
    2
  ),
  "%\n"
)


cat(
  "Genes FDR < 0.05 in BOTH:",
  nrow(
    both_sig
  ),
  "\n"
)


cat(
  "Both-significant direction agreement:",
  round(
    100 *
      both_sig_sign_agreement,
    2
  ),
  "%\n"
)


# ============================================================
# 17. edgeR HALLMARK GSEA
# ============================================================

rank_edgeR <-
  de_edgeR$QL_signed_stat


names(
  rank_edgeR
) <- de_edgeR$gene


rank_edgeR <- rank_edgeR[
  is.finite(
    rank_edgeR
  )
]


rank_edgeR <- sort(
  rank_edgeR,
  decreasing = TRUE
)


edgeR_tie_fraction <-
  1 -
  length(
    unique(
      rank_edgeR
    )
  ) /
  length(
    rank_edgeR
  )


cat(
  "\nedgeR ranking genes:",
  length(
    rank_edgeR
  ),
  "\n"
)


cat(
  "edgeR exact-tie fraction:",
  round(
    edgeR_tie_fraction,
    6
  ),
  "\n"
)


gsea_edgeR <- run_hallmark_gsea(
  ranked_vector = rank_edgeR,
  label = "edgeR_QL"
)


write.csv(
  gsea_edgeR,
  file.path(
    analysis_dir,
    "D8_edgeR_QL_Hallmark_GSEA.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. HALLMARK SENSITIVITY COMPARISON
# ============================================================

gsea_compare <- merge(
  
  gsea_limma[
    ,
    c(
      "ID",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ],
  
  gsea_edgeR[
    ,
    c(
      "ID",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ],
  
  by = "ID",
  
  suffixes = c(
    "_limma",
    "_edgeR"
  ),
  
  all = TRUE
)


complete_gsea <- gsea_compare[
  is.finite(
    gsea_compare$NES_limma
  ) &
    is.finite(
      gsea_compare$NES_edgeR
    ),
  ,
  drop = FALSE
]


hallmark_rho <- cor(
  complete_gsea$NES_limma,
  complete_gsea$NES_edgeR,
  method = "spearman"
)


hallmark_same_direction <- mean(
  sign(
    complete_gsea$NES_limma
  ) ==
    sign(
      complete_gsea$NES_edgeR
    )
)


gsea_compare$same_direction <-
  sign(
    gsea_compare$NES_limma
  ) ==
  sign(
    gsea_compare$NES_edgeR
  )


gsea_compare$significant_both <-
  gsea_compare$p.adjust_limma < 0.05 &
  gsea_compare$p.adjust_edgeR < 0.05


gsea_compare$robust_significant <-
  gsea_compare$same_direction &
  gsea_compare$significant_both


write.csv(
  gsea_compare,
  file.path(
    analysis_dir,
    "D8_Hallmark_limma_vs_edgeR.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("HALLMARK SENSITIVITY\n")
cat("============================================================\n\n")


cat(
  "Hallmarks compared:",
  nrow(
    complete_gsea
  ),
  "\n"
)


cat(
  "NES Spearman rho:",
  round(
    hallmark_rho,
    5
  ),
  "\n"
)


cat(
  "Same-direction Hallmarks:",
  sum(
    complete_gsea$NES_limma *
      complete_gsea$NES_edgeR > 0
  ),
  "/",
  nrow(
    complete_gsea
  ),
  "\n"
)


cat(
  "Hallmarks significant in both with same direction:",
  sum(
    gsea_compare$robust_significant,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 19. PRINT ROBUST HALLMARKS
# ============================================================

robust_hallmarks <- gsea_compare[
  gsea_compare$robust_significant,
  ,
  drop = FALSE
]


robust_hallmarks <- robust_hallmarks[
  order(
    robust_hallmarks$p.adjust_limma,
    -abs(
      robust_hallmarks$NES_limma
    )
  ),
  ,
  drop = FALSE
]


cat(
  "\nRobust Hallmarks:\n\n"
)


print(
  robust_hallmarks[
    ,
    c(
      "ID",
      "NES_limma",
      "p.adjust_limma",
      "NES_edgeR",
      "p.adjust_edgeR"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 20. D1 FROZEN VALIDATION HYPOTHESES
# ============================================================

d1_frozen <- data.frame(
  
  pathway = c(
    
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
  
  expected_sign = c(
    1,
    1,
    -1,
    -1,
    -1,
    -1,
    -1,
    -1,
    -1,
    -1
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 21. MERGE D8 GSEA WITH D1 FROZEN HYPOTHESES
# ============================================================

d1_validation <- merge(
  
  d1_frozen,
  
  gsea_compare,
  
  by.x = "pathway",
  
  by.y = "ID",
  
  all.x = TRUE,
  
  sort = FALSE
)


# Restore original hypothesis order.

d1_validation$order_index <- match(
  d1_validation$pathway,
  d1_frozen$pathway
)


d1_validation <- d1_validation[
  order(
    d1_validation$order_index
  ),
  ,
  drop = FALSE
]


d1_validation$primary_direction <-
  ifelse(
    d1_validation$NES_limma > 0,
    "UP",
    "DOWN"
  )


d1_validation$primary_direction_match <-
  sign(
    d1_validation$NES_limma
  ) ==
  d1_validation$expected_sign


d1_validation$primary_significant_match <-
  d1_validation$primary_direction_match &
  d1_validation$p.adjust_limma < 0.05


d1_validation$sensitivity_direction_match <-
  sign(
    d1_validation$NES_edgeR
  ) ==
  d1_validation$expected_sign


d1_validation$robust_significant_match <-
  d1_validation$primary_direction_match &
  d1_validation$sensitivity_direction_match &
  d1_validation$p.adjust_limma < 0.05 &
  d1_validation$p.adjust_edgeR < 0.05


write.csv(
  d1_validation,
  file.path(
    analysis_dir,
    "D8_vs_D1_frozen_validation.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("D1 FROZEN HYPOTHESIS VALIDATION IN D8\n")
cat("============================================================\n\n")


print(
  d1_validation[
    ,
    c(
      "pathway",
      "expected_direction",
      "NES_limma",
      "p.adjust_limma",
      "NES_edgeR",
      "p.adjust_edgeR",
      "primary_direction_match",
      "robust_significant_match"
    )
  ],
  row.names = FALSE
)


cat(
  "\nD1 directional matches:",
  sum(
    d1_validation$primary_direction_match,
    na.rm = TRUE
  ),
  "/10\n"
)


cat(
  "D1 primary significant matches:",
  sum(
    d1_validation$primary_significant_match,
    na.rm = TRUE
  ),
  "/10\n"
)


cat(
  "D1 robust significant matches:",
  sum(
    d1_validation$robust_significant_match,
    na.rm = TRUE
  ),
  "/10\n"
)


# ============================================================
# 22. ACCESSION-LEVEL HALLMARK TABLE FOR FUTURE INTEGRATION
#
# D8 contains one frozen primary resistance contrast,
# so primary limma GSEA is the accession-level vote.
# ============================================================

D8_hallmark_accession <- gsea_limma[
  ,
  c(
    "ID",
    "NES",
    "pvalue",
    "p.adjust",
    "setSize"
  )
]


colnames(
  D8_hallmark_accession
) <- c(
  "pathway",
  "D8_NES",
  "D8_pvalue",
  "D8_FDR",
  "setSize"
)


D8_hallmark_accession$direction <-
  ifelse(
    D8_hallmark_accession$D8_NES > 0,
    "UP",
    "DOWN"
  )


D8_hallmark_accession$significant <-
  D8_hallmark_accession$D8_FDR < 0.05


write.csv(
  D8_hallmark_accession,
  file.path(
    analysis_dir,
    "D8_accession_Hallmark_vote.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. FINAL ANALYSIS SUMMARY
# ============================================================

analysis_summary <- data.frame(
  
  metric = c(
    
    "Genes_tested",
    
    "Limma_FDR_lt_0.05",
    
    "Limma_FDR_lt_0.05_abs_logFC_ge_1",
    
    "Large_DEG_UP",
    
    "Large_DEG_DOWN",
    
    "Primary_Hallmark_FDR_lt_0.05",
    
    "Gene_logFC_limma_edgeR_rho",
    
    "Gene_genomewide_direction_agreement",
    
    "Large_DEG_direction_agreement",
    
    "Hallmark_NES_limma_edgeR_rho",
    
    "Hallmark_same_direction_fraction",
    
    "Robust_Hallmark_count",
    
    "D1_directional_matches",
    
    "D1_primary_significant_matches",
    
    "D1_robust_significant_matches"
    
  ),
  
  value = c(
    
    nrow(
      de_limma
    ),
    
    n_fdr,
    
    n_large,
    
    n_up,
    
    n_down,
    
    sum(
      gsea_limma$p.adjust < 0.05
    ),
    
    logfc_rho,
    
    genomewide_sign_agreement,
    
    large_sign_agreement,
    
    hallmark_rho,
    
    hallmark_same_direction,
    
    sum(
      gsea_compare$robust_significant,
      na.rm = TRUE
    ),
    
    sum(
      d1_validation$primary_direction_match,
      na.rm = TRUE
    ),
    
    sum(
      d1_validation$primary_significant_match,
      na.rm = TRUE
    ),
    
    sum(
      d1_validation$robust_significant_match,
      na.rm = TRUE
    )
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  analysis_summary,
  file.path(
    analysis_dir,
    "D8_analysis_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("FINAL D8 ANALYSIS SUMMARY\n")
cat("============================================================\n\n")


print(
  analysis_summary,
  row.names = FALSE
)


# ============================================================
# 24. SESSION INFO
# ============================================================

sink(
  file.path(
    analysis_dir,
    "D8_analysis_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 25. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D8 PRIMARY DE + GSEA COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Frozen primary comparison:\n",
  "  Sensitive: P2A1, P2A2, P2A3\n",
  "  Resistant: F451, F452, F453\n\n"
)


cat(
  "Primary pipeline:\n",
  "  TMM -> voom -> limma -> moderated-t Hallmark GSEA\n\n"
)


cat(
  "Sensitivity pipeline:\n",
  "  edgeR quasi-likelihood -> signed sqrt(QLF) Hallmark GSEA\n\n"
)


cat(
  "No samples, genes, pathways, or thresholds were changed ",
  "after viewing the results.\n"
)


cat(
  "\n============================================================\n"
)