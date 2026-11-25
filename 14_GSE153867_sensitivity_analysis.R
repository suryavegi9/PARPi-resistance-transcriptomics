# ============================================================
# GSE153867 - SENSITIVITY ANALYSIS
#
# Pre-specified before viewing Dataset-3 results:
#
# Primary filter:
#   FPKM >= 1 in >= 4/16 samples
#
# Sensitivity filter:
#   FPKM >= 1 in >= 8/16 samples
#
# Everything else remains unchanged.
# ============================================================


# ============================================================
# 0. SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)
library(clusterProfiler)
library(msigdbr)

dir.create(
  "results/GSE153867/sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. READ FPKM MATRIX
# ============================================================

raw <- read.delim(
  "data_raw/GSE153867_fpkm.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 2. SELECT ONLY C vs O
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

expression <- raw[
  ,
  c(
    "GeneID",
    selected_samples
  )
]


# ============================================================
# 3. COLLAPSE DUPLICATE SYMBOLS
# ============================================================

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


# ============================================================
# 4. STRICTER PRE-SPECIFIED FILTER
# ============================================================
#
# Keep genes with FPKM >= 1
# in at least 8 of 16 samples.
# ============================================================

keep <- rowSums(
  fpkm >= 1
) >= 8

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
  "Genes after stricter filter:",
  nrow(fpkm_filtered),
  "\n"
)


# ============================================================
# 5. LOG2 TRANSFORMATION
# ============================================================

log_expr <- log2(
  fpkm_filtered + 1
)


# ============================================================
# 6. METADATA
# ============================================================

metadata <- data.frame(
  
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
  
  row.names = colnames(log_expr)
)


# ============================================================
# 7. LIMMA
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


cat("\n========================================\n")
cat("SENSITIVITY LIMMA RESULTS\n")
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
# 8. RANK BY MODERATED t
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


# ============================================================
# 9. HALLMARK DATABASE
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


# ============================================================
# 10. HALLMARK GSEA
# ============================================================

set.seed(123)

gsea_sensitivity <- GSEA(
  
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


sens_results <- as.data.frame(
  gsea_sensitivity
)


write.csv(
  sens_results,
  "results/GSE153867/sensitivity/GSE153867_Hallmark_GSEA_strict_filter.csv",
  row.names = FALSE
)


# ============================================================
# 11. LOAD PRIMARY GSEA
# ============================================================

primary <- read.csv(
  "results/GSE153867/GSE153867_Hallmark_GSEA_all.csv"
)


# ============================================================
# 12. COMPARE PRIMARY VS SENSITIVITY
# ============================================================

comparison <- merge(
  
  primary[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  
  sens_results[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  
  by = "ID",
  
  suffixes = c(
    "_primary",
    "_strict"
  )
)


comparison$direction_same <-
  sign(
    comparison$NES_primary
  ) ==
  sign(
    comparison$NES_strict
  )


write.csv(
  comparison,
  "results/GSE153867/sensitivity/GSE153867_primary_vs_strict_filter.csv",
  row.names = FALSE
)


# ============================================================
# 13. GLOBAL ROBUSTNESS
# ============================================================

rho <- cor(
  comparison$NES_primary,
  comparison$NES_strict,
  method = "spearman",
  use = "complete.obs"
)


same_direction <- sum(
  comparison$direction_same,
  na.rm = TRUE
)


cat("\n========================================\n")
cat("GLOBAL HALLMARK ROBUSTNESS\n")
cat("========================================\n")


cat(
  "\nSpearman correlation of Hallmark NES:",
  round(
    rho,
    3
  ),
  "\n"
)


cat(
  "Hallmarks retaining NES direction:",
  same_direction,
  "/",
  nrow(comparison),
  "\n"
)


# ============================================================
# 14. SIGNIFICANT PATHWAYS IN SENSITIVITY ANALYSIS
# ============================================================

sens_sig <- subset(
  sens_results,
  p.adjust < 0.05
)


cat("\n========================================\n")
cat("SIGNIFICANT HALLMARKS - STRICT FILTER\n")
cat("========================================\n")


cat(
  "\nNumber significant:",
  nrow(sens_sig),
  "\n\n"
)


if (nrow(sens_sig) > 0) {
  
  print(
    sens_sig[
      order(
        sens_sig$p.adjust
      ),
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    row.names = FALSE
  )
}


# ============================================================
# 15. FROZEN PATHWAYS
# ============================================================

frozen <- c(
  
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
)


frozen_comparison <- comparison[
  comparison$ID %in% frozen,
]


cat("\n========================================\n")
cat("FROZEN PATHWAY ROBUSTNESS\n")
cat("========================================\n\n")


print(
  frozen_comparison[
    ,
    c(
      "ID",
      "NES_primary",
      "NES_strict",
      "p.adjust_primary",
      "p.adjust_strict",
      "direction_same"
    )
  ],
  row.names = FALSE
)


cat(
  "\nFrozen pathways retaining direction:",
  sum(
    frozen_comparison$direction_same,
    na.rm = TRUE
  ),
  "/",
  nrow(frozen_comparison),
  "\n"
)


# ============================================================
# 16. IDENTIFY ANY DIRECTION FLIPS
# ============================================================

flipped <- comparison[
  comparison$direction_same == FALSE,
]


cat("\n========================================\n")
cat("DIRECTIONALLY UNSTABLE HALLMARKS\n")
cat("========================================\n")


if (nrow(flipped) == 0) {
  
  cat(
    "\nNo Hallmark pathways changed NES direction.\n"
  )
  
} else {
  
  print(
    flipped[
      ,
      c(
        "ID",
        "NES_primary",
        "NES_strict",
        "p.adjust_primary",
        "p.adjust_strict"
      )
    ],
    row.names = FALSE
  )
}


# ============================================================
# COMPLETE
# ============================================================

cat("\n========================================\n")
cat("GSE153867 SENSITIVITY ANALYSIS COMPLETE\n")
cat("========================================\n")