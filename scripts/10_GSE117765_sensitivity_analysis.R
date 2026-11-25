# ============================================================
# GSE117765 - SENSITIVITY ANALYSIS
#
# Question:
# Do the main Hallmark conclusions depend on treating
# missing TPM entries as zero?
#
# Sensitivity strategy:
# retain only genes with complete TPM observations
# across all 6 samples.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)
library(clusterProfiler)
library(msigdbr)

dir.create(
  "results/GSE117765/sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. READ RAW TPM MATRIX
# ============================================================

candidate_files <- c(
  "data_raw/GSE117765/GSE117765_matrix.txt.gz",
  "data_raw/GSE117765/GSE117765_matrix.txt"
)

existing_files <- candidate_files[file.exists(candidate_files)]

if (length(existing_files) == 0) {
  stop("GSE117765 matrix not found.")
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


# ============================================================
# 2. COMPLETE-CASE GENES ONLY
# ============================================================

complete_keep <- complete.cases(
  raw[, -1]
)

complete_raw <- raw[
  complete_keep,
]

cat(
  "\nGenes with complete TPM observations:",
  nrow(complete_raw),
  "\n"
)


# ============================================================
# 3. COLLAPSE DUPLICATE SYMBOLS
# ============================================================

complete_collapsed <- aggregate(
  . ~ gene,
  data = complete_raw,
  FUN = sum
)

rownames(complete_collapsed) <-
  complete_collapsed$gene

tpm <- as.matrix(
  complete_collapsed[, -1]
)

storage.mode(tpm) <- "numeric"


# ============================================================
# 4. SAME EXPRESSION FILTER
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
  "Genes after complete-case + TPM filter:",
  nrow(tpm_filtered),
  "\n"
)


# ============================================================
# 5. LOG TRANSFORMATION
# ============================================================

log_expr <- log2(
  tpm_filtered + 1
)


# ============================================================
# 6. METADATA
# ============================================================

metadata <- data.frame(
  
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


# ============================================================
# 8. RANK BY MODERATED t
# ============================================================

ranked_genes <- fit2$t[, 1]

names(ranked_genes) <-
  rownames(fit2$t)

ranked_genes <- ranked_genes[
  is.finite(ranked_genes)
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

colnames(hallmark_t2g) <- c(
  "term",
  "gene"
)


# ============================================================
# 10. GSEA
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
  "results/GSE117765/sensitivity/Hallmark_complete_case_GSEA.csv",
  row.names = FALSE
)


# ============================================================
# 11. COMPARE PRIMARY VS SENSITIVITY
# ============================================================

primary <- read.csv(
  "results/GSE117765/GSE117765_Hallmark_GSEA_all.csv"
)

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
    "_complete_case"
  )
)


comparison$direction_same <-
  sign(comparison$NES_primary) ==
  sign(comparison$NES_complete_case)


write.csv(
  comparison,
  "results/GSE117765/sensitivity/Hallmark_primary_vs_complete_case.csv",
  row.names = FALSE
)


# ============================================================
# 12. PRINT MAIN RESULTS
# ============================================================

cat("\n========================================\n")
cat("COMPLETE-CASE SENSITIVITY ANALYSIS\n")
cat("========================================\n")

cat(
  "\nHallmarks significant in complete-case analysis:",
  sum(
    sens_results$p.adjust < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


significant_sensitivity <- subset(
  sens_results,
  p.adjust < 0.05
)

cat("\nSignificant pathways:\n")

print(
  significant_sensitivity[
    order(significant_sensitivity$p.adjust),
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  row.names = FALSE
)


cat("\n========================================\n")
cat("FROZEN PATHWAY COMPARISON\n")
cat("========================================\n")

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

print(
  frozen_comparison[
    ,
    c(
      "ID",
      "NES_primary",
      "NES_complete_case",
      "p.adjust_primary",
      "p.adjust_complete_case",
      "direction_same"
    )
  ],
  row.names = FALSE
)


cat("\n========================================\n")
cat("GLOBAL DIRECTIONAL STABILITY\n")
cat("========================================\n")

cat(
  "Hallmark pathways with same NES direction:",
  sum(
    comparison$direction_same,
    na.rm = TRUE
  ),
  "/",
  nrow(comparison),
  "\n"
)

cat("\nSENSITIVITY ANALYSIS COMPLETE\n")