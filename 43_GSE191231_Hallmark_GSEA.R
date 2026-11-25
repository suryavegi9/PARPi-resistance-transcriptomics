# ============================================================
# GSE191231 - DATASET 6
# CLEAN HALLMARK GSEA SCRIPT
#
# PRIMARY:
#   TPM / limma-trend DE results from script 42
#
# SENSITIVITY:
#   RSEM expected-count / voom DE results from script 42
#
# Rank statistic:
#   moderated t-statistic
#
# Positive NES:
#   higher in olaparib-selected resistant PDX
#
# Negative NES:
#   higher in vehicle-control PDX
#
# This script is completely self-contained.
# Run it from the TOP using source().
# ============================================================


# ============================================================
# 0. PROJECT DIRECTORY
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 1. PACKAGES
# ============================================================

library(readxl)
library(clusterProfiler)
library(msigdbr)


# ============================================================
# 2. FILES
# ============================================================

primary_file <- paste0(
  "results/GSE191231/",
  "GSE191231_primary_limma_TPM_DE.csv"
)

sensitivity_file <- paste0(
  "results/GSE191231/",
  "GSE191231_sensitivity_RSEM_voom_DE.csv"
)

workbook_file <- paste0(
  "data_raw/",
  "GSE191231_Hyun_Analysis_tpm_edger.xlsx"
)


stopifnot(
  file.exists(primary_file),
  file.exists(sensitivity_file),
  file.exists(workbook_file)
)


# ============================================================
# 3. LOAD DE RESULTS
# ============================================================

primary <- read.csv(
  primary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


sensitivity <- read.csv(
  sensitivity_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("INPUT DE TABLES\n")
cat("========================================\n\n")

cat(
  "Primary rows:",
  nrow(primary),
  "\n"
)

cat(
  "Sensitivity rows:",
  nrow(sensitivity),
  "\n"
)


stopifnot(
  "ensembl_id" %in% colnames(primary),
  "symbol" %in% colnames(primary),
  "t" %in% colnames(primary),
  "ensembl_id" %in% colnames(sensitivity),
  "sensitivity_t" %in% colnames(sensitivity)
)


# ============================================================
# 4. REBUILD CLEAN ENSEMBL -> SYMBOL MAPPING
# ============================================================

annotation <- read_excel(
  workbook_file,
  sheet = 1
)

annotation <- as.data.frame(
  annotation,
  check.names = FALSE
)


stopifnot(
  "gene" %in% colnames(annotation),
  "id" %in% colnames(annotation)
)


annotation$ensembl_id <- as.character(
  annotation$gene
)

annotation$symbol <- as.character(
  annotation$id
)


# Prefer informative annotation over "PAR"

annotation$priority <- ifelse(
  annotation$symbol == "PAR",
  1,
  0
)


annotation <- annotation[
  order(
    annotation$ensembl_id,
    annotation$priority
  ),
  ,
  drop = FALSE
]


annotation <- annotation[
  !duplicated(
    annotation$ensembl_id
  ),
  ,
  drop = FALSE
]


mapping <- annotation[
  ,
  c(
    "ensembl_id",
    "symbol"
  ),
  drop = FALSE
]


# ============================================================
# 5. ADD SYMBOLS TO SENSITIVITY TABLE
# ============================================================

sensitivity <- merge(
  sensitivity,
  mapping,
  by = "ensembl_id",
  all.x = TRUE,
  sort = FALSE
)


cat("\n========================================\n")
cat("SYMBOL MAPPING QC\n")
cat("========================================\n\n")

cat(
  "Sensitivity rows:",
  nrow(sensitivity),
  "\n"
)

cat(
  "Missing symbols:",
  sum(
    is.na(sensitivity$symbol) |
      sensitivity$symbol == ""
  ),
  "\n"
)


# ============================================================
# 6. FUNCTION TO BUILD GSEA RANK
#
# Duplicate gene symbols are collapsed by keeping the
# Ensembl feature with the largest absolute moderated t.
# ============================================================

make_rank <- function(
    dat,
    stat_col,
    label
) {
  
  stopifnot(
    "ensembl_id" %in% colnames(dat),
    "symbol" %in% colnames(dat),
    stat_col %in% colnames(dat)
  )
  
  
  # Explicit vectors
  symbol_vec <- as.character(
    dat$symbol
  )
  
  stat_vec <- as.numeric(
    dat[, stat_col]
  )
  
  
  keep <- !is.na(symbol_vec) &
    symbol_vec != "" &
    !is.na(stat_vec) &
    is.finite(stat_vec)
  
  
  tmp <- data.frame(
    
    ensembl_id = as.character(
      dat$ensembl_id[keep]
    ),
    
    symbol = symbol_vec[keep],
    
    statistic = stat_vec[keep],
    
    stringsAsFactors = FALSE
  )
  
  
  cat("\n----------------------------------------\n")
  cat(label, "\n")
  cat("----------------------------------------\n")
  
  cat(
    "Rows before symbol collapse:",
    nrow(tmp),
    "\n"
  )
  
  cat(
    "Unique symbols:",
    length(
      unique(
        tmp$symbol
      )
    ),
    "\n"
  )
  
  
  dup_symbols <- unique(
    tmp$symbol[
      duplicated(tmp$symbol) |
        duplicated(
          tmp$symbol,
          fromLast = TRUE
        )
    ]
  )
  
  
  cat(
    "Duplicated symbols:",
    length(
      dup_symbols
    ),
    "\n"
  )
  
  
  # Sort:
  # symbol alphabetically
  # strongest |t| first inside each symbol
  
  tmp <- tmp[
    order(
      tmp$symbol,
      -abs(tmp$statistic),
      tmp$ensembl_id
    ),
    ,
    drop = FALSE
  ]
  
  
  # Keep strongest feature for each symbol
  
  tmp <- tmp[
    !duplicated(
      tmp$symbol
    ),
    ,
    drop = FALSE
  ]
  
  
  rank_vector <- tmp$statistic
  
  names(rank_vector) <- tmp$symbol
  
  
  # Strong QC checks
  
  stopifnot(
    length(rank_vector) == nrow(tmp)
  )
  
  stopifnot(
    length(rank_vector) ==
      length(
        names(rank_vector)
      )
  )
  
  stopifnot(
    !any(
      duplicated(
        names(rank_vector)
      )
    )
  )
  
  stopifnot(
    all(
      is.finite(
        rank_vector
      )
    )
  )
  
  
  rank_vector <- sort(
    rank_vector,
    decreasing = TRUE
  )
  
  
  cat(
    "Final ranked genes:",
    length(rank_vector),
    "\n"
  )
  
  cat(
    "Maximum statistic:",
    round(
      max(rank_vector),
      4
    ),
    "\n"
  )
  
  cat(
    "Minimum statistic:",
    round(
      min(rank_vector),
      4
    ),
    "\n"
  )
  
  
  return(
    rank_vector
  )
}


# ============================================================
# 7. BUILD PRIMARY AND SENSITIVITY RANKS
# ============================================================

primary_rank <- make_rank(
  primary,
  "t",
  "PRIMARY TPM / LIMMA RANK"
)


sensitivity_rank <- make_rank(
  sensitivity,
  "sensitivity_t",
  "SENSITIVITY RSEM / VOOM RANK"
)


cat("\n========================================\n")
cat("RANK VECTOR QC\n")
cat("========================================\n\n")

cat(
  "Primary ranked genes:",
  length(primary_rank),
  "\n"
)

cat(
  "Sensitivity ranked genes:",
  length(sensitivity_rank),
  "\n"
)


# ============================================================
# 8. LOAD HALLMARK GENE SETS
# ============================================================

hallmark_raw <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)


cat("\n========================================\n")
cat("MSIGDB COLUMN QC\n")
cat("========================================\n\n")

print(
  colnames(hallmark_raw)
)


# Detect gene-symbol column automatically

if (
  "gene_symbol" %in%
  colnames(hallmark_raw)
) {
  
  symbol_column <- "gene_symbol"
  
} else if (
  "human_gene_symbol" %in%
  colnames(hallmark_raw)
) {
  
  symbol_column <- "human_gene_symbol"
  
} else {
  
  stop(
    "Could not identify gene-symbol column in msigdbr output."
  )
}


stopifnot(
  "gs_name" %in%
    colnames(hallmark_raw)
)


hallmark_sets <- data.frame(
  
  pathway =
    hallmark_raw$gs_name,
  
  gene =
    hallmark_raw[
      ,
      symbol_column
    ],
  
  stringsAsFactors = FALSE
)


hallmark_sets <- hallmark_sets[
  !is.na(hallmark_sets$gene) &
    hallmark_sets$gene != "",
  ,
  drop = FALSE
]


hallmark_sets <- unique(
  hallmark_sets
)


n_hallmarks <- length(
  unique(
    hallmark_sets$pathway
  )
)


cat("\n========================================\n")
cat("HALLMARK GENE SET QC\n")
cat("========================================\n\n")

cat(
  "Hallmark pathways:",
  n_hallmarks,
  "\n"
)

cat(
  "Unique Hallmark genes:",
  length(
    unique(
      hallmark_sets$gene
    )
  ),
  "\n"
)


stopifnot(
  n_hallmarks == 50
)


# ============================================================
# 9. GSEA FUNCTION
# ============================================================

run_gsea <- function(
    rank_vector,
    label
) {
  
  set.seed(123)
  
  
  result_object <- clusterProfiler::GSEA(
    
    geneList =
      rank_vector,
    
    TERM2GENE =
      hallmark_sets,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    verbose =
      FALSE,
    
    by =
      "fgsea"
  )
  
  
  result <- as.data.frame(
    result_object
  )
  
  
  cat("\n========================================\n")
  cat(label, "\n")
  cat("========================================\n\n")
  
  
  cat(
    "Hallmarks returned:",
    nrow(result),
    "\n"
  )
  
  
  cat(
    "FDR < 0.05:",
    sum(
      result$p.adjust < 0.05
    ),
    "\n"
  )
  
  
  return(
    result
  )
}


# ============================================================
# 10. PRIMARY GSEA
# ============================================================

primary_gsea <- run_gsea(
  primary_rank,
  "PRIMARY HALLMARK GSEA"
)


# ============================================================
# 11. SENSITIVITY GSEA
# ============================================================

sensitivity_gsea <- run_gsea(
  sensitivity_rank,
  "SENSITIVITY HALLMARK GSEA"
)


stopifnot(
  nrow(primary_gsea) == 50,
  nrow(sensitivity_gsea) == 50
)


# ============================================================
# 12. PREPARE COMPARISON TABLE
# ============================================================

primary_clean <- data.frame(
  
  pathway =
    primary_gsea$ID,
  
  primary_NES =
    primary_gsea$NES,
  
  primary_pvalue =
    primary_gsea$pvalue,
  
  primary_padj =
    primary_gsea$p.adjust,
  
  stringsAsFactors = FALSE
)


sensitivity_clean <- data.frame(
  
  pathway =
    sensitivity_gsea$ID,
  
  sensitivity_NES =
    sensitivity_gsea$NES,
  
  sensitivity_pvalue =
    sensitivity_gsea$pvalue,
  
  sensitivity_padj =
    sensitivity_gsea$p.adjust,
  
  stringsAsFactors = FALSE
)


comparison <- merge(
  
  primary_clean,
  
  sensitivity_clean,
  
  by = "pathway",
  
  all = TRUE,
  
  sort = FALSE
)


stopifnot(
  nrow(comparison) == 50
)


# ============================================================
# 13. GLOBAL ROBUSTNESS
# ============================================================

comparison$same_direction <-
  sign(
    comparison$primary_NES
  ) ==
  sign(
    comparison$sensitivity_NES
  )


comparison$primary_significant <-
  comparison$primary_padj < 0.05


comparison$sensitivity_significant <-
  comparison$sensitivity_padj < 0.05


comparison$robust_significant <-
  
  comparison$same_direction &
  
  comparison$primary_significant &
  
  comparison$sensitivity_significant


nes_rho <- cor(
  
  comparison$primary_NES,
  
  comparison$sensitivity_NES,
  
  method = "spearman"
)


cat("\n========================================\n")
cat("PRIMARY vs SENSITIVITY HALLMARK CONCORDANCE\n")
cat("========================================\n\n")


cat(
  "Hallmarks compared:",
  nrow(comparison),
  "\n"
)


cat(
  "Spearman NES rho:",
  round(
    nes_rho,
    4
  ),
  "\n"
)


cat(
  "Same direction:",
  sum(
    comparison$same_direction
  ),
  "/ 50\n"
)


cat(
  "Direction agreement:",
  round(
    100 *
      mean(
        comparison$same_direction
      ),
    1
  ),
  "%\n"
)


cat(
  "Primary significant:",
  sum(
    comparison$primary_significant
  ),
  "/ 50\n"
)


cat(
  "Sensitivity significant:",
  sum(
    comparison$sensitivity_significant
  ),
  "/ 50\n"
)


cat(
  "Significant in both + same direction:",
  sum(
    comparison$robust_significant
  ),
  "/ 50\n"
)


# ============================================================
# 14. PRIMARY SIGNIFICANT PATHWAYS
# ============================================================

primary_sig <- comparison[
  comparison$primary_significant,
  ,
  drop = FALSE
]


primary_sig <- primary_sig[
  order(
    -abs(
      primary_sig$primary_NES
    )
  ),
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("PRIMARY SIGNIFICANT HALLMARKS\n")
cat("========================================\n\n")


print(
  primary_sig[
    ,
    c(
      "pathway",
      "primary_NES",
      "primary_padj",
      "sensitivity_NES",
      "sensitivity_padj",
      "same_direction",
      "robust_significant"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. ROBUST SIGNIFICANT PATHWAYS
# ============================================================

robust_sig <- comparison[
  comparison$robust_significant,
  ,
  drop = FALSE
]


robust_sig <- robust_sig[
  order(
    -abs(
      robust_sig$primary_NES
    )
  ),
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("ROBUST HALLMARKS\n")
cat("Significant + same direction in both analyses\n")
cat("========================================\n\n")


print(
  robust_sig[
    ,
    c(
      "pathway",
      "primary_NES",
      "primary_padj",
      "sensitivity_NES",
      "sensitivity_padj"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 16. PRE-SPECIFIED D1 VALIDATION HYPOTHESIS
# ============================================================

frozen <- data.frame(
  
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


frozen_result <- merge(
  
  frozen,
  
  comparison,
  
  by = "pathway",
  
  all.x = TRUE,
  
  sort = FALSE
)


frozen_result$primary_direction <-
  ifelse(
    frozen_result$primary_NES > 0,
    "UP",
    "DOWN"
  )


frozen_result$sensitivity_direction <-
  ifelse(
    frozen_result$sensitivity_NES > 0,
    "UP",
    "DOWN"
  )


frozen_result$primary_direction_match <-
  
  sign(
    frozen_result$primary_NES
  ) ==
  frozen_result$expected_sign


frozen_result$sensitivity_direction_match <-
  
  sign(
    frozen_result$sensitivity_NES
  ) ==
  frozen_result$expected_sign


frozen_result$primary_match_significant <-
  
  frozen_result$primary_direction_match &
  
  frozen_result$primary_padj < 0.05


frozen_result$robust_match_significant <-
  
  frozen_result$primary_direction_match &
  
  frozen_result$sensitivity_direction_match &
  
  frozen_result$primary_padj < 0.05 &
  
  frozen_result$sensitivity_padj < 0.05


cat("\n========================================\n")
cat("PRE-SPECIFIED D1 VALIDATION HYPOTHESIS\n")
cat("========================================\n\n")


cat(
  "Primary directional matches:",
  sum(
    frozen_result$primary_direction_match
  ),
  "/ 10\n"
)


cat(
  "Primary directional + significant:",
  sum(
    frozen_result$primary_match_significant
  ),
  "/ 10\n"
)


cat(
  "Robust directional + significant:",
  sum(
    frozen_result$robust_match_significant
  ),
  "/ 10\n\n"
)


print(
  frozen_result[
    ,
    c(
      "pathway",
      "expected_direction",
      "primary_NES",
      "primary_padj",
      "sensitivity_NES",
      "sensitivity_padj",
      "primary_direction_match",
      "robust_match_significant"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 17. DIRECTION-DISCORDANT PATHWAYS
# ============================================================

discordant <- comparison[
  !comparison$same_direction,
  ,
  drop = FALSE
]


if (
  nrow(discordant) > 0
) {
  
  discordant$NES_difference <-
    abs(
      discordant$primary_NES -
        discordant$sensitivity_NES
    )
  
  
  discordant <- discordant[
    order(
      -discordant$NES_difference
    ),
    ,
    drop = FALSE
  ]
}


cat("\n========================================\n")
cat("DIRECTION-DISCORDANT HALLMARKS\n")
cat("========================================\n\n")


if (
  nrow(discordant) == 0
) {
  
  cat(
    "None.\n"
  )
  
} else {
  
  print(
    discordant[
      ,
      c(
        "pathway",
        "primary_NES",
        "primary_padj",
        "sensitivity_NES",
        "sensitivity_padj"
      )
    ],
    row.names = FALSE
  )
}


# ============================================================
# 18. SAVE OUTPUTS
# ============================================================

outdir <- paste0(
  "results/GSE191231/",
  "Hallmark_GSEA"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  
  primary_gsea,
  
  file.path(
    outdir,
    "GSE191231_Hallmark_GSEA_primary_all.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  sensitivity_gsea,
  
  file.path(
    outdir,
    "GSE191231_Hallmark_GSEA_sensitivity_all.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  comparison,
  
  file.path(
    outdir,
    "GSE191231_Hallmark_GSEA_primary_vs_sensitivity.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  primary_sig,
  
  file.path(
    outdir,
    "GSE191231_Hallmark_GSEA_primary_significant.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  robust_sig,
  
  file.path(
    outdir,
    "GSE191231_Hallmark_GSEA_robust_significant.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  frozen_result,
  
  file.path(
    outdir,
    "GSE191231_D1_frozen_validation.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  data.frame(
    symbol = names(primary_rank),
    statistic = as.numeric(primary_rank)
  ),
  
  file.path(
    outdir,
    "GSE191231_primary_GSEA_rank.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  data.frame(
    symbol = names(sensitivity_rank),
    statistic = as.numeric(sensitivity_rank)
  ),
  
  file.path(
    outdir,
    "GSE191231_sensitivity_GSEA_rank.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 19. FINISH
# ============================================================

cat("\n========================================\n")
cat("GSE191231 HALLMARK GSEA COMPLETE\n")
cat("========================================\n\n")

cat(
  "Output directory:",
  outdir,
  "\n"
)