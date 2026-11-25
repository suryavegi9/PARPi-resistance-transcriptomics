# ============================================================
# GSE165914 - D7
# DIFFERENTIAL EXPRESSION + HALLMARK GSEA
#
# D7A:
#   MDA-MB-468 OlaR vs parental
#
# D7B:
#   SUM1315 OlaR vs parental
#
# Primary input:
#   Author-supplied continuous processed expression matrix
#
# Analysis:
#   - remove ambiguous duplicated identifiers 1-Mar / 2-Mar
#   - filter expression >= 1 in >= 3/6 samples
#   - log2(x + 1)
#   - limma + empirical-Bayes trend
#   - DEG: BH FDR < 0.05 and |log2FC| >= 1
#   - Hallmark GSEA ranked by moderated t-statistic
#
# D7 accession-level collapse:
#   - same-direction support only when D7A and D7B NES signs agree
#   - significant concordant support only when BOTH FDR < 0.05
#
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

raw_dir <- "data_raw/GSE165914"
qc_dir <- "results/GSE165914/QC"
out_dir <- "results/GSE165914/D7_analysis"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "limma",
  "clusterProfiler"
)

cran_packages <- c(
  "msigdbr",
  "dplyr",
  "ggplot2"
)


for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
}


for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}


library(limma)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(ggplot2)


set.seed(20261110)


cat("\n")
cat("============================================================\n")
cat("GSE165914 - D7 DE + HALLMARK GSEA\n")
cat("============================================================\n\n")


# ============================================================
# 2. LOAD AUTHOR MATRIX
# ============================================================

input_file <- file.path(
  raw_dir,
  "GSE165914_author_counts_raw.rds"
)


if (!file.exists(input_file)) {
  stop("Author expression matrix not found.")
}


dat <- readRDS(
  input_file
)


gene_symbol <- trimws(
  as.character(
    dat[[1]]
  )
)


expr_all <- as.data.frame(
  dat[
    ,
    -1,
    drop = FALSE
  ],
  check.names = FALSE
)


cat(
  "Original matrix:",
  nrow(expr_all),
  "genes x",
  ncol(expr_all),
  "samples\n"
)


# ============================================================
# 3. REMOVE INVALID / AMBIGUOUS IDENTIFIERS
# ============================================================

# These appeared as duplicated spreadsheet-like identifiers.
# We do NOT attempt retrospective biological renaming.

ambiguous_ids <- c(
  "1-Mar",
  "2-Mar"
)


valid_gene <- (
  !is.na(gene_symbol) &
    gene_symbol != "" &
    !(gene_symbol %in% ambiguous_ids)
)


gene_symbol <- gene_symbol[
  valid_gene
]


expr_all <- expr_all[
  valid_gene,
  ,
  drop = FALSE
]


cat(
  "Rows after removing blank/ambiguous identifiers:",
  nrow(expr_all),
  "\n"
)


# Verify uniqueness

if (anyDuplicated(gene_symbol) > 0) {
  
  stop(
    "Unexpected duplicated gene symbols remain after removing ambiguous rows."
  )
  
}


rownames(expr_all) <- gene_symbol


# ============================================================
# 4. DEFINE PRIMARY CONTRASTS
# ============================================================

m468_parental <- c(
  "M468PCtrl1",
  "M468PCtrl2",
  "M468PCtrl3"
)


m468_resistant <- c(
  "M468OR15Ola1",
  "M468OR15Ola2",
  "M468OR15Ola3"
)


sum_parental <- c(
  "SUMPCtrl1",
  "SUMPCtrl2",
  "SUMPCtrl3"
)


sum_resistant <- c(
  "SUMOR25Ola1",
  "SUMOR25Ola2",
  "SUMOR25Ola3"
)


# ============================================================
# 5. HALLMARK GENE SETS
# ============================================================

hallmark_raw <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)


gene_col <- if (
  "gene_symbol" %in% colnames(hallmark_raw)
) {
  "gene_symbol"
} else if (
  "human_gene_symbol" %in% colnames(hallmark_raw)
) {
  "human_gene_symbol"
} else {
  stop("Could not identify Hallmark gene-symbol column in msigdbr.")
}


hallmark_term2gene <- hallmark_raw |>
  dplyr::select(
    gs_name,
    all_of(gene_col)
  ) |>
  dplyr::rename(
    gene = all_of(gene_col)
  ) |>
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) |>
  dplyr::distinct()


cat(
  "Hallmark pathways retrieved:",
  length(
    unique(
      hallmark_term2gene$gs_name
    )
  ),
  "\n\n"
)


# ============================================================
# 6. FUNCTION: ONE MODEL
# ============================================================

analyse_model <- function(
    model_name,
    parental_samples,
    resistant_samples
) {
  
  
  cat("\n")
  cat("============================================================\n")
  cat(model_name, "\n")
  cat("============================================================\n\n")
  
  
  samples <- c(
    parental_samples,
    resistant_samples
  )
  
  
  missing_samples <- setdiff(
    samples,
    colnames(expr_all)
  )
  
  
  if (length(missing_samples) > 0) {
    
    stop(
      paste(
        "Missing samples:",
        paste(
          missing_samples,
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Extract
  # ----------------------------------------------------------
  
  expr <- as.matrix(
    expr_all[
      ,
      samples,
      drop = FALSE
    ]
  )
  
  
  storage.mode(expr) <- "numeric"
  
  
  # ----------------------------------------------------------
  # Filter
  #
  # Frozen before DE:
  # expression >= 1 in >= 3/6 samples
  # ----------------------------------------------------------
  
  keep <- rowSums(
    expr >= 1
  ) >= 3
  
  
  expr_filtered <- expr[
    keep,
    ,
    drop = FALSE
  ]
  
  
  cat(
    "Genes before filtering:",
    nrow(expr),
    "\n"
  )
  
  
  cat(
    "Genes retained:",
    nrow(expr_filtered),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # log2 transformation
  # ----------------------------------------------------------
  
  log_expr <- log2(
    expr_filtered + 1
  )
  
  
  # ----------------------------------------------------------
  # Design
  # ----------------------------------------------------------
  
  state <- factor(
    c(
      rep("Parental", length(parental_samples)),
      rep("OlaR", length(resistant_samples))
    ),
    levels = c(
      "Parental",
      "OlaR"
    )
  )
  
  
  design <- model.matrix(
    ~ state
  )
  
  
  rownames(design) <- samples
  
  
  cat(
    "\nDesign matrix:\n"
  )
  
  print(
    design
  )
  
  
  # ----------------------------------------------------------
  # limma
  #
  # Coefficient stateOlaR:
  # positive = higher in resistant
  # negative = lower in resistant
  # ----------------------------------------------------------
  
  fit <- limma::lmFit(
    log_expr,
    design
  )
  
  
  fit <- limma::eBayes(
    fit,
    trend = TRUE
  )
  
  
  de <- limma::topTable(
    fit,
    coef = "stateOlaR",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  
  de$gene <- rownames(de)
  
  
  de <- de[
    ,
    c(
      "gene",
      setdiff(
        colnames(de),
        "gene"
      )
    )
  ]
  
  
  # ----------------------------------------------------------
  # DEG classification
  # ----------------------------------------------------------
  
  de$DEG <- (
    de$adj.P.Val < 0.05 &
      abs(de$logFC) >= 1
  )
  
  
  de$direction <- "NS"
  
  
  de$direction[
    de$DEG &
      de$logFC > 0
  ] <- "UP"
  
  
  de$direction[
    de$DEG &
      de$logFC < 0
  ] <- "DOWN"
  
  
  n_fdr <- sum(
    de$adj.P.Val < 0.05,
    na.rm = TRUE
  )
  
  
  n_deg <- sum(
    de$DEG,
    na.rm = TRUE
  )
  
  
  n_up <- sum(
    de$direction == "UP",
    na.rm = TRUE
  )
  
  
  n_down <- sum(
    de$direction == "DOWN",
    na.rm = TRUE
  )
  
  
  cat(
    "\nFDR < 0.05 genes:",
    n_fdr,
    "\n"
  )
  
  
  cat(
    "FDR < 0.05 & |log2FC| >= 1:",
    n_deg,
    "\n"
  )
  
  
  cat(
    "UP:",
    n_up,
    "\n"
  )
  
  
  cat(
    "DOWN:",
    n_down,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Save DE
  # ----------------------------------------------------------
  
  write.csv(
    de,
    file.path(
      out_dir,
      paste0(
        model_name,
        "_DE_limma.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Ranked list for GSEA
  #
  # Moderated t-statistic
  # ----------------------------------------------------------
  
  rank_df <- data.frame(
    
    gene =
      rownames(fit$t),
    
    t =
      fit$t[
        ,
        "stateOlaR"
      ],
    
    stringsAsFactors = FALSE
    
  )
  
  
  rank_df <- rank_df[
    is.finite(rank_df$t),
    ,
    drop = FALSE
  ]
  
  
  # There should be one row per gene,
  # but keep strongest |t| defensively if duplicates appear.
  
  rank_df <- rank_df[
    order(
      abs(rank_df$t),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  rank_df <- rank_df[
    !duplicated(rank_df$gene),
    ,
    drop = FALSE
  ]
  
  
  gene_list <- rank_df$t
  
  names(gene_list) <-
    rank_df$gene
  
  
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
  
  
  cat(
    "\nRanked genes for GSEA:",
    length(gene_list),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Hallmark GSEA
  # ----------------------------------------------------------
  
  gsea <- clusterProfiler::GSEA(
    
    geneList =
      gene_list,
    
    TERM2GENE =
      hallmark_term2gene,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    eps =
      0,
    
    seed =
      TRUE,
    
    verbose =
      FALSE
    
  )
  
  
  gsea_df <- as.data.frame(
    gsea
  )
  
  
  # ----------------------------------------------------------
  # Ensure all 50 Hallmarks represented
  # ----------------------------------------------------------
  
  all_hallmarks <- sort(
    unique(
      hallmark_term2gene$gs_name
    )
  )
  
  
  if (nrow(gsea_df) > 0) {
    
    gsea_keep <- gsea_df[
      ,
      c(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue"
      ),
      drop = FALSE
    ]
    
    
    names(gsea_keep)[
      names(gsea_keep) == "ID"
    ] <- "pathway"
    
    
  } else {
    
    gsea_keep <- data.frame(
      pathway = character(),
      Description = character(),
      setSize = numeric(),
      enrichmentScore = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric()
    )
    
  }
  
  
  hallmark_complete <- data.frame(
    pathway = all_hallmarks,
    stringsAsFactors = FALSE
  ) |>
    left_join(
      gsea_keep,
      by = "pathway"
    )
  
  
  hallmark_complete$significant <-
    hallmark_complete$p.adjust < 0.05
  
  
  hallmark_complete$direction <-
    ifelse(
      hallmark_complete$NES > 0,
      "UP",
      ifelse(
        hallmark_complete$NES < 0,
        "DOWN",
        NA
      )
    )
  
  
  hallmark_complete <- hallmark_complete[
    order(
      hallmark_complete$p.adjust,
      -abs(
        hallmark_complete$NES
      )
    ),
    ,
    drop = FALSE
  ]
  
  
  write.csv(
    hallmark_complete,
    file.path(
      out_dir,
      paste0(
        model_name,
        "_Hallmark_GSEA.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Console GSEA summary
  # ----------------------------------------------------------
  
  sig_gsea <- hallmark_complete[
    !is.na(hallmark_complete$p.adjust) &
      hallmark_complete$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  cat(
    "\nSignificant Hallmarks:",
    nrow(sig_gsea),
    "/ 50\n\n"
  )
  
  
  if (nrow(sig_gsea) > 0) {
    
    print(
      sig_gsea[
        order(
          -abs(sig_gsea$NES)
        ),
        c(
          "pathway",
          "NES",
          "p.adjust",
          "direction"
        )
      ],
      row.names = FALSE
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  
  list(
    
    filtered_genes =
      nrow(expr_filtered),
    
    de =
      de,
    
    n_fdr =
      n_fdr,
    
    n_deg =
      n_deg,
    
    n_up =
      n_up,
    
    n_down =
      n_down,
    
    rank =
      gene_list,
    
    gsea =
      hallmark_complete
    
  )
  
}


# ============================================================
# 7. D7A — MDA-MB-468
# ============================================================

D7A <- analyse_model(
  
  model_name =
    "D7A_M468",
  
  parental_samples =
    m468_parental,
  
  resistant_samples =
    m468_resistant
  
)


# ============================================================
# 8. D7B — SUM1315
# ============================================================

D7B <- analyse_model(
  
  model_name =
    "D7B_SUM1315",
  
  parental_samples =
    sum_parental,
  
  resistant_samples =
    sum_resistant
  
)


# ============================================================
# 9. MERGE D7A + D7B HALLMARK RESULTS
# ============================================================

gsea_A <- D7A$gsea |>
  dplyr::select(
    pathway,
    NES,
    p.adjust
  ) |>
  dplyr::rename(
    M468_NES = NES,
    M468_FDR = p.adjust
  )


gsea_B <- D7B$gsea |>
  dplyr::select(
    pathway,
    NES,
    p.adjust
  ) |>
  dplyr::rename(
    SUM1315_NES = NES,
    SUM1315_FDR = p.adjust
  )


D7_integrated <- full_join(
  gsea_A,
  gsea_B,
  by = "pathway"
)


# ============================================================
# 10. CONSERVATIVE ACCESSION-LEVEL COLLAPSE
# ============================================================

D7_integrated$same_direction <-
  with(
    D7_integrated,
    is.finite(M468_NES) &
      is.finite(SUM1315_NES) &
      sign(M468_NES) ==
      sign(SUM1315_NES)
  )


D7_integrated$both_significant <-
  with(
    D7_integrated,
    !is.na(M468_FDR) &
      !is.na(SUM1315_FDR) &
      M468_FDR < 0.05 &
      SUM1315_FDR < 0.05
  )


D7_integrated$concordant_significant <-
  D7_integrated$same_direction &
  D7_integrated$both_significant


# Conservative accession NES:
# only defined when both models agree in direction.
#
# Mean is used only as a descriptive summary within accession;
# model-level NES values remain preserved separately.

D7_integrated$D7_NES <-
  ifelse(
    D7_integrated$same_direction,
    rowMeans(
      D7_integrated[
        ,
        c(
          "M468_NES",
          "SUM1315_NES"
        )
      ],
      na.rm = FALSE
    ),
    NA_real_
  )


D7_integrated$D7_direction <-
  ifelse(
    is.na(D7_integrated$D7_NES),
    "Discordant",
    ifelse(
      D7_integrated$D7_NES > 0,
      "UP",
      "DOWN"
    )
  )


D7_integrated$D7_significant <-
  D7_integrated$concordant_significant


D7_integrated <- D7_integrated[
  order(
    !D7_integrated$D7_significant,
    -abs(
      D7_integrated$D7_NES
    )
  ),
  ,
  drop = FALSE
]


write.csv(
  D7_integrated,
  file.path(
    out_dir,
    "D7_accession_Hallmark_integration.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. ACCESSION-LEVEL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D7 ACCESSION-LEVEL HALLMARK SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Same-direction Hallmarks:",
  sum(
    D7_integrated$same_direction,
    na.rm = TRUE
  ),
  "/ 50\n"
)


cat(
  "Opposite-direction Hallmarks:",
  sum(
    !D7_integrated$same_direction,
    na.rm = TRUE
  ),
  "/ 50\n"
)


cat(
  "Both models significant + same direction:",
  sum(
    D7_integrated$concordant_significant,
    na.rm = TRUE
  ),
  "/ 50\n\n"
)


cat(
  "Concordant significant Hallmarks:\n\n"
)


print(
  D7_integrated[
    D7_integrated$concordant_significant,
    c(
      "pathway",
      "M468_NES",
      "M468_FDR",
      "SUM1315_NES",
      "SUM1315_FDR",
      "D7_NES",
      "D7_direction"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 12. ORIGINAL D1 FROZEN VALIDATION PATHWAYS
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
  
  D1_expected_direction = c(
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


D7_frozen <- frozen |>
  left_join(
    D7_integrated,
    by = "pathway"
  )


D7_frozen$direction_match_D1 <-
  with(
    D7_frozen,
    D7_direction ==
      D1_expected_direction
  )


D7_frozen$significant_direction_match_D1 <-
  with(
    D7_frozen,
    direction_match_D1 &
      D7_significant
  )


write.csv(
  D7_frozen,
  file.path(
    out_dir,
    "D7_vs_D1_frozen_validation.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("D7 vs ORIGINAL D1 FROZEN HYPOTHESES\n")
cat("============================================================\n\n")


print(
  D7_frozen[
    ,
    c(
      "pathway",
      "D1_expected_direction",
      "M468_NES",
      "M468_FDR",
      "SUM1315_NES",
      "SUM1315_FDR",
      "D7_direction",
      "D7_significant",
      "direction_match_D1",
      "significant_direction_match_D1"
    )
  ],
  row.names = FALSE
)


cat(
  "\nD1-direction matches:",
  sum(
    D7_frozen$direction_match_D1,
    na.rm = TRUE
  ),
  "/ 10\n"
)


cat(
  "D1-direction matches with conservative D7 significance:",
  sum(
    D7_frozen$significant_direction_match_D1,
    na.rm = TRUE
  ),
  "/ 10\n"
)


# ============================================================
# 13. MODEL-TO-MODEL HALLMARK AGREEMENT
# ============================================================

complete_models <- D7_integrated[
  is.finite(D7_integrated$M468_NES) &
    is.finite(D7_integrated$SUM1315_NES),
  ,
  drop = FALSE
]


nes_rho <- cor(
  complete_models$M468_NES,
  complete_models$SUM1315_NES,
  method = "spearman"
)


same_direction_fraction <- mean(
  sign(
    complete_models$M468_NES
  ) ==
    sign(
      complete_models$SUM1315_NES
    )
)


cat("\n")
cat("============================================================\n")
cat("D7 MODEL-TO-MODEL AGREEMENT\n")
cat("============================================================\n\n")


cat(
  "Spearman rho of Hallmark NES:",
  round(
    nes_rho,
    4
  ),
  "\n"
)


cat(
  "Same-direction Hallmarks:",
  sum(
    sign(
      complete_models$M468_NES
    ) ==
      sign(
        complete_models$SUM1315_NES
      )
  ),
  "/",
  nrow(complete_models),
  "\n"
)


cat(
  "Same-direction fraction:",
  round(
    same_direction_fraction,
    4
  ),
  "\n"
)


# ============================================================
# 14. SAVE ANALYSIS SUMMARY
# ============================================================

summary_df <- data.frame(
  
  metric = c(
    
    "D7A_filtered_genes",
    "D7A_FDR_lt_0.05",
    "D7A_large_DEGs",
    "D7A_UP_large_DEGs",
    "D7A_DOWN_large_DEGs",
    
    "D7B_filtered_genes",
    "D7B_FDR_lt_0.05",
    "D7B_large_DEGs",
    "D7B_UP_large_DEGs",
    "D7B_DOWN_large_DEGs",
    
    "Hallmark_NES_Spearman_rho",
    "Hallmark_same_direction_n",
    "Hallmark_same_direction_fraction",
    "D7_concordant_significant_Hallmarks",
    
    "D1_frozen_direction_matches",
    "D1_frozen_significant_direction_matches"
    
  ),
  
  value = c(
    
    D7A$filtered_genes,
    D7A$n_fdr,
    D7A$n_deg,
    D7A$n_up,
    D7A$n_down,
    
    D7B$filtered_genes,
    D7B$n_fdr,
    D7B$n_deg,
    D7B$n_up,
    D7B$n_down,
    
    nes_rho,
    
    sum(
      D7_integrated$same_direction,
      na.rm = TRUE
    ),
    
    same_direction_fraction,
    
    sum(
      D7_integrated$concordant_significant,
      na.rm = TRUE
    ),
    
    sum(
      D7_frozen$direction_match_D1,
      na.rm = TRUE
    ),
    
    sum(
      D7_frozen$significant_direction_match_D1,
      na.rm = TRUE
    )
    
  )
  
)


write.csv(
  summary_df,
  file.path(
    out_dir,
    "D7_analysis_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. SAVE OBJECTS
# ============================================================

saveRDS(
  
  list(
    
    D7A =
      D7A,
    
    D7B =
      D7B,
    
    D7_accession =
      D7_integrated,
    
    D7_frozen =
      D7_frozen,
    
    summary =
      summary_df
    
  ),
  
  file.path(
    out_dir,
    "D7_analysis_objects.rds"
  )
  
)


# ============================================================
# 16. SESSION INFO
# ============================================================

sink(
  file.path(
    out_dir,
    "D7_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 17. FINAL OUTPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D7 ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


print(
  summary_df,
  row.names = FALSE
)


cat("\n")
cat(
  "Interpretation must be based on the actual outputs above.\n",
  "Do NOT modify filtering/GSEA rules after seeing results.\n"
)

cat("============================================================\n")