# ============================================================
# GSE165914 - D7
# AUTHOR DE VALIDATION
# FINAL CLEAN VERSION
#
# Purpose:
#   Validate our independent limma results against the
#   differential-expression results deposited by the authors.
#
# D7A:
#   MDA-MB-468 OlaR vs parental
#   Author contrast: M468OR15Ola/Ctrl
#
# D7B:
#   SUM1315 OlaR vs parental
#   Author contrast: SUMOR25Ola/Ctrl
#
# Our direction:
#   positive logFC = OlaR - Parental
#
# IMPORTANT:
#   Script 53 is NOT modified.
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

raw_dir <- "data_raw/GSE165914"

analysis_dir <- "results/GSE165914/D7_analysis"

validation_dir <- "results/GSE165914/D7_validation"

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(data.table)
library(ggplot2)


cat("\n")
cat("============================================================\n")
cat("GSE165914 - D7 AUTHOR DE VALIDATION\n")
cat("============================================================\n\n")


# ============================================================
# 2. FILE PATHS
# ============================================================

m468_author_file <- file.path(
  raw_dir,
  "GSE165914_M468P_OR_Diff_Expression_mshekhar.csv.gz"
)

sum_author_file <- file.path(
  raw_dir,
  "GSE165914_SUMP_OR_Diff_Expression_mshekhar.csv.gz"
)

m468_our_file <- file.path(
  analysis_dir,
  "D7A_M468_DE_limma.csv"
)

sum_our_file <- file.path(
  analysis_dir,
  "D7B_SUM1315_DE_limma.csv"
)


# ============================================================
# 3. VERIFY REQUIRED FILES
# ============================================================

required_files <- c(
  m468_author_file,
  sum_author_file,
  m468_our_file,
  sum_our_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat(
    "Missing files:\n"
  )
  
  print(
    missing_files
  )
  
  stop(
    "Required files are missing."
  )
}


cat(
  "All required files found.\n"
)


# ============================================================
# 4. SAFE COLUMN-NAME NORMALIZATION
#
# IMPORTANT:
# Lowercase FIRST.
# Then remove punctuation.
#
# Example:
#   log2FC -> log2fc
#   PValue -> pvalue
#   Pvalue -> pvalue
#   FDR    -> fdr
# ============================================================

normalize_colname <- function(x) {
  
  x <- tolower(x)
  
  x <- gsub(
    "[^a-z0-9]",
    "",
    x
  )
  
  return(x)
}


# Sanity check

test_names <- normalize_colname(
  c(
    "log2FC",
    "PValue",
    "Pvalue",
    "FDR"
  )
)

cat(
  "\nColumn-name normalization sanity check:\n"
)

print(
  test_names
)


if (!identical(
  test_names,
  c(
    "log2fc",
    "pvalue",
    "pvalue",
    "fdr"
  )
)) {
  
  stop(
    "Column normalization sanity check failed."
  )
}


# ============================================================
# 5. EXTRACT EXACT AUTHOR CONTRAST
# ============================================================

extract_author_contrast <- function(
    author_file,
    target_header
) {
  
  author <- data.table::fread(
    author_file,
    data.table = FALSE,
    check.names = FALSE
  )
  
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat(
    "Author file:",
    basename(author_file),
    "\n"
  )
  
  cat(
    "Dimensions:",
    nrow(author),
    "rows x",
    ncol(author),
    "columns\n"
  )
  
  
  # ----------------------------------------------------------
  # Author tables are repeated blocks:
  #
  # Gene | log2FC | PValue/Pvalue | FDR
  #
  # Block starts therefore occur at:
  # 1, 5, 9, 13, ...
  # ----------------------------------------------------------
  
  block_starts <- seq(
    from = 1,
    to = ncol(author),
    by = 4
  )
  
  block_headers <- colnames(author)[
    block_starts
  ]
  
  
  cat(
    "\nAvailable contrasts:\n\n"
  )
  
  print(
    data.frame(
      start_column = block_starts,
      contrast = block_headers,
      stringsAsFactors = FALSE
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Find EXACT target header
  # ----------------------------------------------------------
  
  matches <- which(
    block_headers == target_header
  )
  
  
  if (length(matches) == 0) {
    
    stop(
      paste0(
        "Exact author contrast not found: ",
        target_header
      )
    )
  }
  
  
  if (length(matches) > 1) {
    
    stop(
      paste0(
        "Multiple exact matches found for: ",
        target_header
      )
    )
  }
  
  
  start_col <- block_starts[
    matches
  ]
  
  
  if (start_col + 3 > ncol(author)) {
    
    stop(
      "Selected contrast does not contain four complete columns."
    )
  }
  
  
  # ----------------------------------------------------------
  # Inspect following three columns
  # ----------------------------------------------------------
  
  following_columns <- colnames(author)[
    (start_col + 1):(start_col + 3)
  ]
  
  
  cat(
    "\nSelected contrast:",
    target_header,
    "\n"
  )
  
  cat(
    "Start column:",
    start_col,
    "\n"
  )
  
  cat(
    "Following columns:\n"
  )
  
  print(
    following_columns
  )
  
  
  following_normalized <- normalize_colname(
    following_columns
  )
  
  
  expected_normalized <- c(
    "log2fc",
    "pvalue",
    "fdr"
  )
  
  
  cat(
    "\nNormalized following columns:\n"
  )
  
  print(
    following_normalized
  )
  
  
  if (!identical(
    following_normalized,
    expected_normalized
  )) {
    
    stop(
      paste0(
        "Unexpected block structure for ",
        target_header
      )
    )
  }
  
  
  cat(
    "Author block structure verified.\n"
  )
  
  
  # ----------------------------------------------------------
  # Extract exact 4-column block
  # ----------------------------------------------------------
  
  block <- author[
    ,
    start_col:(start_col + 3),
    drop = FALSE
  ]
  
  
  colnames(block) <- c(
    "gene",
    "author_logFC",
    "author_PValue",
    "author_FDR"
  )
  
  
  # ----------------------------------------------------------
  # Convert types
  # ----------------------------------------------------------
  
  block$gene <- trimws(
    as.character(
      block$gene
    )
  )
  
  
  block$author_logFC <- suppressWarnings(
    as.numeric(
      block$author_logFC
    )
  )
  
  
  block$author_PValue <- suppressWarnings(
    as.numeric(
      block$author_PValue
    )
  )
  
  
  block$author_FDR <- suppressWarnings(
    as.numeric(
      block$author_FDR
    )
  )
  
  
  # ----------------------------------------------------------
  # Remove unusable rows
  # ----------------------------------------------------------
  
  block <- block[
    !is.na(block$gene) &
      block$gene != "" &
      is.finite(block$author_logFC),
    ,
    drop = FALSE
  ]
  
  
  cat(
    "Valid rows:",
    nrow(block),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Duplicate handling
  # ----------------------------------------------------------
  
  duplicate_n <- sum(
    duplicated(
      block$gene
    )
  )
  
  
  cat(
    "Duplicated gene symbols:",
    duplicate_n,
    "\n"
  )
  
  
  if (duplicate_n > 0) {
    
    block <- block[
      order(
        abs(
          block$author_logFC
        ),
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    
    block <- block[
      !duplicated(
        block$gene
      ),
      ,
      drop = FALSE
    ]
  }
  
  
  cat(
    "Unique genes retained:",
    nrow(block),
    "\n"
  )
  
  
  return(block)
}


# ============================================================
# 6. VALIDATION FUNCTION
# ============================================================

validate_model <- function(
    model_name,
    author_file,
    our_file,
    target_header
) {
  
  
  cat("\n\n")
  cat("============================================================\n")
  cat(model_name, "\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # AUTHOR RESULTS
  # ----------------------------------------------------------
  
  author_block <- extract_author_contrast(
    author_file = author_file,
    target_header = target_header
  )
  
  
  cat(
    "\nFirst six rows of selected author contrast:\n\n"
  )
  
  print(
    head(
      author_block,
      6
    )
  )
  
  
  # ----------------------------------------------------------
  # OUR LIMMA RESULTS
  # ----------------------------------------------------------
  
  ours <- read.csv(
    our_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  required_columns <- c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
  
  
  if (!all(
    required_columns %in%
    colnames(ours)
  )) {
    
    stop(
      paste0(
        model_name,
        ": expected columns missing from our DE file."
      )
    )
  }
  
  
  ours_small <- data.frame(
    
    gene = trimws(
      as.character(
        ours$gene
      )
    ),
    
    our_logFC = as.numeric(
      ours$logFC
    ),
    
    our_FDR = as.numeric(
      ours$adj.P.Val
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # MERGE
  # ----------------------------------------------------------
  
  merged <- merge(
    ours_small,
    author_block,
    by = "gene"
  )
  
  
  merged <- merged[
    is.finite(
      merged$our_logFC
    ) &
      is.finite(
        merged$author_logFC
      ),
    ,
    drop = FALSE
  ]
  
  
  cat(
    "\nMatched genes:",
    nrow(merged),
    "\n"
  )
  
  
  if (nrow(merged) < 1000) {
    
    stop(
      paste0(
        model_name,
        ": too few overlapping genes."
      )
    )
  }
  
  
  # ==========================================================
  # 7. RAW CORRELATION
  # ==========================================================
  
  raw_spearman <- cor(
    merged$our_logFC,
    merged$author_logFC,
    method = "spearman"
  )
  
  
  raw_pearson <- cor(
    merged$our_logFC,
    merged$author_logFC,
    method = "pearson"
  )
  
  
  cat(
    "\nRaw Spearman rho:",
    round(
      raw_spearman,
      4
    ),
    "\n"
  )
  
  
  cat(
    "Raw Pearson r:",
    round(
      raw_pearson,
      4
    ),
    "\n"
  )
  
  
  # ==========================================================
  # 8. DETERMINE ORIENTATION
  #
  # If raw correlation is negative, author logFC is opposite
  # to our OlaR - Parental convention.
  # ==========================================================
  
  if (raw_spearman >= 0) {
    
    multiplier <- 1
    
    orientation <-
      "Author orientation matches OlaR - Parental"
    
  } else {
    
    multiplier <- -1
    
    orientation <-
      paste0(
        "Author orientation is reversed; ",
        "logFC multiplied by -1 for validation"
      )
  }
  
  
  merged$author_logFC_aligned <-
    merged$author_logFC *
    multiplier
  
  
  cat(
    "\nOrientation:",
    orientation,
    "\n"
  )
  
  
  # ==========================================================
  # 9. ALIGNED CORRELATIONS
  # ==========================================================
  
  aligned_spearman <- cor(
    merged$our_logFC,
    merged$author_logFC_aligned,
    method = "spearman"
  )
  
  
  aligned_pearson <- cor(
    merged$our_logFC,
    merged$author_logFC_aligned,
    method = "pearson"
  )
  
  
  cat(
    "\nAligned Spearman rho:",
    round(
      aligned_spearman,
      4
    ),
    "\n"
  )
  
  
  cat(
    "Aligned Pearson r:",
    round(
      aligned_pearson,
      4
    ),
    "\n"
  )
  
  
  # ==========================================================
  # 10. GENOME-WIDE DIRECTION AGREEMENT
  # ==========================================================
  
  nonzero <- merged[
    merged$our_logFC != 0 &
      merged$author_logFC_aligned != 0,
    ,
    drop = FALSE
  ]
  
  
  genomewide_sign_agreement <- mean(
    sign(
      nonzero$our_logFC
    ) ==
      sign(
        nonzero$author_logFC_aligned
      )
  )
  
  
  cat(
    "\nGenome-wide non-zero genes:",
    nrow(nonzero),
    "\n"
  )
  
  
  cat(
    "Genome-wide sign agreement:",
    round(
      100 *
        genomewide_sign_agreement,
      2
    ),
    "%\n"
  )
  
  
  # ==========================================================
  # 11. OUR LARGE DEGs
  # ==========================================================
  
  our_large <- merged[
    merged$our_FDR < 0.05 &
      abs(
        merged$our_logFC
      ) >= 1,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(our_large) >= 2) {
    
    large_spearman <- cor(
      our_large$our_logFC,
      our_large$author_logFC_aligned,
      method = "spearman"
    )
    
    
    large_sign_agreement <- mean(
      sign(
        our_large$our_logFC
      ) ==
        sign(
          our_large$author_logFC_aligned
        )
    )
    
  } else {
    
    large_spearman <- NA_real_
    large_sign_agreement <- NA_real_
    
  }
  
  
  cat(
    "\nOur large DEGs:",
    nrow(our_large),
    "\n"
  )
  
  
  cat(
    "Large-DEG Spearman rho:",
    round(
      large_spearman,
      4
    ),
    "\n"
  )
  
  
  cat(
    "Large-DEG sign agreement:",
    round(
      100 *
        large_sign_agreement,
      2
    ),
    "%\n"
  )
  
  
  # ==========================================================
  # 12. AUTHOR-SIGNIFICANT GENES
  # ==========================================================
  
  author_sig <- merged[
    is.finite(
      merged$author_FDR
    ) &
      merged$author_FDR < 0.05,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(author_sig) > 0) {
    
    author_sig_sign_agreement <- mean(
      sign(
        author_sig$our_logFC
      ) ==
        sign(
          author_sig$author_logFC_aligned
        )
    )
    
  } else {
    
    author_sig_sign_agreement <- NA_real_
    
  }
  
  
  cat(
    "\nAuthor FDR < 0.05 genes:",
    nrow(author_sig),
    "\n"
  )
  
  
  cat(
    "Direction agreement among author-significant genes:",
    round(
      100 *
        author_sig_sign_agreement,
      2
    ),
    "%\n"
  )
  
  
  # ==========================================================
  # 13. SIGNIFICANT IN BOTH
  # ==========================================================
  
  both_sig <- merged[
    merged$our_FDR < 0.05 &
      is.finite(
        merged$author_FDR
      ) &
      merged$author_FDR < 0.05,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(both_sig) > 0) {
    
    both_sig_sign_agreement <- mean(
      sign(
        both_sig$our_logFC
      ) ==
        sign(
          both_sig$author_logFC_aligned
        )
    )
    
  } else {
    
    both_sig_sign_agreement <- NA_real_
    
  }
  
  
  cat(
    "\nGenes significant in BOTH analyses:",
    nrow(both_sig),
    "\n"
  )
  
  
  cat(
    "Direction agreement among both-significant genes:",
    round(
      100 *
        both_sig_sign_agreement,
      2
    ),
    "%\n"
  )
  
  
  # ==========================================================
  # 14. TOP 20 EFFECTS
  # ==========================================================
  
  top_effects <- merged[
    order(
      abs(
        merged$our_logFC
      ),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  top_effects <- head(
    top_effects,
    20
  )
  
  
  cat(
    "\nTop 20 genes by |our logFC|:\n\n"
  )
  
  
  print(
    top_effects[
      ,
      c(
        "gene",
        "our_logFC",
        "author_logFC",
        "author_logFC_aligned",
        "our_FDR",
        "author_FDR"
      )
    ],
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 15. SAVE MERGED TABLE
  # ==========================================================
  
  write.csv(
    merged,
    file.path(
      validation_dir,
      paste0(
        model_name,
        "_author_vs_our_DE.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 16. SCATTER PLOT
  # ==========================================================
  
  p <- ggplot(
    merged,
    aes(
      x = author_logFC_aligned,
      y = our_logFC
    )
  ) +
    
    geom_point(
      alpha = 0.20,
      size = 0.75
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    
    labs(
      
      title = paste0(
        model_name,
        ": author vs independent DE"
      ),
      
      subtitle = paste0(
        "Spearman rho = ",
        round(
          aligned_spearman,
          3
        ),
        "; sign agreement = ",
        round(
          100 *
            genomewide_sign_agreement,
          1
        ),
        "%"
      ),
      
      x =
        "Author log2FC aligned to OlaR - Parental",
      
      y =
        "Our limma logFC (OlaR - Parental)"
    ) +
    
    theme_classic(
      base_size = 12
    )
  
  
  ggsave(
    filename = file.path(
      validation_dir,
      paste0(
        model_name,
        "_author_vs_our_logFC.pdf"
      )
    ),
    plot = p,
    width = 7,
    height = 6
  )
  
  
  # ==========================================================
  # 17. RETURN SUMMARY
  # ==========================================================
  
  result <- data.frame(
    
    model =
      model_name,
    
    author_contrast =
      target_header,
    
    matched_genes =
      nrow(merged),
    
    raw_spearman_rho =
      raw_spearman,
    
    raw_pearson_r =
      raw_pearson,
    
    orientation_multiplier =
      multiplier,
    
    aligned_spearman_rho =
      aligned_spearman,
    
    aligned_pearson_r =
      aligned_pearson,
    
    genomewide_sign_agreement =
      genomewide_sign_agreement,
    
    our_large_DEGs =
      nrow(our_large),
    
    large_DEG_spearman_rho =
      large_spearman,
    
    large_DEG_sign_agreement =
      large_sign_agreement,
    
    author_FDR_lt_0.05 =
      nrow(author_sig),
    
    author_sig_sign_agreement =
      author_sig_sign_agreement,
    
    both_FDR_lt_0.05 =
      nrow(both_sig),
    
    both_sig_sign_agreement =
      both_sig_sign_agreement,
    
    stringsAsFactors = FALSE
  )
  
  
  return(result)
}


# ============================================================
# 18. D7A — MDA-MB-468
# ============================================================

D7A_validation <- validate_model(
  
  model_name =
    "D7A_M468",
  
  author_file =
    m468_author_file,
  
  our_file =
    m468_our_file,
  
  target_header =
    "M468OR15Ola/Ctrl"
  
)


# ============================================================
# 19. D7B — SUM1315
# ============================================================

D7B_validation <- validate_model(
  
  model_name =
    "D7B_SUM1315",
  
  author_file =
    sum_author_file,
  
  our_file =
    sum_our_file,
  
  target_header =
    "SUMOR25Ola/Ctrl"
  
)


# ============================================================
# 20. COMBINE RESULTS
# ============================================================

validation_summary <- rbind(
  D7A_validation,
  D7B_validation
)


rownames(
  validation_summary
) <- NULL


write.csv(
  validation_summary,
  file.path(
    validation_dir,
    "D7_author_DE_validation_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. FINAL SUMMARY
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("FINAL D7 VALIDATION SUMMARY\n")
cat("============================================================\n\n")


print(
  validation_summary,
  row.names = FALSE
)


# ============================================================
# 22. DIAGNOSTIC FLAGS
#
# These are NOT biological thresholds.
# They are simply validation diagnostics.
# ============================================================

validation_summary$effect_rho_ge_0.8 <-
  validation_summary$aligned_spearman_rho >= 0.8


validation_summary$genomewide_sign_ge_0.8 <-
  validation_summary$genomewide_sign_agreement >= 0.8


validation_summary$large_DEG_sign_ge_0.9 <-
  validation_summary$large_DEG_sign_agreement >= 0.9


cat("\n")
cat("============================================================\n")
cat("VALIDATION DIAGNOSTICS\n")
cat("============================================================\n\n")


print(
  validation_summary[
    ,
    c(
      "model",
      "aligned_spearman_rho",
      "aligned_pearson_r",
      "genomewide_sign_agreement",
      "large_DEG_spearman_rho",
      "large_DEG_sign_agreement",
      "effect_rho_ge_0.8",
      "genomewide_sign_ge_0.8",
      "large_DEG_sign_ge_0.9"
    )
  ],
  row.names = FALSE
)


write.csv(
  validation_summary,
  file.path(
    validation_dir,
    "D7_author_DE_validation_summary_with_flags.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. SESSION INFO
# ============================================================

sink(
  file.path(
    validation_dir,
    "D7_validation_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D7 AUTHOR VALIDATION COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Script 53 was NOT modified.\n",
  "No biological thresholds were changed.\n",
  "This analysis only validates our effect estimates against ",
  "the authors' deposited DE results.\n"
)

cat("\n============================================================\n")