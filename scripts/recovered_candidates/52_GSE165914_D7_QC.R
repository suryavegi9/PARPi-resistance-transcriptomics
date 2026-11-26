# ============================================================
# GSE165914 - D7
# PRIMARY RESISTANCE CONTRAST QC
# CORRECTED VERSION
#
# D7A:
#   MDA-MB-468 parental vs OlaR
#
# D7B:
#   SUM1315 parental vs OlaR
#
# Author-supplied matrix contains processed continuous
# expression values, NOT integer raw counts.
#
# NO differential-expression testing is performed here.
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

raw_dir <- "data_raw/GSE165914"
results_dir <- "results/GSE165914/QC"

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}

library(ggplot2)
library(ggrepel)
library(pheatmap)


# ============================================================
# 2. LOAD AUTHOR-SUPPLIED MATRIX
# ============================================================

input_file <- file.path(
  raw_dir,
  "GSE165914_author_counts_raw.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Cannot find: ",
    input_file
  )
}

dat <- readRDS(
  input_file
)


cat("\n")
cat("============================================================\n")
cat("GSE165914 D7 QC\n")
cat("============================================================\n\n")


cat(
  "Loaded matrix:",
  nrow(dat),
  "rows x",
  ncol(dat),
  "columns\n"
)


# ============================================================
# 3. GENE IDENTIFIER QC
# ============================================================

gene_symbol <- as.character(
  dat[[1]]
)


expr_all <- dat[
  ,
  -1,
  drop = FALSE
]


cat(
  "Gene rows:",
  length(gene_symbol),
  "\n"
)


cat(
  "Unique gene symbols:",
  length(
    unique(gene_symbol)
  ),
  "\n"
)


cat(
  "Duplicated rows:",
  sum(
    duplicated(gene_symbol)
  ),
  "\n"
)


duplicate_symbols <- unique(
  gene_symbol[
    duplicated(gene_symbol) |
      duplicated(
        gene_symbol,
        fromLast = TRUE
      )
  ]
)


cat(
  "\nDuplicated gene identifiers:\n"
)

print(
  duplicate_symbols
)


if (length(duplicate_symbols) > 0) {
  
  duplicate_rows <- dat[
    gene_symbol %in% duplicate_symbols,
    ,
    drop = FALSE
  ]
  
  
  write.csv(
    duplicate_rows,
    file.path(
      results_dir,
      "D7_duplicate_gene_rows.csv"
    ),
    row.names = FALSE
  )
  
}


# ============================================================
# 4. DEFINE/FREEZE PRIMARY SAMPLE SET
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


selected_samples <- c(
  m468_parental,
  m468_resistant,
  sum_parental,
  sum_resistant
)


# ============================================================
# 5. VERIFY ALL REQUIRED SAMPLES EXIST
# ============================================================

missing_samples <- setdiff(
  selected_samples,
  colnames(expr_all)
)


if (length(missing_samples) > 0) {
  
  stop(
    paste(
      "Required samples missing:",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "\nAll 12 required D7 samples found.\n"
)


# ============================================================
# 6. SAMPLE METADATA
# ============================================================

sample_info <- data.frame(
  
  sample = selected_samples,
  
  model = c(
    rep("MDA-MB-468", 6),
    rep("SUM1315", 6)
  ),
  
  state = c(
    rep("Parental", 3),
    rep("OlaR", 3),
    rep("Parental", 3),
    rep("OlaR", 3)
  ),
  
  replicate = c(
    1:3,
    1:3,
    1:3,
    1:3
  ),
  
  treatment_at_collection = c(
    rep("No_PARPi", 3),
    rep("Olaparib_15uM", 3),
    rep("No_PARPi", 3),
    rep("Olaparib_25uM", 3)
  ),
  
  stringsAsFactors = FALSE
  
)


rownames(sample_info) <-
  sample_info$sample


write.csv(
  sample_info,
  file.path(
    results_dir,
    "D7_primary_sample_metadata.csv"
  ),
  row.names = FALSE
)


cat(
  "\nPrimary D7 sample design:\n\n"
)

print(
  sample_info,
  row.names = FALSE
)


# ============================================================
# 7. EXTRACT PRIMARY 12-SAMPLE MATRIX
# ============================================================

expr <- as.matrix(
  expr_all[
    ,
    selected_samples,
    drop = FALSE
  ]
)


storage.mode(expr) <- "numeric"


# ------------------------------------------------------------
# IMPORTANT FIX
#
# Give every matrix row a technical unique identifier.
#
# This does NOT biologically collapse or rename duplicate genes.
# It only gives R stable row labels for PCA/QC operations.
# ------------------------------------------------------------

rownames(expr) <- make.unique(
  gene_symbol,
  sep = "__duplicate_"
)


# ============================================================
# 8. BASIC MATRIX QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRIMARY MATRIX QC\n")
cat("============================================================\n\n")


cat(
  "Dimensions:",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n"
)


cat(
  "NA:",
  sum(
    is.na(expr)
  ),
  "\n"
)


cat(
  "Negative values:",
  sum(
    expr < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Exact zeros:",
  sum(
    expr == 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nSample column sums:\n"
)


print(
  colSums(
    expr
  )
)


# ============================================================
# 9. EXPRESSION FILTER FOR QC
#
# This is NOT yet the final DE filtering rule.
#
# For visualization:
# retain genes with expression >= 1 in >= 3 samples.
# ============================================================

keep_qc <- rowSums(
  expr >= 1
) >= 3


expr_qc <- expr[
  keep_qc,
  ,
  drop = FALSE
]


cat(
  "\nGenes retained for QC:",
  nrow(expr_qc),
  "/",
  nrow(expr),
  "\n"
)


# ============================================================
# 10. LOG2 TRANSFORMATION
# ============================================================

log_expr <- log2(
  expr_qc + 1
)


# ============================================================
# 11. REMOVE ZERO-VARIANCE GENES
# ============================================================

gene_var <- apply(
  log_expr,
  1,
  var
)


keep_var <- is.finite(gene_var) &
  gene_var > 0


log_expr <- log_expr[
  keep_var,
  ,
  drop = FALSE
]


cat(
  "Genes after variance check:",
  nrow(log_expr),
  "\n"
)


if (nrow(log_expr) == 0) {
  
  stop(
    "No genes remain after variance filtering."
  )
  
}


if (ncol(log_expr) != 12) {
  
  stop(
    "Unexpected number of samples after preprocessing."
  )
  
}


# ============================================================
# 12. SAMPLE CORRELATION MATRIX
# ============================================================

sample_cor <- cor(
  log_expr,
  method = "pearson"
)


write.csv(
  sample_cor,
  file.path(
    results_dir,
    "D7_sample_correlation_matrix.csv"
  )
)


pdf(
  file.path(
    results_dir,
    "D7_sample_correlation_heatmap.pdf"
  ),
  width = 9,
  height = 8
)


pheatmap::pheatmap(
  
  sample_cor,
  
  annotation_col =
    sample_info[
      colnames(sample_cor),
      c(
        "model",
        "state"
      ),
      drop = FALSE
    ],
  
  annotation_row =
    sample_info[
      rownames(sample_cor),
      c(
        "model",
        "state"
      ),
      drop = FALSE
    ],
  
  main =
    "GSE165914 — sample correlations"
  
)


dev.off()


# ============================================================
# 13. ROBUST PCA FUNCTION
#
# IMPORTANT FIX:
# PCA genes are selected using integer row indices.
# We no longer depend on names(vars).
# ============================================================

run_pca <- function(
    mat,
    info,
    title,
    filename
) {
  
  
  # ----------------------------------------------------------
  # Basic safety checks
  # ----------------------------------------------------------
  
  if (nrow(mat) < 2) {
    
    stop(
      "PCA matrix contains fewer than 2 genes."
    )
    
  }
  
  
  if (ncol(mat) < 3) {
    
    stop(
      "PCA matrix contains fewer than 3 samples."
    )
    
  }
  
  
  if (any(!is.finite(mat))) {
    
    stop(
      "Non-finite values detected in PCA matrix."
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Calculate gene variances
  # ----------------------------------------------------------
  
  vars <- apply(
    mat,
    1,
    var
  )
  
  
  valid <- is.finite(vars) &
    vars > 0
  
  
  vars <- vars[
    valid
  ]
  
  
  mat2 <- mat[
    valid,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(mat2) < 2) {
    
    stop(
      "Too few variable genes remain for PCA."
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Select top variable genes BY ROW INDEX
  #
  # This fixes the previous zero-dimension error.
  # ----------------------------------------------------------
  
  order_idx <- order(
    vars,
    decreasing = TRUE
  )
  
  
  n_top <- min(
    5000,
    length(order_idx)
  )
  
  
  top_idx <- order_idx[
    seq_len(n_top)
  ]
  
  
  pca_matrix <- mat2[
    top_idx,
    ,
    drop = FALSE
  ]
  
  
  cat(
    "\nPCA:",
    title,
    "\n"
  )
  
  
  cat(
    "Input genes:",
    nrow(mat),
    "\n"
  )
  
  
  cat(
    "Variable genes:",
    nrow(mat2),
    "\n"
  )
  
  
  cat(
    "Top genes used for PCA:",
    nrow(pca_matrix),
    "\n"
  )
  
  
  cat(
    "Samples used for PCA:",
    ncol(pca_matrix),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # PCA
  # ----------------------------------------------------------
  
  pca <- prcomp(
    t(pca_matrix),
    center = TRUE,
    scale. = FALSE
  )
  
  
  variance_percent <-
    100 *
    pca$sdev^2 /
    sum(
      pca$sdev^2
    )
  
  
  # ----------------------------------------------------------
  # Coordinates
  # ----------------------------------------------------------
  
  pca_df <- data.frame(
    
    sample =
      rownames(
        pca$x
      ),
    
    PC1 =
      pca$x[
        ,
        1
      ],
    
    PC2 =
      pca$x[
        ,
        2
      ],
    
    stringsAsFactors = FALSE
    
  )
  
  
  # Add metadata without merge/reordering problems
  
  meta_match <- match(
    pca_df$sample,
    info$sample
  )
  
  
  if (any(is.na(meta_match))) {
    
    stop(
      "Could not match all PCA samples to sample metadata."
    )
    
  }
  
  
  pca_df$model <-
    info$model[
      meta_match
    ]
  
  
  pca_df$state <-
    info$state[
      meta_match
    ]
  
  
  pca_df$replicate <-
    info$replicate[
      meta_match
    ]
  
  
  pca_df$treatment_at_collection <-
    info$treatment_at_collection[
      meta_match
    ]
  
  
  # ----------------------------------------------------------
  # Save coordinates
  # ----------------------------------------------------------
  
  coordinate_file <- sub(
    "\\.pdf$",
    "_coordinates.csv",
    filename
  )
  
  
  write.csv(
    pca_df,
    coordinate_file,
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      shape = state,
      label = sample
    )
  ) +
    
    geom_point(
      size = 4
    ) +
    
    ggrepel::geom_text_repel(
      size = 3,
      max.overlaps = Inf
    ) +
    
    labs(
      
      title = title,
      
      x = paste0(
        "PC1 (",
        round(
          variance_percent[1],
          1
        ),
        "%)"
      ),
      
      y = paste0(
        "PC2 (",
        round(
          variance_percent[2],
          1
        ),
        "%)"
      )
      
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    filename = filename,
    plot = p,
    width = 8,
    height = 6
  )
  
  
  # ----------------------------------------------------------
  # Console output
  # ----------------------------------------------------------
  
  cat(
    "\n",
    title,
    "\n",
    sep = ""
  )
  
  
  cat(
    "PC1 variance:",
    round(
      variance_percent[1],
      2
    ),
    "%\n"
  )
  
  
  cat(
    "PC2 variance:",
    round(
      variance_percent[2],
      2
    ),
    "%\n\n"
  )
  
  
  print(
    pca_df[
      ,
      c(
        "sample",
        "model",
        "state",
        "PC1",
        "PC2"
      )
    ],
    row.names = FALSE
  )
  
  
  invisible(
    
    list(
      
      pca =
        pca,
      
      coordinates =
        pca_df,
      
      variance =
        variance_percent,
      
      n_genes_used =
        nrow(pca_matrix)
      
    )
    
  )
  
}


# ============================================================
# 14. COMBINED 12-SAMPLE PCA
# ============================================================

combined_pca <- run_pca(
  
  mat = log_expr,
  
  info = sample_info,
  
  title =
    "GSE165914 — D7 primary samples",
  
  filename =
    file.path(
      results_dir,
      "D7_combined_PCA.pdf"
    )
  
)


# ============================================================
# 15. D7A — MDA-MB-468 PCA
# ============================================================

m468_samples <- c(
  m468_parental,
  m468_resistant
)


m468_pca <- run_pca(
  
  mat =
    log_expr[
      ,
      m468_samples,
      drop = FALSE
    ],
  
  info =
    sample_info[
      m468_samples,
      ,
      drop = FALSE
    ],
  
  title =
    "GSE165914 — MDA-MB-468 parental vs OlaR",
  
  filename =
    file.path(
      results_dir,
      "D7A_M468_PCA.pdf"
    )
  
)


# ============================================================
# 16. D7B — SUM1315 PCA
# ============================================================

sum_samples <- c(
  sum_parental,
  sum_resistant
)


sum_pca <- run_pca(
  
  mat =
    log_expr[
      ,
      sum_samples,
      drop = FALSE
    ],
  
  info =
    sample_info[
      sum_samples,
      ,
      drop = FALSE
    ],
  
  title =
    "GSE165914 — SUM1315 parental vs OlaR",
  
  filename =
    file.path(
      results_dir,
      "D7B_SUM1315_PCA.pdf"
    )
  
)


# ============================================================
# 17. WITHIN-GROUP REPLICATE CORRELATIONS
# ============================================================

group_cor_summary <- function(samples) {
  
  
  m <- sample_cor[
    samples,
    samples,
    drop = FALSE
  ]
  
  
  vals <- m[
    upper.tri(m)
  ]
  
  
  c(
    
    mean =
      mean(vals),
    
    min =
      min(vals),
    
    max =
      max(vals)
    
  )
  
}


replicate_summary <- rbind(
  
  M468_Parental =
    group_cor_summary(
      m468_parental
    ),
  
  M468_OlaR =
    group_cor_summary(
      m468_resistant
    ),
  
  SUM_Parental =
    group_cor_summary(
      sum_parental
    ),
  
  SUM_OlaR =
    group_cor_summary(
      sum_resistant
    )
  
)


cat("\n")
cat("============================================================\n")
cat("WITHIN-GROUP REPLICATE CORRELATIONS\n")
cat("============================================================\n\n")


print(
  replicate_summary
)


write.csv(
  replicate_summary,
  file.path(
    results_dir,
    "D7_within_group_replicate_correlations.csv"
  )
)


# ============================================================
# 18. SAVE PROCESSED QC MATRICES
# ============================================================

write.csv(
  
  data.frame(
    gene_id =
      rownames(log_expr),
    
    log_expr,
    
    check.names = FALSE
  ),
  
  file.path(
    results_dir,
    "D7_QC_log2_expression_matrix.csv"
  ),
  
  row.names = FALSE
  
)


# ============================================================
# 19. SAVE QC OBJECT
# ============================================================

saveRDS(
  
  list(
    
    original_gene_symbols =
      gene_symbol,
    
    technical_row_ids =
      rownames(expr),
    
    expression_processed =
      expr,
    
    expression_log2_filtered =
      log_expr,
    
    sample_info =
      sample_info,
    
    sample_correlations =
      sample_cor,
    
    replicate_correlations =
      replicate_summary,
    
    combined_pca =
      combined_pca,
    
    m468_pca =
      m468_pca,
    
    sum_pca =
      sum_pca
    
  ),
  
  file.path(
    results_dir,
    "D7_QC_objects.rds"
  )
  
)


# ============================================================
# 20. SESSION INFO
# ============================================================

sink(
  file.path(
    results_dir,
    "D7_QC_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 21. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D7 QC COMPLETE\n")
cat("============================================================\n\n")


cat(
  "NO differential-expression testing was performed.\n\n"
)


cat(
  "Report these outputs:\n\n",
  "1. Genes retained for QC\n",
  "2. Combined PCA PC1/PC2 variance + coordinates\n",
  "3. MDA-MB-468 PCA PC1/PC2 variance + coordinates\n",
  "4. SUM1315 PCA PC1/PC2 variance + coordinates\n",
  "5. Within-group replicate correlations\n"
)


cat("\n============================================================\n")