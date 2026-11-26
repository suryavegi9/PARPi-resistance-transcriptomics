# ============================================================
# GSE153867 - DATASET 3 INSPECTION
# Publication-standard workflow
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

dir.create(
  "results/GSE153867",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. LOCATE FILE
# ============================================================

candidate_files <- c(
  "data_raw/GSE153867_fpkm.txt.gz",
  "data_raw/GSE153867_fpkm.txt"
)

existing_files <- candidate_files[
  file.exists(candidate_files)
]

if (length(existing_files) == 0) {
  stop(
    "GSE153867 FPKM file not found in data_raw/GSE153867/"
  )
}

file_path <- existing_files[1]

cat("\nUsing file:\n")
cat(file_path, "\n\n")


# ============================================================
# 2. READ MATRIX
# ============================================================

if (grepl("\\.gz$", file_path)) {
  
  dat <- read.delim(
    gzfile(file_path),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  dat <- read.delim(
    file_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 3. BASIC STRUCTURE
# ============================================================

cat("========================================\n")
cat("GSE153867 MATRIX INSPECTION\n")
cat("========================================\n\n")

cat("Dimensions:\n")
cat(
  nrow(dat),
  "rows x",
  ncol(dat),
  "columns\n\n"
)

cat("Column names:\n")
print(colnames(dat))

cat("\nFirst 6 rows:\n")
print(head(dat))

cat("\nColumn classes:\n")
print(
  sapply(dat, class)
)

cat("\nMissing values per column:\n")
print(
  colSums(is.na(dat))
)


# ============================================================
# 4. IDENTIFIER COLUMN INSPECTION
# ============================================================

cat("\n========================================\n")
cat("FIRST COLUMN INSPECTION\n")
cat("========================================\n")

cat("\nFirst column name:\n")
print(
  colnames(dat)[1]
)

cat("\nFirst 20 values:\n")
print(
  head(dat[[1]], 20)
)

cat("\nUnique identifiers:\n")
print(
  length(unique(dat[[1]]))
)

cat("\nDuplicated identifiers:\n")
print(
  sum(duplicated(dat[[1]]))
)


# ============================================================
# 5. NUMERIC EXPRESSION COLUMNS
# ============================================================

numeric_columns <- sapply(
  dat,
  is.numeric
)

expr_numeric <- dat[
  ,
  numeric_columns,
  drop = FALSE
]

cat("\n========================================\n")
cat("NUMERIC EXPRESSION INSPECTION\n")
cat("========================================\n")

cat("\nNumeric columns:\n")
print(
  colnames(expr_numeric)
)

cat("\nNumber of numeric columns:\n")
print(
  ncol(expr_numeric)
)

cat("\nExpression summary:\n")
print(
  summary(
    as.vector(
      as.matrix(expr_numeric)
    )
  )
)

cat("\nMinimum:\n")
print(
  min(
    as.matrix(expr_numeric),
    na.rm = TRUE
  )
)

cat("\nMaximum:\n")
print(
  max(
    as.matrix(expr_numeric),
    na.rm = TRUE
  )
)

cat("\nNegative values:\n")
print(
  sum(
    as.matrix(expr_numeric) < 0,
    na.rm = TRUE
  )
)

cat("\nExact zeros:\n")
print(
  sum(
    as.matrix(expr_numeric) == 0,
    na.rm = TRUE
  )
)


# ============================================================
# 6. INTEGER-LIKE CHECK
# ============================================================

vals <- as.vector(
  as.matrix(expr_numeric)
)

vals <- vals[
  is.finite(vals)
]

integer_like <- abs(
  vals - round(vals)
) < 1e-8

cat("\nProportion integer-like:\n")
print(
  mean(integer_like)
)


# ============================================================
# 7. PER-SAMPLE SUMMARIES
# ============================================================

sample_summary <- data.frame(
  
  sample = colnames(expr_numeric),
  
  minimum = sapply(
    expr_numeric,
    min,
    na.rm = TRUE
  ),
  
  median = sapply(
    expr_numeric,
    median,
    na.rm = TRUE
  ),
  
  mean = sapply(
    expr_numeric,
    mean,
    na.rm = TRUE
  ),
  
  maximum = sapply(
    expr_numeric,
    max,
    na.rm = TRUE
  ),
  
  zeros = sapply(
    expr_numeric,
    function(x)
      sum(x == 0, na.rm = TRUE)
  ),
  
  missing = sapply(
    expr_numeric,
    function(x)
      sum(is.na(x))
  )
)

cat("\n========================================\n")
cat("PER-SAMPLE SUMMARY\n")
cat("========================================\n\n")

print(sample_summary)


# ============================================================
# 8. SIMPLE SAMPLE CORRELATION
# ============================================================

if (ncol(expr_numeric) > 1) {
  
  log_expr_preview <- log2(
    as.matrix(expr_numeric) + 1
  )
  
  sample_cor <- cor(
    log_expr_preview,
    use = "pairwise.complete.obs",
    method = "spearman"
  )
  
  cat("\n========================================\n")
  cat("SAMPLE SPEARMAN CORRELATION\n")
  cat("========================================\n\n")
  
  print(
    round(sample_cor, 3)
  )
  
  write.csv(
    sample_cor,
    "results/GSE153867/sample_correlation_preview.csv"
  )
}


# ============================================================
# 9. SAVE REPORT
# ============================================================

write.csv(
  sample_summary,
  "results/GSE153867/sample_expression_summary.csv",
  row.names = FALSE
)

inspection_report <- capture.output({
  
  cat("GSE153867 INSPECTION REPORT\n\n")
  
  cat(
    "Dimensions:",
    nrow(dat),
    "rows x",
    ncol(dat),
    "columns\n\n"
  )
  
  cat("Column names:\n")
  print(colnames(dat))
  
  cat("\nColumn classes:\n")
  print(sapply(dat, class))
  
  cat("\nIdentifier column:\n")
  print(colnames(dat)[1])
  
  cat("\nUnique identifiers:\n")
  print(length(unique(dat[[1]])))
  
  cat("\nDuplicated identifiers:\n")
  print(sum(duplicated(dat[[1]])))
  
  cat("\nNumeric columns:\n")
  print(colnames(expr_numeric))
  
  cat("\nExpression range:\n")
  print(
    range(
      as.matrix(expr_numeric),
      na.rm = TRUE
    )
  )
  
  cat("\nExact zeros:\n")
  print(
    sum(
      as.matrix(expr_numeric) == 0,
      na.rm = TRUE
    )
  )
  
  cat("\nMissing values:\n")
  print(
    sum(
      is.na(
        as.matrix(expr_numeric)
      )
    )
  )
  
  cat("\nInteger-like proportion:\n")
  print(
    mean(integer_like)
  )
  
  cat("\nPer-sample summary:\n")
  print(sample_summary)
})

writeLines(
  inspection_report,
  "results/GSE153867/GSE153867_inspection_report.txt"
)


cat("\n========================================\n")
cat("GSE153867 INSPECTION COMPLETE\n")
cat("========================================\n")