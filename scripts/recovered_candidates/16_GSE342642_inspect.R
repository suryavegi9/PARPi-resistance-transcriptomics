# ============================================================
# GSE342642 - DATASET 4 INSPECTION
# No biological analysis yet
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

dir.create(
  "results/GSE342642",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Select downloaded GEO processed file manually
# ------------------------------------------------------------

file_path <- file.choose()

cat("\nSelected file:\n")
cat(file_path, "\n\n")


# ------------------------------------------------------------
# 2. Read file
# ------------------------------------------------------------

if (grepl("\\.csv\\.gz$|\\.csv$", file_path, ignore.case = TRUE)) {
  
  dat <- read.csv(
    file_path,
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


# ------------------------------------------------------------
# 3. Basic structure
# ------------------------------------------------------------

cat("========================================\n")
cat("GSE342642 MATRIX INSPECTION\n")
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
print(sapply(dat, class))

cat("\nMissing values per column:\n")
print(colSums(is.na(dat)))


# ------------------------------------------------------------
# 4. First column / gene ID inspection
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FIRST COLUMN INSPECTION\n")
cat("========================================\n")

cat("\nFirst column name:\n")
print(colnames(dat)[1])

cat("\nFirst 20 identifiers:\n")
print(head(dat[[1]], 20))

cat("\nUnique identifiers:\n")
print(length(unique(dat[[1]])))

cat("\nDuplicated identifiers:\n")
print(sum(duplicated(dat[[1]])))


# ------------------------------------------------------------
# 5. Numeric expression columns
# ------------------------------------------------------------

numeric_columns <- sapply(
  dat,
  is.numeric
)

expr <- dat[
  ,
  numeric_columns,
  drop = FALSE
]

cat("\n========================================\n")
cat("NUMERIC DATA INSPECTION\n")
cat("========================================\n")

cat("\nNumeric columns:\n")
print(colnames(expr))

cat("\nNumber of numeric columns:\n")
print(ncol(expr))

vals <- as.vector(
  as.matrix(expr)
)

vals <- vals[
  is.finite(vals)
]

cat("\nSummary of numeric values:\n")
print(summary(vals))

cat("\nMinimum:\n")
print(min(vals))

cat("\nMaximum:\n")
print(max(vals))

cat("\nNegative values:\n")
print(sum(vals < 0))

cat("\nExact zeros:\n")
print(sum(vals == 0))

integer_like <- abs(
  vals - round(vals)
) < 1e-8

cat("\nProportion integer-like:\n")
print(mean(integer_like))


# ------------------------------------------------------------
# 6. Per-sample summary
# ------------------------------------------------------------

if (ncol(expr) > 0) {
  
  sample_summary <- data.frame(
    
    sample = colnames(expr),
    
    minimum = sapply(
      expr,
      min,
      na.rm = TRUE
    ),
    
    median = sapply(
      expr,
      median,
      na.rm = TRUE
    ),
    
    mean = sapply(
      expr,
      mean,
      na.rm = TRUE
    ),
    
    maximum = sapply(
      expr,
      max,
      na.rm = TRUE
    ),
    
    zeros = sapply(
      expr,
      function(x)
        sum(x == 0, na.rm = TRUE)
    ),
    
    missing = sapply(
      expr,
      function(x)
        sum(is.na(x))
    )
  )
  
  cat("\n========================================\n")
  cat("PER-SAMPLE SUMMARY\n")
  cat("========================================\n\n")
  
  print(sample_summary)
  
  write.csv(
    sample_summary,
    "results/GSE342642/sample_expression_summary.csv",
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 7. Save inspection report
# ------------------------------------------------------------

report <- capture.output({
  
  cat("GSE342642 INSPECTION REPORT\n\n")
  
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
  
  cat("\nFirst column:\n")
  print(colnames(dat)[1])
  
  cat("\nDuplicated identifiers:\n")
  print(sum(duplicated(dat[[1]])))
  
  cat("\nNumeric columns:\n")
  print(colnames(expr))
  
  cat("\nNumeric value range:\n")
  print(range(vals))
  
  cat("\nExact zeros:\n")
  print(sum(vals == 0))
  
  cat("\nInteger-like proportion:\n")
  print(mean(integer_like))
})

writeLines(
  report,
  "results/GSE342642/GSE342642_inspection_report.txt"
)


cat("\n========================================\n")
cat("GSE342642 INSPECTION COMPLETE\n")
cat("========================================\n")