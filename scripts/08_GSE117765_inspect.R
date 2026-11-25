# ============================================================
# GSE117765 - DATA INSPECTION
# Publication workflow: inspect before choosing statistics
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

dir.create(
  "results/GSE117765",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Find downloaded matrix
# ------------------------------------------------------------

candidate_files <- c(
  "data_raw/GSE117765/GSE117765_matrix.txt.gz",
  "data_raw/GSE117765/GSE117765_matrix.txt"
)

existing_files <- candidate_files[file.exists(candidate_files)]

if (length(existing_files) == 0) {
  stop(
    "GSE117765 matrix file not found in data_raw/GSE117765/"
  )
}

file_path <- existing_files[1]

cat("\nUsing file:\n")
cat(file_path, "\n\n")


# ------------------------------------------------------------
# 2. Read matrix
# ------------------------------------------------------------

if (grepl("\\.gz$", file_path)) {
  
  expression_raw <- read.delim(
    gzfile(file_path),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  expression_raw <- read.delim(
    file_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 3. Basic inspection
# ------------------------------------------------------------

cat("========================================\n")
cat("GSE117765 MATRIX INSPECTION\n")
cat("========================================\n\n")

cat("Dimensions:\n")
cat(
  nrow(expression_raw),
  "rows x",
  ncol(expression_raw),
  "columns\n\n"
)

cat("Column names:\n")
print(colnames(expression_raw))

cat("\nFirst 6 rows:\n")
print(head(expression_raw))

cat("\nColumn classes:\n")
print(sapply(expression_raw, class))

cat("\nMissing values per column:\n")
print(colSums(is.na(expression_raw)))


# ------------------------------------------------------------
# 4. Inspect first column
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FIRST COLUMN INSPECTION\n")
cat("========================================\n")

cat("\nFirst column name:\n")
print(colnames(expression_raw)[1])

cat("\nFirst 20 values:\n")
print(head(expression_raw[[1]], 20))

cat("\nNumber of unique first-column values:\n")
print(length(unique(expression_raw[[1]])))

cat("\nDuplicated first-column values:\n")
print(sum(duplicated(expression_raw[[1]])))


# ------------------------------------------------------------
# 5. Inspect numeric expression columns
# ------------------------------------------------------------

numeric_columns <- sapply(
  expression_raw,
  is.numeric
)

cat("\n========================================\n")
cat("NUMERIC COLUMN INSPECTION\n")
cat("========================================\n")

cat("\nNumeric columns detected:\n")
print(names(expression_raw)[numeric_columns])

expr_numeric <- expression_raw[, numeric_columns, drop = FALSE]

cat("\nExpression summary:\n")
print(summary(as.vector(as.matrix(expr_numeric))))

cat("\nMinimum value:\n")
print(min(as.matrix(expr_numeric), na.rm = TRUE))

cat("\nMaximum value:\n")
print(max(as.matrix(expr_numeric), na.rm = TRUE))

cat("\nNumber of negative values:\n")
print(
  sum(
    as.matrix(expr_numeric) < 0,
    na.rm = TRUE
  )
)

cat("\nNumber of exact zero values:\n")
print(
  sum(
    as.matrix(expr_numeric) == 0,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 6. Are values integer-like?
# ------------------------------------------------------------

all_values <- as.vector(
  as.matrix(expr_numeric)
)

all_values <- all_values[
  is.finite(all_values)
]

integer_like <- abs(
  all_values - round(all_values)
) < 1e-8

cat("\nProportion of values that are integer-like:\n")
print(mean(integer_like))

cat("\nExample expression values:\n")
print(head(all_values, 30))


# ------------------------------------------------------------
# 7. Per-sample summaries
# ------------------------------------------------------------

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
    function(x) sum(x == 0, na.rm = TRUE)
  )
)

cat("\n========================================\n")
cat("PER-SAMPLE SUMMARY\n")
cat("========================================\n\n")

print(sample_summary)


# ------------------------------------------------------------
# 8. Save inspection information
# ------------------------------------------------------------

write.csv(
  sample_summary,
  "results/GSE117765/sample_expression_summary.csv",
  row.names = FALSE
)

inspection_lines <- capture.output({
  
  cat("GSE117765 DATA INSPECTION REPORT\n\n")
  
  cat(
    "Dimensions:",
    nrow(expression_raw),
    "rows x",
    ncol(expression_raw),
    "columns\n\n"
  )
  
  cat("Column names:\n")
  print(colnames(expression_raw))
  
  cat("\nColumn classes:\n")
  print(sapply(expression_raw, class))
  
  cat("\nFirst column:\n")
  print(colnames(expression_raw)[1])
  
  cat("\nDuplicated identifiers:\n")
  print(sum(duplicated(expression_raw[[1]])))
  
  cat("\nNumeric expression columns:\n")
  print(colnames(expr_numeric))
  
  cat("\nExpression range:\n")
  print(range(as.matrix(expr_numeric), na.rm = TRUE))
  
  cat("\nInteger-like proportion:\n")
  print(mean(integer_like))
  
  cat("\nNegative values:\n")
  print(sum(as.matrix(expr_numeric) < 0, na.rm = TRUE))
  
  cat("\nExact zero values:\n")
  print(sum(as.matrix(expr_numeric) == 0, na.rm = TRUE))
  
  cat("\nPer-sample summary:\n")
  print(sample_summary)
})

writeLines(
  inspection_lines,
  "results/GSE117765/GSE117765_inspection_report.txt"
)


cat("\n========================================\n")
cat("INSPECTION COMPLETE\n")
cat("========================================\n")

cat(
  "\nSaved report:\n",
  "results/GSE117765/GSE117765_inspection_report.txt\n"
)

cat(
  "\nSaved sample summary:\n",
  "results/GSE117765/sample_expression_summary.csv\n"
)