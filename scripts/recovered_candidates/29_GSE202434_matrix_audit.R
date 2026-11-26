# ============================================================
# GSE202434 - DATASET 5
# STEP 1: PROCESSED MATRIX AUDIT
#
# Goal:
# Determine exactly what the downloaded "normalized counts"
# contain BEFORE choosing a differential-expression model.
#
# No DE analysis is performed here.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 1. FILE PATH
# ============================================================
file_path <- file.path(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()), "data_raw/GSE202434_Exp43_Tala_Normalized_Counts_ExpressedGenesONLY_ProteinCoding_ONLY.csv.gz")


if (!file.exists(file_path)) {
  
  stop(
    paste0(
      "\nFile not found:\n",
      file_path,
      "\n\nCheck the filename inside data/GSE202434/."
    )
  )
}


# ============================================================
# 2. READ MATRIX
# ============================================================

dat <- read.csv(
  gzfile(file_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("\n========================================\n")
cat("GSE202434 MATRIX AUDIT\n")
cat("========================================\n\n")


# ============================================================
# 3. BASIC DIMENSIONS
# ============================================================

cat("Rows:", nrow(dat), "\n")
cat("Columns:", ncol(dat), "\n\n")


cat("Column names:\n\n")

print(
  colnames(dat)
)


# ============================================================
# 4. FIRST FEW ROWS
# ============================================================

cat("\n========================================\n")
cat("FIRST 6 ROWS / FIRST 10 COLUMNS\n")
cat("========================================\n\n")

print(
  dat[
    1:min(6, nrow(dat)),
    1:min(10, ncol(dat))
  ]
)


# ============================================================
# 5. COLUMN CLASSES
# ============================================================

classes <- sapply(
  dat,
  class
)


class_table <- data.frame(
  column = names(classes),
  class = as.character(classes),
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("COLUMN CLASSES\n")
cat("========================================\n\n")

print(
  class_table,
  row.names = FALSE
)


# ============================================================
# 6. IDENTIFY NUMERIC COLUMNS
# ============================================================

numeric_cols <- which(
  sapply(
    dat,
    is.numeric
  )
)


non_numeric_cols <- which(
  !sapply(
    dat,
    is.numeric
  )
)


cat("\nNumeric columns:", length(numeric_cols), "\n")

cat(
  "Non-numeric columns:",
  length(non_numeric_cols),
  "\n"
)


if (
  length(non_numeric_cols) > 0
) {
  
  cat("\nNon-numeric column names:\n")
  
  print(
    colnames(dat)[
      non_numeric_cols
    ]
  )
}


# ============================================================
# 7. NUMERIC MATRIX QC
# ============================================================

num <- as.matrix(
  dat[
    ,
    numeric_cols,
    drop = FALSE
  ]
)


cat("\n========================================\n")
cat("NUMERIC MATRIX QC\n")
cat("========================================\n\n")


cat(
  "Total numeric values:",
  length(num),
  "\n"
)


cat(
  "NA values:",
  sum(is.na(num)),
  "\n"
)


cat(
  "Zero values:",
  sum(num == 0, na.rm = TRUE),
  "\n"
)


cat(
  "Negative values:",
  sum(num < 0, na.rm = TRUE),
  "\n"
)


# ============================================================
# 8. INTEGER VS CONTINUOUS CHECK
#
# Raw RNA-seq counts should normally be integer-valued.
#
# Normalized counts often contain decimals.
# ============================================================

is_integer_like <-
  abs(
    num - round(num)
  ) < 1e-8


integer_fraction <- mean(
  is_integer_like,
  na.rm = TRUE
)


cat(
  "\nFraction of numeric values that are integer-like:",
  round(
    integer_fraction,
    6
  ),
  "\n"
)


cat(
  "Fraction containing decimals:",
  round(
    1 - integer_fraction,
    6
  ),
  "\n"
)


# ============================================================
# 9. GLOBAL VALUE DISTRIBUTION
# ============================================================

all_values <- as.vector(
  num
)


cat("\n========================================\n")
cat("GLOBAL VALUE DISTRIBUTION\n")
cat("========================================\n\n")


print(
  quantile(
    all_values,
    probs = c(
      0,
      0.01,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


# ============================================================
# 10. PER-COLUMN SUMMARY
#
# Useful for identifying sample columns, depth differences,
# and suspicious annotation-like numeric columns.
# ============================================================

column_summary <- data.frame(
  
  column =
    colnames(dat)[numeric_cols],
  
  min =
    apply(
      num,
      2,
      min,
      na.rm = TRUE
    ),
  
  median =
    apply(
      num,
      2,
      median,
      na.rm = TRUE
    ),
  
  mean =
    colMeans(
      num,
      na.rm = TRUE
    ),
  
  max =
    apply(
      num,
      2,
      max,
      na.rm = TRUE
    ),
  
  sum =
    colSums(
      num,
      na.rm = TRUE
    ),
  
  zero_fraction =
    colMeans(
      num == 0,
      na.rm = TRUE
    ),
  
  integer_fraction =
    apply(
      is_integer_like,
      2,
      mean,
      na.rm = TRUE
    ),
  
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("PER-COLUMN SUMMARY\n")
cat("========================================\n\n")

print(
  column_summary,
  row.names = FALSE
)


# ============================================================
# 11. CHECK FOR DUPLICATED IDENTIFIERS
#
# We do not yet know which annotation column is the gene ID,
# so inspect every non-numeric column.
# ============================================================

if (
  length(non_numeric_cols) > 0
) {
  
  cat("\n========================================\n")
  cat("NON-NUMERIC IDENTIFIER AUDIT\n")
  cat("========================================\n\n")
  
  
  for (
    j in non_numeric_cols
  ) {
    
    x <- dat[[j]]
    
    cat(
      "\nColumn:",
      colnames(dat)[j],
      "\n"
    )
    
    cat(
      "Unique values:",
      length(unique(x)),
      "\n"
    )
    
    cat(
      "Duplicated values:",
      sum(duplicated(x)),
      "\n"
    )
    
    cat(
      "Missing values:",
      sum(is.na(x) | x == ""),
      "\n"
    )
  }
}


# ============================================================
# 12. LOOK FOR EXPECTED D5 SAMPLE LABELS
# ============================================================

expected_terms <- c(
  "PSN1",
  "HCC1806",
  "Tala",
  "Parental",
  "DF",
  "Maint",
  "DMSO"
)


cat("\n========================================\n")
cat("SAMPLE-NAME KEYWORD SEARCH\n")
cat("========================================\n\n")


for (
  term in expected_terms
) {
  
  hits <- grep(
    term,
    colnames(dat),
    ignore.case = TRUE,
    value = TRUE
  )
  
  cat(
    term,
    ":",
    length(hits),
    "matches\n"
  )
  
  if (
    length(hits) > 0
  ) {
    
    print(
      hits
    )
  }
}


# ============================================================
# 13. SAVE AUDIT TABLES
# ============================================================

outdir <- "results/GSE202434/audit"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  class_table,
  file.path(
    outdir,
    "GSE202434_column_classes.csv"
  ),
  row.names = FALSE
)


write.csv(
  column_summary,
  file.path(
    outdir,
    "GSE202434_numeric_column_summary.csv"
  ),
  row.names = FALSE
)


writeLines(
  colnames(dat),
  file.path(
    outdir,
    "GSE202434_column_names.txt"
  )
)


# ============================================================
# 14. FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("AUDIT SUMMARY\n")
cat("========================================\n\n")

cat(
  "Dimensions:",
  nrow(dat),
  "rows x",
  ncol(dat),
  "columns\n"
)


cat(
  "Numeric columns:",
  length(numeric_cols),
  "\n"
)


cat(
  "Non-numeric columns:",
  length(non_numeric_cols),
  "\n"
)


cat(
  "NA numeric values:",
  sum(is.na(num)),
  "\n"
)


cat(
  "Zero numeric values:",
  sum(num == 0, na.rm = TRUE),
  "\n"
)


cat(
  "Integer-like fraction:",
  round(
    integer_fraction,
    6
  ),
  "\n"
)


cat("\n========================================\n")
cat("GSE202434 MATRIX AUDIT COMPLETE\n")
cat("========================================\n")
