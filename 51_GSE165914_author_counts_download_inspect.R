# ============================================================
# GSE165914 - D7
# DOWNLOAD + INSPECT AUTHOR-SUPPLIED ORIGINAL COUNT MATRIX
#
# NO differential-expression analysis yet.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

raw_dir <- "data_raw/GSE165914"

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}

library(data.table)
library(R.utils)


# ============================================================
# 1. FILE INFORMATION
# ============================================================

filename <-
  "GSE165914_SUM_M468_Parental_OR_counts_mshekhar.csv.gz"

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE165nnn/GSE165914/suppl/",
  filename
)

gz_file <- file.path(
  raw_dir,
  filename
)

csv_file <- sub(
  "\\.gz$",
  "",
  gz_file
)


cat("\n")
cat("============================================================\n")
cat("GSE165914 AUTHOR COUNT MATRIX\n")
cat("============================================================\n\n")


# ============================================================
# 2. DOWNLOAD
# ============================================================

if (!file.exists(gz_file)) {
  
  cat(
    "Downloading author-supplied count matrix...\n\n"
  )
  
  download.file(
    url = url,
    destfile = gz_file,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  
} else {
  
  cat(
    "Compressed count file already exists.\n"
  )
  
}


cat(
  "Downloaded file size:",
  round(
    file.info(gz_file)$size / 1024^2,
    3
  ),
  "MB\n"
)


# ============================================================
# 3. DECOMPRESS
# ============================================================

if (!file.exists(csv_file)) {
  
  R.utils::gunzip(
    filename = gz_file,
    destname = csv_file,
    remove = FALSE,
    overwrite = TRUE
  )
  
}


# ============================================================
# 4. READ DATA
# ============================================================

dat <- data.table::fread(
  csv_file,
  data.table = FALSE,
  check.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("MATRIX STRUCTURE\n")
cat("============================================================\n\n")


cat(
  "Dimensions:",
  nrow(dat),
  "rows x",
  ncol(dat),
  "columns\n\n"
)


cat(
  "Column names:\n\n"
)

print(
  colnames(dat)
)


cat(
  "\nFirst 10 rows:\n\n"
)

print(
  head(dat, 10)
)


cat(
  "\nColumn classes:\n\n"
)

print(
  sapply(
    dat,
    class
  )
)


# ============================================================
# 5. IDENTIFIER QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("IDENTIFIER QC\n")
cat("============================================================\n\n")


cat(
  "First column:",
  colnames(dat)[1],
  "\n"
)


cat(
  "Unique identifiers:",
  length(
    unique(dat[[1]])
  ),
  "\n"
)


cat(
  "Duplicated identifiers:",
  sum(
    duplicated(dat[[1]])
  ),
  "\n"
)


# ============================================================
# 6. EXPRESSION COLUMN QC
#
# IMPORTANT:
# Explicitly EXCLUDE first identifier column.
# ============================================================

expr <- dat[
  ,
  -1,
  drop = FALSE
]


numeric_cols <- vapply(
  expr,
  is.numeric,
  logical(1)
)


cat("\n")
cat("============================================================\n")
cat("EXPRESSION MATRIX QC\n")
cat("============================================================\n\n")


cat(
  "Expression columns:",
  ncol(expr),
  "\n"
)


cat(
  "Numeric expression columns:",
  sum(numeric_cols),
  "\n"
)


if (!all(numeric_cols)) {
  
  cat(
    "\nNon-numeric expression columns:\n"
  )
  
  print(
    colnames(expr)[
      !numeric_cols
    ]
  )
  
}


x <- as.matrix(
  expr[
    ,
    numeric_cols,
    drop = FALSE
  ]
)


cat(
  "\nNA values:",
  sum(is.na(x)),
  "\n"
)


cat(
  "Negative values:",
  sum(
    x < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Exact zeros:",
  sum(
    x == 0,
    na.rm = TRUE
  ),
  "\n"
)


integer_like <- mean(
  abs(
    x - round(x)
  ) < 1e-8,
  na.rm = TRUE
)


cat(
  "Integer-like fraction:",
  round(
    integer_like,
    6
  ),
  "\n"
)


cat(
  "\nValue quantiles:\n\n"
)


print(
  
  quantile(
    as.numeric(x),
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


cat(
  "\nColumn sums:\n\n"
)


print(
  sort(
    colSums(
      x,
      na.rm = TRUE
    )
  )
)


# ============================================================
# 7. SEARCH COLUMN NAMES FOR KEY GROUPS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("KEY SAMPLE/GROUP SEARCH\n")
cat("============================================================\n\n")


patterns <- c(
  "MDA",
  "M468",
  "SUM",
  "Parent",
  "Control",
  "OlaR",
  "OR",
  "olaparib"
)


for (pattern in patterns) {
  
  matches <- grep(
    pattern,
    colnames(dat),
    value = TRUE,
    ignore.case = TRUE
  )
  
  
  cat(
    "\nPattern:",
    pattern,
    "\n"
  )
  
  print(
    matches
  )
  
}


# ============================================================
# 8. SAVE INSPECTION OBJECT
# ============================================================

saveRDS(
  dat,
  file.path(
    raw_dir,
    "GSE165914_author_counts_raw.rds"
  )
)


# ============================================================
# 9. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE165914 AUTHOR COUNT INSPECTION COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No filtering, normalization, DE testing, or GSEA was performed.\n"
)

cat(
  "Next step = identify exact D7A and D7B parental/OlaR columns ",
  "from the author matrix and perform replicate QC.\n"
)

cat("============================================================\n")