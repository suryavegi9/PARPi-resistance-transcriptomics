# 30_GSE202434_HCC1806_audit.R

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(readxl)

file_path <- paste0(
  "data_raw/",
  "GSE202434_HCC1806_Talazoparib_TPM_ExpressedGenesONLY_ProteinCoding_ONLY.xlsx"
)

stopifnot(file.exists(file_path))

cat("\nSHEETS\n")
print(excel_sheets(file_path))

dat <- read_excel(
  file_path,
  sheet = 1
)

dat <- as.data.frame(
  dat,
  check.names = FALSE
)

cat("\n========================================\n")
cat("HCC1806 MATRIX AUDIT\n")
cat("========================================\n\n")

cat(
  "Dimensions:",
  nrow(dat),
  "rows x",
  ncol(dat),
  "columns\n\n"
)

cat("Column names:\n")
print(colnames(dat))

cat("\nFIRST 6 ROWS:\n")
print(
  head(dat)
)

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

num <- as.matrix(
  dat[
    ,
    numeric_cols,
    drop = FALSE
  ]
)

cat("\n========================================\n")
cat("NUMERIC QC\n")
cat("========================================\n\n")

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

integer_fraction <- mean(
  abs(
    num - round(num)
  ) < 1e-8,
  na.rm = TRUE
)

cat(
  "Integer-like fraction:",
  round(
    integer_fraction,
    6
  ),
  "\n"
)

cat("\nVALUE DISTRIBUTION:\n")

print(
  quantile(
    num,
    probs = c(
      0,
      .01,
      .05,
      .25,
      .5,
      .75,
      .95,
      .99,
      1
    ),
    na.rm = TRUE
  )
)

if (
  length(non_numeric_cols) > 0
) {
  
  cat("\nIDENTIFIER COLUMNS:\n")
  
  for (j in non_numeric_cols) {
    
    cat(
      colnames(dat)[j],
      ": unique =",
      length(unique(dat[[j]])),
      "; duplicates =",
      sum(duplicated(dat[[j]])),
      "\n"
    )
  }
}

cat("\n========================================\n")
cat("HCC1806 AUDIT COMPLETE\n")
cat("========================================\n")