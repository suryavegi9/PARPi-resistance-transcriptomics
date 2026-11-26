# ============================================================
# GSE191231 - DATASET 6
# MATRIX / WORKBOOK AUDIT
#
# Goal:
# Determine workbook structure, sample columns,
# whether raw integer counts are available,
# TPM structure, identifiers and duplicate genes.
#
# NO differential expression yet.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(readxl)


file_path <- paste0(
  "data_raw/",
  "GSE191231_Hyun_Analysis_tpm_edger.xlsx"
)

stopifnot(
  file.exists(file_path)
)


# ============================================================
# 1. WORKBOOK SHEETS
# ============================================================

sheets <- excel_sheets(
  file_path
)

cat("\n========================================\n")
cat("WORKBOOK SHEETS\n")
cat("========================================\n\n")

print(sheets)


# ============================================================
# 2. AUDIT EVERY SHEET
# ============================================================

for (sheet_name in sheets) {
  
  cat("\n\n")
  cat("########################################\n")
  cat("SHEET:", sheet_name, "\n")
  cat("########################################\n\n")
  
  
  dat <- read_excel(
    file_path,
    sheet = sheet_name
  )
  
  dat <- as.data.frame(
    dat,
    check.names = FALSE
  )
  
  
  cat(
    "Dimensions:",
    nrow(dat),
    "rows x",
    ncol(dat),
    "columns\n\n"
  )
  
  
  cat("Column names:\n")
  
  print(
    colnames(dat)
  )
  
  
  cat("\nFIRST 6 ROWS:\n")
  
  print(
    head(dat)
  )
  
  
  # ----------------------------------------------------------
  # Column types
  # ----------------------------------------------------------
  
  classes <- sapply(
    dat,
    class
  )
  
  
  cat("\nCOLUMN CLASSES:\n")
  
  print(classes)
  
  
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
  
  
  cat(
    "\nNumeric columns:",
    length(numeric_cols),
    "\n"
  )
  
  cat(
    "Non-numeric columns:",
    length(non_numeric_cols),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Numeric QC
  # ----------------------------------------------------------
  
  if (length(numeric_cols) > 0) {
    
    num <- as.matrix(
      dat[
        ,
        numeric_cols,
        drop = FALSE
      ]
    )
    
    
    integer_fraction <- mean(
      abs(
        num - round(num)
      ) < 1e-8,
      na.rm = TRUE
    )
    
    
    cat("\nNUMERIC QC:\n")
    
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
    
    cat(
      "Integer-like fraction:",
      round(integer_fraction, 6),
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
          .50,
          .75,
          .95,
          .99,
          1
        ),
        na.rm = TRUE
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Identifier QC
  # ----------------------------------------------------------
  
  if (length(non_numeric_cols) > 0) {
    
    cat("\nIDENTIFIER AUDIT:\n")
    
    for (j in non_numeric_cols) {
      
      x <- dat[[j]]
      
      cat(
        "\nColumn:",
        colnames(dat)[j],
        "\n"
      )
      
      cat(
        "Unique:",
        length(unique(x)),
        "\n"
      )
      
      cat(
        "Duplicates:",
        sum(duplicated(x)),
        "\n"
      )
      
      cat(
        "Missing:",
        sum(is.na(x) | x == ""),
        "\n"
      )
    }
  }
}


cat("\n========================================\n")
cat("GSE191231 WORKBOOK AUDIT COMPLETE\n")
cat("========================================\n")