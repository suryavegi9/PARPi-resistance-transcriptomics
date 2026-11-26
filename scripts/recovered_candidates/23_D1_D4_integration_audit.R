# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 1: RESULT-FILE AUDIT
#
# Goal:
# Find all Hallmark/GSEA result tables already produced
# for Datasets 1-4 and inspect their structure before merging.
#
# NO cross-dataset conclusions are generated here.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ------------------------------------------------------------
# 1. Find candidate CSV files
# ------------------------------------------------------------

all_csv <- list.files(
  path = ".",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

candidate_files <- all_csv[
  grepl(
    "hallmark|gsea|frozen",
    all_csv,
    ignore.case = TRUE
  )
]


cat("\n========================================\n")
cat("CANDIDATE HALLMARK / GSEA FILES\n")
cat("========================================\n\n")

for (i in seq_along(candidate_files)) {
  
  cat(
    sprintf(
      "[%02d] %s\n",
      i,
      candidate_files[i]
    )
  )
}


# ------------------------------------------------------------
# 2. Inspect each candidate file
# ------------------------------------------------------------

audit_list <- list()

for (i in seq_along(candidate_files)) {
  
  f <- candidate_files[i]
  
  cat("\n\n========================================\n")
  cat("FILE", i, "\n")
  cat("========================================\n")
  
  cat(
    "Path:",
    f,
    "\n\n"
  )
  
  x <- tryCatch(
    
    read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(x)) {
    
    cat("Could not read file.\n")
    
    next
  }
  
  
  cat(
    "Dimensions:",
    nrow(x),
    "rows x",
    ncol(x),
    "columns\n"
  )
  
  
  cat("\nColumn names:\n")
  
  print(
    colnames(x)
  )
  
  
  cat("\nFirst 3 rows:\n")
  
  print(
    head(
      x,
      3
    )
  )
  
  
  # ----------------------------------------------------------
  # Try to detect useful integration columns
  # ----------------------------------------------------------
  
  pathway_col <- intersect(
    c(
      "ID",
      "pathway",
      "Pathway",
      "Description"
    ),
    colnames(x)
  )
  
  nes_col <- grep(
    "^NES$|NES_",
    colnames(x),
    value = TRUE
  )
  
  padj_col <- grep(
    "p.adjust|padj|adj.P",
    colnames(x),
    value = TRUE,
    ignore.case = TRUE
  )
  
  
  audit_list[[length(audit_list) + 1]] <- data.frame(
    
    file = f,
    
    rows = nrow(x),
    
    columns = ncol(x),
    
    pathway_column = ifelse(
      length(pathway_col) > 0,
      paste(pathway_col, collapse = "; "),
      NA
    ),
    
    NES_columns = ifelse(
      length(nes_col) > 0,
      paste(nes_col, collapse = "; "),
      NA
    ),
    
    adjusted_p_columns = ifelse(
      length(padj_col) > 0,
      paste(padj_col, collapse = "; "),
      NA
    ),
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 3. Save audit summary
# ------------------------------------------------------------

audit_summary <- do.call(
  rbind,
  audit_list
)

dir.create(
  "results/integration_D1_D4",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  audit_summary,
  "results/integration_D1_D4/Hallmark_file_audit.csv",
  row.names = FALSE
)


cat("\n\n========================================\n")
cat("AUDIT SUMMARY\n")
cat("========================================\n\n")

print(
  audit_summary,
  row.names = FALSE
)


cat("\n========================================\n")
cat("INTEGRATION AUDIT COMPLETE\n")
cat("========================================\n")