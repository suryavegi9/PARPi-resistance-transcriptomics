# ============================================================
# GSE163854 - D8
# DOWNLOAD + INSPECT RAW COUNTS AND METADATA
#
# HGSOC PH039 PDX
# Serial niraparib-resistance experiment
#
# THIS SCRIPT ONLY:
#   1. Downloads author-supplied raw counts
#   2. Downloads GEO SOFT metadata
#   3. Inspects gene identifiers and count values
#   4. Prints sample structure and treatment information
#
# NO FILTERING
# NO NORMALIZATION
# NO PCA
# NO DIFFERENTIAL EXPRESSION
# NO GSEA
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

accession <- "GSE163854"

raw_dir <- file.path(
  "data_raw",
  accession
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("GSE163854 - D8 DOWNLOAD + INSPECTION\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGE
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

library(data.table)


# ============================================================
# 2. SAFE DOWNLOAD FUNCTION
# ============================================================

download_if_missing <- function(
    url,
    destination,
    max_attempts = 5
) {
  
  if (
    file.exists(destination) &&
    file.info(destination)$size > 0
  ) {
    
    cat(
      "Already available:",
      basename(destination),
      "\n"
    )
    
    return(TRUE)
  }
  
  
  for (attempt in seq_len(max_attempts)) {
    
    cat(
      "\nDownloading:",
      basename(destination),
      "\nAttempt",
      attempt,
      "of",
      max_attempts,
      "\n"
    )
    
    
    success <- tryCatch(
      {
        
        download.file(
          url = url,
          destfile = destination,
          method = "libcurl",
          mode = "wb",
          quiet = FALSE
        )
        
        TRUE
      },
      
      error = function(e) {
        
        cat(
          "\nDownload error:",
          conditionMessage(e),
          "\n"
        )
        
        FALSE
      }
    )
    
    
    if (
      success &&
      file.exists(destination) &&
      file.info(destination)$size > 0
    ) {
      
      cat(
        "Download complete:",
        round(
          file.info(destination)$size / 1024^2,
          3
        ),
        "MB\n"
      )
      
      return(TRUE)
    }
    
    
    if (file.exists(destination)) {
      unlink(destination)
    }
    
    
    if (attempt < max_attempts) {
      Sys.sleep(3)
    }
  }
  
  
  stop(
    paste(
      "Failed to download:",
      basename(destination)
    )
  )
}


# ============================================================
# 3. DOWNLOAD AUTHOR RAW COUNT MATRIX
# ============================================================

counts_filename <- "GSE163854_counts.txt.gz"

counts_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE163nnn/GSE163854/suppl/",
  counts_filename
)

counts_file <- file.path(
  raw_dir,
  counts_filename
)


download_if_missing(
  url = counts_url,
  destination = counts_file
)


# ============================================================
# 4. DOWNLOAD GEO FAMILY SOFT METADATA
# ============================================================

soft_filename <- "GSE163854_family.soft.gz"

soft_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE163nnn/GSE163854/soft/",
  soft_filename
)

soft_file <- file.path(
  raw_dir,
  soft_filename
)


download_if_missing(
  url = soft_url,
  destination = soft_file
)


# ============================================================
# 5. READ RAW COUNTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RAW COUNT MATRIX\n")
cat("============================================================\n\n")


counts <- data.table::fread(
  counts_file,
  data.table = FALSE,
  check.names = FALSE
)


cat(
  "Dimensions:",
  nrow(counts),
  "rows x",
  ncol(counts),
  "columns\n\n"
)


cat(
  "Column names:\n\n"
)

print(
  colnames(counts)
)


cat(
  "\nFirst 10 rows:\n\n"
)

print(
  head(
    counts,
    10
  )
)


cat(
  "\nColumn classes:\n\n"
)

print(
  vapply(
    counts,
    class,
    character(1)
  )
)


# ============================================================
# 6. GENE IDENTIFIER QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GENE IDENTIFIER QC\n")
cat("============================================================\n\n")


gene_ids <- counts[[1]]


cat(
  "Identifier column:",
  colnames(counts)[1],
  "\n"
)


cat(
  "Total gene rows:",
  length(gene_ids),
  "\n"
)


cat(
  "Unique identifiers:",
  length(
    unique(gene_ids)
  ),
  "\n"
)


cat(
  "Duplicated identifiers:",
  sum(
    duplicated(gene_ids)
  ),
  "\n"
)


cat(
  "Missing identifiers:",
  sum(
    is.na(gene_ids)
  ),
  "\n"
)


cat(
  "Blank identifiers:",
  sum(
    trimws(
      as.character(gene_ids)
    ) == ""
  ),
  "\n"
)


# ============================================================
# 7. EXPRESSION MATRIX QC
#
# Explicitly exclude identifier column.
# ============================================================

expr <- counts[
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
cat("COUNT VALUE QC\n")
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

storage.mode(x) <- "numeric"


cat(
  "\nNA values:",
  sum(
    is.na(x)
  ),
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


integer_fraction <- mean(
  abs(
    x - round(x)
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


cat(
  "\nCOUNT QUANTILES:\n\n"
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
  "\nLIBRARY SIZES:\n\n"
)


library_sizes <- colSums(
  x,
  na.rm = TRUE
)


print(
  sort(
    library_sizes
  )
)


# ============================================================
# 8. CHECK FOR VERY LOW-DEPTH SAMPLES
#
# Inspection only - nothing is excluded.
# ============================================================

median_library <- median(
  library_sizes
)


library_fraction <- library_sizes /
  median_library


depth_table <- data.frame(
  
  sample =
    names(library_sizes),
  
  library_size =
    as.numeric(
      library_sizes
    ),
  
  fraction_of_median =
    as.numeric(
      library_fraction
    ),
  
  stringsAsFactors = FALSE
)


depth_table <- depth_table[
  order(
    depth_table$library_size
  ),
  ,
  drop = FALSE
]


cat(
  "\nLibrary depth relative to median:\n\n"
)

print(
  depth_table,
  row.names = FALSE
)


write.csv(
  depth_table,
  file.path(
    raw_dir,
    "GSE163854_library_size_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. PARSE GEO SOFT METADATA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GEO SAMPLE METADATA\n")
cat("============================================================\n\n")


connection <- gzfile(
  soft_file,
  open = "rt"
)


soft_lines <- readLines(
  connection,
  warn = FALSE
)


close(
  connection
)


sample_starts <- grep(
  "^\\^SAMPLE = ",
  soft_lines
)


cat(
  "Samples detected in SOFT:",
  length(sample_starts),
  "\n"
)


if (length(sample_starts) == 0) {
  
  stop(
    "No SAMPLE blocks detected in GEO SOFT metadata."
  )
}


sample_ends <- c(
  sample_starts[-1] - 1,
  length(soft_lines)
)


# ============================================================
# 10. SOFT FIELD EXTRACTION FUNCTIONS
# ============================================================

extract_first <- function(
    block,
    prefix
) {
  
  matches <- grep(
    prefix,
    block,
    value = TRUE
  )
  
  
  if (length(matches) == 0) {
    return(NA_character_)
  }
  
  
  sub(
    prefix,
    "",
    matches[1]
  )
}


extract_all <- function(
    block,
    prefix
) {
  
  matches <- grep(
    prefix,
    block,
    value = TRUE
  )
  
  
  if (length(matches) == 0) {
    return(NA_character_)
  }
  
  
  values <- sub(
    prefix,
    "",
    matches
  )
  
  
  paste(
    values,
    collapse = " | "
  )
}


# ============================================================
# 11. EXTRACT SAMPLE METADATA
# ============================================================

metadata_rows <- vector(
  mode = "list",
  length = length(sample_starts)
)


for (i in seq_along(sample_starts)) {
  
  block <- soft_lines[
    sample_starts[i]:
      sample_ends[i]
  ]
  
  
  sample_id <- sub(
    "^\\^SAMPLE =\\s*",
    "",
    block[1]
  )
  
  
  metadata_rows[[i]] <- data.frame(
    
    sample_id =
      sample_id,
    
    title =
      extract_first(
        block,
        "^!Sample_title = "
      ),
    
    source =
      extract_first(
        block,
        "^!Sample_source_name_ch1 = "
      ),
    
    characteristics =
      extract_all(
        block,
        "^!Sample_characteristics_ch1 = "
      ),
    
    treatment_protocol =
      extract_all(
        block,
        "^!Sample_treatment_protocol_ch1 = "
      ),
    
    growth_protocol =
      extract_all(
        block,
        "^!Sample_growth_protocol_ch1 = "
      ),
    
    description =
      extract_all(
        block,
        "^!Sample_description = "
      ),
    
    stringsAsFactors = FALSE
  )
}


sample_metadata <- do.call(
  rbind,
  metadata_rows
)


write.csv(
  sample_metadata,
  file.path(
    raw_dir,
    "GSE163854_sample_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. PRINT SAMPLE TITLES
# ============================================================

cat(
  "\nSAMPLE TITLES:\n\n"
)


print(
  sample_metadata[
    ,
    c(
      "sample_id",
      "title",
      "source"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 13. PRINT CHARACTERISTICS
# ============================================================

cat(
  "\nSAMPLE CHARACTERISTICS:\n\n"
)


print(
  sample_metadata[
    ,
    c(
      "sample_id",
      "title",
      "characteristics"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14. PRINT TREATMENT PROTOCOL
# ============================================================

cat(
  "\nTREATMENT PROTOCOLS:\n\n"
)


print(
  sample_metadata[
    ,
    c(
      "sample_id",
      "title",
      "treatment_protocol"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. COUNT COLUMN / SAMPLE MAPPING
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COUNT COLUMN / METADATA MAPPING\n")
cat("============================================================\n\n")


count_sample_names <- colnames(
  expr
)


cat(
  "Count expression columns:",
  length(
    count_sample_names
  ),
  "\n"
)


cat(
  "\nCount column names:\n\n"
)

print(
  count_sample_names
)


# ------------------------------------------------------------
# Direct GSM mapping
# ------------------------------------------------------------

gsm_columns <- grep(
  "^GSM[0-9]+$",
  count_sample_names,
  value = TRUE
)


cat(
  "\nGSM-formatted count columns:",
  length(
    gsm_columns
  ),
  "\n"
)


if (length(gsm_columns) > 0) {
  
  matched_gsm <- intersect(
    gsm_columns,
    sample_metadata$sample_id
  )
  
  
  cat(
    "Direct GSM matches:",
    length(matched_gsm),
    "/",
    length(gsm_columns),
    "\n"
  )
  
  
  cat(
    "\nUnmatched GSM columns:\n"
  )
  
  print(
    setdiff(
      gsm_columns,
      sample_metadata$sample_id
    )
  )
}


# ------------------------------------------------------------
# Exact title mapping
# ------------------------------------------------------------

title_matches <- count_sample_names[
  count_sample_names %in%
    sample_metadata$title
]


cat(
  "\nExact count-column/title matches:",
  length(title_matches),
  "/",
  length(count_sample_names),
  "\n"
)


# ============================================================
# 16. PRELIMINARY BIOLOGICAL GROUP LIST
#
# This section ONLY reports GEO titles.
# It does NOT freeze the final contrast.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRELIMINARY GEO PHENOTYPE GROUPS\n")
cat("============================================================\n\n")


sensitive_rows <- grepl(
  "Sensitive",
  sample_metadata$title,
  ignore.case = TRUE
)


mixed_rows <- grepl(
  "Mixed",
  sample_metadata$title,
  ignore.case = TRUE
)


resistant_rows <- grepl(
  "Resistant",
  sample_metadata$title,
  ignore.case = TRUE
)


cat(
  "Sensitive-labelled samples:\n"
)

print(
  sample_metadata[
    sensitive_rows,
    c(
      "sample_id",
      "title"
    )
  ],
  row.names = FALSE
)


cat(
  "\nMixed-phenotype samples:\n"
)

print(
  sample_metadata[
    mixed_rows,
    c(
      "sample_id",
      "title"
    )
  ],
  row.names = FALSE
)


cat(
  "\nResistant-labelled samples:\n"
)

print(
  sample_metadata[
    resistant_rows,
    c(
      "sample_id",
      "title"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 17. SAVE RAW OBJECT
# ============================================================

saveRDS(
  counts,
  file.path(
    raw_dir,
    "GSE163854_raw_counts.rds"
  )
)


# ============================================================
# 18. FILE INVENTORY
# ============================================================

files <- list.files(
  raw_dir,
  full.names = TRUE
)


file_inventory <- data.frame(
  
  filename =
    basename(files),
  
  size_MB =
    round(
      file.info(files)$size /
        1024^2,
      3
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  file_inventory,
  file.path(
    raw_dir,
    "GSE163854_file_inventory.csv"
  ),
  row.names = FALSE
)


cat(
  "\nLOCAL FILE INVENTORY:\n\n"
)

print(
  file_inventory,
  row.names = FALSE
)


# ============================================================
# 19. SESSION INFO
# ============================================================

sink(
  file.path(
    raw_dir,
    "GSE163854_download_inspect_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 20. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE163854 D8 DOWNLOAD + INSPECTION COMPLETE\n")
cat("============================================================\n\n")


cat(
  "NO filtering, normalization, PCA, differential expression, ",
  "or GSEA was performed.\n\n"
)


cat(
  "Next decisions AFTER inspection:\n",
  "  1. exact sensitive-vs-resistant sample mapping\n",
  "  2. whether passage itself is a confounder\n",
  "  3. treatment state at RNA collection\n",
  "  4. primary count filtering rule\n",
  "  5. normalization / DE model\n",
  "  6. sensitivity analysis design\n"
)


cat("\n============================================================\n")