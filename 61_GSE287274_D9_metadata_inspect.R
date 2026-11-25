# ============================================================
# GSE287274 - POSSIBLE D9
# GEO METADATA + FILE-STRUCTURE INSPECTION
#
# PURPOSE:
#
# Determine whether GSE287274 satisfies the frozen eligibility
# criteria for the PARPi resistance cross-dataset study.
#
# THIS SCRIPT DOES NOT:
#   - select a final contrast
#   - perform filtering
#   - normalize expression
#   - perform PCA
#   - perform DE
#   - perform GSEA
#
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

accession <- "GSE287274"

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
cat("GSE287274 - POSSIBLE D9 METADATA AUDIT\n")
cat("============================================================\n\n")


# ============================================================
# 1. DOWNLOAD GEO FAMILY SOFT
# ============================================================

soft_filename <-
  "GSE287274_family.soft.gz"


soft_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE287nnn/GSE287274/soft/",
  soft_filename
)


soft_file <- file.path(
  raw_dir,
  soft_filename
)


if (
  !file.exists(soft_file) ||
  file.info(soft_file)$size == 0
) {
  
  cat(
    "Downloading GEO SOFT metadata...\n"
  )
  
  
  download.file(
    url = soft_url,
    destfile = soft_file,
    method = "libcurl",
    mode = "wb",
    quiet = FALSE
  )
  
} else {
  
  cat(
    "SOFT metadata already downloaded.\n"
  )
}


if (
  !file.exists(soft_file) ||
  file.info(soft_file)$size == 0
) {
  
  stop(
    "Failed to obtain GSE287274 family SOFT file."
  )
}


cat(
  "SOFT file size:",
  round(
    file.info(soft_file)$size / 1024^2,
    3
  ),
  "MB\n\n"
)


# ============================================================
# 2. READ SOFT FILE
# ============================================================

con <- gzfile(
  soft_file,
  open = "rt"
)


soft_lines <- readLines(
  con,
  warn = FALSE
)


close(con)


cat(
  "SOFT lines:",
  length(soft_lines),
  "\n"
)


# ============================================================
# 3. SERIES-LEVEL INFORMATION
# ============================================================

series_title <- grep(
  "^!Series_title = ",
  soft_lines,
  value = TRUE
)


series_summary <- grep(
  "^!Series_summary = ",
  soft_lines,
  value = TRUE
)


series_overall_design <- grep(
  "^!Series_overall_design = ",
  soft_lines,
  value = TRUE
)


series_supplementary <- grep(
  "^!Series_supplementary_file",
  soft_lines,
  value = TRUE
)


series_relation <- grep(
  "^!Series_relation",
  soft_lines,
  value = TRUE
)


cat("\n")
cat("============================================================\n")
cat("SERIES INFORMATION\n")
cat("============================================================\n\n")


cat(
  "TITLE:\n"
)

print(
  series_title
)


cat(
  "\nSUMMARY:\n"
)

print(
  series_summary
)


cat(
  "\nOVERALL DESIGN:\n"
)

print(
  series_overall_design
)


cat(
  "\nSUPPLEMENTARY FILES:\n"
)

print(
  series_supplementary
)


cat(
  "\nSERIES RELATIONS:\n"
)

print(
  series_relation
)


# ============================================================
# 4. SAVE SERIES-LEVEL TEXT
# ============================================================

series_output <- c(
  
  "===== TITLE =====",
  series_title,
  
  "",
  "===== SUMMARY =====",
  series_summary,
  
  "",
  "===== OVERALL DESIGN =====",
  series_overall_design,
  
  "",
  "===== SUPPLEMENTARY FILES =====",
  series_supplementary,
  
  "",
  "===== RELATIONS =====",
  series_relation
)


writeLines(
  series_output,
  file.path(
    raw_dir,
    "GSE287274_series_information.txt"
  )
)


# ============================================================
# 5. IDENTIFY SAMPLE BLOCKS
# ============================================================

sample_starts <- grep(
  "^\\^SAMPLE = ",
  soft_lines
)


if (length(sample_starts) == 0) {
  
  stop(
    "No GEO sample blocks detected."
  )
}


sample_ends <- c(
  sample_starts[-1] - 1,
  length(soft_lines)
)


cat(
  "\nNumber of GEO samples:",
  length(sample_starts),
  "\n"
)


# ============================================================
# 6. FIELD EXTRACTION FUNCTIONS
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
    
    return(
      NA_character_
    )
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
    
    return(
      NA_character_
    )
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
# 7. PARSE EACH SAMPLE
# ============================================================

sample_rows <- vector(
  mode = "list",
  length = length(sample_starts)
)


for (i in seq_along(sample_starts)) {
  
  block <- soft_lines[
    sample_starts[i]:
      sample_ends[i]
  ]
  
  
  gsm <- sub(
    "^\\^SAMPLE =\\s*",
    "",
    block[1]
  )
  
  
  sample_rows[[i]] <- data.frame(
    
    sample_id =
      gsm,
    
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
    
    organism =
      extract_first(
        block,
        "^!Sample_organism_ch1 = "
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
    
    extract_protocol =
      extract_all(
        block,
        "^!Sample_extract_protocol_ch1 = "
      ),
    
    description =
      extract_all(
        block,
        "^!Sample_description = "
      ),
    
    library_strategy =
      extract_first(
        block,
        "^!Sample_library_strategy = "
      ),
    
    library_source =
      extract_first(
        block,
        "^!Sample_library_source = "
      ),
    
    supplementary_file =
      extract_all(
        block,
        "^!Sample_supplementary_file"
      ),
    
    relation =
      extract_all(
        block,
        "^!Sample_relation"
      ),
    
    stringsAsFactors = FALSE
  )
}


sample_metadata <- do.call(
  rbind,
  sample_rows
)


rownames(
  sample_metadata
) <- NULL


# ============================================================
# 8. SAVE COMPLETE METADATA
# ============================================================

write.csv(
  sample_metadata,
  file.path(
    raw_dir,
    "GSE287274_complete_sample_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. PRINT BASIC SAMPLE TABLE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SAMPLE TITLES\n")
cat("============================================================\n\n")


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
# 10. PRINT CHARACTERISTICS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SAMPLE CHARACTERISTICS\n")
cat("============================================================\n\n")


for (i in seq_len(nrow(sample_metadata))) {
  
  cat(
    "\n--- ",
    sample_metadata$sample_id[i],
    " | ",
    sample_metadata$title[i],
    " ---\n",
    sep = ""
  )
  
  
  cat(
    "Source:\n",
    sample_metadata$source[i],
    "\n\n",
    sep = ""
  )
  
  
  cat(
    "Characteristics:\n",
    sample_metadata$characteristics[i],
    "\n\n",
    sep = ""
  )
  
  
  cat(
    "Treatment:\n",
    sample_metadata$treatment_protocol[i],
    "\n\n",
    sep = ""
  )
  
  
  cat(
    "Description:\n",
    sample_metadata$description[i],
    "\n",
    sep = ""
  )
}


# ============================================================
# 11. KEYWORD SCREEN
#
# This does NOT classify samples biologically.
# It only helps us find candidate resistant/control samples.
# ============================================================

combined_text <- paste(
  
  sample_metadata$title,
  
  sample_metadata$source,
  
  sample_metadata$characteristics,
  
  sample_metadata$treatment_protocol,
  
  sample_metadata$description,
  
  sep = " | "
)


combined_lower <- tolower(
  combined_text
)


keyword_table <- data.frame(
  
  sample_id =
    sample_metadata$sample_id,
  
  title =
    sample_metadata$title,
  
  resistant =
    grepl(
      "resist",
      combined_lower
    ),
  
  parental =
    grepl(
      "parent",
      combined_lower
    ),
  
  control =
    grepl(
      "control|ctrl",
      combined_lower
    ),
  
  olaparib =
    grepl(
      "olaparib",
      combined_lower
    ),
  
  parpi =
    grepl(
      "parp",
      combined_lower
    ),
  
  yap1 =
    grepl(
      "yap",
      combined_lower
    ),
  
  taz =
    grepl(
      "taz",
      combined_lower
    ),
  
  brd4 =
    grepl(
      "brd4",
      combined_lower
    ),
  
  knockdown =
    grepl(
      "knockdown|kd|shr|sirna",
      combined_lower
    ),
  
  overexpression =
    grepl(
      "overexpress|oe|transfect",
      combined_lower
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  keyword_table,
  file.path(
    raw_dir,
    "GSE287274_keyword_screen.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("KEYWORD SCREEN\n")
cat("============================================================\n\n")


print(
  keyword_table,
  row.names = FALSE
)


# ============================================================
# 12. CANDIDATE RESISTANT / CONTROL SAMPLES
# ============================================================

candidate_resistant <- keyword_table[
  keyword_table$resistant |
    keyword_table$olaparib |
    keyword_table$parpi,
  ,
  drop = FALSE
]


candidate_control <- keyword_table[
  keyword_table$parental |
    keyword_table$control,
  ,
  drop = FALSE
]


cat("\n")
cat("============================================================\n")
cat("CANDIDATE PARPI / RESISTANCE-RELATED SAMPLES\n")
cat("============================================================\n\n")


print(
  candidate_resistant,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CANDIDATE CONTROL / PARENTAL SAMPLES\n")
cat("============================================================\n\n")


print(
  candidate_control,
  row.names = FALSE
)


# ============================================================
# 13. LIBRARY STRATEGY COUNTS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("LIBRARY STRATEGIES\n")
cat("============================================================\n\n")


print(
  table(
    sample_metadata$library_strategy,
    useNA = "ifany"
  )
)


cat("\nLibrary sources:\n\n")


print(
  table(
    sample_metadata$library_source,
    useNA = "ifany"
  )
)


# ============================================================
# 14. SAMPLE SUPPLEMENTARY FILES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SAMPLE SUPPLEMENTARY FILES\n")
cat("============================================================\n\n")


print(
  sample_metadata[
    ,
    c(
      "sample_id",
      "title",
      "supplementary_file"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. SAVE RAW SOFT SAMPLE BLOCKS
#
# Useful if any metadata field was not parsed as expected.
# ============================================================

sample_block_dir <- file.path(
  raw_dir,
  "SOFT_sample_blocks"
)


dir.create(
  sample_block_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


for (i in seq_along(sample_starts)) {
  
  block <- soft_lines[
    sample_starts[i]:
      sample_ends[i]
  ]
  
  
  gsm <- sample_metadata$sample_id[i]
  
  
  writeLines(
    block,
    file.path(
      sample_block_dir,
      paste0(
        gsm,
        ".txt"
      )
    )
  )
}


# ============================================================
# 16. AUDIT SUMMARY
# ============================================================

audit_summary <- data.frame(
  
  metric = c(
    "Total_GEO_samples",
    "RNA_seq_samples",
    "Resistance_keyword_samples",
    "Parental_keyword_samples",
    "Control_keyword_samples",
    "Olaparib_keyword_samples",
    "PARP_keyword_samples",
    "YAP_keyword_samples",
    "TAZ_keyword_samples",
    "BRD4_keyword_samples"
  ),
  
  value = c(
    nrow(sample_metadata),
    
    sum(
      sample_metadata$library_strategy == "RNA-Seq",
      na.rm = TRUE
    ),
    
    sum(
      keyword_table$resistant
    ),
    
    sum(
      keyword_table$parental
    ),
    
    sum(
      keyword_table$control
    ),
    
    sum(
      keyword_table$olaparib
    ),
    
    sum(
      keyword_table$parpi
    ),
    
    sum(
      keyword_table$yap1
    ),
    
    sum(
      keyword_table$taz
    ),
    
    sum(
      keyword_table$brd4
    )
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  audit_summary,
  file.path(
    raw_dir,
    "GSE287274_metadata_audit_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("AUDIT SUMMARY\n")
cat("============================================================\n\n")


print(
  audit_summary,
  row.names = FALSE
)


# ============================================================
# 17. SESSION INFO
# ============================================================

sink(
  file.path(
    raw_dir,
    "GSE287274_metadata_inspect_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 18. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE287274 METADATA AUDIT COMPLETE\n")
cat("============================================================\n\n")


cat(
  "NO eligibility decision has been made yet.\n\n"
)


cat(
  "We will include GSE287274 only if GEO confirms:\n",
  "  1. experimentally acquired PARPi resistance\n",
  "  2. appropriate sensitive/parental comparator\n",
  "  3. adequate biological replication\n",
  "  4. unperturbed resistance contrast suitable for DE\n",
  "  5. expression data suitable for independent reanalysis\n"
)


cat(
  "\n============================================================\n"
)