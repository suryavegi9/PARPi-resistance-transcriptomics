# ============================================================
# GSE165914 - D7
# CLEAN DOWNLOAD + INSPECTION SCRIPT
#
# 1. Download NCBI standardized RNA-seq raw counts
# 2. Retry automatically if download is interrupted
# 3. Validate and decompress counts
# 4. Download GEO family SOFT metadata
# 5. Extract sample titles / characteristics
# 6. Perform basic count QC
#
# NO DIFFERENTIAL EXPRESSION IS PERFORMED HERE
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

accession <- "GSE165914"

raw_dir <- file.path(
  "data_raw",
  accession
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. PACKAGES
# ============================================================

packages <- c(
  "curl",
  "data.table",
  "R.utils"
)

for (pkg in packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}

library(curl)
library(data.table)
library(R.utils)


cat("\n")
cat("============================================================\n")
cat("GSE165914 - DOWNLOAD + INSPECTION\n")
cat("============================================================\n\n")


# ============================================================
# 2. DOWNLOAD FUNCTION WITH RETRIES
# ============================================================

download_with_retry <- function(
    url,
    dest,
    max_attempts = 8
) {
  
  for (attempt in seq_len(max_attempts)) {
    
    cat(
      "\nDownload attempt",
      attempt,
      "of",
      max_attempts,
      "\n"
    )
    
    result <- tryCatch(
      
      {
        
        h <- curl::new_handle()
        
        curl::handle_setopt(
          h,
          connecttimeout = 30,
          timeout = 0
        )
        
        curl::curl_download(
          url = url,
          destfile = dest,
          quiet = FALSE,
          mode = "wb",
          handle = h
        )
        
        TRUE
      },
      
      error = function(e) {
        
        cat(
          "\nDownload error:\n",
          conditionMessage(e),
          "\n"
        )
        
        FALSE
      }
      
    )
    
    
    if (isTRUE(result) && file.exists(dest)) {
      
      cat(
        "\nDownload finished.\n"
      )
      
      cat(
        "File size:",
        round(
          file.info(dest)$size / 1024^2,
          2
        ),
        "MB\n"
      )
      
      return(TRUE)
    }
    
    
    if (file.exists(dest)) {
      
      cat(
        "Removing incomplete file before retry.\n"
      )
      
      unlink(dest)
    }
    
    
    if (attempt < max_attempts) {
      
      cat(
        "Retrying...\n"
      )
      
      Sys.sleep(3)
    }
    
  }
  
  
  FALSE
}


# ============================================================
# 3. GZIP VALIDATION FUNCTION
# ============================================================

validate_gzip <- function(path) {
  
  if (!file.exists(path)) {
    return(FALSE)
  }
  
  
  temp_out <- tempfile(
    fileext = ".tsv"
  )
  
  
  valid <- tryCatch(
    
    {
      
      R.utils::gunzip(
        filename = path,
        destname = temp_out,
        remove = FALSE,
        overwrite = TRUE
      )
      
      TRUE
    },
    
    error = function(e) {
      
      cat(
        "\nGZIP validation failed:\n",
        conditionMessage(e),
        "\n"
      )
      
      FALSE
    }
    
  )
  
  
  if (file.exists(temp_out)) {
    unlink(temp_out)
  }
  
  
  valid
}


# ============================================================
# 4. DOWNLOAD NCBI RAW COUNTS
# ============================================================

counts_filename <- paste0(
  accession,
  "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)


counts_gz <- file.path(
  raw_dir,
  counts_filename
)


counts_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/",
  "?type=rnaseq_counts",
  "&acc=", accession,
  "&format=file",
  "&file=", counts_filename
)


cat("\n")
cat("============================================================\n")
cat("DOWNLOADING NCBI RAW COUNTS\n")
cat("============================================================\n\n")


already_valid <- validate_gzip(
  counts_gz
)


if (already_valid) {
  
  cat(
    "Valid count file already exists.\n"
  )
  
} else {
  
  if (file.exists(counts_gz)) {
    unlink(counts_gz)
  }
  
  
  success <- download_with_retry(
    url = counts_url,
    dest = counts_gz,
    max_attempts = 8
  )
  
  
  if (!success) {
    
    stop(
      "Count matrix could not be downloaded after repeated attempts."
    )
    
  }
  
  
  if (!validate_gzip(counts_gz)) {
    
    stop(
      "Count file downloaded but gzip validation failed."
    )
    
  }
  
}


cat(
  "\nCount-file gzip validation: PASSED\n"
)


# ============================================================
# 5. DECOMPRESS COUNT MATRIX
# ============================================================

counts_tsv <- file.path(
  raw_dir,
  paste0(
    accession,
    "_raw_counts_GRCh38.p13_NCBI.tsv"
  )
)


if (!file.exists(counts_tsv)) {
  
  cat(
    "\nDecompressing count matrix...\n"
  )
  
  
  R.utils::gunzip(
    filename = counts_gz,
    destname = counts_tsv,
    remove = FALSE,
    overwrite = TRUE
  )
  
}


# ============================================================
# 6. READ COUNT MATRIX
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RAW COUNT MATRIX QC\n")
cat("============================================================\n\n")


counts <- data.table::fread(
  counts_tsv,
  data.table = FALSE
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
  "\nFirst 6 rows:\n\n"
)

print(
  head(
    counts,
    6
  )
)


# ============================================================
# 7. IDENTIFY EXPRESSION COLUMNS
# ============================================================

numeric_columns <- vapply(
  counts,
  is.numeric,
  logical(1)
)


cat(
  "\nNumeric columns:",
  sum(numeric_columns),
  "\n"
)


cat(
  "Non-numeric columns:",
  sum(!numeric_columns),
  "\n\n"
)


if (any(!numeric_columns)) {
  
  cat(
    "Non-numeric columns:\n"
  )
  
  print(
    colnames(counts)[
      !numeric_columns
    ]
  )
  
}


count_matrix <- as.matrix(
  counts[
    ,
    numeric_columns,
    drop = FALSE
  ]
)


# ============================================================
# 8. COUNT QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COUNT VALUE QC\n")
cat("============================================================\n\n")


cat(
  "NA values:",
  sum(
    is.na(count_matrix)
  ),
  "\n"
)


cat(
  "Negative values:",
  sum(
    count_matrix < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Exact zeros:",
  sum(
    count_matrix == 0,
    na.rm = TRUE
  ),
  "\n"
)


integer_fraction <- mean(
  abs(
    count_matrix -
      round(count_matrix)
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
    as.numeric(count_matrix),
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
  count_matrix,
  na.rm = TRUE
)


print(
  sort(
    library_sizes
  )
)


# ============================================================
# 9. SAVE COUNT MATRIX AS RDS
# ============================================================

saveRDS(
  counts,
  file.path(
    raw_dir,
    "GSE165914_NCBI_raw_counts.rds"
  )
)


# ============================================================
# 10. DOWNLOAD GEO FAMILY SOFT METADATA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("DOWNLOADING GEO SAMPLE METADATA\n")
cat("============================================================\n\n")


soft_filename <- paste0(
  accession,
  "_family.soft.gz"
)


soft_gz <- file.path(
  raw_dir,
  soft_filename
)


soft_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE165nnn/",
  accession,
  "/soft/",
  soft_filename
)


if (!file.exists(soft_gz)) {
  
  soft_success <- download_with_retry(
    url = soft_url,
    dest = soft_gz,
    max_attempts = 8
  )
  
} else {
  
  soft_success <- TRUE
  
  cat(
    "SOFT metadata file already exists.\n"
  )
  
}


# ============================================================
# 11. PARSE SOFT METADATA MANUALLY
# ============================================================

sample_metadata <- NULL


if (
  isTRUE(soft_success) &&
  file.exists(soft_gz)
) {
  
  cat(
    "\nParsing GEO SOFT metadata...\n"
  )
  
  
  con <- gzfile(
    soft_gz,
    open = "rt"
  )
  
  
  soft_lines <- readLines(
    con,
    warn = FALSE
  )
  
  
  close(con)
  
  
  sample_starts <- grep(
    "^\\^SAMPLE = ",
    soft_lines
  )
  
  
  cat(
    "Samples detected in SOFT:",
    length(sample_starts),
    "\n"
  )
  
  
  if (length(sample_starts) > 0) {
    
    sample_ends <- c(
      sample_starts[-1] - 1,
      length(soft_lines)
    )
    
    
    extract_one <- function(
    block,
    prefix
    ) {
      
      hits <- grep(
        prefix,
        block,
        value = TRUE
      )
      
      
      if (length(hits) == 0) {
        return(NA_character_)
      }
      
      
      sub(
        prefix,
        "",
        hits[1]
      )
    }
    
    
    extract_many <- function(
    block,
    prefix
    ) {
      
      hits <- grep(
        prefix,
        block,
        value = TRUE
      )
      
      
      if (length(hits) == 0) {
        return(NA_character_)
      }
      
      
      paste(
        sub(
          prefix,
          "",
          hits
        ),
        collapse = " | "
      )
    }
    
    
    metadata_list <- vector(
      "list",
      length(sample_starts)
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
      
      
      metadata_list[[i]] <- data.frame(
        
        sample_id = sample_id,
        
        title = extract_one(
          block,
          "^!Sample_title = "
        ),
        
        source_name = extract_one(
          block,
          "^!Sample_source_name_ch1 = "
        ),
        
        characteristics = extract_many(
          block,
          "^!Sample_characteristics_ch1 = "
        ),
        
        growth_protocol = extract_many(
          block,
          "^!Sample_growth_protocol_ch1 = "
        ),
        
        treatment_protocol = extract_many(
          block,
          "^!Sample_treatment_protocol_ch1 = "
        ),
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
    
    sample_metadata <- do.call(
      rbind,
      metadata_list
    )
    
    
    write.csv(
      sample_metadata,
      file.path(
        raw_dir,
        "GSE165914_sample_metadata.csv"
      ),
      row.names = FALSE
    )
    
    
    cat("\n")
    cat("========================================\n")
    cat("GEO SAMPLE TITLES\n")
    cat("========================================\n\n")
    
    
    print(
      sample_metadata[
        ,
        c(
          "sample_id",
          "title",
          "source_name"
        )
      ],
      row.names = FALSE
    )
    
    
    cat("\n")
    cat("========================================\n")
    cat("SAMPLE CHARACTERISTICS\n")
    cat("========================================\n\n")
    
    
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
    
  }
  
}


# ============================================================
# 12. COUNT COLUMN / GSM MATCHING
# ============================================================

cat("\n")
cat("============================================================\n")
cat("COUNT / METADATA MATCHING\n")
cat("============================================================\n\n")


count_columns <- colnames(
  count_matrix
)


cat(
  "Numeric count columns:",
  length(count_columns),
  "\n"
)


if (!is.null(sample_metadata)) {
  
  gsm_in_counts <- grep(
    "^GSM[0-9]+$",
    count_columns,
    value = TRUE
  )
  
  
  matched <- intersect(
    gsm_in_counts,
    sample_metadata$sample_id
  )
  
  
  cat(
    "GSM-style count columns:",
    length(gsm_in_counts),
    "\n"
  )
  
  
  cat(
    "Matched to metadata:",
    length(matched),
    "\n"
  )
  
  
  if (length(gsm_in_counts) > 0) {
    
    cat(
      "\nGSM columns without metadata match:\n"
    )
    
    print(
      setdiff(
        gsm_in_counts,
        sample_metadata$sample_id
      )
    )
    
  }
  
}


# ============================================================
# 13. SAVE FILE INVENTORY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("LOCAL FILE INVENTORY\n")
cat("============================================================\n\n")


files <- list.files(
  raw_dir,
  full.names = TRUE
)


inventory <- data.frame(
  
  filename = basename(files),
  
  size_MB = round(
    file.info(files)$size /
      1024^2,
    3
  ),
  
  stringsAsFactors = FALSE
  
)


print(
  inventory,
  row.names = FALSE
)


write.csv(
  inventory,
  file.path(
    raw_dir,
    "GSE165914_file_inventory.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. SAVE SESSION INFO
# ============================================================

sink(
  file.path(
    raw_dir,
    "GSE165914_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 15. FINISH
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE165914 DOWNLOAD + INSPECTION COMPLETE\n")
cat("============================================================\n\n")

cat(
  "No differential-expression analysis was performed.\n\n"
)

cat(
  "Next step:\n",
  "1. Identify MDA-MB-468 parental and OlaR triplicates.\n",
  "2. Identify SUM1315 parental and OlaR triplicates.\n",
  "3. Exclude acute-treatment / cisplatin / combination samples.\n",
  "4. Perform PCA and replicate QC.\n",
  "5. Freeze D7A and D7B before DE analysis.\n"
)

cat("\n============================================================\n")