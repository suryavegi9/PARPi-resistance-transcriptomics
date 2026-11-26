# ============================================================
# GSE287274 - POSSIBLE D9
# ACCESS RECOVERY + AVAILABILITY AUDIT
#
# Purpose:
#   Determine whether GSE287274 can currently be retrieved
#   reproducibly from official NCBI/GEO infrastructure.
#
# NO biological analysis is performed.
# NO inclusion/exclusion decision is made inside this script.
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
cat("GSE287274 - D9 ACCESS RECOVERY AUDIT\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

if (!requireNamespace(
  "curl",
  quietly = TRUE
)) {
  
  install.packages("curl")
}


if (!requireNamespace(
  "jsonlite",
  quietly = TRUE
)) {
  
  install.packages("jsonlite")
}


if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages("BiocManager")
}


if (!requireNamespace(
  "GEOquery",
  quietly = TRUE
)) {
  
  BiocManager::install(
    "GEOquery",
    ask = FALSE,
    update = FALSE
  )
}


library(curl)
library(jsonlite)
library(GEOquery)


# ============================================================
# 2. REMOVE FAILED ZERO-BYTE FILE FROM OLD SCRIPT
# ============================================================

old_soft <- file.path(
  raw_dir,
  "GSE287274_family.soft.gz"
)


if (
  file.exists(old_soft) &&
  file.info(old_soft)$size == 0
) {
  
  unlink(old_soft)
  
  cat(
    "Removed zero-byte SOFT file from previous failed attempt.\n"
  )
}


# ============================================================
# 3. HELPER: SAFE TEXT REQUEST
# ============================================================

safe_text_request <- function(
    url,
    label
) {
  
  cat(
    "\n------------------------------------------------------------\n"
  )
  
  cat(
    label,
    "\n"
  )
  
  cat(
    url,
    "\n"
  )
  
  
  result <- tryCatch(
    
    {
      
      response <- curl::curl_fetch_memory(
        url
      )
      
      
      text <- rawToChar(
        response$content
      )
      
      
      cat(
        "HTTP status:",
        response$status_code,
        "\n"
      )
      
      
      list(
        success =
          response$status_code >= 200 &&
          response$status_code < 300,
        
        status =
          response$status_code,
        
        text =
          text
      )
      
    },
    
    error = function(e) {
      
      cat(
        "Request failed:",
        conditionMessage(e),
        "\n"
      )
      
      
      list(
        success = FALSE,
        status = NA_integer_,
        text = NA_character_
      )
    }
  )
  
  
  return(result)
}


# ============================================================
# 4. NCBI GDS E-UTILS SEARCH
#
# This independently asks NCBI whether the accession is
# indexed in the GEO DataSets database.
# ============================================================

gds_search_url <- paste0(
  
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  "esearch.fcgi?db=gds&term=",
  accession,
  "%5BACCN%5D&retmode=json"
)


gds_search <- safe_text_request(
  gds_search_url,
  "NCBI GDS accession search"
)


gds_ids <- character(0)


if (gds_search$success) {
  
  writeLines(
    gds_search$text,
    file.path(
      raw_dir,
      "GSE287274_NCBI_GDS_esearch.json"
    )
  )
  
  
  parsed <- tryCatch(
    
    jsonlite::fromJSON(
      gds_search$text
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  if (!is.null(parsed)) {
    
    gds_ids <- as.character(
      parsed$esearchresult$idlist
    )
    
    
    cat(
      "NCBI GDS records found:",
      length(gds_ids),
      "\n"
    )
    
    
    if (length(gds_ids) > 0) {
      
      cat(
        "GDS internal IDs:",
        paste(
          gds_ids,
          collapse = ", "
        ),
        "\n"
      )
    }
  }
}


# ============================================================
# 5. NCBI GDS SUMMARY
# ============================================================

gds_summary_success <- FALSE


if (length(gds_ids) > 0) {
  
  gds_summary_url <- paste0(
    
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "esummary.fcgi?db=gds&id=",
    
    paste(
      gds_ids,
      collapse = ","
    ),
    
    "&retmode=json"
  )
  
  
  gds_summary <- safe_text_request(
    gds_summary_url,
    "NCBI GDS record summary"
  )
  
  
  if (gds_summary$success) {
    
    gds_summary_success <- TRUE
    
    
    writeLines(
      gds_summary$text,
      file.path(
        raw_dir,
        "GSE287274_NCBI_GDS_esummary.json"
      )
    )
    
    
    parsed_summary <- tryCatch(
      
      jsonlite::fromJSON(
        gds_summary$text,
        simplifyVector = FALSE
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (!is.null(parsed_summary)) {
      
      cat(
        "\nGDS summary successfully retrieved.\n"
      )
      
      
      cat(
        "Saved complete JSON for inspection.\n"
      )
    }
  }
}


# ============================================================
# 6. SEARCH NCBI SRA FOR THE GEO ACCESSION
#
# Useful if GEO family metadata disappeared but sequencing
# runs remain archived in SRA.
# ============================================================

sra_search_url <- paste0(
  
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  "esearch.fcgi?db=sra&term=",
  accession,
  "%5BAll%20Fields%5D&retmax=500&retmode=json"
)


sra_search <- safe_text_request(
  sra_search_url,
  "NCBI SRA search"
)


sra_ids <- character(0)


if (sra_search$success) {
  
  writeLines(
    sra_search$text,
    file.path(
      raw_dir,
      "GSE287274_NCBI_SRA_esearch.json"
    )
  )
  
  
  parsed_sra <- tryCatch(
    
    jsonlite::fromJSON(
      sra_search$text
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  if (!is.null(parsed_sra)) {
    
    sra_ids <- as.character(
      parsed_sra$esearchresult$idlist
    )
    
    
    cat(
      "SRA records linked to accession:",
      length(sra_ids),
      "\n"
    )
  }
}


# ============================================================
# 7. SRA SUMMARY
# ============================================================

sra_summary_success <- FALSE


if (length(sra_ids) > 0) {
  
  sra_summary_url <- paste0(
    
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "esummary.fcgi?db=sra&id=",
    
    paste(
      sra_ids,
      collapse = ","
    ),
    
    "&retmode=json"
  )
  
  
  sra_summary <- safe_text_request(
    sra_summary_url,
    "NCBI SRA summaries"
  )
  
  
  if (sra_summary$success) {
    
    sra_summary_success <- TRUE
    
    
    writeLines(
      sra_summary$text,
      file.path(
        raw_dir,
        "GSE287274_NCBI_SRA_esummary.json"
      )
    )
  }
}


# ============================================================
# 8. GEO ACCESSION VIEWER - TEXT MODE
#
# This route does not depend on the GEO FTP SOFT directory.
# ============================================================

geo_text_url <- paste0(
  
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
  "acc=",
  accession,
  "&targ=self&form=text&view=full"
)


geo_text <- safe_text_request(
  geo_text_url,
  "GEO accession viewer - text/full"
)


geo_text_success <- FALSE


if (geo_text$success) {
  
  blocked <- grepl(
    "Checking your browser|reCAPTCHA",
    geo_text$text,
    ignore.case = TRUE
  )
  
  
  if (!blocked) {
    
    geo_text_success <- TRUE
    
    
    writeLines(
      geo_text$text,
      file.path(
        raw_dir,
        "GSE287274_GEO_accession_text.txt"
      )
    )
    
    
    geo_lines <- strsplit(
      geo_text$text,
      "\n",
      fixed = TRUE
    )[[1]]
    
    
    interesting <- geo_lines[
      grepl(
        "^\\^SERIES|^!Series_|^\\^SAMPLE|^!Sample_",
        geo_lines
      )
    ]
    
    
    cat(
      "\nRelevant GEO text lines:\n\n"
    )
    
    
    print(
      interesting
    )
    
  } else {
    
    cat(
      "GEO viewer returned browser/reCAPTCHA page rather than metadata.\n"
    )
  }
}


# ============================================================
# 9. GEOQUERY: SERIES MATRIX
#
# GEOquery discovers available matrix files itself.
# We are NOT hardcoding the matrix filename.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GEOQUERY SERIES MATRIX ATTEMPT\n")
cat("============================================================\n\n")


geo_matrix_success <- FALSE


geo_matrix <- tryCatch(
  
  {
    
    obj <- GEOquery::getGEO(
      
      GEO =
        accession,
      
      GSEMatrix =
        TRUE,
      
      AnnotGPL =
        FALSE,
      
      getGPL =
        FALSE,
      
      destdir =
        raw_dir
    )
    
    
    geo_matrix_success <- TRUE
    
    obj
  },
  
  error = function(e) {
    
    cat(
      "GEOquery series-matrix retrieval failed:\n",
      conditionMessage(e),
      "\n"
    )
    
    NULL
  }
)


if (geo_matrix_success) {
  
  cat(
    "GEOquery series matrix retrieval succeeded.\n"
  )
  
  
  cat(
    "Number of ExpressionSet objects:",
    length(geo_matrix),
    "\n"
  )
  
  
  for (i in seq_along(geo_matrix)) {
    
    eset <- geo_matrix[[i]]
    
    
    metadata <- Biobase::pData(
      eset
    )
    
    
    write.csv(
      
      metadata,
      
      file.path(
        raw_dir,
        paste0(
          "GSE287274_GEOquery_metadata_",
          i,
          ".csv"
        )
      ),
      
      row.names = TRUE
    )
    
    
    cat(
      "\nExpressionSet",
      i,
      "\n"
    )
    
    
    cat(
      "Features:",
      nrow(eset),
      "\n"
    )
    
    
    cat(
      "Samples:",
      ncol(eset),
      "\n"
    )
    
    
    cat(
      "Sample names:\n"
    )
    
    
    print(
      sampleNames(eset)
    )
  }
}


# ============================================================
# 10. GEOQUERY SUPPLEMENTARY FILE DISCOVERY
#
# This asks GEOquery to discover the current supplementary
# directory instead of us guessing filenames.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GEOQUERY SUPPLEMENTARY FILE ATTEMPT\n")
cat("============================================================\n\n")


supp_success <- FALSE


supp_result <- tryCatch(
  
  {
    
    result <- GEOquery::getGEOSuppFiles(
      
      GEO =
        accession,
      
      makeDirectory =
        FALSE,
      
      baseDir =
        raw_dir,
      
      fetch_files =
        FALSE
    )
    
    
    supp_success <- TRUE
    
    result
  },
  
  error = function(e) {
    
    cat(
      "Supplementary file discovery failed:\n",
      conditionMessage(e),
      "\n"
    )
    
    NULL
  }
)


if (supp_success) {
  
  cat(
    "Supplementary file listing succeeded.\n\n"
  )
  
  
  print(
    supp_result
  )
  
  
  write.csv(
    
    as.data.frame(
      supp_result
    ),
    
    file.path(
      raw_dir,
      "GSE287274_GEOquery_supplementary_listing.csv"
    ),
    
    row.names = TRUE
  )
}


# ============================================================
# 11. CHECK ORIGINAL SOFT ENDPOINT AGAIN
#
# Diagnostic only.
# We already know this returned 404.
# ============================================================

soft_url <- paste0(
  
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE287nnn/",
  accession,
  "/soft/",
  accession,
  "_family.soft.gz"
)


soft_test <- safe_text_request(
  soft_url,
  "Canonical GEO SOFT endpoint"
)


# ============================================================
# 12. FINAL ACCESS SUMMARY
# ============================================================

access_summary <- data.frame(
  
  route = c(
    
    "NCBI_GDS_esearch",
    
    "NCBI_GDS_esummary",
    
    "NCBI_SRA_esearch",
    
    "NCBI_SRA_esummary",
    
    "GEO_accession_text",
    
    "GEOquery_series_matrix",
    
    "GEOquery_supplementary_listing",
    
    "Canonical_SOFT_endpoint"
  ),
  
  success = c(
    
    length(gds_ids) > 0,
    
    gds_summary_success,
    
    length(sra_ids) > 0,
    
    sra_summary_success,
    
    geo_text_success,
    
    geo_matrix_success,
    
    supp_success,
    
    soft_test$success
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  access_summary,
  file.path(
    raw_dir,
    "GSE287274_access_recovery_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("ACCESS RECOVERY SUMMARY\n")
cat("============================================================\n\n")


print(
  access_summary,
  row.names = FALSE
)


# ============================================================
# 13. FILE INVENTORY
# ============================================================

local_files <- list.files(
  raw_dir,
  full.names = TRUE
)


inventory <- data.frame(
  
  filename =
    basename(local_files),
  
  size_MB =
    round(
      file.info(local_files)$size /
        1024^2,
      4
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  inventory,
  file.path(
    raw_dir,
    "GSE287274_access_file_inventory.csv"
  ),
  row.names = FALSE
)


cat(
  "\nLocal files:\n\n"
)


print(
  inventory,
  row.names = FALSE
)


# ============================================================
# 14. SESSION INFO
# ============================================================

sink(
  file.path(
    raw_dir,
    "GSE287274_access_recovery_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 15. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE287274 ACCESS AUDIT COMPLETE\n")
cat("============================================================\n\n")


cat(
  "No biological analysis was performed.\n"
)


cat(
  "No inclusion/exclusion decision was made automatically.\n\n"
)


cat(
  "Interpretation:\n",
  "  - If GEO metadata/matrix/supplementary files are recoverable,\n",
  "    we inspect the exact resistance contrast next.\n",
  "  - If only bibliographic indexing remains but no reproducible\n",
  "    GEO/SRA expression data can be retrieved, D9 should be\n",
  "    excluded from the primary analysis with that reason documented.\n"
)


cat("\n============================================================\n")