
# ============================================================
# VERIFY_FINAL_62_63.R
# FINAL SYNTHESIS REPRODUCIBILITY AUDIT — CORRECTED VERSION
#
# - Repairs 13 malformed indexing expressions in script 62
# - Executes scripts 62 and 63 in isolated directories
# - Compares outputs against frozen publication results
# - Preserves all original files
# - No differential expression or GSEA execution
# ============================================================

(function() {
  
  # ==========================================================
  # 1. LOCATE ORIGINAL PROJECT
  # ==========================================================
  
  root <- Sys.getenv(
    "PARPI_PROJECT_ROOT",
    unset = getwd()
  )
  
  old_rel <- paste0(
    "results/FINAL_D1_D6_integration/",
    "FINAL_D1_D6_Hallmark_master.csv"
  )
  
  if (!file.exists(file.path(root, old_rel))) {
    
    message(
      "Select the ORIGINAL FINAL_D1_D6_Hallmark_master.csv"
    )
    
    selected <- file.choose()
    
    stopifnot(
      basename(selected) == basename(old_rel)
    )
    
    root <- dirname(dirname(dirname(selected)))
  }
  
  root <- normalizePath(
    root,
    winslash = "/",
    mustWork = TRUE
  )
  
  previous_wd <- getwd()
  
  on.exit(
    setwd(previous_wd),
    add = TRUE
  )
  
  # ==========================================================
  # 2. VERIFY ORIGINAL INPUT FINGERPRINTS
  # ==========================================================
  
  expected <- c(
    
    "62_FINAL_D1_D8_accession_integration.R" =
      "d351f4247fe8cb2e79ff0b0c40ab5d84",
    
    "63_FINAL_FIGURE8_D1_D8.R" =
      "6000be909dbd63135a2e1a6e075ea67a",
    
    "results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv" =
      "bc0c690745af3d35ba38a94edbe046a8",
    
    "results/GSE165914/D7_analysis/D7_accession_Hallmark_integration.csv" =
      "e9db9a3ba741833509e979359b112f04",
    
    "results/GSE165914/D7_analysis/D7A_M468_Hallmark_GSEA.csv" =
      "63ce51c7d04e788e24a15a08007e50f3",
    
    "results/GSE165914/D7_analysis/D7B_SUM1315_Hallmark_GSEA.csv" =
      "f076b4ba25f43fad068dcaa40a937c58",
    
    "results/GSE163854/D8_analysis/D8_primary_Hallmark_GSEA.csv" =
      "1a9aa04c12a9dc1fec24440f8b21aa35",
    
    "results/GSE202434/D5_integration/D5_H3_P3_Hallmark_master_comparison.csv" =
      "b98ab4154b4dff63f843c8c6ddfb5072"
  )
  
  fig_rel <- paste0(
    "results/FINAL_D1_D8_integration/",
    "Figure8_FINAL"
  )
  
  ref_hash <- c(
    
    "Figure8_accession_source_table.csv" =
      "f55f438c94215d3491382d70a7b0d92d",
    
    "Figure8_FINAL_checkpoint.csv" =
      "859498a42398247531adc7b60a6fff0d",
    
    "Figure8_frozen_hypothesis_heatmap_source.csv" =
      "d06adf51e987d577e931a85256d5fe1b",
    
    "Figure8_model_level_source_table.csv" =
      "d7f1e70d7dacf67504a9b9df0310de9e",
    
    "Figure8_treatment_state_source_table.csv" =
      "75d8a1ddb9aeb75f9655c97851bb21ed"
  )
  
  ref_expected <- setNames(
    unname(ref_hash),
    file.path(fig_rel, names(ref_hash))
  )
  
  all_expected <- c(expected, ref_expected)
  
  orig <- file.path(
    root,
    names(all_expected)
  )
  
  before <- tools::md5sum(orig)
  
  ok <- !is.na(before) &
    unname(before) == unname(all_expected)
  
  if (!all(ok)) {
    
    print(
      data.frame(
        path = names(all_expected)[!ok],
        expected = unname(all_expected)[!ok],
        actual = unname(before)[!ok]
      )
    )
    
    stop(
      paste(
        "Missing or non-identical archived files.",
        "Restore exact archive copies before testing."
      )
    )
  }
  
  message(
    "All ", length(all_expected),
    " original file fingerprints verified."
  )
  
  # Check original files even if execution fails
  
  on.exit({
    
    unchanged_on_exit <- identical(
      before,
      tools::md5sum(orig)
    )
    
    message(
      "Original files unchanged: ",
      unchanged_on_exit
    )
    
  }, add = TRUE)
  
  # ==========================================================
  # 3. CHECK REQUIRED PACKAGES
  # ==========================================================
  
  required_packages <- c(
    "ggplot2",
    "patchwork"
  )
  
  available <- vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  
  if (!all(available)) {
    
    stop(
      "Missing packages: ",
      paste(
        required_packages[!available],
        collapse = ", "
      )
    )
  }
  
  # ==========================================================
  # 4. CREATE ISOLATED VALIDATION DIRECTORY
  # ==========================================================
  
  out <- file.path(
    root,
    "reproducibility_checks",
    "FINAL_62_63_execution",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    )
  )
  
  stopifnot(!dir.exists(out))
  
  dir.create(
    out,
    recursive = TRUE
  )
  
  write.csv(
    
    data.frame(
      path = names(all_expected),
      MD5 = unname(before)
    ),
    
    file.path(
      out,
      "input_fingerprints.csv"
    ),
    
    row.names = FALSE
  )
  
  message(
    "Isolated test directory: ",
    out
  )
  
  # ==========================================================
  # 5. FUNCTION TO EXECUTE ISOLATED SCRIPTS
  # ==========================================================
  
  run_script <- function(
    script,
    stage,
    input_paths
  ) {
    
    dir.create(
      stage,
      recursive = TRUE
    )
    
    # Copy only required original inputs
    
    for (p in input_paths) {
      
      target <- file.path(
        stage,
        p
      )
      
      dir.create(
        dirname(target),
        recursive = TRUE,
        showWarnings = FALSE
      )
      
      stopifnot(
        
        file.copy(
          file.path(root, p),
          target
        ),
        
        unname(
          tools::md5sum(target)
        ) == expected[[p]]
      )
    }
    
    # Read source without altering original file
    
    script_text <- readLines(
      file.path(root, script),
      warn = FALSE
    )
    
    script_text <- paste(
      script_text,
      collapse = "\n"
    )
    
    # --------------------------------------------------------
    # SYNTAX REPAIR — SCRIPT 62 ONLY
    # --------------------------------------------------------
    
    if (startsWith(script, "62_")) {
      
      pattern <- paste0(
        "(old_raw|d7_raw|d8_raw)",
        "\\[\\s*\\[\\s*",
        "((?:old|d7|d8)_cols\\[\"[^\"]+\"\\])",
        "\\s*\\]\\s*\\]"
      )
      
      matches <- gregexpr(
        pattern,
        script_text,
        perl = TRUE
      )
      
      original_matches <- regmatches(
        script_text,
        matches
      )[[1]]
      
      n_repairs <- length(
        original_matches
      )
      
      # Exactly 13 malformed indexing expressions
      # were identified in the archived script.
      
      stopifnot(n_repairs == 13L)
      
      repaired_text <- gsub(
        pattern,
        "\\1[[\\2]]",
        script_text,
        perl = TRUE
      )
      
      # Preserve a record of every syntax repair
      
      repaired_matches <- gsub(
        pattern,
        "\\1[[\\2]]",
        original_matches,
        perl = TRUE
      )
      
      write.csv(
        
        data.frame(
          repair_number = seq_len(n_repairs),
          original = original_matches,
          corrected = repaired_matches
        ),
        
        file.path(
          out,
          "script62_syntax_repair_log.csv"
        ),
        
        row.names = FALSE
      )
      
      script_text <- repaired_text
      
      writeLines(
        
        script_text,
        
        file.path(
          out,
          "62_SYNTAX_REPAIRED.R"
        )
      )
      
      message(
        "Script 62: repaired ",
        n_repairs,
        " malformed indexing expressions."
      )
    }
    
    # --------------------------------------------------------
    # PARSE SCRIPT
    # --------------------------------------------------------
    
    expr <- parse(
      text = script_text,
      keep.source = FALSE
    )
    
    message(
      script,
      ": syntax validation passed."
    )
    
    # --------------------------------------------------------
    # REDIRECT WORKING DIRECTORY
    # --------------------------------------------------------
    
    idx <- which(
      
      vapply(
        
        expr,
        
        function(e) {
          
          is.call(e) &&
            identical(
              e[[1]],
              as.name("setwd")
            )
        },
        
        logical(1)
      )
    )
    
    stopifnot(
      length(idx) == 1L,
      idx == 1L
    )
    
    expr[[idx]] <- substitute(
      setwd(P),
      list(P = stage)
    )
    
    env <- new.env(
      parent = globalenv()
    )
    
    setwd(stage)
    
    on.exit(
      setwd(root),
      add = TRUE
    )
    
    message(
      "\nExecuting: ",
      script,
      "\n"
    )
    
    eval(
      expr,
      envir = env
    )
    
    stopifnot(
      
      identical(
        
        normalizePath(
          getwd(),
          winslash = "/"
        ),
        
        normalizePath(
          stage,
          winslash = "/"
        )
      )
    )
    
    message(
      script,
      ": execution completed."
    )
    
    invisible(NULL)
  }
  
  # ==========================================================
  # 6. EXECUTE SCRIPTS 62 AND 63
  # ==========================================================
  
  inputs <- names(expected)[-(1:2)]
  
  a <- file.path(
    out,
    "script62"
  )
  
  b <- file.path(
    out,
    "script63"
  )
  
  # Script 62: final accession integration
  
  run_script(
    
    names(expected)[1],
    
    a,
    
    inputs[c(1, 2, 5)]
  )
  
  # Script 63: final synthesis figure
  
  run_script(
    
    names(expected)[2],
    
    b,
    
    inputs[c(1, 3, 4, 5, 6)]
  )
  
  # ==========================================================
  # 7. TABLE COMPARISON FUNCTIONS
  # ==========================================================
  
  read_tab <- function(p) {
    
    read.csv(
      p,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  compare <- function(
    x,
    y,
    keys,
    label
  ) {
    
    stopifnot(
      
      identical(
        names(x),
        names(y)
      ),
      
      nrow(x) == nrow(y),
      
      all(keys %in% names(x))
    )
    
    make_key <- function(d) {
      
      do.call(
        
        paste,
        
        c(
          unname(d[keys]),
          sep = "|"
        )
      )
    }
    
    kx <- make_key(x)
    ky <- make_key(y)
    
    stopifnot(
      
      !anyDuplicated(kx),
      
      !anyDuplicated(ky),
      
      setequal(kx, ky)
    )
    
    # Align rows by pathway or composite key
    
    y <- y[
      match(kx, ky),
      ,
      drop = FALSE
    ]
    
    do.call(
      
      rbind,
      
      lapply(
        
        names(x),
        
        function(n) {
          
          u <- x[[n]]
          v <- y[[n]]
          
          same_missing <- identical(
            is.na(u),
            is.na(v)
          )
          
          both <- !is.na(u) &
            !is.na(v)
          
          if (
            is.numeric(u) &&
            is.numeric(v)
          ) {
            
            err <- abs(
              u[both] - v[both]
            )
            
            # Integer vote counts require exact equality.
            
            tol <- if (
              grepl(
                "^(final_n_|model12_n_)",
                n
              )
            ) {
              0
            } else {
              1e-12
            }
            
            pass <- same_missing &&
              all(is.finite(err)) &&
              all(err <= tol)
            
            delta <- if (length(err)) {
              
              max(err)
              
            } else {
              
              0
            }
            
          } else {
            
            pass <- same_missing &&
              
              identical(
                as.character(u[both]),
                as.character(v[both])
              )
            
            delta <- NA_real_
          }
          
          data.frame(
            
            table = label,
            
            column = n,
            
            max_abs_difference = delta,
            
            pass = pass
          )
        }
      )
    )
  }
  
  # ==========================================================
  # 8. VERIFY ALL FIVE SCRIPT 63 SOURCE TABLES
  # ==========================================================
  
  metrics <- list()
  
  for (f in names(ref_hash)) {
    
    keys <- if (
      grepl("checkpoint", f)
    ) {
      
      "metric"
      
    } else if (
      grepl("heatmap", f)
    ) {
      
      c(
        "pathway",
        "dataset"
      )
      
    } else {
      
      "pathway"
    }
    
    regenerated <- read_tab(
      file.path(
        b,
        fig_rel,
        f
      )
    )
    
    original <- read_tab(
      file.path(
        root,
        fig_rel,
        f
      )
    )
    
    metrics[[f]] <- compare(
      
      regenerated,
      
      original,
      
      keys,
      
      f
    )
    
    message(
      
      f,
      
      ": ",
      
      sum(metrics[[f]]$pass),
      
      " / ",
      
      nrow(metrics[[f]]),
      
      " columns passed."
    )
  }
  
  # ==========================================================
  # 9. VERIFY SCRIPT 62 AGAINST FROZEN ACCESSION TABLE
  # ==========================================================
  
  # The complete historical script 62 output is unavailable.
  # Compare its 12 overlapping fields (including pathway ID)
  # against the independently archived script 63 table.
  
  m <- read_tab(
    
    file.path(
      
      a,
      
      paste0(
        "results/FINAL_D1_D8_integration/",
        "accession_level/",
        "FINAL_D1_D8_accession_Hallmark_master.csv"
      )
    )
  )
  
  ref <- read_tab(
    
    file.path(
      
      root,
      
      fig_rel,
      
      "Figure8_accession_source_table.csv"
    )
  )
  
  # Correct case-sensitive column mapping
  
  map <- c(
    
    pathway = "pathway",
    
    old_class = "old_class",
    
    final_class = "final_class",
    
    final_n_UP = "nUP",
    
    final_n_DOWN = "nDOWN",
    
    final_n_sig_UP = "nSigUP",
    
    final_n_sig_DOWN = "nSigDown",
    
    D7_NES = "D7_NES",
    
    D7_significant = "D7_significant",
    
    D8_NES = "D8_NES",
    
    D8_significant = "D8_significant"
  )
  
  stopifnot(
    
    all(unname(map) %in% names(m)),
    
    identical(
      names(map),
      names(ref)
    )
  )
  
  sub <- m[
    ,
    unname(map),
    drop = FALSE
  ]
  
  names(sub) <- names(map)
  
  metrics[["62_vs_frozen_63"]] <- compare(
    
    sub,
    
    ref,
    
    "pathway",
    
    "script62_vs_archived_accession_table"
  )
  
  # ==========================================================
  # 10. FINAL COMPARISON RESULTS
  # ==========================================================
  
  all_metrics <- do.call(
    rbind,
    metrics
  )
  
  rownames(all_metrics) <- NULL
  
  unchanged <- identical(
    before,
    tools::md5sum(orig)
  )
  
  write.csv(
    
    all_metrics,
    
    file.path(
      out,
      "final_comparison.csv"
    ),
    
    row.names = FALSE
  )
  
  capture.output(
    
    sessionInfo(),
    
    file = file.path(
      out,
      "sessionInfo.txt"
    )
  )
  
  passed <- all(
    all_metrics$pass
  ) && unchanged
  
  # ==========================================================
  # 11. GENERATE FINAL AUDIT REPORT
  # ==========================================================
  
  report <- c(
    
    "FINAL SYNTHESIS EXECUTION AUDIT",
    
    paste(
      "Script 62:",
      "13 syntax repairs applied;",
      "original archived script unchanged."
    ),
    
    paste(
      "Scripts 62 and 63 executed",
      "in separate minimal project copies."
    ),
    
    "Script 63: five archived CSV outputs compared.",
    
    paste(
      "Script 62:",
      "12 overlapping accession fields compared",
      "with archived script 63 source table."
    ),
    
    paste(
      "Column comparisons passed:",
      sum(all_metrics$pass),
      "/",
      nrow(all_metrics)
    ),
    
    paste(
      "Frozen input and reference MD5 unchanged:",
      unchanged
    ),
    
    paste(
      "Audit result:",
      if (passed) "PASS" else "NOT VERIFIED"
    ),
    
    paste(
      "No DE or GSEA executed.",
      "Figure appearance still requires visual review."
    )
  )
  
  writeLines(
    
    report,
    
    file.path(
      out,
      "FINAL_execution_report.txt"
    )
  )
  
  cat(
    "\n",
    paste(
      report,
      collapse = "\n"
    ),
    "\n"
  )
  
  cat(
    "\nReport directory:",
    out,
    "\n"
  )
  
  if (!passed) {
    
    print(
      all_metrics[
        !all_metrics$pass,
        ,
        drop = FALSE
      ]
    )
    
    stop(
      paste(
        "Inspect final_comparison.csv.",
        "Do not change frozen results."
      )
    )
  }
  
})()
