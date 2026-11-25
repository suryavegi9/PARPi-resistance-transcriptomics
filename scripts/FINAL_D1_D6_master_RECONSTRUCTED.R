# Historical D1-D6 Hallmark master RECONSTRUCTION, not recovered source.
# Rebuild from eleven saved Hallmark tables. No DE or GSEA execution.
# Run only in a separate validation directory; frozen files are read-only.

local({
  root <- Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd())
  rel_ref <- "results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv"
  ref <- file.path(root, rel_ref)
  if (!file.exists(ref)) {
    message("Select the ORIGINAL FINAL_D1_D6_Hallmark_master.csv")
    ref <- file.choose()
    stopifnot(basename(ref) == "FINAL_D1_D6_Hallmark_master.csv")
    root <- dirname(dirname(dirname(ref)))
    stopifnot(file.exists(file.path(root, rel_ref)))
  }
  source_paths <- c(
    D1_def = "results/pathway_analysis/Hallmark_GSEA_BRCA1_deficient_all.csv",
    D1_prof = "results/pathway_analysis/Hallmark_GSEA_BRCA1_proficient_all.csv",
    D2 = "results/GSE117765/GSE117765_Hallmark_GSEA_all.csv",
    D3 = "results/GSE153867/GSE153867_Hallmark_GSEA_all.csv",
    D4 = "results/GSE342642/Hallmark_GSEA/GSE342642_Hallmark_GSEA_all.csv",
    H1 = "results/GSE202434/HCC1806/Hallmark_GSEA/HCC1806_H1_Hallmark_GSEA_all.csv",
    H2 = "results/GSE202434/HCC1806/Hallmark_GSEA/HCC1806_H2_Hallmark_GSEA_all.csv",
    P1 = "results/GSE202434/PSN1/Hallmark_GSEA/PSN1_P1_Hallmark_GSEA_all.csv",
    P1_sens = "results/GSE202434/PSN1/Hallmark_GSEA/PSN1_P1_Hallmark_GSEA_without_maintain2.csv",
    D6 = "results/GSE191231/Hallmark_GSEA/GSE191231_Hallmark_GSEA_primary_all.csv",
    D6_sens = "results/GSE191231/Hallmark_GSEA/GSE191231_Hallmark_GSEA_sensitivity_all.csv"
  )
  files <- setNames(file.path(root, source_paths), names(source_paths))
  stopifnot(all(file.exists(c(files, ref))))
  before <- tools::md5sum(c(files, ref))
  out <- file.path(root, "reproducibility_checks", "D1_D6_master_reconstruction",
                   format(Sys.time(), "%Y%m%d_%H%M%S"))
  stopifnot(!dir.exists(out))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  
  # 1. Load and align original pathway-level inputs by Hallmark identifier.
  tabs <- lapply(files, function(p) {
    d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("ID", "NES", "p.adjust") %in% names(d)),
              nrow(d) == 50L, !anyNA(d$ID), !anyDuplicated(d$ID),
              all(grepl("^HALLMARK_", d$ID)),
              all(is.finite(d$NES)), all(is.finite(d$p.adjust)),
              all(d$p.adjust >= 0 & d$p.adjust <= 1))
    d
  })
  keys <- sort(tabs[[1]]$ID)
  rec <- data.frame(pathway = keys, stringsAsFactors = FALSE)
  for (k in names(tabs)) {
    d <- tabs[[k]]
    stopifnot(setequal(d$ID, keys))
    d <- d[match(keys, d$ID), , drop = FALSE]
    rec[[paste0(k, "_NES")]] <- d$NES
    rec[[paste0(k, "_padj")]] <- d$p.adjust
  }
  
  # 2. Conservative votes, with sensitivity annotations kept separate.
  vote <- function(keys) {
    s <- sign(as.matrix(rec[, paste0(keys, "_NES"), drop = FALSE]))
    p <- as.matrix(rec[, paste0(keys, "_padj"), drop = FALSE])
    direction <- ifelse(rowSums(s == s[, 1]) == length(keys), s[, 1], 0)
    significant <- ifelse(rowSums(p < 0.05) == length(keys), direction, 0)
    list(direction = direction, significant = significant)
  }
  add_vote <- function(d, key, components) {
    v <- vote(components)
    d[[paste0(key, "_direction")]] <- v$direction
    d[[paste0(key, "_sig_direction")]] <- v$significant
    d
  }
  
  rec$D1_NES <- rowMeans(rec[, c("D1_def_NES", "D1_prof_NES")])
  rec <- add_vote(rec, "D1", c("D1_def", "D1_prof"))
  for (k in c("D2", "D3", "D4")) rec <- add_vote(rec, k, k)
  
  # Historical D5 descriptive mean includes the three primary contrasts.
  # Formal vote additionally requires the P1 sensitivity result to agree.
  rec$D5_NES <- rowMeans(rec[, c("H1_NES", "H2_NES", "P1_NES")])
  rec <- add_vote(rec, "D5", c("H1", "H2", "P1", "P1_sens"))
  
  # D6 primary vote is not replaced by its sensitivity significance.
  rec <- add_vote(rec, "D6", "D6")
  rec$D6_sensitivity_same_direction <- sign(rec$D6_NES) == sign(rec$D6_sens_NES)
  rec$D6_robust_significant <- rec$D6_sensitivity_same_direction &
    rec$D6_padj < 0.05 & rec$D6_sens_padj < 0.05
  
  # 3. Historical six-accession classification, not the final D1-D8 result.
  dirs <- as.matrix(rec[, paste0("D", 1:6, "_direction")])
  sigs <- as.matrix(rec[, paste0("D", 1:6, "_sig_direction")])
  rec$n_UP <- rowSums(dirs == 1)
  rec$n_DOWN <- rowSums(dirs == -1)
  rec$n_no_formal_direction <- rowSums(dirs == 0)
  rec$n_sig_UP <- rowSums(sigs == 1)
  rec$n_sig_DOWN <- rowSums(sigs == -1)
  rec$n_sig_total <- rec$n_sig_UP + rec$n_sig_DOWN
  rec$max_direction_support <- pmax(rec$n_UP, rec$n_DOWN)
  rec$max_sig_support <- pmax(rec$n_sig_UP, rec$n_sig_DOWN)
  rec$majority_direction <- ifelse(rec$n_UP > rec$n_DOWN, "UP",
                                   ifelse(rec$n_DOWN > rec$n_UP, "DOWN", "TIE"))
  
  context <- rec$n_sig_UP > 0 & rec$n_sig_DOWN > 0
  strong <- (rec$n_UP >= 4 & rec$n_sig_UP >= 3 & rec$n_sig_DOWN == 0) |
    (rec$n_DOWN >= 4 & rec$n_sig_DOWN >= 3 & rec$n_sig_UP == 0)
  directional <- (rec$n_UP >= 4 & rec$n_sig_DOWN == 0) |
    (rec$n_DOWN >= 4 & rec$n_sig_UP == 0)
  
  rec$classification <- ifelse(context, "Context-dependent",
                               ifelse(strong, "Strong recurrent",
                                      ifelse(directional, "Directionally recurrent",
                                             ifelse(rec$n_sig_total == 1, "Dataset-specific significant", "Inconsistent/weak"))))
  
  # 4. Only now read the archived master for comparison.
  archived <- read.csv(ref, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(nrow(archived) == 50L, ncol(archived) == 49L,
            !anyDuplicated(archived$pathway), setequal(archived$pathway, keys),
            identical(names(rec), names(archived)))
  archived <- archived[match(rec$pathway, archived$pathway), , drop = FALSE]
  tol <- 1e-12
  metrics <- do.call(rbind, lapply(names(rec), function(k) {
    a <- archived[[k]]
    b <- rec[[k]]
    if (is.numeric(a) && is.numeric(b)) {
      stopifnot(all(is.finite(a)), all(is.finite(b)))
      err <- abs(a - b)
      # Integer vote counts must be exactly equal.
      exact_required <- !grepl("(_NES|_padj)$", k)
      ok <- if (exact_required) a == b else err <= tol
      data.frame(field = k, type = if (exact_required) "exact_numeric" else "numeric",
                 max_abs_difference = max(err),
                 n_mismatch = sum(!ok), pass = all(ok))
    } else {
      stopifnot(!anyNA(a), !anyNA(b))
      ok <- as.character(a) == as.character(b)
      data.frame(field = k, type = "exact", max_abs_difference = NA_real_,
                 n_mismatch = sum(!ok), pass = all(ok))
    }
  }))
  classes <- c("Strong recurrent", "Directionally recurrent", "Context-dependent",
               "Dataset-specific significant", "Inconsistent/weak")
  counts <- table(factor(rec$classification, levels = classes))
  original_unchanged <- identical(before, tools::md5sum(names(before)))
  passed <- all(metrics$pass) && identical(as.integer(counts), c(1L,12L,15L,12L,10L)) &&
    original_unchanged
  
  write.csv(rec, file.path(out, "reconstructed_Hallmark_master.csv"), row.names = FALSE)
  write.csv(metrics, file.path(out, "column_comparison.csv"), row.names = FALSE)
  write.csv(data.frame(classification = classes, n = as.integer(counts)),
            file.path(out, "classification_counts.csv"), row.names = FALSE)
  # Store relative paths rather than workstation-specific paths.
  write.csv(data.frame(path = c(unname(source_paths), rel_ref), MD5 = unname(before)),
            file.path(out, "input_fingerprints.csv"), row.names = FALSE)
  report <- c(
    "HISTORICAL D1-D6 MASTER RECONSTRUCTION (not historical source)",
    paste("Source Hallmark tables:", length(tabs)),
    paste("Pathways:", nrow(rec)),
    paste("Columns:", ncol(rec)),
    paste("Matching columns:", sum(metrics$pass), "/", nrow(metrics)),
    paste("Maximum absolute numerical difference:", format(max(metrics$max_abs_difference, na.rm = TRUE), digits = 12)),
    paste("Absolute numerical tolerance for NES/FDR:", tol),
    paste("Original input and reference files unchanged (MD5):", original_unchanged),
    "Classification counts:", capture.output(print(counts)),
    paste("Audit result:", if (passed) "PASS within stated numerical tolerance" else "NOT VERIFIED"),
    "No differential expression or GSEA executed. No frozen output overwritten."
  )
  writeLines(report, file.path(out, "D1_D6_reconstruction_report.txt"))
  cat(paste(report, collapse = "\n"), "\n")
  cat("Report directory:", out, "\n")
  if (!passed) stop("Reconstruction did not pass. Inspect column_comparison.csv.")
})
