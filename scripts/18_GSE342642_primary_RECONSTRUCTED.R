# D4 primary blocked DE model: RECONSTRUCTION, not recovered historical code.
# Derived from archived 21_GSE342642_Hallmark_sensitivity.R, Part A,
# with all 12 samples retained. Compare against frozen primary files.
# Run only in validation mode. Never overwrites original outputs.

local({
  library(edgeR)
  library(limma)
  
  # Find project root; if needed select the original archived raw file.
  root <- Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd())
  raw_path <- file.path(root, "data_raw", "GSE342642_Data.txt.gz")
  if (!file.exists(raw_path)) {
    message("Select data_raw/GSE342642_Data.txt.gz in the ORIGINAL project")
    raw_path <- file.choose()
    root <- dirname(dirname(raw_path))
  }
  ref_dir <- file.path(root, "results/GSE342642/primary_DEG")
  ref_paths <- file.path(ref_dir, c("GSE342642_primary_DEG_results.csv",
                                    "tumor_aware_filter.csv",
                                    "GSE342642_moderated_t_ranking.csv"))
  meta_path <- file.path(root, "results/GSE342642/QC/sample_metadata.csv")
  stopifnot(all(file.exists(c(raw_path, ref_paths, meta_path))))
  stopifnot(unname(tools::md5sum(raw_path)) == "3e2797609fda01938848901111616414")
  before <- tools::md5sum(c(ref_paths, raw_path, meta_path))
  out <- file.path(root, "reproducibility_checks", "D4_primary_reconstruction",
                   format(Sys.time(), "%Y%m%d_%H%M%S"))
  stopifnot(!dir.exists(out))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  message("Validation outputs: ", out)
  
  # A. Load original count matrix and independently check sample identity.
  dat <- read.delim(raw_path, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE)
  count_cols <- grep("_Read_Count$", colnames(dat), value = TRUE)
  sample_names <- sub("_Read_Count$", "", count_cols)
  counts <- as.matrix(dat[, count_cols, drop = FALSE])
  stopifnot(nrow(dat) == 15623L, ncol(counts) == 12L,
            !anyNA(counts), all(is.finite(counts)),
            all(counts >= 0), all(counts == floor(counts)),
            !anyNA(dat$Gene_ID), !anyDuplicated(dat$Gene_ID))
  storage.mode(counts) <- "integer"
  rownames(counts) <- as.character(dat$Gene_ID)
  colnames(counts) <- sample_names
  metadata <- data.frame(
    sample = sample_names,
    condition = factor(ifelse(grepl("^Cont_", sample_names), "Control", "Resistant"),
                       levels = c("Control", "Resistant")),
    tumor = factor(sub("^(Cont|Resi)_([0-9]+).*", "\\2", sample_names)),
    row.names = sample_names
  )
  saved_meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(saved_meta$sample),
            setequal(saved_meta$sample, sample_names))
  saved_meta <- saved_meta[match(sample_names, saved_meta$sample), ]
  stopifnot(identical(as.character(metadata$condition), as.character(saved_meta$condition)),
            identical(as.character(metadata$tumor), as.character(saved_meta$tumor)),
            all(table(metadata$tumor) == 3L),
            all(table(metadata$condition) == 6L),
            "Cont_440_P3_2-3" %in% sample_names)
  
  # B. Tumor-aware expression filter, identical rule to archived sensitivity code.
  # normLibSizes is called as in script 21; no alternative fitting variants tried.
  dge <- DGEList(counts = counts)
  dge <- normLibSizes(dge, method = "TMM")
  cpm_mat <- cpm(dge, log = FALSE)
  med <- sapply(levels(metadata$tumor), function(tumor_id) {
    apply(cpm_mat[, metadata$tumor == tumor_id, drop = FALSE], 1, median)
  })
  colnames(med) <- paste0("Tumor_", levels(metadata$tumor))
  keep <- rowSums(med >= 1) >= 2
  filt <- data.frame(Gene_ID = rownames(counts), med, retained = keep,
                     check.names = FALSE)
  write.csv(filt, file.path(out, "reconstructed_filter.csv"), row.names = FALSE)
  saved_filter <- read.csv(ref_paths[2], stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(saved_filter$Gene_ID),
            setequal(as.character(saved_filter$Gene_ID), filt$Gene_ID))
  saved_filter <- saved_filter[match(filt$Gene_ID, as.character(saved_filter$Gene_ID)), ]
  median_error <- max(abs(as.matrix(filt[, colnames(med)]) -
                            as.matrix(saved_filter[, colnames(med)])))
  filter_match <- identical(filt$retained, saved_filter$retained)
  message("Retained genes: ", sum(keep), "; frozen filter match: ", filter_match)
  message("Max absolute median CPM difference: ", format(median_error, digits = 10))
  if (!filter_match || sum(keep) != 14277L) {
    stop("Filtering differs from frozen result. Inspect reconstructed_filter.csv; do not fit a different model.")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- normLibSizes(dge, method = "TMM")
  
  # C. Reconstruct candidate two-pass blocked voom model from script 21.
  # The exact original primary script is missing; this is a testable candidate.
  design <- model.matrix(~ condition, data = metadata)
  stopifnot(identical(rownames(design), colnames(dge)))
  v1 <- voom(dge, design, plot = FALSE)
  cor1 <- duplicateCorrelation(v1, design, block = metadata$tumor)
  v2 <- voom(dge, design, plot = FALSE, block = metadata$tumor,
             correlation = cor1$consensus.correlation)
  cor2 <- duplicateCorrelation(v2, design, block = metadata$tumor)
  fit <- lmFit(v2, design, block = metadata$tumor,
               correlation = cor2$consensus.correlation)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "conditionResistant", number = Inf,
                 sort.by = "none", adjust.method = "BH")
  rec <- data.frame(Gene_ID = rownames(tt), tt, check.names = FALSE)
  rec$status <- ifelse(rec$adj.P.Val < .05 & rec$logFC >= 1, "Up",
                       ifelse(rec$adj.P.Val < .05 & rec$logFC <= -1, "Down", "Not_significant"))
  write.csv(rec, file.path(out, "reconstructed_DE.csv"), row.names = FALSE)
  
  # D. Compare the full gene-level output, not merely DEG counts.
  old <- read.csv(ref_paths[1], stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(old$Gene_ID), !anyDuplicated(rec$Gene_ID),
            setequal(as.character(old$Gene_ID), rec$Gene_ID))
  old <- old[match(rec$Gene_ID, as.character(old$Gene_ID)), ]
  fields <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  stopifnot(all(is.finite(as.matrix(rec[, fields]))),
            all(is.finite(as.matrix(old[, fields]))))
  # Numerical tolerance is an audit criterion, not a biological threshold.
  tol <- 1e-7
  metrics <- do.call(rbind, lapply(fields, function(f) {
    a <- old[[f]]; b <- rec[[f]]
    data.frame(field = f, max_abs_difference = max(abs(a-b)),
               spearman = cor(a,b,method="spearman"),
               all_within_tolerance = all(abs(a-b) <= tol * pmax(1,abs(a))))
  }))
  old_rank <- read.csv(ref_paths[3], stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(old_rank$Gene_ID),
            setequal(as.character(old_rank$Gene_ID), rec$Gene_ID))
  old_rank <- old_rank[match(rec$Gene_ID, as.character(old_rank$Gene_ID)), ]
  ranking_match <- all(abs(old_rank$t-rec$t) <= tol * pmax(1,abs(old_rank$t)))
  rho_target <- 0.39477637734058
  rho_match <- abs(cor2$consensus.correlation - rho_target) <= tol
  fdr_match <- identical(rec$adj.P.Val < .05, old$adj.P.Val < .05)
  status_match <- identical(rec$status, old$status)
  sign_match <- all(sign(rec$logFC) == sign(old$logFC))
  n_fdr <- sum(rec$adj.P.Val < .05)
  n_up <- sum(rec$status == "Up")
  n_down <- sum(rec$status == "Down")
  pass <- filter_match && median_error <= tol && rho_match &&
    all(metrics$all_within_tolerance) && ranking_match && fdr_match &&
    status_match && sign_match && n_fdr == 1984L && n_up == 588L && n_down == 453L
  stopifnot(identical(before, tools::md5sum(names(before))))
  write.csv(metrics, file.path(out, "numeric_comparison.csv"), row.names = FALSE)
  report <- c(
    "D4 PRIMARY MODEL RECONSTRUCTION (not historical source)",
    paste("Genes retained:", sum(keep), "(expected 14277)"),
    paste("Exact retained gene membership:", filter_match),
    paste("Max absolute median CPM difference:", format(median_error, digits = 12)),
    paste("Pass 1 correlation:", format(cor1$consensus.correlation, digits = 15)),
    paste("Pass 2 correlation:", format(cor2$consensus.correlation, digits = 15)),
    paste("Archived correlation:", rho_target),
    paste("FDR-significant:", n_fdr, "(expected 1984)"),
    paste("Large-effect up/down:", n_up, "/", n_down, "(expected 588 / 453)"),
    paste("Identical significance membership:", fdr_match),
    paste("Identical large-effect classification:", status_match),
    paste("Identical logFC direction:", sign_match),
    paste("Ranking statistics within tolerance:", ranking_match),
    "Numeric comparison:", capture.output(print(metrics, row.names = FALSE)),
    paste("Audit result:", if (pass) "PASS within stated numerical tolerance" else "NOT VERIFIED — investigate differences"),
    "Original files unchanged (MD5 verified). No GSEA executed."
  )
  writeLines(report, file.path(out, "D4_reconstruction_report.txt"))
  cat(paste(report, collapse = "\n"), "\n")
  cat("Report directory:", out, "\n")
})
