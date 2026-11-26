# Portable release validation, distinct from the exact author-executed wrapper.
# Executes only frozen-input integration/plotting. Never performs DE or GSEA.
local({
  args <- commandArgs(trailingOnly = TRUE)
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_arg) != 1L) stop("Invoke with Rscript --vanilla.")
  script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
  root <- dirname(dirname(dirname(script_path)))
  option <- function(key, default) {
    x <- args[startsWith(args, paste0("--", key, "="))]
    if (length(x) > 1L) stop("Duplicate option: ", key)
    if (length(x)) substring(x, nchar(key) + 4L) else default
  }
  if (any(!grepl("^--output=", args))) stop("Only --output=NEW_DIRECTORY is supported.")
  output <- option("output", tempfile("PARPi_synthesis_", tmpdir = dirname(root)))
  parent <- dirname(output)
  if (!dir.exists(parent)) stop("Output parent must already exist.")
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  output <- file.path(parent, basename(output))
  if (tolower(parent) == tolower(root) || startsWith(tolower(parent), paste0(tolower(root), "/")))
    stop("Output must be outside the repository, to protect frozen results.")
  if (file.exists(output)) stop("Output already exists; refusing to overwrite.")
  required <- c("ggplot2", "patchwork")
  if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
    stop("Install ggplot2 and patchwork separately; see environment/README.md.")
  manifest <- read.csv(file.path(root, "provenance/portable_synthesis/INPUT_MANIFEST.csv"),
                       stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(manifest$path), all(!grepl("^/|[.][.]", manifest$path)))
  source_paths <- file.path(root, manifest$path)
  stopifnot(all(file.exists(source_paths)))
  before <- unname(tools::md5sum(source_paths))
  if (!identical(before, manifest$md5)) stop("Source fingerprint mismatch; do not change references.")
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  dir.create(output)
  # Capture at the end as well, after the plotting libraries have loaded.
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo_before.txt"))
  run_one <- function(name, selected) {
    stage <- file.path(output, name)
    dir.create(stage)
    for (relative in selected) {
      destination <- file.path(stage, relative)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(file.path(root, relative), destination, overwrite = FALSE)) stop("Copy failed")
      stopifnot(unname(tools::md5sum(destination)) == unname(tools::md5sum(file.path(root, relative))))
    }
    filename <- if (name == "script62") "62_FINAL_D1_D8_accession_integration_SYNTAX_REPAIRED.R" else "63_FINAL_FIGURE8_D1_D8.R"
    expr <- parse(file.path(root, "scripts", filename), keep.source = FALSE)
    is_setwd <- vapply(expr, function(e) is.call(e) && identical(e[[1]], as.name("setwd")), logical(1))
    stopifnot(sum(is_setwd) == 1L, which(is_setwd) == 1L)
    # Same isolation mechanism as the historical wrapper; no scientific expressions changed.
    expr[[1]] <- substitute(setwd(P), list(P = stage))
    setwd(stage)
    eval(expr, envir = new.env(parent = globalenv()))
    setwd(previous_wd)
    stage
  }
  old <- "results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv"
  d8 <- "results/GSE163854/D8_analysis/D8_primary_Hallmark_GSEA.csv"
  d7 <- "results/GSE165914/D7_analysis/"
  s62 <- run_one("script62", c(old, d8, paste0(d7, "D7_accession_Hallmark_integration.csv")))
  s63 <- run_one("script63", c(old, d8, paste0(d7, "D7A_M468_Hallmark_GSEA.csv"),
                              paste0(d7, "D7B_SUM1315_Hallmark_GSEA.csv"),
                              "results/GSE202434/D5_integration/D5_H3_P3_Hallmark_master_comparison.csv"))
  fig <- "results/FINAL_D1_D8_integration/Figure8_FINAL"
  tables <- c("Figure8_accession_source_table.csv", "Figure8_FINAL_checkpoint.csv",
              "Figure8_frozen_hypothesis_heatmap_source.csv", "Figure8_model_level_source_table.csv",
              "Figure8_treatment_state_source_table.csv")
  metrics <- list()
  align <- function(a, b, keys) {
    key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
    ka <- key(a); kb <- key(b)
    stopifnot(!anyDuplicated(ka), !anyDuplicated(kb), setequal(ka, kb))
    b[match(ka, kb), , drop = FALSE]
  }
  compare <- function(a, b, name, mapping) {
    for (column in names(mapping)) {
      x <- a[[column]]; y <- b[[mapping[[column]]]]
      stopifnot(length(x) == length(y))
      same_na <- identical(is.na(x), is.na(y))
      numeric <- is.numeric(x) && is.numeric(y)
      difference <- if (numeric && any(!is.na(x))) max(abs(x[!is.na(x)] - y[!is.na(x)])) else NA_real_
      passed <- same_na && if (numeric) (is.na(difference) || difference <= 1e-12) else
        identical(as.character(x), as.character(y))
      metrics[[length(metrics) + 1L]] <<- data.frame(table = name, column = column,
            max_abs_difference = difference, pass = passed, stringsAsFactors = FALSE)
    }
  }
  byte_records <- list()
  for (name in tables) {
    ref <- file.path(root, fig, name); generated <- file.path(s63, fig, name)
    a <- read.csv(ref, stringsAsFactors = FALSE, check.names = FALSE)
    b <- read.csv(generated, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(identical(names(a), names(b)))
    keys <- if (name == "Figure8_FINAL_checkpoint.csv") "metric" else if
      (name == "Figure8_frozen_hypothesis_heatmap_source.csv") c("pathway", "dataset") else "pathway"
    b <- align(a, b, keys)
    compare(a, b, name, setNames(names(a), names(a)))
    byte_records[[name]] <- data.frame(table = name,
      byte_identical = unname(tools::md5sum(ref)) == unname(tools::md5sum(generated)))
  }
  a <- read.csv(file.path(root, fig, tables[1]), stringsAsFactors = FALSE)
  b <- read.csv(file.path(s62, "results/FINAL_D1_D8_integration/accession_level/FINAL_D1_D8_accession_Hallmark_master.csv"), stringsAsFactors = FALSE)
  b <- align(a, b, "pathway")
  mapping <- c(pathway="pathway", old_class="old_class", final_class="final_class",
    final_n_UP="nUP", final_n_DOWN="nDOWN", final_n_sig_UP="nSigUP", final_n_sig_DOWN="nSigDown",
    D7_NES="D7_NES", D7_significant="D7_significant", D8_NES="D8_NES", D8_significant="D8_significant")
  compare(a, b, "script62_vs_archived_accession_table", mapping)
  comparison <- do.call(rbind, metrics)
  unchanged <- identical(before, unname(tools::md5sum(source_paths)))
  stopifnot(nrow(comparison) == 46L)
  write.csv(comparison, file.path(output, "final_comparison.csv"), row.names = FALSE)
  write.csv(do.call(rbind, byte_records), file.path(output, "script63_byte_comparison.csv"), row.names = FALSE)
  write.csv(data.frame(path=manifest$path, md5_before=before,
            md5_after=unname(tools::md5sum(source_paths))),file.path(output,"input_integrity.csv"),row.names=FALSE)
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  report <- c("Portable final synthesis verification", paste("Columns passed:",sum(comparison$pass),"/ 46"),
    "Script 63: 35 columns in five tables. Script 62: 11 overlapping columns.",
    paste("Frozen inputs and portable scripts unchanged:",unchanged),
    "Tolerance: absolute 1e-12 for numeric fields; exact IDs/categories/logicals and missingness.",
    "No DE or GSEA executed. Generated Figure 8 is validation artwork, not publication Figure 5.",
    "Complete historical script-62 master equivalence and raw-to-results reproduction are not claimed.")
  writeLines(report, file.path(output,"VALIDATION_REPORT.txt"))
  cat(paste(report,collapse="\n"),"\n")
  if (!all(comparison$pass) || !unchanged) stop("Validation failed; references must remain unchanged.")
})
