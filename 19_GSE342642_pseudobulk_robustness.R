# ============================================================
# GSE342642 - DATASET 4
# TUMOR-LEVEL PSEUDOBULK ROBUSTNESS ANALYSIS
#
# Goal:
# Collapse 3 spatial regions within each tumor by SUMMING
# raw counts, then compare 2 control vs 2 resistant tumors.
#
# This is a SENSITIVITY ANALYSIS, not the primary model.
# ============================================================


library(edgeR)
library(limma)


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/GSE342642/pseudobulk"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load original GEO file
# ------------------------------------------------------------

geo_file <- file.choose()

dat <- read.delim(
  geo_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Extract read counts
# ------------------------------------------------------------

count_cols <- grep(
  "_Read_Count$",
  colnames(dat),
  value = TRUE
)

sample_names <- sub(
  "_Read_Count$",
  "",
  count_cols
)

counts <- as.matrix(
  dat[, count_cols]
)

storage.mode(counts) <- "integer"

colnames(counts) <- sample_names
rownames(counts) <- as.character(dat$Gene_ID)


# ------------------------------------------------------------
# 3. Recover tumor identity
# ------------------------------------------------------------

tumor <- sub(
  "^(Cont|Resi)_([0-9]+).*",
  "\\2",
  sample_names
)

condition <- ifelse(
  tumor %in% c("437", "440"),
  "Control",
  "Resistant"
)


# ------------------------------------------------------------
# 4. Build tumor-level pseudobulk counts
#
# SUM counts from all three regions belonging
# to the same biological tumor.
# ------------------------------------------------------------

tumor_ids <- c(
  "437",
  "440",
  "433",
  "438"
)

pb_counts <- sapply(
  tumor_ids,
  function(t) {
    
    idx <- tumor == t
    
    rowSums(
      counts[, idx, drop = FALSE]
    )
  }
)

colnames(pb_counts) <- paste0(
  "Tumor_",
  tumor_ids
)

rownames(pb_counts) <- rownames(counts)


cat("\n========================================\n")
cat("PSEUDOBULK LIBRARIES\n")
cat("========================================\n\n")

print(
  colSums(pb_counts) / 1e6
)


# ------------------------------------------------------------
# 5. Pseudobulk metadata
# ------------------------------------------------------------

pb_meta <- data.frame(
  tumor = tumor_ids,
  
  condition = factor(
    c(
      "Control",
      "Control",
      "Resistant",
      "Resistant"
    ),
    levels = c(
      "Control",
      "Resistant"
    )
  )
)

rownames(pb_meta) <- colnames(pb_counts)

print(pb_meta)


# ------------------------------------------------------------
# 6. Load PRIMARY DEG result
#
# We intentionally use the SAME genes tested in
# the primary blocked analysis so that the two
# logFC estimates are directly comparable.
# ------------------------------------------------------------

primary_file <- paste0(
  "results/GSE342642/primary_DEG/",
  "GSE342642_primary_DEG_results.csv"
)

primary <- read.csv(
  primary_file,
  stringsAsFactors = FALSE
)

tested_genes <- as.character(
  primary$Gene_ID
)

common_genes <- intersect(
  tested_genes,
  rownames(pb_counts)
)

pb_counts_f <- pb_counts[
  common_genes,
  ,
  drop = FALSE
]

cat("\nGenes used for pseudobulk comparison:",
    nrow(pb_counts_f),
    "\n")


# ------------------------------------------------------------
# 7. edgeR object + TMM normalization
# ------------------------------------------------------------

dge_pb <- DGEList(
  counts = pb_counts_f,
  group = pb_meta$condition
)

dge_pb <- normLibSizes(
  dge_pb,
  method = "TMM"
)


cat("\n========================================\n")
cat("PSEUDOBULK TMM NORMALIZATION\n")
cat("========================================\n\n")

print(
  dge_pb$samples
)


# ------------------------------------------------------------
# 8. Design matrix
# ------------------------------------------------------------

design_pb <- model.matrix(
  ~ condition,
  data = pb_meta
)

print(design_pb)


# ------------------------------------------------------------
# 9. voom
# ------------------------------------------------------------

pdf(
  file.path(
    outdir,
    "pseudobulk_voom.pdf"
  ),
  width = 7,
  height = 6
)

v_pb <- voom(
  dge_pb,
  design_pb,
  plot = TRUE
)

dev.off()


# ------------------------------------------------------------
# 10. Fit tumor-level model
# ------------------------------------------------------------

fit_pb <- lmFit(
  v_pb,
  design_pb
)

fit_pb <- eBayes(
  fit_pb
)


pb_res <- topTable(
  fit_pb,
  coef = "conditionResistant",
  number = Inf,
  sort.by = "none"
)

pb_res$Gene_ID <- rownames(pb_res)


# ------------------------------------------------------------
# 11. Merge primary and pseudobulk results
# ------------------------------------------------------------

compare <- merge(
  primary[
    ,
    c(
      "Gene_ID",
      "Gene_Symbol",
      "logFC",
      "t",
      "adj.P.Val",
      "status"
    )
  ],
  pb_res[
    ,
    c(
      "Gene_ID",
      "logFC",
      "t",
      "P.Value",
      "adj.P.Val"
    )
  ],
  by = "Gene_ID",
  suffixes = c(
    "_primary",
    "_pseudobulk"
  )
)


# ------------------------------------------------------------
# 12. Direction concordance
# ------------------------------------------------------------

compare$same_direction <-
  sign(compare$logFC_primary) ==
  sign(compare$logFC_pseudobulk)


all_direction_concordance <- mean(
  compare$same_direction
)


# Primary DEG subset
primary_deg <- compare[
  compare$status != "Not_significant",
]

deg_direction_concordance <- mean(
  primary_deg$same_direction
)


# ------------------------------------------------------------
# 13. Effect-size correlation
# ------------------------------------------------------------

rho_all <- cor(
  compare$logFC_primary,
  compare$logFC_pseudobulk,
  method = "spearman"
)

rho_deg <- cor(
  primary_deg$logFC_primary,
  primary_deg$logFC_pseudobulk,
  method = "spearman"
)


# ------------------------------------------------------------
# 14. Individual resistant-tumor consistency
#
# Convert pseudobulk counts to normalized logCPM.
#
# For each gene:
# compare tumor 433 and tumor 438 individually
# against the mean of control tumors 437 and 440.
# ------------------------------------------------------------

logCPM_pb <- cpm(
  dge_pb,
  log = TRUE,
  prior.count = 2
)

control_mean <- rowMeans(
  logCPM_pb[
    ,
    c(
      "Tumor_437",
      "Tumor_440"
    ),
    drop = FALSE
  ]
)

diff_433 <-
  logCPM_pb[, "Tumor_433"] -
  control_mean

diff_438 <-
  logCPM_pb[, "Tumor_438"] -
  control_mean


tumor_consistency <- data.frame(
  Gene_ID = rownames(logCPM_pb),
  diff_433 = diff_433,
  diff_438 = diff_438,
  same_resistant_direction =
    sign(diff_433) ==
    sign(diff_438)
)


compare <- merge(
  compare,
  tumor_consistency,
  by = "Gene_ID"
)


# ------------------------------------------------------------
# 15. Primary DEG consistency across BOTH
# resistant tumors
# ------------------------------------------------------------

primary_deg2 <- compare[
  compare$status != "Not_significant",
]

resistant_tumor_consistency <- mean(
  primary_deg2$same_resistant_direction
)


# ------------------------------------------------------------
# 16. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("PSEUDOBULK ROBUSTNESS SUMMARY\n")
cat("========================================\n\n")

cat(
  "Genes compared:",
  nrow(compare),
  "\n"
)

cat(
  "Spearman rho - all genes:",
  round(rho_all, 4),
  "\n"
)

cat(
  "Same direction - all genes:",
  round(
    100 * all_direction_concordance,
    2
  ),
  "%\n"
)

cat(
  "\nPrimary DEGs:",
  nrow(primary_deg),
  "\n"
)

cat(
  "Spearman rho - primary DEGs:",
  round(rho_deg, 4),
  "\n"
)

cat(
  "Same direction - primary DEGs:",
  round(
    100 * deg_direction_concordance,
    2
  ),
  "%\n"
)

cat(
  "Primary DEGs with same direction in both",
  "resistant tumors:",
  round(
    100 * resistant_tumor_consistency,
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 17. Save results
# ------------------------------------------------------------

write.csv(
  compare,
  file.path(
    outdir,
    "primary_vs_pseudobulk_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  logCPM_pb,
  file.path(
    outdir,
    "tumor_pseudobulk_logCPM.csv"
  )
)


# ------------------------------------------------------------
# 18. Effect-size comparison plot
# ------------------------------------------------------------

pdf(
  file.path(
    outdir,
    "primary_vs_pseudobulk_logFC.pdf"
  ),
  width = 7,
  height = 7
)

plot(
  compare$logFC_primary,
  compare$logFC_pseudobulk,
  pch = 16,
  cex = 0.35,
  xlab = "Primary blocked-model log2FC",
  ylab = "Tumor-pseudobulk log2FC",
  main = paste0(
    "Primary vs pseudobulk effect sizes\nSpearman rho = ",
    round(rho_all, 3)
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

abline(
  h = 0,
  v = 0,
  lty = 3
)

dev.off()


# ------------------------------------------------------------
# 19. Save summary text
# ------------------------------------------------------------

summary_lines <- c(
  
  "GSE342642 PSEUDOBULK ROBUSTNESS ANALYSIS",
  
  "",
  
  paste(
    "Genes compared:",
    nrow(compare)
  ),
  
  paste(
    "Spearman rho all genes:",
    rho_all
  ),
  
  paste(
    "Direction concordance all genes:",
    all_direction_concordance
  ),
  
  paste(
    "Primary DEG count:",
    nrow(primary_deg)
  ),
  
  paste(
    "Spearman rho primary DEGs:",
    rho_deg
  ),
  
  paste(
    "Direction concordance primary DEGs:",
    deg_direction_concordance
  ),
  
  paste(
    "Primary DEG same-direction consistency",
    "across resistant tumors:",
    resistant_tumor_consistency
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "pseudobulk_robustness_summary.txt"
  )
)


cat("\n========================================\n")
cat("PSEUDOBULK ANALYSIS COMPLETE\n")
cat("========================================\n")