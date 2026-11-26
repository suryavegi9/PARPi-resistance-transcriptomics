# ============================================================
# GSE342642 - DATASET 4
# PRIMARY DIFFERENTIAL EXPRESSION ANALYSIS
#
# Biological question:
# Which genes differ between acquired olaparib-resistant
# and control PDX tumors while accounting for repeated
# spatial samples from the same tumor?
#
# PRIMARY ANALYSIS ONLY
# ============================================================


library(edgeR)
library(limma)


# ------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/GSE342642/primary_DEG"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load GEO file
# ------------------------------------------------------------

file_path <- file.choose()

dat <- read.delim(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Extract RAW READ COUNTS only
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
# 3. Gene annotation
# ------------------------------------------------------------

gene_annotation <- data.frame(
  Gene_ID = as.character(dat$Gene_ID),
  Gene_Symbol = dat$Gene_Symbol,
  Description = dat$Description,
  gene_biotype = dat$gene_biotype,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Build sample metadata
# ------------------------------------------------------------

condition <- ifelse(
  grepl("^Cont_", sample_names),
  "Control",
  "Resistant"
)

tumor <- sub(
  "^(Cont|Resi)_([0-9]+).*",
  "\\2",
  sample_names
)

region <- sub(
  ".*-",
  "",
  sample_names
)

metadata <- data.frame(
  sample = sample_names,
  
  condition = factor(
    condition,
    levels = c("Control", "Resistant")
  ),
  
  tumor = factor(tumor),
  
  region = region,
  
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample


# Safety check
stopifnot(
  identical(
    colnames(counts),
    rownames(metadata)
  )
)

cat("\n========================================\n")
cat("SAMPLE DESIGN\n")
cat("========================================\n\n")

print(metadata)


# ------------------------------------------------------------
# 5. Create edgeR object
# ------------------------------------------------------------

dge <- DGEList(
  counts = counts,
  genes = gene_annotation
)

dge <- normLibSizes(
  dge,
  method = "TMM"
)


# ------------------------------------------------------------
# 6. Tumor-aware expression filter
# ------------------------------------------------------------

# Important:
# We do NOT count three regions from one tumor as
# three independent confirmations of expression.
#
# For every tumor:
# calculate median CPM across its three regions.
#
# Keep gene if median CPM >= 1
# in at least 2 of the 4 independent tumors.

cpm_all <- cpm(
  dge,
  log = FALSE
)

tumor_levels <- levels(
  metadata$tumor
)

tumor_median_cpm <- sapply(
  tumor_levels,
  function(t) {
    
    sample_index <- metadata$tumor == t
    
    apply(
      cpm_all[, sample_index, drop = FALSE],
      1,
      median
    )
  }
)

colnames(tumor_median_cpm) <- paste0(
  "Tumor_",
  tumor_levels
)

keep <- rowSums(
  tumor_median_cpm >= 1
) >= 2


cat("\n========================================\n")
cat("TUMOR-AWARE FILTER\n")
cat("========================================\n\n")

cat(
  "Genes before filtering:",
  nrow(dge),
  "\n"
)

cat(
  "Genes retained:",
  sum(keep),
  "\n"
)

cat(
  "Genes removed:",
  sum(!keep),
  "\n"
)


filter_table <- data.frame(
  Gene_ID = rownames(counts),
  tumor_median_cpm,
  retained = keep
)

write.csv(
  filter_table,
  file.path(
    outdir,
    "tumor_aware_filter.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Apply filter
# ------------------------------------------------------------

dge_f <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

dge_f <- normLibSizes(
  dge_f,
  method = "TMM"
)


# ------------------------------------------------------------
# 8. Design matrix
# ------------------------------------------------------------

# Model:
#
# expression ~ condition
#
# Tumor is NOT included as a fixed-effect coefficient
# because each tumor belongs to only one condition.
#
# Tumor will instead enter as a BLOCKING variable.

design <- model.matrix(
  ~ condition,
  data = metadata
)

cat("\n========================================\n")
cat("DESIGN MATRIX\n")
cat("========================================\n\n")

print(design)

cat(
  "\nDesign rank:",
  qr(design)$rank,
  "of",
  ncol(design),
  "columns\n"
)

if (
  qr(design)$rank != ncol(design)
) {
  stop(
    "Design matrix is not full rank."
  )
}


# ------------------------------------------------------------
# 9. First voom pass
# ------------------------------------------------------------

pdf(
  file.path(
    outdir,
    "voom_pass1.pdf"
  ),
  width = 7,
  height = 6
)

v1 <- voom(
  dge_f,
  design = design,
  plot = TRUE
)

dev.off()


# ------------------------------------------------------------
# 10. Estimate within-tumor correlation
# ------------------------------------------------------------

corfit1 <- duplicateCorrelation(
  v1,
  design = design,
  block = metadata$tumor
)

cat("\n========================================\n")
cat("DUPLICATE CORRELATION - PASS 1\n")
cat("========================================\n\n")

cat(
  "Consensus within-tumor correlation:",
  corfit1$consensus.correlation,
  "\n"
)


# ------------------------------------------------------------
# 11. Second voom pass
#
# Now voom knows that samples from the same tumor
# are correlated.
# ------------------------------------------------------------

pdf(
  file.path(
    outdir,
    "voom_pass2_blocked.pdf"
  ),
  width = 7,
  height = 6
)

v2 <- voom(
  dge_f,
  design = design,
  plot = TRUE,
  block = metadata$tumor,
  correlation = corfit1$consensus.correlation
)

dev.off()


# ------------------------------------------------------------
# 12. Re-estimate correlation after updated voom weights
# ------------------------------------------------------------

corfit2 <- duplicateCorrelation(
  v2,
  design = design,
  block = metadata$tumor
)

cat("\n========================================\n")
cat("DUPLICATE CORRELATION - PASS 2\n")
cat("========================================\n\n")

cat(
  "Final consensus within-tumor correlation:",
  corfit2$consensus.correlation,
  "\n"
)


# ------------------------------------------------------------
# 13. Fit blocked linear model
# ------------------------------------------------------------

fit <- lmFit(
  v2,
  design = design,
  block = metadata$tumor,
  correlation = corfit2$consensus.correlation
)

fit <- eBayes(
  fit
)


# ------------------------------------------------------------
# 14. Extract Resistant vs Control result
# ------------------------------------------------------------

coef_name <- "conditionResistant"

res <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "P"
)

res$Gene_ID <- rownames(res)


# Add gene annotations back
idx <- match(
  res$Gene_ID,
  gene_annotation$Gene_ID
)

res$Gene_Symbol <- gene_annotation$Gene_Symbol[idx]
res$Description <- gene_annotation$Description[idx]
res$gene_biotype <- gene_annotation$gene_biotype[idx]


# Reorder columns
res <- res[
  ,
  c(
    "Gene_ID",
    "Gene_Symbol",
    "Description",
    "gene_biotype",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]


# ------------------------------------------------------------
# 15. Apply LOCKED DEG definition
#
# BH FDR < 0.05
# AND
# |log2FC| >= 1
# ------------------------------------------------------------

res$status <- "Not_significant"

res$status[
  res$adj.P.Val < 0.05 &
    res$logFC >= 1
] <- "Up"

res$status[
  res$adj.P.Val < 0.05 &
    res$logFC <= -1
] <- "Down"


n_up <- sum(
  res$status == "Up"
)

n_down <- sum(
  res$status == "Down"
)

n_sig_fdr <- sum(
  res$adj.P.Val < 0.05
)


cat("\n========================================\n")
cat("PRIMARY DEG RESULTS\n")
cat("========================================\n\n")

cat(
  "Genes tested:",
  nrow(res),
  "\n"
)

cat(
  "BH FDR < 0.05:",
  n_sig_fdr,
  "\n"
)

cat(
  "Up: FDR <0.05 and log2FC >=1:",
  n_up,
  "\n"
)

cat(
  "Down: FDR <0.05 and log2FC <=-1:",
  n_down,
  "\n"
)


# ------------------------------------------------------------
# 16. Save complete DEG table
# ------------------------------------------------------------

write.csv(
  res,
  file.path(
    outdir,
    "GSE342642_primary_DEG_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Save ranked statistics for FUTURE GSEA
#
# We are NOT running GSEA yet.
# ------------------------------------------------------------

rank_table <- res[
  order(
    res$t,
    decreasing = TRUE
  ),
  c(
    "Gene_ID",
    "Gene_Symbol",
    "t",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
]

write.csv(
  rank_table,
  file.path(
    outdir,
    "GSE342642_moderated_t_ranking.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Display top genes
# ------------------------------------------------------------

cat("\n========================================\n")
cat("TOP 15 GENES BY P-VALUE\n")
cat("========================================\n\n")

print(
  head(
    res[
      ,
      c(
        "Gene_Symbol",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
        "status"
      )
    ],
    15
  )
)


# ------------------------------------------------------------
# 19. Save analysis summary
# ------------------------------------------------------------

summary_lines <- c(
  "GSE342642 PRIMARY ANALYSIS",
  "",
  paste(
    "Genes before filtering:",
    nrow(dge)
  ),
  paste(
    "Genes after tumor-aware filtering:",
    nrow(dge_f)
  ),
  paste(
    "Consensus within-tumor correlation:",
    corfit2$consensus.correlation
  ),
  paste(
    "Genes BH FDR < 0.05:",
    n_sig_fdr
  ),
  paste(
    "Up DEGs FDR<0.05 and logFC>=1:",
    n_up
  ),
  paste(
    "Down DEGs FDR<0.05 and logFC<=-1:",
    n_down
  )
)

writeLines(
  summary_lines,
  file.path(
    outdir,
    "primary_analysis_summary.txt"
  )
)


cat("\n========================================\n")
cat("PRIMARY ANALYSIS COMPLETE\n")
cat("========================================\n")