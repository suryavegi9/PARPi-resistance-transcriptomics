# ============================================================
# GSE342642 - DATASET 4 RAW COUNT QC
# No differential-expression testing yet
# ============================================================

library(edgeR)
library(limma)

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

dir.create(
  "results/GSE342642/QC",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load GEO table
# ------------------------------------------------------------

file_path <- file.choose()

dat <- read.delim(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Extract ONLY raw read-count columns
# ------------------------------------------------------------

count_cols <- grep(
  "_Read_Count$",
  colnames(dat),
  value = TRUE
)

counts <- as.matrix(
  dat[, count_cols]
)

storage.mode(counts) <- "integer"

rownames(counts) <- dat$Gene_ID

cat("\nNumber of genes:", nrow(counts), "\n")
cat("Number of samples:", ncol(counts), "\n\n")


# ------------------------------------------------------------
# 3. Build sample metadata
# ------------------------------------------------------------

sample_names <- sub(
  "_Read_Count$",
  "",
  count_cols
)

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
  region = region
)

rownames(metadata) <- sample_names

print(metadata)

write.csv(
  metadata,
  "results/GSE342642/QC/sample_metadata.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Raw library sizes
# ------------------------------------------------------------

library_sizes <- colSums(counts)

lib_table <- data.frame(
  sample = sample_names,
  condition = condition,
  tumor = tumor,
  region = region,
  library_size = library_sizes,
  library_size_millions = library_sizes / 1e6
)

cat("\n========================================\n")
cat("RAW LIBRARY SIZES\n")
cat("========================================\n")

print(lib_table)

write.csv(
  lib_table,
  "results/GSE342642/QC/library_sizes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Inspect suspicious lack of zeros
# ------------------------------------------------------------

cat("\n========================================\n")
cat("LOW-COUNT INSPECTION\n")
cat("========================================\n")

cat("Total zero counts:", sum(counts == 0), "\n")
cat("Total counts equal to 1:", sum(counts == 1), "\n")
cat("Total counts equal to 2:", sum(counts == 2), "\n")
cat("Total counts equal to 3:", sum(counts == 3), "\n")

cat(
  "Genes positive in all 12 samples:",
  sum(rowSums(counts > 0) == 12),
  "\n"
)

cat(
  "Genes with at least one count <= 5:",
  sum(apply(counts, 1, min) <= 5),
  "\n"
)


# ------------------------------------------------------------
# 6. edgeR object + TMM normalization
# ------------------------------------------------------------

dge <- DGEList(
  counts = counts,
  group = metadata$condition
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

norm_table <- data.frame(
  sample = sample_names,
  library_size = dge$samples$lib.size,
  norm_factor = dge$samples$norm.factors,
  effective_library_size =
    dge$samples$lib.size *
    dge$samples$norm.factors
)

cat("\n========================================\n")
cat("TMM NORMALIZATION\n")
cat("========================================\n")

print(norm_table)

write.csv(
  norm_table,
  "results/GSE342642/QC/TMM_normalization.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Expression filtering FOR QC VISUALIZATION ONLY
# ------------------------------------------------------------

# This is NOT yet our final manuscript DEG filter.
# We only remove extremely weak genes so that PCA/MDS
# are not dominated by near-zero noise.

cpm_matrix <- cpm(dge)

keep_qc <- rowSums(
  cpm_matrix >= 1
) >= 3

cat("\nGenes retained for QC:", sum(keep_qc), "\n")

dge_qc <- dge[
  keep_qc,
  ,
  keep.lib.sizes = FALSE
]

dge_qc <- calcNormFactors(
  dge_qc,
  method = "TMM"
)


# ------------------------------------------------------------
# 8. Log-CPM for exploratory QC
# ------------------------------------------------------------

logCPM <- cpm(
  dge_qc,
  log = TRUE,
  prior.count = 2
)


# ------------------------------------------------------------
# 9. PCA
# ------------------------------------------------------------

pca <- prcomp(
  t(logCPM),
  center = TRUE,
  scale. = FALSE
)

pca_var <- 100 *
  pca$sdev^2 /
  sum(pca$sdev^2)

pca_df <- data.frame(
  sample = sample_names,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = metadata$condition,
  tumor = metadata$tumor,
  region = metadata$region
)

cat("\n========================================\n")
cat("PCA\n")
cat("========================================\n")

cat(
  "PC1 variance:",
  round(pca_var[1], 2),
  "%\n"
)

cat(
  "PC2 variance:",
  round(pca_var[2], 2),
  "%\n"
)

print(pca_df)

write.csv(
  pca_df,
  "results/GSE342642/QC/PCA_coordinates.csv",
  row.names = FALSE
)

pdf(
  "results/GSE342642/QC/PCA.pdf",
  width = 8,
  height = 7
)

plot(
  pca_df$PC1,
  pca_df$PC2,
  pch = 19,
  xlab = paste0(
    "PC1 (",
    round(pca_var[1], 1),
    "%)"
  ),
  ylab = paste0(
    "PC2 (",
    round(pca_var[2], 1),
    "%)"
  ),
  main = "GSE342642 count-level QC PCA"
)

text(
  pca_df$PC1,
  pca_df$PC2,
  labels = paste0(
    pca_df$tumor,
    "-R",
    pca_df$region
  ),
  pos = 3,
  cex = 0.8
)

dev.off()


# ------------------------------------------------------------
# 10. Sample correlations
# ------------------------------------------------------------

sample_cor <- cor(
  logCPM,
  method = "pearson"
)

cat("\n========================================\n")
cat("SAMPLE CORRELATIONS\n")
cat("========================================\n")

print(
  round(
    sample_cor,
    3
  )
)

write.csv(
  sample_cor,
  "results/GSE342642/QC/sample_correlations.csv"
)

pdf(
  "results/GSE342642/QC/correlation_heatmap.pdf",
  width = 9,
  height = 9
)

heatmap(
  sample_cor,
  symm = TRUE,
  margins = c(10, 10),
  main = "GSE342642 sample correlations"
)

dev.off()


# ------------------------------------------------------------
# 11. MDS
# ------------------------------------------------------------

pdf(
  "results/GSE342642/QC/MDS.pdf",
  width = 8,
  height = 7
)

plotMDS(
  dge_qc,
  labels = paste0(
    metadata$tumor,
    "-R",
    metadata$region
  ),
  main = "GSE342642 MDS"
)

dev.off()


cat("\n========================================\n")
cat("COUNT QC COMPLETE\n")
cat("========================================\n")