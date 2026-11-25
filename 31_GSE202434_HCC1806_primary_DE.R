# ============================================================
# GSE202434 - DATASET 5
# HCC1806 PRIMARY DIFFERENTIAL EXPRESSION ANALYSIS
#
# Input:
# supplied log2(TPM + 1) matrix
#
# H1 = TalaR-M  vs Parental
# H2 = TalaR-DF vs Parental
# H3 = TalaR-M  vs TalaR-DF
#
# H3 is the PRIMARY treatment-state contrast.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(readxl)
library(limma)
library(ggplot2)


# ============================================================
# 1. LOAD DATA
# ============================================================

file_path <- paste0(
  "data_raw/",
  "GSE202434_HCC1806_Talazoparib_TPM_ExpressedGenesONLY_ProteinCoding_ONLY.xlsx"
)

dat <- read_excel(
  file_path,
  sheet = 1
)

dat <- as.data.frame(
  dat,
  check.names = FALSE
)

colnames(dat)[1] <- "gene"


# ============================================================
# 2. REMOVE AMBIGUOUS EXCEL-CORRUPTED GENE SYMBOLS
#
# 43891 / 43892 correspond to Excel date serial artefacts.
# Their original gene identities cannot be uniquely recovered.
# ============================================================

bad_symbols <- c(
  "43891",
  "43892"
)

n_before <- nrow(dat)

dat <- dat[
  !(dat$gene %in% bad_symbols),
]


cat("\n========================================\n")
cat("GENE SYMBOL QC\n")
cat("========================================\n\n")

cat(
  "Rows before removing ambiguous symbols:",
  n_before,
  "\n"
)

cat(
  "Rows after removing ambiguous symbols:",
  nrow(dat),
  "\n"
)

cat(
  "Rows removed:",
  n_before - nrow(dat),
  "\n"
)

cat(
  "Remaining duplicated gene symbols:",
  sum(duplicated(dat$gene)),
  "\n"
)

stopifnot(
  sum(duplicated(dat$gene)) == 0
)


# ============================================================
# 3. EXPRESSION MATRIX
# ============================================================

expr <- as.matrix(
  dat[, -1]
)

rownames(expr) <- dat$gene

storage.mode(expr) <- "numeric"


# ============================================================
# 4. SAMPLE METADATA
# ============================================================

sample_names <- colnames(expr)

group <- ifelse(
  grepl(
    "parental",
    sample_names,
    ignore.case = TRUE
  ),
  "Parental",
  ifelse(
    grepl(
      "50nM",
      sample_names,
      ignore.case = TRUE
    ),
    "TalaR_M",
    ifelse(
      grepl(
        "TalaR",
        sample_names,
        ignore.case = TRUE
      ),
      "TalaR_DF",
      NA
    )
  )
)


metadata <- data.frame(
  sample = sample_names,
  group = group,
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("SAMPLE METADATA\n")
cat("========================================\n\n")

print(
  metadata,
  row.names = FALSE
)


stopifnot(
  !any(is.na(metadata$group))
)

stopifnot(
  all(
    table(metadata$group) == 3
  )
)


# ============================================================
# 5. BACK-TRANSFORM FOR FILTERING ONLY
#
# Expression supplied by GEO:
# log2(TPM + 1)
#
# TPM = 2^x - 1
# ============================================================

tpm <- 2^expr - 1


# Keep genes with TPM >= 1 in at least 3 samples.
#
# This allows a gene expressed consistently in any one
# biological condition to remain in the analysis.

keep <- rowSums(
  tpm >= 1
) >= 3


cat("\n========================================\n")
cat("EXPRESSION FILTER\n")
cat("========================================\n\n")

cat(
  "Genes before filtering:",
  nrow(expr),
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


expr_f <- expr[
  keep,
  ,
  drop = FALSE
]


# ============================================================
# 6. PCA QC
# ============================================================

pca <- prcomp(
  t(expr_f),
  center = TRUE,
  scale. = FALSE
)


variance <- 100 *
  pca$sdev^2 /
  sum(pca$sdev^2)


pca_data <- data.frame(
  sample = rownames(pca$x),
  group = metadata$group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)


cat("\n========================================\n")
cat("PCA SUMMARY\n")
cat("========================================\n\n")

cat(
  "PC1 variance:",
  round(variance[1], 2),
  "%\n"
)

cat(
  "PC2 variance:",
  round(variance[2], 2),
  "%\n\n"
)

print(
  pca_data,
  row.names = FALSE
)


# ============================================================
# 7. DESIGN MATRIX
# ============================================================

group_factor <- factor(
  metadata$group,
  levels = c(
    "Parental",
    "TalaR_DF",
    "TalaR_M"
  )
)


design <- model.matrix(
  ~ 0 + group_factor
)


colnames(design) <- c(
  "Parental",
  "TalaR_DF",
  "TalaR_M"
)


cat("\n========================================\n")
cat("DESIGN MATRIX\n")
cat("========================================\n\n")

print(design)


# ============================================================
# 8. LIMMA MODEL
# ============================================================

fit <- lmFit(
  expr_f,
  design
)


contrasts <- makeContrasts(
  
  H1_M_vs_Parental =
    TalaR_M - Parental,
  
  H2_DF_vs_Parental =
    TalaR_DF - Parental,
  
  H3_M_vs_DF =
    TalaR_M - TalaR_DF,
  
  levels = design
)


fit2 <- contrasts.fit(
  fit,
  contrasts
)

fit2 <- eBayes(
  fit2,
  trend = TRUE
)


# ============================================================
# 9. EXTRACT COMPLETE RESULTS
# ============================================================

get_result <- function(coef_name) {
  
  res <- topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  
  res$gene <- rownames(res)
  
  res <- res[
    ,
    c(
      "gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]
  
  res
}


H1 <- get_result(
  "H1_M_vs_Parental"
)

H2 <- get_result(
  "H2_DF_vs_Parental"
)

H3 <- get_result(
  "H3_M_vs_DF"
)


# ============================================================
# 10. DEG FLAGS
# ============================================================

add_deg_flag <- function(x) {
  
  x$DEG <- (
    x$adj.P.Val < 0.05 &
      abs(x$logFC) >= 1
  )
  
  x$direction <- ifelse(
    x$DEG & x$logFC >= 1,
    "UP",
    ifelse(
      x$DEG & x$logFC <= -1,
      "DOWN",
      "NS"
    )
  )
  
  x
}


H1 <- add_deg_flag(H1)
H2 <- add_deg_flag(H2)
H3 <- add_deg_flag(H3)


# ============================================================
# 11. DEG SUMMARY
# ============================================================

summary_table <- data.frame(
  
  contrast = c(
    "H1_M_vs_Parental",
    "H2_DF_vs_Parental",
    "H3_M_vs_DF"
  ),
  
  FDR_lt_0.05 = c(
    sum(H1$adj.P.Val < 0.05),
    sum(H2$adj.P.Val < 0.05),
    sum(H3$adj.P.Val < 0.05)
  ),
  
  DEG_FDR_and_absLogFC1 = c(
    sum(H1$DEG),
    sum(H2$DEG),
    sum(H3$DEG)
  ),
  
  UP = c(
    sum(H1$direction == "UP"),
    sum(H2$direction == "UP"),
    sum(H3$direction == "UP")
  ),
  
  DOWN = c(
    sum(H1$direction == "DOWN"),
    sum(H2$direction == "DOWN"),
    sum(H3$direction == "DOWN")
  )
)


cat("\n========================================\n")
cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("========================================\n\n")

print(
  summary_table,
  row.names = FALSE
)


# ============================================================
# 12. H3 TOP GENES
#
# H3 is our primary treatment-state comparison.
# ============================================================

H3_top <- H3[
  order(
    H3$adj.P.Val,
    -abs(H3$logFC)
  ),
]


cat("\n========================================\n")
cat("H3: TalaR-M vs TalaR-DF - TOP 20\n")
cat("========================================\n\n")

print(
  head(
    H3_top,
    20
  ),
  row.names = FALSE
)


# ============================================================
# 13. OUTPUTS
# ============================================================

outdir <- "results/GSE202434/HCC1806"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  metadata,
  file.path(
    outdir,
    "HCC1806_metadata.csv"
  ),
  row.names = FALSE
)


write.csv(
  pca_data,
  file.path(
    outdir,
    "HCC1806_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


write.csv(
  summary_table,
  file.path(
    outdir,
    "HCC1806_DE_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  H1,
  file.path(
    outdir,
    "HCC1806_H1_M_vs_Parental_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  H2,
  file.path(
    outdir,
    "HCC1806_H2_DF_vs_Parental_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  H3,
  file.path(
    outdir,
    "HCC1806_H3_M_vs_DF_all.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. PCA FIGURE
# ============================================================

p <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = group,
    label = sample
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    nudge_y = 0.7,
    size = 3
  ) +
  
  labs(
    title = "GSE202434 HCC1806 expression PCA",
    x = paste0(
      "PC1 (",
      round(variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance[2], 1),
      "%)"
    ),
    colour = "Condition"
  ) +
  
  theme_classic()


ggsave(
  file.path(
    outdir,
    "HCC1806_PCA.png"
  ),
  p,
  width = 8,
  height = 6,
  dpi = 300
)


cat("\n========================================\n")
cat("HCC1806 PRIMARY DE COMPLETE\n")
cat("========================================\n")