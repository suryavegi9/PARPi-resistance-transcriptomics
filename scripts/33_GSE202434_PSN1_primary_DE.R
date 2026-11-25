# ============================================================
# GSE202434 - DATASET 5
# PSN1 PRIMARY DIFFERENTIAL EXPRESSION ANALYSIS
#
# Supplied matrix behaves as log2(x + 1)
#
# P1 = Maintain  vs Parental
# P2 = Withdraw  vs Parental
# P3 = Maintain  vs Withdraw
#
# IMPORTANT:
# P3 represents drug withdrawal + resistance reversion.
# It is NOT a pure exposure-state contrast like HCC1806 H3.
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)
library(ggplot2)


# ============================================================
# 1. LOAD DATA
# ============================================================

file_path <- paste0(
  "data_raw/",
  "GSE202434_Exp43_Tala_Normalized_Counts_ExpressedGenesONLY_ProteinCoding_ONLY.csv.gz"
)

dat <- read.csv(
  gzfile(file_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(dat)[1] <- "gene"


# ============================================================
# 2. IDENTIFIER QC
# ============================================================

cat("\n========================================\n")
cat("GENE SYMBOL QC\n")
cat("========================================\n\n")

cat(
  "Genes:",
  nrow(dat),
  "\n"
)

cat(
  "Duplicated gene symbols:",
  sum(duplicated(dat$gene)),
  "\n"
)

cat(
  "Missing gene symbols:",
  sum(is.na(dat$gene) | dat$gene == ""),
  "\n"
)

stopifnot(
  sum(duplicated(dat$gene)) == 0
)

stopifnot(
  sum(is.na(dat$gene) | dat$gene == "") == 0
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
      "maintain",
      sample_names,
      ignore.case = TRUE
    ),
    
    "TalaR_M",
    
    ifelse(
      
      grepl(
        "withdraw",
        sample_names,
        ignore.case = TRUE
      ),
      
      "Withdraw",
      
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
# 5. BACK-TRANSFORM FOR FILTERING
#
# Supplied values behave as:
#
# log2(x + 1)
#
# Therefore:
#
# x = 2^expression - 1
#
# We use x only for expression filtering.
# ============================================================

expr_back <- 2^expr - 1


# Keep genes with back-transformed expression >= 1
# in at least 3 of the 9 samples.
#
# This allows expression confined to one biological group
# to remain in the analysis.

keep <- rowSums(
  expr_back >= 1
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
    "Withdraw",
    "TalaR_M"
  )
)


design <- model.matrix(
  ~ 0 + group_factor
)


colnames(design) <- c(
  "Parental",
  "Withdraw",
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
  
  P1_M_vs_Parental =
    TalaR_M - Parental,
  
  P2_WD_vs_Parental =
    Withdraw - Parental,
  
  P3_M_vs_WD =
    TalaR_M - Withdraw,
  
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


P1 <- get_result(
  "P1_M_vs_Parental"
)

P2 <- get_result(
  "P2_WD_vs_Parental"
)

P3 <- get_result(
  "P3_M_vs_WD"
)


# ============================================================
# 10. DEG FLAGS
# ============================================================

add_deg_flag <- function(x) {
  
  x$DEG <-
    x$adj.P.Val < 0.05 &
    abs(x$logFC) >= 1
  
  x$direction <- ifelse(
    
    x$DEG &
      x$logFC >= 1,
    
    "UP",
    
    ifelse(
      
      x$DEG &
        x$logFC <= -1,
      
      "DOWN",
      
      "NS"
    )
  )
  
  x
}


P1 <- add_deg_flag(P1)
P2 <- add_deg_flag(P2)
P3 <- add_deg_flag(P3)


# ============================================================
# 11. DEG SUMMARY
# ============================================================

summary_table <- data.frame(
  
  contrast = c(
    "P1_M_vs_Parental",
    "P2_WD_vs_Parental",
    "P3_M_vs_WD"
  ),
  
  FDR_lt_0.05 = c(
    sum(P1$adj.P.Val < 0.05),
    sum(P2$adj.P.Val < 0.05),
    sum(P3$adj.P.Val < 0.05)
  ),
  
  DEG_FDR_and_absLogFC1 = c(
    sum(P1$DEG),
    sum(P2$DEG),
    sum(P3$DEG)
  ),
  
  UP = c(
    sum(P1$direction == "UP"),
    sum(P2$direction == "UP"),
    sum(P3$direction == "UP")
  ),
  
  DOWN = c(
    sum(P1$direction == "DOWN"),
    sum(P2$direction == "DOWN"),
    sum(P3$direction == "DOWN")
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
# 12. P3 TOP GENES
#
# P3 = maintained resistant vs withdrawn/reversion state
# ============================================================

P3_top <- P3[
  order(
    P3$adj.P.Val,
    -abs(P3$logFC)
  ),
]


cat("\n========================================\n")
cat("P3: Maintain vs Withdraw - TOP 20\n")
cat("========================================\n\n")

print(
  head(
    P3_top,
    20
  ),
  row.names = FALSE
)


# ============================================================
# 13. OUTPUTS
# ============================================================

outdir <- "results/GSE202434/PSN1"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  metadata,
  file.path(
    outdir,
    "PSN1_metadata.csv"
  ),
  row.names = FALSE
)


write.csv(
  pca_data,
  file.path(
    outdir,
    "PSN1_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


write.csv(
  summary_table,
  file.path(
    outdir,
    "PSN1_DE_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  P1,
  file.path(
    outdir,
    "PSN1_P1_M_vs_Parental_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  P2,
  file.path(
    outdir,
    "PSN1_P2_WD_vs_Parental_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  P3,
  file.path(
    outdir,
    "PSN1_P3_M_vs_WD_all.csv"
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
    title = "GSE202434 PSN1 expression PCA",
    
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
    "PSN1_PCA.png"
  ),
  p,
  width = 8,
  height = 6,
  dpi = 300
)


cat("\n========================================\n")
cat("PSN1 PRIMARY DE COMPLETE\n")
cat("========================================\n")