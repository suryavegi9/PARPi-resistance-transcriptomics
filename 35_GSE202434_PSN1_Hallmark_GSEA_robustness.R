# ============================================================
# GSE202434 - DATASET 5
# PSN1 HALLMARK GSEA + MAINTAIN-2 SENSITIVITY
#
# PRIMARY:
# all 9 samples
#
# SENSITIVITY:
# exclude PSN1-maintain-2
#
# P1 = Maintain vs Parental
# P2 = Withdraw vs Parental
# P3 = Maintain vs Withdraw
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(limma)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(tidyr)


# ============================================================
# 1. LOAD PRIMARY DE RESULTS
# ============================================================

indir <- "results/GSE202434/PSN1"

P1 <- read.csv(
  file.path(
    indir,
    "PSN1_P1_M_vs_Parental_all.csv"
  ),
  stringsAsFactors = FALSE
)

P2 <- read.csv(
  file.path(
    indir,
    "PSN1_P2_WD_vs_Parental_all.csv"
  ),
  stringsAsFactors = FALSE
)

P3 <- read.csv(
  file.path(
    indir,
    "PSN1_P3_M_vs_WD_all.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. HALLMARK DATABASE
# ============================================================

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

term2gene <- hallmark %>%
  select(
    gs_name,
    gene_symbol
  ) %>%
  distinct()


# ============================================================
# 3. GSEA FUNCTION
# ============================================================

run_gsea <- function(
    table,
    label
) {
  
  ranks <- table$t
  
  names(ranks) <- table$gene
  
  valid <-
    is.finite(ranks) &
    !is.na(names(ranks)) &
    names(ranks) != ""
  
  ranks <- ranks[valid]
  
  stopifnot(
    sum(duplicated(names(ranks))) == 0
  )
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  set.seed(123)
  
  x <- GSEA(
    geneList = ranks,
    TERM2GENE = term2gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    exponent = 1,
    eps = 1e-10,
    verbose = FALSE,
    by = "fgsea"
  )
  
  out <- as.data.frame(x)
  
  out$contrast <- label
  
  out
}


# ============================================================
# 4. PRIMARY GSEA
# ============================================================

gsea_P1 <- run_gsea(
  P1,
  "P1_M_vs_Parental"
)

gsea_P2 <- run_gsea(
  P2,
  "P2_WD_vs_Parental"
)

gsea_P3 <- run_gsea(
  P3,
  "P3_M_vs_WD"
)


cat("\n========================================\n")
cat("PRIMARY PSN1 GSEA\n")
cat("========================================\n\n")

cat(
  "P1 significant Hallmarks:",
  sum(gsea_P1$p.adjust < 0.05),
  "/ 50\n"
)

cat(
  "P2 significant Hallmarks:",
  sum(gsea_P2$p.adjust < 0.05),
  "/ 50\n"
)

cat(
  "P3 significant Hallmarks:",
  sum(gsea_P3$p.adjust < 0.05),
  "/ 50\n"
)


# ============================================================
# 5. REBUILD NO-MAINTAIN-2 EXPRESSION MODEL
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

expr <- as.matrix(
  dat[, -1]
)

rownames(expr) <- dat$gene

storage.mode(expr) <- "numeric"


# Same primary filter

back <- 2^expr - 1

keep <- rowSums(
  back >= 1
) >= 3

expr <- expr[
  keep,
  ,
  drop = FALSE
]


# Remove maintain-2 only for sensitivity analysis

expr_sens <- expr[
  ,
  colnames(expr) != "PSN1-maintain-2",
  drop = FALSE
]


samples <- colnames(expr_sens)


group <- ifelse(
  
  grepl(
    "parental",
    samples,
    ignore.case = TRUE
  ),
  
  "Parental",
  
  ifelse(
    
    grepl(
      "maintain",
      samples,
      ignore.case = TRUE
    ),
    
    "TalaR_M",
    
    "Withdraw"
  )
)


group <- factor(
  group,
  levels = c(
    "Parental",
    "Withdraw",
    "TalaR_M"
  )
)


design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- c(
  "Parental",
  "Withdraw",
  "TalaR_M"
)


fit <- lmFit(
  expr_sens,
  design
)


contrasts <- makeContrasts(
  
  P1_M_vs_Parental =
    TalaR_M - Parental,
  
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
# 6. EXTRACT SENSITIVITY RESULTS
# ============================================================

extract_result <- function(coef_name) {
  
  x <- topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  
  x$gene <- rownames(x)
  
  x
}


P1_sens <- extract_result(
  "P1_M_vs_Parental"
)

P3_sens <- extract_result(
  "P3_M_vs_WD"
)


# ============================================================
# 7. SENSITIVITY GSEA
# ============================================================

gsea_P1_sens <- run_gsea(
  P1_sens,
  "P1_without_maintain2"
)

gsea_P3_sens <- run_gsea(
  P3_sens,
  "P3_without_maintain2"
)


# ============================================================
# 8. GSEA ROBUSTNESS FUNCTION
# ============================================================

compare_gsea <- function(
    primary,
    sensitivity,
    contrast_name
) {
  
  x <- merge(
    
    primary[
      ,
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    
    sensitivity[
      ,
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    
    by = "ID",
    
    suffixes = c(
      "_primary",
      "_sensitivity"
    )
  )
  
  
  x$same_direction <-
    sign(x$NES_primary) ==
    sign(x$NES_sensitivity)
  
  
  x$sig_primary <-
    x$p.adjust_primary < 0.05
  
  
  x$sig_sensitivity <-
    x$p.adjust_sensitivity < 0.05
  
  
  rho <- cor(
    x$NES_primary,
    x$NES_sensitivity,
    method = "spearman"
  )
  
  
  cat("\n========================================\n")
  cat(
    contrast_name,
    "GSEA ROBUSTNESS\n"
  )
  cat("========================================\n\n")
  
  
  cat(
    "Spearman NES rho:",
    round(rho, 4),
    "\n"
  )
  
  
  cat(
    "Same direction:",
    sum(x$same_direction),
    "/ 50\n"
  )
  
  
  cat(
    "Significant primary:",
    sum(x$sig_primary),
    "/ 50\n"
  )
  
  
  cat(
    "Significant sensitivity:",
    sum(x$sig_sensitivity),
    "/ 50\n"
  )
  
  
  cat(
    "Significant in both:",
    sum(
      x$sig_primary &
        x$sig_sensitivity
    ),
    "/ 50\n"
  )
  
  
  return(x)
}


P1_compare <- compare_gsea(
  gsea_P1,
  gsea_P1_sens,
  "P1"
)


P3_compare <- compare_gsea(
  gsea_P3,
  gsea_P3_sens,
  "P3"
)


# ============================================================
# 9. PROGRAMS OF INTEREST
# ============================================================

poi <- c(
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_GLYCOLYSIS"
)


P3_poi <- P3_compare[
  P3_compare$ID %in% poi,
]


P3_poi <- P3_poi[
  order(
    P3_poi$NES_primary
  ),
]


cat("\n========================================\n")
cat("P3 PATHWAYS OF INTEREST ROBUSTNESS\n")
cat("========================================\n\n")

print(
  P3_poi,
  row.names = FALSE
)


# ============================================================
# 10. PRIMARY P1 / P2 / P3 POI TABLE
# ============================================================

extract_poi <- function(
    gsea,
    label
) {
  
  x <- gsea[
    gsea$ID %in% poi,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ]
  
  x$contrast <- label
  
  x
}


poi_primary <- bind_rows(
  
  extract_poi(
    gsea_P1,
    "P1_M_vs_Parental"
  ),
  
  extract_poi(
    gsea_P2,
    "P2_WD_vs_Parental"
  ),
  
  extract_poi(
    gsea_P3,
    "P3_M_vs_WD"
  )
)


poi_wide <- poi_primary %>%
  
  pivot_wider(
    names_from = contrast,
    values_from = c(
      NES,
      p.adjust
    )
  )


cat("\n========================================\n")
cat("PRIMARY P1-P3 PATHWAYS OF INTEREST\n")
cat("========================================\n\n")

print(
  poi_wide,
  row.names = FALSE
)


# ============================================================
# 11. SAVE OUTPUTS
# ============================================================

outdir <- file.path(
  indir,
  "Hallmark_GSEA"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  gsea_P1,
  file.path(
    outdir,
    "PSN1_P1_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_P2,
  file.path(
    outdir,
    "PSN1_P2_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_P3,
  file.path(
    outdir,
    "PSN1_P3_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_P1_sens,
  file.path(
    outdir,
    "PSN1_P1_Hallmark_GSEA_without_maintain2.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_P3_sens,
  file.path(
    outdir,
    "PSN1_P3_Hallmark_GSEA_without_maintain2.csv"
  ),
  row.names = FALSE
)


write.csv(
  P1_compare,
  file.path(
    outdir,
    "PSN1_P1_GSEA_robustness.csv"
  ),
  row.names = FALSE
)


write.csv(
  P3_compare,
  file.path(
    outdir,
    "PSN1_P3_GSEA_robustness.csv"
  ),
  row.names = FALSE
)


write.csv(
  poi_wide,
  file.path(
    outdir,
    "PSN1_primary_pathways_of_interest.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("PSN1 HALLMARK ROBUSTNESS COMPLETE\n")
cat("========================================\n")