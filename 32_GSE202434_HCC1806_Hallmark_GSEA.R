# ============================================================
# GSE202434 - DATASET 5
# HCC1806 HALLMARK GSEA
#
# H1 = TalaR-M  vs Parental
# H2 = TalaR-DF vs Parental
# H3 = TalaR-M  vs TalaR-DF
#
# H3 is the PRIMARY exposure-state contrast.
#
# Ranking metric:
# limma moderated t-statistic
#
# IMPORTANT:
# GSEA uses ALL filtered genes, not only DEGs.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))
library(tidyr)
library(clusterProfiler)
library(msigdbr)
library(dplyr)


# ============================================================
# 1. LOAD COMPLETE DE RESULTS
# ============================================================

indir <- "results/GSE202434/HCC1806"

H1 <- read.csv(
  file.path(
    indir,
    "HCC1806_H1_M_vs_Parental_all.csv"
  ),
  stringsAsFactors = FALSE
)

H2 <- read.csv(
  file.path(
    indir,
    "HCC1806_H2_DF_vs_Parental_all.csv"
  ),
  stringsAsFactors = FALSE
)

H3 <- read.csv(
  file.path(
    indir,
    "HCC1806_H3_M_vs_DF_all.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. GET HALLMARK GENE SETS
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


cat("\n========================================\n")
cat("HALLMARK DATABASE\n")
cat("========================================\n\n")

cat(
  "MSigDB version:",
  unique(hallmark$db_version),
  "\n"
)

cat(
  "Hallmark pathways:",
  length(unique(term2gene$gs_name)),
  "\n"
)


# ============================================================
# 3. FUNCTION TO RUN GSEA
# ============================================================

run_hallmark_gsea <- function(
    de_table,
    contrast_name
) {
  
  # ----------------------------------------------------------
  # Ranked gene list
  #
  # Positive t:
  # enriched toward first condition in contrast
  #
  # Negative t:
  # enriched toward second condition
  # ----------------------------------------------------------
  
  ranks <- de_table$t
  
  names(ranks) <- de_table$gene
  
  
  # Remove NA / infinite values
  
  valid <- is.finite(ranks) &
    !is.na(names(ranks)) &
    names(ranks) != ""
  
  ranks <- ranks[valid]
  
  
  # Gene symbols should already be unique
  
  stopifnot(
    sum(duplicated(names(ranks))) == 0
  )
  
  
  ranks <- sort(
    ranks,
    decreasing = TRUE
  )
  
  
  cat("\n========================================\n")
  cat(
    "RUNNING:",
    contrast_name,
    "\n"
  )
  cat("========================================\n\n")
  
  cat(
    "Ranked genes:",
    length(ranks),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # GSEA
  # ----------------------------------------------------------
  
  set.seed(123)
  
  gsea <- GSEA(
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
  
  
  res <- as.data.frame(
    gsea
  )
  
  
  # ----------------------------------------------------------
  # Sanity checks
  # ----------------------------------------------------------
  
  cat(
    "Hallmarks returned:",
    nrow(res),
    "\n"
  )
  
  cat(
    "Significant Hallmarks (BH < 0.05):",
    sum(
      res$p.adjust < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Add contrast name
  # ----------------------------------------------------------
  
  res$contrast <- contrast_name
  
  
  # ----------------------------------------------------------
  # Order by adjusted p-value
  # ----------------------------------------------------------
  
  res <- res[
    order(
      res$p.adjust,
      -abs(res$NES)
    ),
  ]
  
  
  return(
    res
  )
}


# ============================================================
# 4. RUN H1-H3
# ============================================================

gsea_H1 <- run_hallmark_gsea(
  H1,
  "H1_M_vs_Parental"
)

gsea_H2 <- run_hallmark_gsea(
  H2,
  "H2_DF_vs_Parental"
)

gsea_H3 <- run_hallmark_gsea(
  H3,
  "H3_M_vs_DF"
)


# ============================================================
# 5. SAVE ALL RESULTS
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
  gsea_H1,
  file.path(
    outdir,
    "HCC1806_H1_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_H2,
  file.path(
    outdir,
    "HCC1806_H2_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_H3,
  file.path(
    outdir,
    "HCC1806_H3_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_H1[
    gsea_H1$p.adjust < 0.05,
  ],
  file.path(
    outdir,
    "HCC1806_H1_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_H2[
    gsea_H2$p.adjust < 0.05,
  ],
  file.path(
    outdir,
    "HCC1806_H2_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_H3[
    gsea_H3$p.adjust < 0.05,
  ],
  file.path(
    outdir,
    "HCC1806_H3_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. PRINT H3 SIGNIFICANT HALLMARKS
# ============================================================

H3_sig <- gsea_H3[
  gsea_H3$p.adjust < 0.05,
  c(
    "ID",
    "NES",
    "pvalue",
    "p.adjust"
  )
]


cat("\n========================================\n")
cat("H3 SIGNIFICANT HALLMARKS\n")
cat("TalaR-M vs TalaR-DF\n")
cat("========================================\n\n")

print(
  H3_sig,
  row.names = FALSE
)


# ============================================================
# 7. D1-D4 PATHWAYS OF SPECIAL INTEREST
#
# These directly test the hypotheses generated by the
# previous cross-dataset integration.
# ============================================================

programs_of_interest <- c(
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_GLYCOLYSIS"
)


extract_programs <- function(
    gsea_result,
    contrast_label
) {
  
  x <- gsea_result[
    gsea_result$ID %in%
      programs_of_interest,
    c(
      "ID",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ]
  
  x$contrast <- contrast_label
  
  x
}


poi <- bind_rows(
  
  extract_programs(
    gsea_H1,
    "H1_M_vs_Parental"
  ),
  
  extract_programs(
    gsea_H2,
    "H2_DF_vs_Parental"
  ),
  
  extract_programs(
    gsea_H3,
    "H3_M_vs_DF"
  )
)


poi$significant <- poi$p.adjust < 0.05


# ============================================================
# 8. WIDE COMPARISON TABLE
# ============================================================

poi_wide <- poi %>%
  
  select(
    ID,
    contrast,
    NES,
    p.adjust
  ) %>%
  
  pivot_wider(
    names_from = contrast,
    values_from = c(
      NES,
      p.adjust
    )
  )


cat("\n========================================\n")
cat("PATHWAYS OF INTEREST ACROSS H1-H3\n")
cat("========================================\n\n")

print(
  poi_wide,
  row.names = FALSE
)


write.csv(
  poi,
  file.path(
    outdir,
    "HCC1806_pathways_of_interest_long.csv"
  ),
  row.names = FALSE
)

write.csv(
  poi_wide,
  file.path(
    outdir,
    "HCC1806_pathways_of_interest_wide.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. H3 EXPOSURE-STATE SUMMARY
# ============================================================

H3_poi <- poi[
  poi$contrast == "H3_M_vs_DF",
]


H3_poi <- H3_poi[
  order(
    H3_poi$NES
  ),
]


cat("\n========================================\n")
cat("H3 EXPOSURE-STATE PROGRAMS OF INTEREST\n")
cat("========================================\n\n")

print(
  H3_poi,
  row.names = FALSE
)


# ============================================================
# 10. FINISHED
# ============================================================

cat("\n========================================\n")
cat("HCC1806 HALLMARK GSEA COMPLETE\n")
cat("========================================\n")