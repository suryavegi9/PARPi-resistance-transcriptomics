# ============================================================
# GSE202434 - DATASET 5
# HCC1806 H3 vs PSN1 P3 HALLMARK COMPARISON
#
# HCC1806 H3:
#   TalaR-M vs TalaR-DF
#   Drug exposure changes while resistance persists
#
# PSN1 P3:
#   Maintain vs Withdraw
#   Drug withdrawal accompanies resistance reversion
#
# Goal:
# Identify shared exposure-associated pathway shifts
# and pathways that diverge between the two models.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 1. LOAD GSEA RESULTS
# ============================================================

h3 <- read.csv(
  paste0(
    "results/GSE202434/HCC1806/Hallmark_GSEA/",
    "HCC1806_H3_Hallmark_GSEA_all.csv"
  ),
  stringsAsFactors = FALSE
)


p3 <- read.csv(
  paste0(
    "results/GSE202434/PSN1/Hallmark_GSEA/",
    "PSN1_P3_Hallmark_GSEA_all.csv"
  ),
  stringsAsFactors = FALSE
)


p3_sens <- read.csv(
  paste0(
    "results/GSE202434/PSN1/Hallmark_GSEA/",
    "PSN1_P3_Hallmark_GSEA_without_maintain2.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. BASIC QC
# ============================================================

cat("\n========================================\n")
cat("SOURCE TABLE QC\n")
cat("========================================\n\n")

cat(
  "HCC1806 H3 Hallmarks:",
  nrow(h3),
  "\n"
)

cat(
  "PSN1 P3 Hallmarks:",
  nrow(p3),
  "\n"
)

cat(
  "PSN1 P3 sensitivity Hallmarks:",
  nrow(p3_sens),
  "\n"
)


stopifnot(
  nrow(h3) == 50,
  nrow(p3) == 50,
  nrow(p3_sens) == 50
)


stopifnot(
  sum(duplicated(h3$ID)) == 0,
  sum(duplicated(p3$ID)) == 0,
  sum(duplicated(p3_sens$ID)) == 0
)


# ============================================================
# 3. BUILD MASTER D5 WITHDRAWAL-COMPARISON TABLE
# ============================================================

H3 <- h3[
  ,
  c(
    "ID",
    "NES",
    "p.adjust"
  )
]

colnames(H3) <- c(
  "pathway",
  "H3_NES",
  "H3_padj"
)


P3 <- p3[
  ,
  c(
    "ID",
    "NES",
    "p.adjust"
  )
]

colnames(P3) <- c(
  "pathway",
  "P3_NES",
  "P3_padj"
)


P3_sens <- p3_sens[
  ,
  c(
    "ID",
    "NES",
    "p.adjust"
  )
]

colnames(P3_sens) <- c(
  "pathway",
  "P3_sensitivity_NES",
  "P3_sensitivity_padj"
)


d5 <- Reduce(
  
  function(x, y) {
    
    merge(
      x,
      y,
      by = "pathway",
      all = TRUE
    )
    
  },
  
  list(
    H3,
    P3,
    P3_sens
  )
)


stopifnot(
  nrow(d5) == 50,
  sum(is.na(d5)) == 0
)


# ============================================================
# 4. DIRECTION + SIGNIFICANCE FLAGS
# ============================================================

d5$H3_direction <- ifelse(
  d5$H3_NES > 0,
  "UP",
  "DOWN"
)


d5$P3_direction <- ifelse(
  d5$P3_NES > 0,
  "UP",
  "DOWN"
)


d5$P3_sensitivity_direction <- ifelse(
  d5$P3_sensitivity_NES > 0,
  "UP",
  "DOWN"
)


d5$H3_significant <-
  d5$H3_padj < 0.05


d5$P3_significant <-
  d5$P3_padj < 0.05


d5$P3_sensitivity_significant <-
  d5$P3_sensitivity_padj < 0.05


# ============================================================
# 5. GLOBAL H3 vs P3 CONCORDANCE
# ============================================================

d5$same_direction <-
  sign(d5$H3_NES) ==
  sign(d5$P3_NES)


d5$P3_robust_direction <-
  sign(d5$P3_NES) ==
  sign(d5$P3_sensitivity_NES)


rho_H3_P3 <- cor(
  d5$H3_NES,
  d5$P3_NES,
  method = "spearman"
)


cat("\n========================================\n")
cat("H3 vs P3 GLOBAL CONCORDANCE\n")
cat("========================================\n\n")

cat(
  "Spearman NES rho:",
  round(rho_H3_P3, 4),
  "\n"
)

cat(
  "Same direction:",
  sum(d5$same_direction),
  "/ 50\n"
)

cat(
  "Direction agreement:",
  round(
    100 * mean(d5$same_direction),
    1
  ),
  "%\n"
)


# ============================================================
# 6. SHARED EXPOSURE-ASSOCIATED CANDIDATES
#
# Same direction in both withdrawal comparisons
# AND significant in both primary analyses.
#
# Reminder:
# Positive NES means higher in maintained/on-drug state.
# Negative NES means higher after withdrawal.
# ============================================================

d5$shared_exposure_candidate <-
  
  d5$same_direction &
  
  d5$H3_significant &
  
  d5$P3_significant


# Stronger version:
#
# PSN1 direction must survive maintain-2 exclusion
# AND remain significant in sensitivity analysis.

d5$robust_shared_exposure_candidate <-
  
  d5$shared_exposure_candidate &
  
  d5$P3_robust_direction &
  
  d5$P3_sensitivity_significant


# ============================================================
# 7. DIVERGENT CANDIDATES
# ============================================================

d5$divergent_candidate <-
  !d5$same_direction


d5$strong_divergent_candidate <-
  
  !d5$same_direction &
  
  d5$H3_significant &
  
  d5$P3_significant


# ============================================================
# 8. DESCRIPTIVE DIFFERENCE IN NES
# ============================================================

d5$NES_difference_H3_minus_P3 <-
  d5$H3_NES -
  d5$P3_NES


d5$abs_NES_difference <-
  abs(
    d5$NES_difference_H3_minus_P3
  )


# ============================================================
# 9. PRINT GLOBAL CLASSIFICATION SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D5 PATHWAY COMPARISON SUMMARY\n")
cat("========================================\n\n")

cat(
  "Shared exposure-associated candidates:",
  sum(d5$shared_exposure_candidate),
  "/ 50\n"
)

cat(
  "Robust shared exposure-associated candidates:",
  sum(d5$robust_shared_exposure_candidate),
  "/ 50\n"
)

cat(
  "Divergent pathways:",
  sum(d5$divergent_candidate),
  "/ 50\n"
)

cat(
  "Strong divergent pathways:",
  sum(d5$strong_divergent_candidate),
  "/ 50\n"
)


# ============================================================
# 10. ROBUST SHARED EXPOSURE CANDIDATES
# ============================================================

shared <- d5[
  d5$robust_shared_exposure_candidate,
  c(
    "pathway",
    "H3_NES",
    "H3_padj",
    "P3_NES",
    "P3_padj",
    "P3_sensitivity_NES",
    "P3_sensitivity_padj"
  )
]


shared$exposure_direction <- ifelse(
  shared$H3_NES > 0,
  "Higher_on_drug",
  "Higher_after_withdrawal"
)


shared <- shared[
  order(
    shared$exposure_direction,
    -abs(shared$H3_NES)
  ),
]


cat("\n========================================\n")
cat("ROBUST SHARED EXPOSURE-ASSOCIATED CANDIDATES\n")
cat("========================================\n\n")

print(
  shared,
  row.names = FALSE
)


# ============================================================
# 11. STRONG DIVERGENT CANDIDATES
# ============================================================

divergent <- d5[
  d5$strong_divergent_candidate,
  c(
    "pathway",
    "H3_NES",
    "H3_padj",
    "P3_NES",
    "P3_padj",
    "P3_sensitivity_NES",
    "P3_sensitivity_padj",
    "abs_NES_difference"
  )
]


divergent <- divergent[
  order(
    -divergent$abs_NES_difference
  ),
]


cat("\n========================================\n")
cat("STRONG DIVERGENT CANDIDATES\n")
cat("========================================\n\n")

print(
  divergent,
  row.names = FALSE
)


# ============================================================
# 12. PRE-SPECIFIED PATHWAYS OF INTEREST
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


poi_table <- d5[
  d5$pathway %in% poi,
  c(
    "pathway",
    "H3_NES",
    "H3_padj",
    "P3_NES",
    "P3_padj",
    "P3_sensitivity_NES",
    "P3_sensitivity_padj",
    "same_direction",
    "robust_shared_exposure_candidate",
    "strong_divergent_candidate"
  )
]


poi_table <- poi_table[
  order(
    poi_table$H3_NES
  ),
]


cat("\n========================================\n")
cat("PRE-SPECIFIED PATHWAYS OF INTEREST\n")
cat("========================================\n\n")

print(
  poi_table,
  row.names = FALSE
)


# ============================================================
# 13. SAVE OUTPUTS
# ============================================================

outdir <- "results/GSE202434/D5_integration"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  d5,
  file.path(
    outdir,
    "D5_H3_P3_Hallmark_master_comparison.csv"
  ),
  row.names = FALSE
)


write.csv(
  shared,
  file.path(
    outdir,
    "D5_robust_shared_exposure_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  divergent,
  file.path(
    outdir,
    "D5_strong_divergent_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  poi_table,
  file.path(
    outdir,
    "D5_predefined_pathways_H3_vs_P3.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("D5 H3-P3 COMPARISON COMPLETE\n")
cat("========================================\n")