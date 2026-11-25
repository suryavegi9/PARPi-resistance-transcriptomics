# ============================================================
# GSE202434 - DATASET 5
# FINAL MECHANISTIC PATHWAY INTEGRATION
#
# HCC1806:
# H1 = Maintained resistant vs Parental
# H2 = Drug-free persistent resistant vs Parental
# H3 = Maintained resistant vs Drug-free resistant
#
# PSN1:
# P1 = Maintained resistant vs Parental
# P2 = Withdrawn / re-sensitized vs Parental
# P3 = Maintained resistant vs Withdrawn
#
# Goal:
# Separate exposure-associated programs from candidate
# resistance-persistence / reversion-associated programs.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 1. FILE LOCATIONS
# ============================================================

hcc_dir <-
  "results/GSE202434/HCC1806/Hallmark_GSEA"

psn_dir <-
  "results/GSE202434/PSN1/Hallmark_GSEA"

outdir <-
  "results/GSE202434/D5_integration"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. HELPER TO LOAD CLEAN TABLE
# ============================================================

load_gsea <- function(
    file,
    prefix
) {
  
  x <- read.csv(
    file,
    stringsAsFactors = FALSE
  )
  
  x <- x[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ]
  
  colnames(x) <- c(
    "pathway",
    paste0(prefix, "_NES"),
    paste0(prefix, "_padj")
  )
  
  x
}


# ============================================================
# 3. LOAD ALL SIX PRIMARY CONTRASTS
# ============================================================

H1 <- load_gsea(
  file.path(
    hcc_dir,
    "HCC1806_H1_Hallmark_GSEA_all.csv"
  ),
  "H1"
)

H2 <- load_gsea(
  file.path(
    hcc_dir,
    "HCC1806_H2_Hallmark_GSEA_all.csv"
  ),
  "H2"
)

H3 <- load_gsea(
  file.path(
    hcc_dir,
    "HCC1806_H3_Hallmark_GSEA_all.csv"
  ),
  "H3"
)


P1 <- load_gsea(
  file.path(
    psn_dir,
    "PSN1_P1_Hallmark_GSEA_all.csv"
  ),
  "P1"
)

P2 <- load_gsea(
  file.path(
    psn_dir,
    "PSN1_P2_Hallmark_GSEA_all.csv"
  ),
  "P2"
)

P3 <- load_gsea(
  file.path(
    psn_dir,
    "PSN1_P3_Hallmark_GSEA_all.csv"
  ),
  "P3"
)


# ============================================================
# 4. LOAD PSN1 OUTLIER SENSITIVITY
# ============================================================

P1_sens <- load_gsea(
  file.path(
    psn_dir,
    "PSN1_P1_Hallmark_GSEA_without_maintain2.csv"
  ),
  "P1_sens"
)

P3_sens <- load_gsea(
  file.path(
    psn_dir,
    "PSN1_P3_Hallmark_GSEA_without_maintain2.csv"
  ),
  "P3_sens"
)


# ============================================================
# 5. MERGE MASTER TABLE
# ============================================================

master <- Reduce(
  
  function(x, y) {
    
    merge(
      x,
      y,
      by = "pathway",
      all = TRUE
    )
    
  },
  
  list(
    H1,
    H2,
    H3,
    P1,
    P2,
    P3,
    P1_sens,
    P3_sens
  )
)


stopifnot(
  nrow(master) == 50
)

stopifnot(
  sum(is.na(master)) == 0
)


cat("\n========================================\n")
cat("MASTER TABLE QC\n")
cat("========================================\n\n")

cat(
  "Hallmarks:",
  nrow(master),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(master)),
  "\n"
)


# ============================================================
# 6. GLOBAL CROSS-MODEL COMPARISONS
# ============================================================

compare_profiles <- function(
    x,
    y,
    label
) {
  
  rho <- cor(
    x,
    y,
    method = "spearman"
  )
  
  same <- mean(
    sign(x) ==
      sign(y)
  )
  
  cat(
    label,
    "\n"
  )
  
  cat(
    "  Spearman rho:",
    round(rho, 4),
    "\n"
  )
  
  cat(
    "  Direction agreement:",
    round(
      100 * same,
      1
    ),
    "%\n\n"
  )
  
  c(
    rho = rho,
    direction_agreement = same
  )
}


cat("\n========================================\n")
cat("CROSS-MODEL STATE COMPARISONS\n")
cat("========================================\n\n")


on_drug_compare <- compare_profiles(
  master$H1_NES,
  master$P1_NES,
  "H1 vs P1: resistant / maintained-on-drug"
)


withdraw_compare <- compare_profiles(
  master$H2_NES,
  master$P2_NES,
  "H2 vs P2: post-withdrawal states"
)


exposure_compare <- compare_profiles(
  master$H3_NES,
  master$P3_NES,
  "H3 vs P3: maintained vs withdrawn"
)


# ============================================================
# 7. ROBUST SHARED EXPOSURE PROGRAMS
#
# Must:
# - agree H3 and P3
# - significant in both
# - PSN1 direction survives maintain-2 exclusion
# - PSN1 remains significant in sensitivity
# ============================================================

master$exposure_same_direction <-
  sign(master$H3_NES) ==
  sign(master$P3_NES)


master$P3_robust_direction <-
  sign(master$P3_NES) ==
  sign(master$P3_sens_NES)


master$robust_exposure_program <-
  
  master$exposure_same_direction &
  
  master$H3_padj < 0.05 &
  
  master$P3_padj < 0.05 &
  
  master$P3_robust_direction &
  
  master$P3_sens_padj < 0.05


# ============================================================
# 8. RESISTANT-STATE CONSENSUS
#
# A pathway must point the same way in:
#
# H1 = HCC resistant on drug
# H2 = HCC resistant off drug
# P1 = PSN resistant on drug
#
# This asks whether the pathway is repeatedly associated
# with the resistant phenotype before PSN1 withdrawal.
# ============================================================

master$resistant_consensus_direction <-
  
  sign(master$H1_NES) ==
  sign(master$H2_NES) &
  
  sign(master$H1_NES) ==
  sign(master$P1_NES)


master$resistant_consensus_significant <-
  
  master$resistant_consensus_direction &
  
  master$H1_padj < 0.05 &
  
  master$H2_padj < 0.05 &
  
  master$P1_padj < 0.05


# ============================================================
# 9. PSN1 P1 OUTLIER ROBUSTNESS
# ============================================================

master$P1_robust <-
  
  sign(master$P1_NES) ==
  sign(master$P1_sens_NES) &
  
  master$P1_sens_padj < 0.05


# ============================================================
# 10. WHAT HAPPENS AFTER PSN1 RESENSITIZATION?
#
# Does P2 still significantly support the SAME direction
# observed in the resistant-state consensus?
# ============================================================

master$P2_same_as_resistant <- ifelse(
  
  master$resistant_consensus_direction,
  
  sign(master$P2_NES) ==
    sign(master$H2_NES),
  
  FALSE
)


master$P2_same_direction_significant <-
  
  master$P2_same_as_resistant &
  
  master$P2_padj < 0.05


# ============================================================
# 11. CLEAN RESISTANCE-PERSISTENCE CANDIDATE
#
# Requirements:
#
# 1. H1/H2/P1 all agree
# 2. all three significant
# 3. P1 survives outlier sensitivity
# 4. P2 no longer significantly supports that same direction
# 5. pathway is NOT already classified as robust exposure-driven
#
# Interpretation:
# altered in resistant states,
# persists in HCC after withdrawal,
# but disappears or reverses after PSN1 resensitization.
#
# This remains ASSOCIATIONAL.
# ============================================================

master$resistance_persistence_candidate <-
  
  master$resistant_consensus_significant &
  
  master$P1_robust &
  
  !master$P2_same_direction_significant &
  
  !master$robust_exposure_program


# ============================================================
# 12. STRONG REVERSION CANDIDATE
#
# Even stricter:
# P2 becomes significantly OPPOSITE to the resistant state.
# ============================================================

master$P2_opposite_significant <-
  
  master$resistant_consensus_direction &
  
  sign(master$P2_NES) !=
  sign(master$H2_NES) &
  
  master$P2_padj < 0.05


master$strong_reversion_candidate <-
  
  master$resistant_consensus_significant &
  
  master$P1_robust &
  
  master$P2_opposite_significant &
  
  !master$robust_exposure_program


# ============================================================
# 13. SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D5 MECHANISTIC CLASSIFICATION SUMMARY\n")
cat("========================================\n\n")

cat(
  "Robust shared exposure programs:",
  sum(master$robust_exposure_program),
  "/ 50\n"
)

cat(
  "Resistant-state consensus pathways:",
  sum(master$resistant_consensus_significant),
  "/ 50\n"
)

cat(
  "Clean resistance-persistence candidates:",
  sum(master$resistance_persistence_candidate),
  "/ 50\n"
)

cat(
  "Strong reversion candidates:",
  sum(master$strong_reversion_candidate),
  "/ 50\n"
)


# ============================================================
# 14. EXPOSURE PROGRAM TABLE
# ============================================================

exposure_table <- master[
  master$robust_exposure_program,
  c(
    "pathway",
    "H3_NES",
    "H3_padj",
    "P3_NES",
    "P3_padj",
    "P3_sens_NES",
    "P3_sens_padj"
  )
]


exposure_table$exposure_state <- ifelse(
  exposure_table$H3_NES > 0,
  "Higher_on_drug",
  "Higher_after_withdrawal"
)


exposure_table <- exposure_table[
  order(
    exposure_table$exposure_state,
    -abs(exposure_table$H3_NES)
  ),
]


cat("\n========================================\n")
cat("ROBUST EXPOSURE-ASSOCIATED PROGRAMS\n")
cat("========================================\n\n")

print(
  exposure_table,
  row.names = FALSE
)


# ============================================================
# 15. RESISTANCE-PERSISTENCE TABLE
# ============================================================

resistance_table <- master[
  master$resistance_persistence_candidate,
  c(
    "pathway",
    
    "H1_NES",
    "H1_padj",
    
    "H2_NES",
    "H2_padj",
    
    "P1_NES",
    "P1_padj",
    
    "P1_sens_NES",
    "P1_sens_padj",
    
    "P2_NES",
    "P2_padj"
  )
]


if (nrow(resistance_table) > 0) {
  
  resistance_table$direction <- ifelse(
    resistance_table$H2_NES > 0,
    "UP_in_resistance",
    "DOWN_in_resistance"
  )
  
}


cat("\n========================================\n")
cat("CLEAN RESISTANCE-PERSISTENCE CANDIDATES\n")
cat("========================================\n\n")

print(
  resistance_table,
  row.names = FALSE
)


# ============================================================
# 16. STRONG REVERSION TABLE
# ============================================================

reversion_table <- master[
  master$strong_reversion_candidate,
  c(
    "pathway",
    
    "H1_NES",
    "H2_NES",
    "P1_NES",
    
    "P2_NES",
    
    "H1_padj",
    "H2_padj",
    "P1_padj",
    "P2_padj"
  )
]


cat("\n========================================\n")
cat("STRONG REVERSION CANDIDATES\n")
cat("========================================\n\n")

print(
  reversion_table,
  row.names = FALSE
)


# ============================================================
# 17. PRE-SPECIFIED PATHWAYS
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


poi_table <- master[
  master$pathway %in% poi,
  c(
    "pathway",
    
    "H1_NES",
    "H2_NES",
    "H3_NES",
    
    "P1_NES",
    "P2_NES",
    "P3_NES",
    
    "robust_exposure_program",
    "resistant_consensus_significant",
    "resistance_persistence_candidate",
    "strong_reversion_candidate"
  )
]


cat("\n========================================\n")
cat("PRE-SPECIFIED PATHWAY MECHANISTIC SUMMARY\n")
cat("========================================\n\n")

print(
  poi_table,
  row.names = FALSE
)


# ============================================================
# 18. SAVE EVERYTHING
# ============================================================

write.csv(
  master,
  file.path(
    outdir,
    "D5_six_contrast_mechanistic_master.csv"
  ),
  row.names = FALSE
)


write.csv(
  exposure_table,
  file.path(
    outdir,
    "D5_robust_exposure_programs.csv"
  ),
  row.names = FALSE
)


write.csv(
  resistance_table,
  file.path(
    outdir,
    "D5_resistance_persistence_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  reversion_table,
  file.path(
    outdir,
    "D5_strong_reversion_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  poi_table,
  file.path(
    outdir,
    "D5_predefined_mechanistic_summary.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("D5 MECHANISTIC INTEGRATION COMPLETE\n")
cat("========================================\n")