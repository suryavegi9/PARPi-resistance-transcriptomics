# ============================================================
# FINAL D1-D6 TREATMENT-STATE INTEGRATION
#
# Goal:
# Quantify whether treatment-at-RNA-collection contributes
# to cross-dataset Hallmark heterogeneity.
#
# BETWEEN-DATASET EXPOSURE GROUPS
# --------------------------------
# OFF treatment:
#   D1 = drug-free resistant cell-line states
#   D4 = untreated acquired-resistant PDX
#   D6 = recurrent PDX collected long after olaparib withdrawal
#
# ON treatment:
#   D2 = resistant cells maintained in olaparib
#   D3 = resistant cells maintained in olaparib
#
# WITHIN-DATASET CONTROLLED EXPOSURE EVIDENCE
# --------------------------------------------
# D5 H3 = HCC1806 maintained vs drug-free resistant
# D5 P3 = PSN1 maintained vs withdrawn
#
# Positive exposure effect = higher ON drug
# Negative exposure effect = higher AFTER withdrawal
#
# IMPORTANT:
# Between-dataset ON/OFF comparisons are confounded by model,
# cancer type, genotype and in-vivo/in-vitro context.
# They are DESCRIPTIVE, not causal.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. FILES
# ============================================================

master_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "FINAL_D1_D6_Hallmark_master.csv"
)

d5_file <- paste0(
  "results/GSE202434/D5_integration/",
  "D5_six_contrast_mechanistic_master.csv"
)

model_context_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "model_level_sensitivity/",
  "FINAL_robust_context_dependent_both_methods.csv"
)


stopifnot(
  file.exists(master_file),
  file.exists(d5_file),
  file.exists(model_context_file)
)


outdir <- paste0(
  "results/FINAL_D1_D6_integration/",
  "treatment_state"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. LOAD ACCESSION-LEVEL MASTER
# ============================================================

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_master <- c(
  "pathway",
  "D1_NES",
  "D2_NES",
  "D3_NES",
  "D4_NES",
  "D6_NES",
  "classification"
)


stopifnot(
  all(
    required_master %in%
      colnames(master)
  )
)


# ============================================================
# 2. LOAD D5 EXPOSURE CONTRASTS
# ============================================================

d5 <- read.csv(
  d5_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_d5 <- c(
  "pathway",
  
  "H3_NES",
  "H3_padj",
  
  "P3_NES",
  "P3_padj",
  
  "P3_sens_NES",
  "P3_sens_padj",
  
  "robust_exposure_program"
)


stopifnot(
  all(
    required_d5 %in%
      colnames(d5)
  )
)


d5 <- d5[
  ,
  required_d5
]


# ============================================================
# 3. MERGE
# ============================================================

x <- merge(
  master,
  d5,
  by = "pathway",
  all = FALSE
)


stopifnot(
  nrow(x) == 50
)


cat("\n========================================\n")
cat("TREATMENT-STATE MASTER QC\n")
cat("========================================\n\n")

cat(
  "Hallmarks:",
  nrow(x),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(x)),
  "\n"
)


# ============================================================
# 4. BETWEEN-DATASET OBSERVATIONAL EXPOSURE EFFECT
#
# OFF:
# D1 + D4 + D6
#
# ON:
# D2 + D3
#
# Positive delta:
# higher in ON-treatment datasets
#
# Negative delta:
# higher in OFF-treatment datasets
# ============================================================

x$OFF_mean_NES <- rowMeans(
  x[
    ,
    c(
      "D1_NES",
      "D4_NES",
      "D6_NES"
    )
  ]
)


x$ON_mean_NES <- rowMeans(
  x[
    ,
    c(
      "D2_NES",
      "D3_NES"
    )
  ]
)


x$observational_exposure_delta <-
  x$ON_mean_NES -
  x$OFF_mean_NES


x$observational_direction <- ifelse(
  x$observational_exposure_delta > 0,
  "Higher_on_drug",
  "Higher_off_drug"
)


# ============================================================
# 5. CLEAN BETWEEN-GROUP SEPARATION
#
# Strong descriptive criterion:
#
# Positive exposure:
# both ON datasets higher than all three OFF datasets
#
# Negative exposure:
# both ON datasets lower than all three OFF datasets
#
# This is still DESCRIPTIVE, not an inferential test.
# ============================================================

on_min <- pmin(
  x$D2_NES,
  x$D3_NES
)


on_max <- pmax(
  x$D2_NES,
  x$D3_NES
)


off_min <- pmin(
  x$D1_NES,
  x$D4_NES,
  x$D6_NES
)


off_max <- pmax(
  x$D1_NES,
  x$D4_NES,
  x$D6_NES
)


x$clean_between_group_separation <-
  
  (
    x$observational_exposure_delta > 0 &
      on_min > off_max
  ) |
  
  (
    x$observational_exposure_delta < 0 &
      on_max < off_min
  )


# ============================================================
# 6. D5 WITHIN-DATASET EXPOSURE EFFECT
#
# H3 and P3:
# positive = higher ON drug
# negative = higher after withdrawal
# ============================================================

x$D5_exposure_mean <- rowMeans(
  x[
    ,
    c(
      "H3_NES",
      "P3_NES"
    )
  ]
)


x$D5_same_direction <-
  
  sign(
    x$H3_NES
  ) ==
  sign(
    x$P3_NES
  )


x$D5_primary_and_sensitivity_direction_stable <-
  
  sign(
    x$P3_NES
  ) ==
  sign(
    x$P3_sens_NES
  )


# ============================================================
# 7. OBSERVATIONAL vs D5 EXPOSURE AGREEMENT
# ============================================================

x$exposure_direction_agreement <-
  
  sign(
    x$observational_exposure_delta
  ) ==
  sign(
    x$D5_exposure_mean
  )


rho <- cor(
  x$observational_exposure_delta,
  x$D5_exposure_mean,
  method = "spearman"
)


direction_agreement <- mean(
  x$exposure_direction_agreement
)


cat("\n========================================\n")
cat("GLOBAL EXPOSURE-SIGNAL CONCORDANCE\n")
cat("========================================\n\n")


cat(
  "Spearman rho:",
  round(
    rho,
    4
  ),
  "\n"
)


cat(
  "Same exposure direction:",
  sum(
    x$exposure_direction_agreement
  ),
  "/ 50\n"
)


cat(
  "Direction agreement:",
  round(
    100 *
      direction_agreement,
    1
  ),
  "%\n"
)


# ============================================================
# 8. D5 ROBUST EXPOSURE PROGRAMS
# ============================================================

d5_robust <- x[
  x$robust_exposure_program,
]


cat("\n========================================\n")
cat("D5 ROBUST EXPOSURE PROGRAMS\n")
cat("========================================\n\n")


cat(
  "D5 robust exposure programs:",
  nrow(d5_robust),
  "/ 50\n"
)


cat(
  "Also agree with observational ON/OFF direction:",
  sum(
    d5_robust$exposure_direction_agreement
  ),
  "/",
  nrow(d5_robust),
  "\n"
)


cat(
  "Also show clean ON/OFF accession separation:",
  sum(
    d5_robust$
      exposure_direction_agreement &
      d5_robust$
      clean_between_group_separation
  ),
  "/",
  nrow(d5_robust),
  "\n"
)


# ============================================================
# 9. EXPOSURE-EVIDENCE CLASSIFICATION
#
# STRONG CONVERGENT:
# - robust D5 exposure program
# - observational direction agrees
# - clean ON/OFF accession separation
#
# CONVERGENT:
# - robust D5 exposure program
# - observational direction agrees
# - no complete separation required
#
# D5-SPECIFIC:
# - robust D5 exposure program
# - observational direction does not agree
#
# OBSERVATIONAL ONLY:
# - clean ON/OFF separation
# - not robust in D5
#
# Otherwise weak/mixed.
# ============================================================

x$exposure_evidence <-
  "Weak/mixed"


strong_convergent <-
  
  x$robust_exposure_program &
  
  x$exposure_direction_agreement &
  
  x$clean_between_group_separation


x$exposure_evidence[
  strong_convergent
] <- "Strong convergent"


convergent <-
  
  x$robust_exposure_program &
  
  x$exposure_direction_agreement &
  
  !x$clean_between_group_separation


x$exposure_evidence[
  convergent
] <- "Convergent"


d5_specific <-
  
  x$robust_exposure_program &
  
  !x$exposure_direction_agreement


x$exposure_evidence[
  d5_specific
] <- "D5-specific exposure"


observational_only <-
  
  !x$robust_exposure_program &
  
  x$clean_between_group_separation


x$exposure_evidence[
  observational_only
] <- "Observational-only"


# ============================================================
# 10. EXPOSURE CLASSIFICATION COUNTS
# ============================================================

exposure_levels <- c(
  "Strong convergent",
  "Convergent",
  "D5-specific exposure",
  "Observational-only",
  "Weak/mixed"
)


exposure_counts <- data.frame(
  
  exposure_evidence =
    exposure_levels,
  
  n = sapply(
    
    exposure_levels,
    
    function(z) {
      
      sum(
        x$exposure_evidence ==
          z
      )
    }
  )
)


cat("\n========================================\n")
cat("EXPOSURE-EVIDENCE CLASSIFICATION\n")
cat("========================================\n\n")


print(
  exposure_counts,
  row.names = FALSE
)


# ============================================================
# 11. STRONG CONVERGENT PATHWAYS
# ============================================================

strong_table <- x[
  x$exposure_evidence ==
    "Strong convergent",
]


strong_table <- strong_table[
  order(
    -abs(
      strong_table$D5_exposure_mean
    )
  ),
]


cat("\n========================================\n")
cat("STRONG CONVERGENT EXPOSURE-ALIGNED PATHWAYS\n")
cat("========================================\n\n")


print(
  strong_table[
    ,
    c(
      "pathway",
      
      "D1_NES",
      "D2_NES",
      "D3_NES",
      "D4_NES",
      "D6_NES",
      
      "OFF_mean_NES",
      "ON_mean_NES",
      
      "observational_exposure_delta",
      
      "H3_NES",
      "P3_NES",
      "D5_exposure_mean"
    )
  ],
  row.names = FALSE
)


# # ============================================================
# 12. ALL D5 ROBUST PROGRAMS
#
# Recreate from x AFTER exposure_evidence has been assigned,
# so all derived classification columns are present.
# ============================================================

d5_robust_final <- x[
  x$robust_exposure_program,
  ,
  drop = FALSE
]


d5_robust_final <- d5_robust_final[
  order(
    -abs(
      d5_robust_final$D5_exposure_mean
    )
  ),
  ,
  drop = FALSE
]


cat("\n========================================\n")
cat("ALL ROBUST D5 EXPOSURE PROGRAMS\n")
cat("========================================\n\n")


print(
  d5_robust_final[
    ,
    c(
      "pathway",
      
      "observational_exposure_delta",
      
      "H3_NES",
      "P3_NES",
      
      "D5_exposure_mean",
      
      "exposure_direction_agreement",
      
      "clean_between_group_separation",
      
      "exposure_evidence"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ============================================================
# 13. ROBUST CONTEXT-DEPENDENT PATHWAYS
#
# These were context-dependent in BOTH accession-level and
# model-level integrations.
# ============================================================

robust_context <- read.csv(
  model_context_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


robust_context_names <- unique(
  robust_context$pathway
)


context_exposure <- x[
  x$pathway %in%
    robust_context_names,
]


context_exposure$exposure_supported <-
  
  context_exposure$robust_exposure_program &
  
  context_exposure$exposure_direction_agreement


cat("\n========================================\n")
cat("ROBUST CONTEXT-DEPENDENT PATHWAYS:\n")
cat("TREATMENT-STATE SUPPORT\n")
cat("========================================\n\n")


cat(
  "Robust context-dependent pathways:",
  nrow(context_exposure),
  "\n"
)


cat(
  "With convergent D5 + observational exposure support:",
  sum(
    context_exposure$exposure_supported
  ),
  "/",
  nrow(context_exposure),
  "\n\n"
)


print(
  context_exposure[
    ,
    c(
      "pathway",
      
      "observational_exposure_delta",
      
      "D5_exposure_mean",
      
      "robust_exposure_program",
      
      "exposure_direction_agreement",
      
      "clean_between_group_separation",
      
      "exposure_evidence"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14. PRE-SPECIFIED / CENTRAL PATHWAYS
# ============================================================

central_pathways <- c(
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_GLYCOLYSIS",
  
  "HALLMARK_INFLAMMATORY_RESPONSE"
)


central <- x[
  x$pathway %in%
    central_pathways,
]


cat("\n========================================\n")
cat("CENTRAL PATHWAYS: TREATMENT-STATE SUMMARY\n")
cat("========================================\n\n")


print(
  central[
    ,
    c(
      "pathway",
      
      "OFF_mean_NES",
      "ON_mean_NES",
      "observational_exposure_delta",
      
      "H3_NES",
      "P3_NES",
      "D5_exposure_mean",
      
      "robust_exposure_program",
      
      "exposure_direction_agreement",
      
      "clean_between_group_separation",
      
      "exposure_evidence"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 15. TNFA/NFKB SPECIFIC AUDIT
# ============================================================

tnfa <- x[
  x$pathway ==
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
]


cat("\n========================================\n")
cat("TNFA / NFKB SPECIFIC CHECK\n")
cat("========================================\n\n")


print(
  tnfa[
    ,
    c(
      "pathway",
      
      "D1_NES",
      "D2_NES",
      "D3_NES",
      "D4_NES",
      "D6_NES",
      
      "OFF_mean_NES",
      "ON_mean_NES",
      
      "observational_exposure_delta",
      
      "H3_NES",
      "P3_NES",
      
      "D5_exposure_mean",
      
      "exposure_evidence"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 16. E2F/G2M/MYC CHECK
# ============================================================

cell_cycle <- x[
  x$pathway %in%
    c(
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_MYC_TARGETS_V1"
    ),
]


cat("\n========================================\n")
cat("E2F / G2M / MYC EXPOSURE CHECK\n")
cat("========================================\n\n")


print(
  cell_cycle[
    ,
    c(
      "pathway",
      
      "D1_NES",
      "D2_NES",
      "D3_NES",
      "D4_NES",
      "D6_NES",
      
      "observational_exposure_delta",
      
      "H3_NES",
      "P3_NES",
      
      "D5_exposure_mean",
      
      "exposure_evidence"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 17. SAVE
# ============================================================

write.csv(
  
  x,
  
  file.path(
    outdir,
    "FINAL_treatment_state_master.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  exposure_counts,
  
  file.path(
    outdir,
    "FINAL_treatment_state_classification_counts.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  strong_table,
  
  file.path(
    outdir,
    "FINAL_strong_convergent_exposure_pathways.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  context_exposure,
  
  file.path(
    outdir,
    "FINAL_context_dependent_treatment_support.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  central,
  
  file.path(
    outdir,
    "FINAL_central_pathway_treatment_summary.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 18. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("FINAL TREATMENT-STATE INTEGRATION COMPLETE\n")
cat("========================================\n\n")

cat(
  "Output directory:",
  outdir,
  "\n"
)