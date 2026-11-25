# ============================================================
# FINAL D1-D6 MODEL-LEVEL SENSITIVITY INTEGRATION
#
# PURPOSE
# -------
# Stress-test the PRIMARY accession-level integration from
# script 45 by representing biologically distinct submodels
# separately.
#
# IMPORTANT:
# These 9 units are NOT all statistically independent.
# This analysis is therefore a SENSITIVITY ANALYSIS only.
#
# Units:
# 1. D1 BRCA1-deficient
# 2. D1 BRCA1-proficient
# 3. D2
# 4. D3
# 5. D4
# 6. D5 HCC1806 maintained resistant vs parental (H1)
# 7. D5 HCC1806 drug-free resistant vs parental (H2)
# 8. D5 PSN1 maintained resistant vs parental (P1)
# 9. D6
#
# H3 and P3 are NOT included because they measure
# maintained-vs-withdrawn exposure state, not resistance
# versus parental state.
#
# P1 sensitivity is used only as robustness evidence.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. OUTPUT
# ============================================================

outdir <- "results/FINAL_D1_D6_integration/model_level_sensitivity"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. HELPER FUNCTION
# ============================================================

read_gsea <- function(file, prefix) {
  
  stopifnot(
    file.exists(file)
  )
  
  x <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  if ("ID" %in% colnames(x)) {
    
    pathway_col <- "ID"
    
  } else if ("pathway" %in% colnames(x)) {
    
    pathway_col <- "pathway"
    
  } else {
    
    stop(
      paste(
        "Could not find pathway column in:",
        file
      )
    )
  }
  
  
  if (!"NES" %in% colnames(x)) {
    
    stop(
      paste(
        "NES column missing in:",
        file
      )
    )
  }
  
  
  if ("p.adjust" %in% colnames(x)) {
    
    padj_col <- "p.adjust"
    
  } else if ("padj" %in% colnames(x)) {
    
    padj_col <- "padj"
    
  } else {
    
    stop(
      paste(
        "Adjusted-P column missing in:",
        file
      )
    )
  }
  
  
  out <- data.frame(
    
    pathway =
      as.character(
        x[[pathway_col]]
      ),
    
    NES =
      as.numeric(
        x$NES
      ),
    
    padj =
      as.numeric(
        x[[padj_col]]
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  colnames(out)[2:3] <- c(
    paste0(prefix, "_NES"),
    paste0(prefix, "_padj")
  )
  
  
  stopifnot(
    !any(
      duplicated(
        out$pathway
      )
    )
  )
  
  
  out
}


# ============================================================
# 2. FILE PATHS
# ============================================================

d1_def_file <- paste0(
  "results/pathway_analysis/",
  "Hallmark_GSEA_BRCA1_deficient_all.csv"
)

d1_prof_file <- paste0(
  "results/pathway_analysis/",
  "Hallmark_GSEA_BRCA1_proficient_all.csv"
)

d2_file <- paste0(
  "results/GSE117765/",
  "GSE117765_Hallmark_GSEA_all.csv"
)

d3_file <- paste0(
  "results/GSE153867/",
  "GSE153867_Hallmark_GSEA_all.csv"
)

d4_file <- paste0(
  "results/GSE342642/Hallmark_GSEA/",
  "GSE342642_Hallmark_GSEA_all.csv"
)

d5_file <- paste0(
  "results/GSE202434/D5_integration/",
  "D5_six_contrast_mechanistic_master.csv"
)

d6_file <- paste0(
  "results/GSE191231/Hallmark_GSEA/",
  "GSE191231_Hallmark_GSEA_primary_all.csv"
)

primary_master_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "FINAL_D1_D6_Hallmark_master.csv"
)


required_files <- c(
  d1_def_file,
  d1_prof_file,
  d2_file,
  d3_file,
  d4_file,
  d5_file,
  d6_file,
  primary_master_file
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (length(missing_files) > 0) {
  
  print(
    missing_files
  )
  
  stop(
    "One or more required files are missing."
  )
}


# ============================================================
# 3. LOAD D1-D4 + D6
# ============================================================

D1_def <- read_gsea(
  d1_def_file,
  "D1_def"
)

D1_prof <- read_gsea(
  d1_prof_file,
  "D1_prof"
)

D2 <- read_gsea(
  d2_file,
  "D2"
)

D3 <- read_gsea(
  d3_file,
  "D3"
)

D4 <- read_gsea(
  d4_file,
  "D4"
)

D6 <- read_gsea(
  d6_file,
  "D6"
)


# ============================================================
# 4. LOAD D5 MODEL-LEVEL CONTRASTS
# ============================================================

D5 <- read.csv(
  d5_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_D5 <- c(
  
  "pathway",
  
  "H1_NES",
  "H1_padj",
  
  "H2_NES",
  "H2_padj",
  
  "P1_NES",
  "P1_padj",
  
  "P1_sens_NES",
  "P1_sens_padj"
)


stopifnot(
  all(
    required_D5 %in%
      colnames(D5)
  )
)


D5 <- D5[
  ,
  required_D5
]


# Rename D5 contrasts clearly

colnames(D5)[
  colnames(D5) == "H1_NES"
] <- "D5_HCC_M_NES"

colnames(D5)[
  colnames(D5) == "H1_padj"
] <- "D5_HCC_M_padj"


colnames(D5)[
  colnames(D5) == "H2_NES"
] <- "D5_HCC_DF_NES"

colnames(D5)[
  colnames(D5) == "H2_padj"
] <- "D5_HCC_DF_padj"


colnames(D5)[
  colnames(D5) == "P1_NES"
] <- "D5_PSN_M_NES"

colnames(D5)[
  colnames(D5) == "P1_padj"
] <- "D5_PSN_M_padj"


colnames(D5)[
  colnames(D5) == "P1_sens_NES"
] <- "D5_PSN_M_sens_NES"

colnames(D5)[
  colnames(D5) == "P1_sens_padj"
] <- "D5_PSN_M_sens_padj"


# ============================================================
# 5. MERGE ALL NINE MODEL-LEVEL UNITS
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
    D1_def,
    D1_prof,
    D2,
    D3,
    D4,
    D5,
    D6
  )
)


cat("\n========================================\n")
cat("MODEL-LEVEL MASTER QC\n")
cat("========================================\n\n")


cat(
  "Hallmarks:",
  nrow(master),
  "\n"
)

cat(
  "Missing values:",
  sum(
    is.na(master)
  ),
  "\n"
)


stopifnot(
  nrow(master) == 50,
  sum(is.na(master)) == 0
)


# ============================================================
# 6. PSN1 P1 ROBUSTNESS
# ============================================================

master$D5_PSN_M_robust_direction <-
  
  sign(
    master$D5_PSN_M_NES
  ) ==
  sign(
    master$D5_PSN_M_sens_NES
  )


master$D5_PSN_M_robust_sig <-
  
  master$D5_PSN_M_robust_direction &
  
  master$D5_PSN_M_padj < 0.05 &
  
  master$D5_PSN_M_sens_padj < 0.05


# ============================================================
# 7. BUILD NINE FORMAL DIRECTION VOTES
#
# For PSN1 P1:
# direction vote only if primary and sensitivity agree.
#
# significance vote only if both are significant and agree.
# ============================================================


units <- c(
  
  "D1_def",
  "D1_prof",
  
  "D2",
  "D3",
  "D4",
  
  "D5_HCC_M",
  "D5_HCC_DF",
  
  "D5_PSN_M",
  
  "D6"
)


# ------------------------------------------------------------
# Ordinary units
# ------------------------------------------------------------

ordinary_units <- c(
  
  "D1_def",
  "D1_prof",
  
  "D2",
  "D3",
  "D4",
  
  "D5_HCC_M",
  "D5_HCC_DF",
  
  "D6"
)


for (unit in ordinary_units) {
  
  nes_col <- paste0(
    unit,
    "_NES"
  )
  
  padj_col <- paste0(
    unit,
    "_padj"
  )
  
  
  master[[paste0(unit, "_direction")]] <-
    sign(
      master[[nes_col]]
    )
  
  
  master[[paste0(unit, "_sig_direction")]] <-
    ifelse(
      master[[padj_col]] < 0.05,
      sign(master[[nes_col]]),
      0
    )
}


# ------------------------------------------------------------
# PSN1 P1 with outlier robustness requirement
# ------------------------------------------------------------

master$D5_PSN_M_direction <- ifelse(
  
  master$D5_PSN_M_robust_direction,
  
  sign(
    master$D5_PSN_M_NES
  ),
  
  0
)


master$D5_PSN_M_sig_direction <- ifelse(
  
  master$D5_PSN_M_robust_sig,
  
  sign(
    master$D5_PSN_M_NES
  ),
  
  0
)


# ============================================================
# 8. DIRECTION MATRICES
# ============================================================

direction_cols <- paste0(
  units,
  "_direction"
)

sig_cols <- paste0(
  units,
  "_sig_direction"
)


direction_matrix <- as.matrix(
  master[
    ,
    direction_cols
  ]
)


sig_matrix <- as.matrix(
  master[
    ,
    sig_cols
  ]
)


# ============================================================
# 9. COUNTS
# ============================================================

master$n_UP <- rowSums(
  direction_matrix == 1
)


master$n_DOWN <- rowSums(
  direction_matrix == -1
)


master$n_abstain <- rowSums(
  direction_matrix == 0
)


master$n_sig_UP <- rowSums(
  sig_matrix == 1
)


master$n_sig_DOWN <- rowSums(
  sig_matrix == -1
)


master$n_sig_total <-
  master$n_sig_UP +
  master$n_sig_DOWN


master$majority_direction <- ifelse(
  
  master$n_UP >
    master$n_DOWN,
  
  "UP",
  
  ifelse(
    
    master$n_DOWN >
      master$n_UP,
    
    "DOWN",
    
    "TIE"
  )
)


# ============================================================
# 10. MODEL-LEVEL CLASSIFICATION
#
# Nine units total.
#
# Strong recurrent:
# >=6/9 same direction
# >=4 significant in that direction
# zero significant opposite
#
# Directionally recurrent:
# >=6/9 same direction
# <4 significant in direction
# zero significant opposite
#
# Context-dependent:
# >=1 significant UP
# AND >=1 significant DOWN
#
# Unit-specific significant:
# exactly one significant unit
#
# Otherwise inconsistent/weak.
# ============================================================

master$model_level_classification <-
  "Inconsistent/weak"


# Context-dependent gets first priority

context_idx <-
  
  master$n_sig_UP >= 1 &
  
  master$n_sig_DOWN >= 1


master$model_level_classification[
  context_idx
] <- "Context-dependent"


# Strong UP

strong_up <-
  
  !context_idx &
  
  master$n_UP >= 6 &
  
  master$n_sig_UP >= 4 &
  
  master$n_sig_DOWN == 0


# Strong DOWN

strong_down <-
  
  !context_idx &
  
  master$n_DOWN >= 6 &
  
  master$n_sig_DOWN >= 4 &
  
  master$n_sig_UP == 0


master$model_level_classification[
  strong_up |
    strong_down
] <- "Strong recurrent"


# Directionally recurrent UP

directional_up <-
  
  !context_idx &
  
  !(strong_up | strong_down) &
  
  master$n_UP >= 6 &
  
  master$n_sig_DOWN == 0


# Directionally recurrent DOWN

directional_down <-
  
  !context_idx &
  
  !(strong_up | strong_down) &
  
  master$n_DOWN >= 6 &
  
  master$n_sig_UP == 0


master$model_level_classification[
  directional_up |
    directional_down
] <- "Directionally recurrent"


# Unit-specific significant

unit_specific <-
  
  master$model_level_classification ==
  "Inconsistent/weak" &
  
  master$n_sig_total == 1


master$model_level_classification[
  unit_specific
] <- "Unit-specific significant"


# ============================================================
# 11. MODEL-LEVEL SUMMARY
# ============================================================

class_order <- c(
  
  "Strong recurrent",
  
  "Directionally recurrent",
  
  "Context-dependent",
  
  "Unit-specific significant",
  
  "Inconsistent/weak"
)


model_counts <- data.frame(
  
  classification =
    class_order,
  
  n = sapply(
    
    class_order,
    
    function(x) {
      
      sum(
        master$model_level_classification ==
          x
      )
    }
  )
)


cat("\n========================================\n")
cat("MODEL-LEVEL CLASSIFICATION COUNTS\n")
cat("========================================\n\n")


print(
  model_counts,
  row.names = FALSE
)


# ============================================================
# 12. STRONG RECURRENT
# ============================================================

strong <- master[
  master$model_level_classification ==
    "Strong recurrent",
]


strong <- strong[
  order(
    -pmax(
      strong$n_UP,
      strong$n_DOWN
    ),
    -strong$n_sig_total
  ),
]


cat("\n========================================\n")
cat("MODEL-LEVEL STRONG RECURRENT PATHWAYS\n")
cat("========================================\n\n")


print(
  strong[
    ,
    c(
      "pathway",
      "majority_direction",
      "n_UP",
      "n_DOWN",
      "n_sig_UP",
      "n_sig_DOWN"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 13. CONTEXT-DEPENDENT
# ============================================================

context <- master[
  master$model_level_classification ==
    "Context-dependent",
]


context <- context[
  order(
    -context$n_sig_total
  ),
]


cat("\n========================================\n")
cat("MODEL-LEVEL CONTEXT-DEPENDENT PATHWAYS\n")
cat("========================================\n\n")


print(
  context[
    ,
    c(
      "pathway",
      "n_UP",
      "n_DOWN",
      "n_sig_UP",
      "n_sig_DOWN"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 14. LOAD PRIMARY ACCESSION-LEVEL CLASSIFICATION
# ============================================================

primary_master <- read.csv(
  primary_master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


primary_class <- primary_master[
  ,
  c(
    "pathway",
    "classification"
  )
]


colnames(primary_class)[2] <-
  "accession_level_classification"


comparison <- merge(
  
  master[
    ,
    c(
      "pathway",
      "model_level_classification"
    )
  ],
  
  primary_class,
  
  by = "pathway",
  all = TRUE
)


comparison$same_classification <-
  
  comparison$model_level_classification ==
  comparison$accession_level_classification


cat("\n========================================\n")
cat("ACCESSION vs MODEL-LEVEL CLASSIFICATION\n")
cat("========================================\n\n")


cat(
  "Same classification:",
  sum(
    comparison$same_classification
  ),
  "/ 50\n"
)


cat(
  "Classification agreement:",
  round(
    100 *
      mean(
        comparison$same_classification
      ),
    1
  ),
  "%\n"
)


# ============================================================
# 15. CROSS-TABULATION
# ============================================================

classification_cross_tab <- table(
  
  Accession =
    comparison$accession_level_classification,
  
  Model =
    comparison$model_level_classification
)


cat("\nClassification cross-tabulation:\n\n")

print(
  classification_cross_tab
)


# ============================================================
# 16. FROZEN D1 PATHWAYS
# ============================================================

frozen_pathways <- c(
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_GLYCOLYSIS",
  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_UV_RESPONSE_DN"
)


frozen <- merge(
  
  master[
    master$pathway %in%
      frozen_pathways,
    c(
      "pathway",
      "n_UP",
      "n_DOWN",
      "n_sig_UP",
      "n_sig_DOWN",
      "model_level_classification"
    )
  ],
  
  primary_class,
  
  by = "pathway"
)


cat("\n========================================\n")
cat("FROZEN D1 PATHWAYS:")
cat("\nACCESSION vs MODEL-LEVEL STATUS\n")
cat("========================================\n\n")


print(
  frozen,
  row.names = FALSE
)


# ============================================================
# 17. SPECIFICALLY TEST THE MAIN CLAIM
#
# How many pathways are context-dependent
# under BOTH integration strategies?
# ============================================================

both_context <- comparison[
  comparison$accession_level_classification ==
    "Context-dependent" &
    comparison$model_level_classification ==
    "Context-dependent",
]


cat("\n========================================\n")
cat("ROBUST CONTEXT-DEPENDENCE\n")
cat("========================================\n\n")


cat(
  "Context-dependent in BOTH approaches:",
  nrow(both_context),
  "/ 50\n\n"
)


print(
  both_context,
  row.names = FALSE
)


# ============================================================
# 18. GLYCOLYSIS SPECIFIC CHECK
# ============================================================

glycolysis <- master[
  master$pathway ==
    "HALLMARK_GLYCOLYSIS",
]


cat("\n========================================\n")
cat("GLYCOLYSIS MODEL-LEVEL CHECK\n")
cat("========================================\n\n")


print(
  glycolysis[
    ,
    c(
      "pathway",
      
      "D1_def_NES",
      "D1_prof_NES",
      
      "D2_NES",
      "D3_NES",
      "D4_NES",
      
      "D5_HCC_M_NES",
      "D5_HCC_DF_NES",
      "D5_PSN_M_NES",
      
      "D6_NES",
      
      "n_UP",
      "n_DOWN",
      "n_sig_UP",
      "n_sig_DOWN",
      
      "model_level_classification"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 19. SAVE
# ============================================================

write.csv(
  
  master,
  
  file.path(
    outdir,
    "FINAL_model_level_Hallmark_master.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  model_counts,
  
  file.path(
    outdir,
    "FINAL_model_level_classification_counts.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  comparison,
  
  file.path(
    outdir,
    "FINAL_accession_vs_model_classification.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  context,
  
  file.path(
    outdir,
    "FINAL_model_level_context_dependent.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  strong,
  
  file.path(
    outdir,
    "FINAL_model_level_strong_recurrent.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  frozen,
  
  file.path(
    outdir,
    "FINAL_frozen_pathways_accession_vs_model.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  both_context,
  
  file.path(
    outdir,
    "FINAL_robust_context_dependent_both_methods.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 20. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("MODEL-LEVEL SENSITIVITY COMPLETE\n")
cat("========================================\n\n")

cat(
  "Output directory:",
  outdir,
  "\n"
)