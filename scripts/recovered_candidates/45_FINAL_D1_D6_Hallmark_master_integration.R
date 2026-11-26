# ============================================================
# FINAL D1-D6 HALLMARK INTEGRATION
#
# PURPOSE
# -------
# Build one accession-level vote per dataset:
#
# D1 = GSE235980
# D2 = GSE117765
# D3 = GSE153867
# D4 = GSE342642
# D5 = GSE202434
# D6 = GSE191231
#
# IMPORTANT:
# D1 and D5 contain multiple biological submodels.
# They are collapsed conservatively so they do NOT receive
# multiple votes in the final cross-dataset integration.
#
# This script performs:
#   1. source-table QC
#   2. accession-level master matrix
#   3. direction/significance voting
#   4. final pathway classification
#
# No biological interpretation is forced by the script.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. OUTPUT
# ============================================================

outdir <- "results/FINAL_D1_D6_integration"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. HELPER: READ STANDARD GSEA FILE
# ============================================================

read_gsea <- function(
    file,
    prefix
) {
  
  stopifnot(
    file.exists(file)
  )
  
  x <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  # ------------------------------
  # Pathway column
  # ------------------------------
  
  if ("ID" %in% colnames(x)) {
    
    pathway_col <- "ID"
    
  } else if ("pathway" %in% colnames(x)) {
    
    pathway_col <- "pathway"
    
  } else {
    
    stop(
      paste(
        "Could not identify pathway column in:",
        file
      )
    )
  }
  
  
  # ------------------------------
  # NES column
  # ------------------------------
  
  if (!"NES" %in% colnames(x)) {
    
    stop(
      paste(
        "NES column missing in:",
        file
      )
    )
  }
  
  
  # ------------------------------
  # adjusted-P column
  # ------------------------------
  
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
  
  
  if (
    any(
      duplicated(
        out$pathway
      )
    )
  ) {
    
    stop(
      paste(
        "Duplicated pathways in:",
        file
      )
    )
  }
  
  
  out
}


# ============================================================
# 2. FILE PATHS
# ============================================================

# ------------------------------
# D1
# ------------------------------

d1_def_file <- paste0(
  "results/pathway_analysis/",
  "Hallmark_GSEA_BRCA1_deficient_all.csv"
)

d1_prof_file <- paste0(
  "results/pathway_analysis/",
  "Hallmark_GSEA_BRCA1_proficient_all.csv"
)


# ------------------------------
# D2
# ------------------------------

d2_file <- paste0(
  "results/GSE117765/",
  "GSE117765_Hallmark_GSEA_all.csv"
)


# ------------------------------
# D3
# ------------------------------

d3_file <- paste0(
  "results/GSE153867/",
  "GSE153867_Hallmark_GSEA_all.csv"
)


# ------------------------------
# D4
# ------------------------------

d4_file <- paste0(
  "results/GSE342642/Hallmark_GSEA/",
  "GSE342642_Hallmark_GSEA_all.csv"
)


# ------------------------------
# D5
# ------------------------------

d5_file <- paste0(
  "results/GSE202434/D5_integration/",
  "D5_six_contrast_mechanistic_master.csv"
)


# ------------------------------
# D6 primary
# ------------------------------

d6_primary_file <- paste0(
  "results/GSE191231/Hallmark_GSEA/",
  "GSE191231_Hallmark_GSEA_primary_all.csv"
)


# ------------------------------
# D6 sensitivity
# ------------------------------

d6_sens_file <- paste0(
  "results/GSE191231/Hallmark_GSEA/",
  "GSE191231_Hallmark_GSEA_sensitivity_all.csv"
)


all_files <- c(
  d1_def_file,
  d1_prof_file,
  d2_file,
  d3_file,
  d4_file,
  d5_file,
  d6_primary_file,
  d6_sens_file
)


missing_files <- all_files[
  !file.exists(
    all_files
  )
]


if (
  length(
    missing_files
  ) > 0
) {
  
  cat(
    "\nMissing files:\n"
  )
  
  print(
    missing_files
  )
  
  stop(
    "One or more input files are missing."
  )
}


# ============================================================
# 3. LOAD D1-D4
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


# ============================================================
# 4. LOAD D6
# ============================================================

D6_primary <- read_gsea(
  d6_primary_file,
  "D6"
)


D6_sens <- read_gsea(
  d6_sens_file,
  "D6_sens"
)


# ============================================================
# 5. LOAD D5 MASTER
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


# ============================================================
# 6. MASTER MERGE
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
    D6_primary,
    D6_sens
  )
)


cat("\n========================================\n")
cat("SOURCE TABLE QC\n")
cat("========================================\n\n")


cat(
  "Master Hallmarks:",
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
  nrow(master) == 50
)


stopifnot(
  sum(
    is.na(master)
  ) == 0
)


# ============================================================
# 7. D1 ACCESSION-LEVEL COLLAPSE
#
# Descriptive NES:
# mean of deficient + proficient
#
# Formal direction:
# only if both BRCA1 contexts agree.
#
# Formal significance:
# both contexts significant AND same direction.
#
# Mean NES is DESCRIPTIVE, not meta-analysis.
# ============================================================

master$D1_NES <- rowMeans(
  master[
    ,
    c(
      "D1_def_NES",
      "D1_prof_NES"
    )
  ]
)


master$D1_direction <- ifelse(
  
  sign(
    master$D1_def_NES
  ) ==
    sign(
      master$D1_prof_NES
    ),
  
  sign(
    master$D1_def_NES
  ),
  
  0
)


master$D1_sig_direction <- ifelse(
  
  master$D1_direction != 0 &
    
    master$D1_def_padj < 0.05 &
    
    master$D1_prof_padj < 0.05,
  
  master$D1_direction,
  
  0
)


# ============================================================
# 8. D2-D4 ACCESSION VOTES
# ============================================================

for (dataset in c("D2", "D3", "D4")) {
  
  nes_col <- paste0(
    dataset,
    "_NES"
  )
  
  padj_col <- paste0(
    dataset,
    "_padj"
  )
  
  
  # Direction irrespective of significance
  master[[paste0(dataset, "_direction")]] <-
    sign(
      master[[nes_col]]
    )
  
  
  # Significant direction only
  master[[paste0(dataset, "_sig_direction")]] <-
    ifelse(
      master[[padj_col]] < 0.05,
      sign(master[[nes_col]]),
      0
    )
}


# ============================================================
# 9. D5 ACCESSION-LEVEL COLLAPSE
#
# Biological resistant-state inputs:
#
# H1 = HCC1806 maintained resistant vs parental
# H2 = HCC1806 drug-free persistent resistant vs parental
# P1 = PSN1 maintained resistant vs parental
#
# P1 sensitivity verifies robustness to maintain-2.
#
# P2 is NOT included because after withdrawal PSN1
# substantially re-sensitizes and therefore does not represent
# the persistent resistant phenotype.
#
# Formal D5 direction requires:
# H1, H2, P1 and P1 sensitivity ALL agree.
#
# Formal D5 significance requires:
# H1, H2, P1 and P1 sensitivity ALL FDR < .05
# and same direction.
#
# This is intentionally conservative.
# ============================================================

master$D5_NES <- rowMeans(
  master[
    ,
    c(
      "H1_NES",
      "H2_NES",
      "P1_NES"
    )
  ]
)


D5_same_direction <-
  
  sign(
    master$H1_NES
  ) ==
  sign(
    master$H2_NES
  ) &
  
  sign(
    master$H1_NES
  ) ==
  sign(
    master$P1_NES
  ) &
  
  sign(
    master$P1_NES
  ) ==
  sign(
    master$P1_sens_NES
  )


master$D5_direction <- ifelse(
  
  D5_same_direction,
  
  sign(
    master$H1_NES
  ),
  
  0
)


master$D5_sig_direction <- ifelse(
  
  D5_same_direction &
    
    master$H1_padj < 0.05 &
    
    master$H2_padj < 0.05 &
    
    master$P1_padj < 0.05 &
    
    master$P1_sens_padj < 0.05,
  
  sign(
    master$H1_NES
  ),
  
  0
)


# ============================================================
# 10. D6 ACCESSION VOTE
#
# Primary analysis determines formal D6 vote.
# Sensitivity is retained separately as robustness evidence.
# ============================================================

master$D6_direction <-
  sign(
    master$D6_NES
  )


master$D6_sig_direction <- ifelse(
  
  master$D6_padj <
    0.05,
  
  master$D6_direction,
  
  0
)


master$D6_sensitivity_same_direction <-
  
  sign(
    master$D6_NES
  ) ==
  sign(
    master$D6_sens_NES
  )


master$D6_robust_significant <-
  
  master$D6_sensitivity_same_direction &
  
  master$D6_padj <
  0.05 &
  
  master$D6_sens_padj <
  0.05


# ============================================================
# 11. ACCESSION-LEVEL DIRECTION MATRIX
# ============================================================

direction_cols <- c(
  
  "D1_direction",
  "D2_direction",
  "D3_direction",
  "D4_direction",
  "D5_direction",
  "D6_direction"
)


sig_cols <- c(
  
  "D1_sig_direction",
  "D2_sig_direction",
  "D3_sig_direction",
  "D4_sig_direction",
  "D5_sig_direction",
  "D6_sig_direction"
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
# 12. FINAL DIRECTION COUNTS
# ============================================================

master$n_UP <- rowSums(
  direction_matrix == 1
)


master$n_DOWN <- rowSums(
  direction_matrix == -1
)


master$n_no_formal_direction <- rowSums(
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


master$max_direction_support <- pmax(
  master$n_UP,
  master$n_DOWN
)


master$max_sig_support <- pmax(
  master$n_sig_UP,
  master$n_sig_DOWN
)


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
# 13. FINAL CLASSIFICATION RULES
#
# CONTEXT-DEPENDENT
#   >=1 significantly UP accession
#   AND
#   >=1 significantly DOWN accession
#
# STRONG RECURRENT
#   >=4/6 accession-level direction support
#   >=3 significant accession-level supports
#   zero significant accessions in opposite direction
#
# DIRECTIONALLY RECURRENT
#   >=4/6 direction support
#   fewer than 3 significant supports
#   zero significant opposite supports
#
# DATASET-SPECIFIC SIGNIFICANT
#   exactly one accession significant
#
# otherwise:
#   inconsistent / weak
#
# These categories describe transcriptomic recurrence,
# NOT causal resistance mechanisms.
# ============================================================

master$classification <- "Inconsistent/weak"


# ------------------------------
# Context-dependent first
# ------------------------------

context_idx <-
  
  master$n_sig_UP >= 1 &
  
  master$n_sig_DOWN >= 1


master$classification[
  context_idx
] <- "Context-dependent"


# ------------------------------
# Strong recurrent
# ------------------------------

strong_up <-
  
  !context_idx &
  
  master$n_UP >= 4 &
  
  master$n_sig_UP >= 3 &
  
  master$n_sig_DOWN == 0


strong_down <-
  
  !context_idx &
  
  master$n_DOWN >= 4 &
  
  master$n_sig_DOWN >= 3 &
  
  master$n_sig_UP == 0


master$classification[
  strong_up |
    strong_down
] <- "Strong recurrent"


# ------------------------------
# Directionally recurrent
# ------------------------------

directional_up <-
  
  !context_idx &
  
  !(strong_up | strong_down) &
  
  master$n_UP >= 4 &
  
  master$n_sig_DOWN == 0


directional_down <-
  
  !context_idx &
  
  !(strong_up | strong_down) &
  
  master$n_DOWN >= 4 &
  
  master$n_sig_UP == 0


master$classification[
  directional_up |
    directional_down
] <- "Directionally recurrent"


# ------------------------------
# Dataset-specific significant
# ------------------------------

dataset_specific <-
  
  master$classification ==
  "Inconsistent/weak" &
  
  master$n_sig_total ==
  1


master$classification[
  dataset_specific
] <- "Dataset-specific significant"


# ============================================================
# 14. CLASSIFICATION SUMMARY
# ============================================================

class_order <- c(
  
  "Strong recurrent",
  
  "Directionally recurrent",
  
  "Context-dependent",
  
  "Dataset-specific significant",
  
  "Inconsistent/weak"
)


classification_counts <- data.frame(
  
  classification =
    class_order,
  
  n = sapply(
    
    class_order,
    
    function(x) {
      
      sum(
        master$classification ==
          x
      )
    }
  )
)


cat("\n========================================\n")
cat("FINAL D1-D6 CLASSIFICATION COUNTS\n")
cat("========================================\n\n")


print(
  classification_counts,
  row.names = FALSE
)


# ============================================================
# 15. STRONG RECURRENT PATHWAYS
# ============================================================

strong <- master[
  master$classification ==
    "Strong recurrent",
]


strong <- strong[
  order(
    -strong$max_direction_support,
    -strong$max_sig_support
  ),
]


cat("\n========================================\n")
cat("STRONG RECURRENT PATHWAYS\n")
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
# 16. DIRECTIONALLY RECURRENT PATHWAYS
# ============================================================

directional <- master[
  master$classification ==
    "Directionally recurrent",
]


directional <- directional[
  order(
    -directional$max_direction_support,
    -directional$max_sig_support
  ),
]


cat("\n========================================\n")
cat("DIRECTIONALLY RECURRENT PATHWAYS\n")
cat("========================================\n\n")


print(
  directional[
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
# 17. CONTEXT-DEPENDENT PATHWAYS
# ============================================================

context <- master[
  master$classification ==
    "Context-dependent",
]


context <- context[
  order(
    -context$n_sig_total,
    -context$max_direction_support
  ),
]


cat("\n========================================\n")
cat("CONTEXT-DEPENDENT PATHWAYS\n")
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
# 18. PRE-SPECIFIED D1 PATHWAYS
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


frozen_summary <- master[
  master$pathway %in%
    frozen_pathways,
]


cat("\n========================================\n")
cat("FINAL STATUS OF FROZEN D1 PATHWAYS\n")
cat("========================================\n\n")


print(
  frozen_summary[
    ,
    c(
      "pathway",
      
      "D1_NES",
      "D2_NES",
      "D3_NES",
      "D4_NES",
      "D5_NES",
      "D6_NES",
      
      "n_UP",
      "n_DOWN",
      
      "n_sig_UP",
      "n_sig_DOWN",
      
      "classification"
    )
  ],
  row.names = FALSE
)


# ============================================================
# 19. D6 ROBUSTNESS SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D6 ROBUSTNESS WITHIN FINAL MASTER\n")
cat("========================================\n\n")


cat(
  "D6 same direction primary vs sensitivity:",
  sum(
    master$D6_sensitivity_same_direction
  ),
  "/ 50\n"
)


cat(
  "D6 robust significant:",
  sum(
    master$D6_robust_significant
  ),
  "/ 50\n"
)


# ============================================================
# 20. D5 CONSENSUS SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D5 ACCESSION-LEVEL CONSENSUS\n")
cat("========================================\n\n")


cat(
  "D5 pathways with formal resistant-state direction:",
  sum(
    master$D5_direction != 0
  ),
  "/ 50\n"
)


cat(
  "D5 pathways significant across H1/H2/P1 + sensitivity:",
  sum(
    master$D5_sig_direction != 0
  ),
  "/ 50\n"
)


# ============================================================
# 21. SAVE MASTER TABLE
# ============================================================

write.csv(
  
  master,
  
  file.path(
    outdir,
    "FINAL_D1_D6_Hallmark_master.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  classification_counts,
  
  file.path(
    outdir,
    "FINAL_D1_D6_classification_counts.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  strong,
  
  file.path(
    outdir,
    "FINAL_D1_D6_strong_recurrent.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  directional,
  
  file.path(
    outdir,
    "FINAL_D1_D6_directionally_recurrent.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  context,
  
  file.path(
    outdir,
    "FINAL_D1_D6_context_dependent.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  frozen_summary,
  
  file.path(
    outdir,
    "FINAL_D1_D6_frozen_D1_pathway_status.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 22. COMPLETE
# ============================================================

cat("\n========================================\n")
cat("FINAL D1-D6 MASTER INTEGRATION COMPLETE\n")
cat("========================================\n\n")

cat(
  "Output directory:",
  outdir,
  "\n"
)