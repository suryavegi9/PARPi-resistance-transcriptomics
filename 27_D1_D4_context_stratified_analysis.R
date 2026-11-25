# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 5: CONTEXT-STRATIFIED EXPLORATORY ANALYSIS
#
# Main exploratory question:
# Do Hallmark pathway patterns align with treatment state
# at RNA collection?
#
# Drug-free / untreated at collection:
#   D1, D4
#
# Ongoing olaparib exposure:
#   D2, D3
#
# IMPORTANT:
# n = 2 datasets per exposure group.
# This is hypothesis-generating only.
# It CANNOT establish causality.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/integration_D1_D4"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. LOAD CLASSIFIED MASTER TABLE
# ============================================================

master <- read.csv(
  file.path(
    outdir,
    "D1_D4_Hallmark_classification.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. LOCK DATASET CONTEXT METADATA
# ============================================================

context <- data.frame(
  
  dataset = c(
    "D1",
    "D2",
    "D3",
    "D4"
  ),
  
  accession = c(
    "GSE235980",
    "GSE117765",
    "GSE153867",
    "GSE342642"
  ),
  
  biological_model = c(
    "UWB1.289 / UWB1.289+BRCA1",
    "PEO1",
    "A2780",
    "HGSOC PDX"
  ),
  
  model_type = c(
    "Cell line",
    "Cell line",
    "Cell line",
    "PDX"
  ),
  
  environment = c(
    "In vitro",
    "In vitro",
    "In vitro",
    "In vivo"
  ),
  
  BRCA_context = c(
    "BRCA1-deficient + BRCA1-proficient isogenic contexts",
    "BRCA2-mutant",
    "Not used as integration stratifier",
    "Germline BRCA1-mutant"
  ),
  
  treatment_at_collection = c(
    "Drug-free",
    "Ongoing olaparib",
    "Ongoing olaparib",
    "Untreated acquired-resistant tumor"
  ),
  
  exposure_group = c(
    "Drug-free_or_untreated",
    "Ongoing_olaparib",
    "Ongoing_olaparib",
    "Drug-free_or_untreated"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  context,
  file.path(
    outdir,
    "D1_D4_context_metadata.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("DATASET CONTEXT METADATA\n")
cat("========================================\n\n")

print(
  context,
  row.names = FALSE
)


# ============================================================
# 3. EXPOSURE-GROUP NES SUMMARIES
#
# Drug-free / untreated:
#   D1 + D4
#
# Ongoing olaparib:
#   D2 + D3
#
# These are DESCRIPTIVE means only.
# ============================================================

master$NES_drugfree_mean <- rowMeans(
  master[
    ,
    c(
      "D1_NES",
      "D4_NES"
    )
  ]
)


master$NES_ongoing_mean <- rowMeans(
  master[
    ,
    c(
      "D2_NES",
      "D3_NES"
    )
  ]
)


master$exposure_NES_difference <-
  master$NES_drugfree_mean -
  master$NES_ongoing_mean


master$abs_exposure_difference <-
  abs(
    master$exposure_NES_difference
  )


# ============================================================
# 4. WITHIN-GROUP DIRECTION CONSISTENCY
# ============================================================

master$drugfree_same_direction <-
  sign(master$D1_NES) ==
  sign(master$D4_NES)


master$ongoing_same_direction <-
  sign(master$D2_NES) ==
  sign(master$D3_NES)


# ============================================================
# 5. EXPOSURE-GROUP DIRECTIONS
# ============================================================

direction_label <- function(x) {
  
  ifelse(
    x > 0,
    "UP",
    ifelse(
      x < 0,
      "DOWN",
      "ZERO"
    )
  )
}


master$drugfree_direction <- direction_label(
  master$NES_drugfree_mean
)

master$ongoing_direction <- direction_label(
  master$NES_ongoing_mean
)


# ============================================================
# 6. CLEAN EXPOSURE-ALIGNED REVERSAL
#
# Strict descriptive criterion:
#
# 1. D1 and D4 agree with one another
# 2. D2 and D3 agree with one another
# 3. The two exposure groups point in opposite directions
#
# Example:
#
# D1 DOWN
# D4 DOWN
#
# D2 UP
# D3 UP
#
# This is the cleanest pattern that COULD be consistent
# with treatment-at-collection effects.
#
# It still does NOT prove exposure caused the reversal.
# ============================================================

master$clean_exposure_reversal <-
  
  master$drugfree_same_direction &
  
  master$ongoing_same_direction &
  
  master$drugfree_direction !=
  master$ongoing_direction


# ============================================================
# 7. SIGNIFICANCE SUPPORT WITHIN EACH EXPOSURE GROUP
#
# D1 significance is conservative:
# both D1 contrasts significant + same direction.
# ============================================================

master$drugfree_sig_n <- rowSums(
  cbind(
    master$D1_formal_significant,
    master$D4_significant
  )
)


master$ongoing_sig_n <- rowSums(
  cbind(
    master$D2_significant,
    master$D3_significant
  )
)


master$total_sig_n <-
  master$drugfree_sig_n +
  master$ongoing_sig_n


# ============================================================
# 8. STRONG EXPOSURE-ALIGNED CANDIDATES
#
# Exploratory definition:
#
# clean 2-vs-2 directional reversal
#
# AND
#
# at least one significant dataset
# from EACH exposure group.
#
# Still hypothesis-generating.
# ============================================================

master$strong_exposure_candidate <-
  
  master$clean_exposure_reversal &
  
  master$drugfree_sig_n >= 1 &
  
  master$ongoing_sig_n >= 1


# ============================================================
# 9. GLOBAL ALIGNMENT SUMMARY
# ============================================================

cat("\n========================================\n")
cat("TREATMENT-STATE ALIGNMENT SUMMARY\n")
cat("========================================\n\n")


cat(
  "D1-D4 same direction (drug-free/untreated pair):",
  sum(master$drugfree_same_direction),
  "/ 50\n"
)


cat(
  "D2-D3 same direction (ongoing-olaparib pair):",
  sum(master$ongoing_same_direction),
  "/ 50\n"
)


cat(
  "Clean 2-vs-2 exposure-aligned reversals:",
  sum(master$clean_exposure_reversal),
  "/ 50\n"
)


cat(
  "Strong exposure-aligned candidates:",
  sum(master$strong_exposure_candidate),
  "/ 50\n"
)


# ============================================================
# 10. CLEAN EXPOSURE-ALIGNED REVERSALS
# ============================================================

reversal_table <- master[
  master$clean_exposure_reversal,
  c(
    "pathway",
    
    "D1_NES",
    "D4_NES",
    
    "D2_NES",
    "D3_NES",
    
    "NES_drugfree_mean",
    "NES_ongoing_mean",
    
    "exposure_NES_difference",
    
    "drugfree_sig_n",
    "ongoing_sig_n",
    
    "strong_exposure_candidate",
    
    "class"
  )
]


 

# Fix ordering safely
reversal_table <- reversal_table[
  order(
    -abs(
      reversal_table$exposure_NES_difference
    )
  ),
]


cat("\n========================================\n")
cat("CLEAN 2-vs-2 EXPOSURE-ALIGNED REVERSALS\n")
cat("========================================\n\n")

print(
  reversal_table,
  row.names = FALSE
)


# ============================================================
# 11. STRONG EXPOSURE-ALIGNED CANDIDATES
# ============================================================

strong_exposure <- master[
  master$strong_exposure_candidate,
  c(
    "pathway",
    
    "D1_NES",
    "D4_NES",
    
    "D2_NES",
    "D3_NES",
    
    "NES_drugfree_mean",
    "NES_ongoing_mean",
    
    "exposure_NES_difference",
    
    "drugfree_sig_n",
    "ongoing_sig_n",
    
    "class"
  )
]


strong_exposure <- strong_exposure[
  order(
    -abs(
      strong_exposure$exposure_NES_difference
    )
  ),
]


cat("\n========================================\n")
cat("STRONG EXPOSURE-ALIGNED CANDIDATES\n")
cat("========================================\n\n")

print(
  strong_exposure,
  row.names = FALSE
)


# ============================================================
# 12. TOP EXPOSURE-GROUP DIFFERENCES
#
# Regardless of whether they satisfy clean 2-vs-2 reversal.
# ============================================================

top_diff <- master[
  order(
    -master$abs_exposure_difference
  ),
  c(
    "pathway",
    
    "D1_NES",
    "D4_NES",
    
    "D2_NES",
    "D3_NES",
    
    "NES_drugfree_mean",
    "NES_ongoing_mean",
    
    "exposure_NES_difference",
    
    "drugfree_same_direction",
    "ongoing_same_direction",
    
    "class"
  )
]


top_diff <- head(
  top_diff,
  15
)


cat("\n========================================\n")
cat("TOP 15 EXPOSURE-GROUP NES DIFFERENCES\n")
cat("========================================\n\n")

print(
  top_diff,
  row.names = FALSE
)


# ============================================================
# 13. SAVE OUTPUTS
# ============================================================

write.csv(
  master,
  file.path(
    outdir,
    "D1_D4_context_stratified_master.csv"
  ),
  row.names = FALSE
)


write.csv(
  reversal_table,
  file.path(
    outdir,
    "D1_D4_clean_exposure_reversals.csv"
  ),
  row.names = FALSE
)


write.csv(
  strong_exposure,
  file.path(
    outdir,
    "D1_D4_strong_exposure_candidates.csv"
  ),
  row.names = FALSE
)


write.csv(
  top_diff,
  file.path(
    outdir,
    "D1_D4_top_exposure_group_differences.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("CONTEXT-STRATIFIED ANALYSIS COMPLETE\n")
cat("========================================\n")