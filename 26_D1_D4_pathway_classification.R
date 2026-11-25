# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 4: FORMAL PATHWAY CLASSIFICATION
#
# Categories:
# 1. Context-dependent
# 2. Strong recurrent
# 3. Directionally recurrent
# 4. Dataset-specific significant
# 5. Inconsistent / weak
#
# IMPORTANT:
# This is a descriptive evidence-classification framework,
# NOT a formal statistical meta-analysis.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/integration_D1_D4"


# ============================================================
# 1. LOAD MASTER MATRIX
# ============================================================

master <- read.csv(
  file.path(
    outdir,
    "D1_D4_Hallmark_master_matrix.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 2. CONSERVATIVE D1 REPRESENTATION
#
# D1 consists of two biological contrasts.
#
# For formal cross-dataset classification:
#
# Direction counts ONLY if both D1 contrasts agree.
#
# Significant D1 support counts ONLY if:
# - both contrasts are significant
# - both point in the same direction
#
# If D1 contrasts oppose each other:
# D1 is labelled MIXED and contributes no directional vote.
# ============================================================

master$D1_formal_direction <- ifelse(
  master$D1_same_direction,
  master$D1_direction,
  "MIXED"
)


master$D1_formal_significant <-
  master$D1_consistent_significant


# ============================================================
# 3. SIGNIFICANT UP / DOWN FLAGS
# ============================================================

master$D1_sig_UP <-
  master$D1_formal_significant &
  master$D1_formal_direction == "UP"

master$D1_sig_DOWN <-
  master$D1_formal_significant &
  master$D1_formal_direction == "DOWN"


master$D2_sig_UP <-
  master$D2_significant &
  master$D2_direction == "UP"

master$D2_sig_DOWN <-
  master$D2_significant &
  master$D2_direction == "DOWN"


master$D3_sig_UP <-
  master$D3_significant &
  master$D3_direction == "UP"

master$D3_sig_DOWN <-
  master$D3_significant &
  master$D3_direction == "DOWN"


master$D4_sig_UP <-
  master$D4_significant &
  master$D4_direction == "UP"

master$D4_sig_DOWN <-
  master$D4_significant &
  master$D4_direction == "DOWN"


# ============================================================
# 4. COUNT SIGNIFICANT SUPPORT
# ============================================================

master$n_sig_UP <- rowSums(
  master[
    ,
    c(
      "D1_sig_UP",
      "D2_sig_UP",
      "D3_sig_UP",
      "D4_sig_UP"
    )
  ]
)

master$n_sig_DOWN <- rowSums(
  master[
    ,
    c(
      "D1_sig_DOWN",
      "D2_sig_DOWN",
      "D3_sig_DOWN",
      "D4_sig_DOWN"
    )
  ]
)

master$n_sig_total <-
  master$n_sig_UP +
  master$n_sig_DOWN


# ============================================================
# 5. CONSERVATIVE DIRECTION COUNTS
#
# D1 contributes only when its two contrasts agree.
# D2-D4 always contribute their NES direction.
# ============================================================

direction_vote <- function(x) {
  
  ifelse(
    x == "UP",
    1,
    ifelse(
      x == "DOWN",
      -1,
      0
    )
  )
}


dir_matrix <- cbind(
  
  direction_vote(
    master$D1_formal_direction
  ),
  
  direction_vote(
    master$D2_direction
  ),
  
  direction_vote(
    master$D3_direction
  ),
  
  direction_vote(
    master$D4_direction
  )
)


colnames(dir_matrix) <- c(
  "D1",
  "D2",
  "D3",
  "D4"
)


master$n_direction_UP <- rowSums(
  dir_matrix == 1
)

master$n_direction_DOWN <- rowSums(
  dir_matrix == -1
)

master$n_direction_mixed_or_missing <- rowSums(
  dir_matrix == 0
)


master$dominant_direction <- ifelse(
  
  master$n_direction_UP >
    master$n_direction_DOWN,
  
  "UP",
  
  ifelse(
    
    master$n_direction_DOWN >
      master$n_direction_UP,
    
    "DOWN",
    
    "TIE"
  )
)


master$dominant_direction_n <- pmax(
  master$n_direction_UP,
  master$n_direction_DOWN
)


# ============================================================
# 6. SUPPORT IN DOMINANT DIRECTION
# ============================================================

master$n_sig_dominant <- ifelse(
  
  master$dominant_direction == "UP",
  
  master$n_sig_UP,
  
  ifelse(
    master$dominant_direction == "DOWN",
    master$n_sig_DOWN,
    0
  )
)


master$n_sig_opposite <- ifelse(
  
  master$dominant_direction == "UP",
  
  master$n_sig_DOWN,
  
  ifelse(
    master$dominant_direction == "DOWN",
    master$n_sig_UP,
    master$n_sig_total
  )
)


# ============================================================
# 7. FORMAL CLASSIFICATION
#
# Order matters.
#
# A pathway with significant evidence in BOTH directions
# is context-dependent even if 3/4 datasets point one way.
# ============================================================

master$class <- "Inconsistent / weak"


# ------------------------------------------------------------
# A. CONTEXT-DEPENDENT
#
# Significant positive enrichment in >=1 dataset
# AND significant negative enrichment in >=1 dataset.
# ------------------------------------------------------------

context_dependent <-
  master$n_sig_UP >= 1 &
  master$n_sig_DOWN >= 1

master$class[
  context_dependent
] <- "Context-dependent"


# ------------------------------------------------------------
# B. STRONG RECURRENT
#
# - At least 3 datasets support same direction
# - At least 2 datasets significant in that direction
# - Zero significant opposite-direction datasets
#
# Do not overwrite context-dependent pathways.
# ------------------------------------------------------------

strong_recurrent <-
  !context_dependent &
  master$dominant_direction_n >= 3 &
  master$n_sig_dominant >= 2 &
  master$n_sig_opposite == 0

master$class[
  strong_recurrent
] <- "Strong recurrent"


# ------------------------------------------------------------
# C. DIRECTIONALLY RECURRENT
#
# - At least 3 datasets support same direction
# - Fewer than 2 significant datasets in that direction
# - No significant opposite-direction evidence
# ------------------------------------------------------------

directional_recurrent <-
  !context_dependent &
  !strong_recurrent &
  master$dominant_direction_n >= 3 &
  master$n_sig_opposite == 0

master$class[
  directional_recurrent
] <- "Directionally recurrent"


# ------------------------------------------------------------
# D. DATASET-SPECIFIC SIGNIFICANT
#
# Exactly one significant dataset overall,
# without enough directional support for recurrent category.
# ------------------------------------------------------------

dataset_specific <-
  master$class == "Inconsistent / weak" &
  master$n_sig_total == 1

master$class[
  dataset_specific
] <- "Dataset-specific significant"


# ============================================================
# 8. ADD HUMAN-READABLE SUPPORT SUMMARY
# ============================================================

master$support_summary <- paste0(
  
  "Direction UP=",
  master$n_direction_UP,
  
  "; DOWN=",
  master$n_direction_DOWN,
  
  "; sig UP=",
  master$n_sig_UP,
  
  "; sig DOWN=",
  master$n_sig_DOWN
)


# ============================================================
# 9. CATEGORY COUNTS
# ============================================================

class_order <- c(
  "Strong recurrent",
  "Directionally recurrent",
  "Context-dependent",
  "Dataset-specific significant",
  "Inconsistent / weak"
)

master$class <- factor(
  master$class,
  levels = class_order
)


cat("\n========================================\n")
cat("PATHWAY CLASSIFICATION COUNTS\n")
cat("========================================\n\n")

print(
  table(
    master$class,
    useNA = "ifany"
  )
)


# ============================================================
# 10. STRONG RECURRENT PATHWAYS
# ============================================================

strong_table <- master[
  master$class == "Strong recurrent",
  c(
    "pathway",
    "dominant_direction",
    "dominant_direction_n",
    "n_sig_dominant",
    "D1_formal_direction",
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES",
    "D1_formal_significant",
    "D2_significant",
    "D3_significant",
    "D4_significant"
  )
]


strong_table <- strong_table[
  order(
    strong_table$dominant_direction,
    -strong_table$n_sig_dominant,
    strong_table$pathway
  ),
]


cat("\n========================================\n")
cat("STRONG RECURRENT PATHWAYS\n")
cat("========================================\n\n")

print(
  strong_table,
  row.names = FALSE
)


# ============================================================
# 11. CONTEXT-DEPENDENT PATHWAYS
# ============================================================

context_table <- master[
  master$class == "Context-dependent",
  c(
    "pathway",
    "n_sig_UP",
    "n_sig_DOWN",
    "D1_formal_direction",
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES",
    "D1_formal_significant",
    "D2_significant",
    "D3_significant",
    "D4_significant"
  )
]


context_table <- context_table[
  order(
    -context_table$n_sig_UP -
      context_table$n_sig_DOWN,
    context_table$pathway
  ),
]


cat("\n========================================\n")
cat("CONTEXT-DEPENDENT PATHWAYS\n")
cat("========================================\n\n")

print(
  context_table,
  row.names = FALSE
)


# ============================================================
# 12. DIRECTIONALLY RECURRENT PATHWAYS
# ============================================================

directional_table <- master[
  master$class == "Directionally recurrent",
  c(
    "pathway",
    "dominant_direction",
    "dominant_direction_n",
    "n_sig_dominant",
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES"
  )
]


cat("\n========================================\n")
cat("DIRECTIONALLY RECURRENT PATHWAYS\n")
cat("========================================\n\n")

print(
  directional_table,
  row.names = FALSE
)


# ============================================================
# 13. SAVE ALL TABLES
# ============================================================

write.csv(
  master,
  file.path(
    outdir,
    "D1_D4_Hallmark_classification.csv"
  ),
  row.names = FALSE
)


write.csv(
  strong_table,
  file.path(
    outdir,
    "D1_D4_strong_recurrent.csv"
  ),
  row.names = FALSE
)


write.csv(
  context_table,
  file.path(
    outdir,
    "D1_D4_context_dependent.csv"
  ),
  row.names = FALSE
)


write.csv(
  directional_table,
  file.path(
    outdir,
    "D1_D4_directionally_recurrent.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. FINISHED
# ============================================================

cat("\n========================================\n")
cat("PATHWAY CLASSIFICATION COMPLETE\n")
cat("========================================\n")