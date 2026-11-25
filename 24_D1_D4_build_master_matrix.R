# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 2: BUILD MASTER 50-HALLMARK MATRIX
#
# D1 = GSE235980
#      BRCA1-deficient + BRCA1-proficient contrasts
#
# D2 = GSE117765
# D3 = GSE153867
# D4 = GSE342642
#
# IMPORTANT:
# This script builds the integration table.
# It does NOT classify pathways as recurrent/context-specific yet.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/integration_D1_D4"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. LOAD PRIMARY HALLMARK TABLES
# ============================================================

d1_def <- read.csv(
  "results/pathway_analysis/Hallmark_GSEA_BRCA1_deficient_all.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

d1_prof <- read.csv(
  "results/pathway_analysis/Hallmark_GSEA_BRCA1_proficient_all.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

d2 <- read.csv(
  "results/GSE117765/GSE117765_Hallmark_GSEA_all.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

d3 <- read.csv(
  "results/GSE153867/GSE153867_Hallmark_GSEA_all.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

d4 <- read.csv(
  paste0(
    "results/GSE342642/Hallmark_GSEA/",
    "GSE342642_Hallmark_GSEA_all.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 2. BASIC QC
# ============================================================

cat("\n========================================\n")
cat("SOURCE TABLE QC\n")
cat("========================================\n\n")

cat("D1 deficient:", nrow(d1_def), "Hallmarks\n")
cat("D1 proficient:", nrow(d1_prof), "Hallmarks\n")
cat("D2:", nrow(d2), "Hallmarks\n")
cat("D3:", nrow(d3), "Hallmarks\n")
cat("D4:", nrow(d4), "Hallmarks\n")


# All primary Hallmark analyses should contain 50 pathways

stopifnot(
  nrow(d1_def) == 50,
  nrow(d1_prof) == 50,
  nrow(d2) == 50,
  nrow(d3) == 50,
  nrow(d4) == 50
)


# Check duplicate pathway IDs

stopifnot(
  sum(duplicated(d1_def$ID)) == 0,
  sum(duplicated(d1_prof$ID)) == 0,
  sum(duplicated(d2$ID)) == 0,
  sum(duplicated(d3$ID)) == 0,
  sum(duplicated(d4$ID)) == 0
)


# ============================================================
# 3. CREATE CLEAN DATASET TABLES
# ============================================================

D1_def <- data.frame(
  
  pathway = d1_def$ID,
  
  D1_def_NES = d1_def$NES,
  
  D1_def_padj = d1_def$p.adjust,
  
  stringsAsFactors = FALSE
)


D1_prof <- data.frame(
  
  pathway = d1_prof$ID,
  
  D1_prof_NES = d1_prof$NES,
  
  D1_prof_padj = d1_prof$p.adjust,
  
  stringsAsFactors = FALSE
)


D2 <- data.frame(
  
  pathway = d2$ID,
  
  D2_NES = d2$NES,
  
  D2_padj = d2$p.adjust,
  
  stringsAsFactors = FALSE
)


D3 <- data.frame(
  
  pathway = d3$ID,
  
  D3_NES = d3$NES,
  
  D3_padj = d3$p.adjust,
  
  stringsAsFactors = FALSE
)


D4 <- data.frame(
  
  pathway = d4$ID,
  
  D4_NES = d4$NES,
  
  D4_padj = d4$p.adjust,
  
  stringsAsFactors = FALSE
)


# ============================================================
# 4. MERGE ALL PRIMARY TABLES
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
    D4
  )
)


cat("\n========================================\n")
cat("MERGED MATRIX QC\n")
cat("========================================\n\n")

cat(
  "Merged Hallmarks:",
  nrow(master),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(master)),
  "\n"
)


# We expect exactly 50 complete Hallmark pathways

stopifnot(
  nrow(master) == 50
)

if (sum(is.na(master)) > 0) {
  
  cat("\nWARNING: Missing values detected.\n")
  
  print(
    master[
      !complete.cases(master),
    ]
  )
  
}


# ============================================================
# 5. D1 DATASET-LEVEL CONSENSUS
# ============================================================

# D1 contains two biological contrasts.
#
# Their NES values are preserved separately.
#
# D1 dataset-level NES = arithmetic mean
# of deficient and proficient NES.
#
# IMPORTANT:
# This is a descriptive integration summary,
# NOT a formal meta-analysis effect size.


master$D1_NES <- rowMeans(
  master[
    ,
    c(
      "D1_def_NES",
      "D1_prof_NES"
    )
  ],
  na.rm = FALSE
)


# Do the two D1 contrasts point in the same direction?

master$D1_same_direction <-
  sign(master$D1_def_NES) ==
  sign(master$D1_prof_NES)


# Were both D1 contrasts statistically significant?

master$D1_both_significant <-
  master$D1_def_padj < 0.05 &
  master$D1_prof_padj < 0.05


# At least one significant?

master$D1_any_significant <-
  master$D1_def_padj < 0.05 |
  master$D1_prof_padj < 0.05


# Strongest conservative D1 support:
# both significant AND same direction

master$D1_consistent_significant <-
  master$D1_same_direction &
  master$D1_both_significant


# ============================================================
# 6. ADD DIRECTION LABELS
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


master$D1_direction <- direction_label(
  master$D1_NES
)

master$D2_direction <- direction_label(
  master$D2_NES
)

master$D3_direction <- direction_label(
  master$D3_NES
)

master$D4_direction <- direction_label(
  master$D4_NES
)


# ============================================================
# 7. COUNT CROSS-DATASET DIRECTION AGREEMENT
# ============================================================

nes_matrix <- as.matrix(
  master[
    ,
    c(
      "D1_NES",
      "D2_NES",
      "D3_NES",
      "D4_NES"
    )
  ]
)


master$n_positive <- rowSums(
  nes_matrix > 0
)

master$n_negative <- rowSums(
  nes_matrix < 0
)


master$majority_direction <- ifelse(
  
  master$n_positive > master$n_negative,
  "UP",
  
  ifelse(
    master$n_negative > master$n_positive,
    "DOWN",
    "TIE"
  )
)


master$direction_agreement_n <- pmax(
  master$n_positive,
  master$n_negative
)


# ============================================================
# 8. SIGNIFICANCE FLAGS FOR D2-D4
# ============================================================

master$D2_significant <-
  master$D2_padj < 0.05

master$D3_significant <-
  master$D3_padj < 0.05

master$D4_significant <-
  master$D4_padj < 0.05


# ============================================================
# 9. ORDER COLUMNS
# ============================================================

master <- master[
  ,
  c(
    "pathway",
    
    "D1_def_NES",
    "D1_def_padj",
    
    "D1_prof_NES",
    "D1_prof_padj",
    
    "D1_NES",
    "D1_same_direction",
    "D1_both_significant",
    "D1_any_significant",
    "D1_consistent_significant",
    "D1_direction",
    
    "D2_NES",
    "D2_padj",
    "D2_significant",
    "D2_direction",
    
    "D3_NES",
    "D3_padj",
    "D3_significant",
    "D3_direction",
    
    "D4_NES",
    "D4_padj",
    "D4_significant",
    "D4_direction",
    
    "n_positive",
    "n_negative",
    "majority_direction",
    "direction_agreement_n"
  )
]


# ============================================================
# 10. SAVE MASTER TABLE
# ============================================================

write.csv(
  master,
  file.path(
    outdir,
    "D1_D4_Hallmark_master_matrix.csv"
  ),
  row.names = FALSE
)


# Also save simple NES-only matrix

nes_output <- master[
  ,
  c(
    "pathway",
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES"
  )
]


write.csv(
  nes_output,
  file.path(
    outdir,
    "D1_D4_Hallmark_NES_matrix.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. D1 INTERNAL CONSISTENCY SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D1 INTERNAL CONSISTENCY\n")
cat("========================================\n\n")

cat(
  "D1 same direction:",
  sum(master$D1_same_direction),
  "/ 50\n"
)

cat(
  "D1 both significant:",
  sum(master$D1_both_significant),
  "/ 50\n"
)

cat(
  "D1 same direction AND both significant:",
  sum(master$D1_consistent_significant),
  "/ 50\n"
)


# ============================================================
# 12. RAW CROSS-DATASET DIRECTION SUMMARY
# ============================================================

cat("\n========================================\n")
cat("D1-D4 RAW DIRECTION AGREEMENT\n")
cat("========================================\n\n")


direction_counts <- table(
  master$direction_agreement_n
)

print(
  direction_counts
)


cat(
  "\n4/4 same direction:",
  sum(master$direction_agreement_n == 4),
  "\n"
)

cat(
  "3/4 same direction:",
  sum(master$direction_agreement_n == 3),
  "\n"
)

cat(
  "2/4 split:",
  sum(master$direction_agreement_n == 2),
  "\n"
)


# ============================================================
# 13. SHOW PERFECTLY CONCORDANT HALLMARKS
# ============================================================

perfect <- master[
  master$direction_agreement_n == 4,
  c(
    "pathway",
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES",
    "majority_direction"
  )
]


perfect <- perfect[
  order(
    perfect$majority_direction,
    perfect$pathway
  ),
]


cat("\n========================================\n")
cat("4/4 DIRECTION-CONCORDANT HALLMARKS\n")
cat("========================================\n\n")

print(
  perfect,
  row.names = FALSE
)


write.csv(
  perfect,
  file.path(
    outdir,
    "D1_D4_4of4_direction_concordant.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. FINISHED
# ============================================================

cat("\n========================================\n")
cat("MASTER MATRIX COMPLETE\n")
cat("========================================\n")