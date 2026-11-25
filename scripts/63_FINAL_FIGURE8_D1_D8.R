# ============================================================
# 63_FINAL_FIGURE8_D1_D8.R
#
# FINAL SYNTHESIS FIGURE
# Cross-dataset transcriptomic programs of acquired
# PARP-inhibitor resistance
#
# PANELS
# A. Accession-level classification: D1-D6 vs D1-D8
# B. Frozen D1 hypothesis heatmap across D1-D8
# C. Model-level sensitivity: 9 units vs 12 units
# D. Treatment-state decomposition
#
# IMPORTANT
# - D1-D8 analyses are NOT rerun here.
# - This script only integrates locked Hallmark outputs.
# - D7 contributes ONE conservative accession vote in A/B.
# - D7A and D7B are separate units only in model sensitivity C.
# - D8 contributes ONE primary S2-vs-C3 vote.
# - D8 trajectory contributes zero additional votes.
# - D8 is excluded from ON/OFF grouping in panel D because
#   treatment-at-collection status is uncertain.
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


required_packages <- c(
  "ggplot2",
  "patchwork"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Missing package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\nInstall them before running Figure 8."
    )
  )
}


library(ggplot2)
library(patchwork)


out_dir <-
  "results/FINAL_D1_D8_integration/Figure8_FINAL"


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("FINAL FIGURE 8 — D1-D8 SYNTHESIS\n")
cat("============================================================\n\n")


# ============================================================
# 1. LOCATE LOCKED SOURCE FILES
# ============================================================

find_one_file <- function(filename) {
  
  hits <- list.files(
    "results",
    pattern = paste0(
      "^",
      filename,
      "$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  
  
  if (length(hits) == 0) {
    
    stop(
      paste0(
        "Could not find required file:\n",
        filename
      )
    )
  }
  
  
  if (length(hits) > 1) {
    
    cat(
      "\nMultiple copies found for:\n",
      filename,
      "\n"
    )
    
    print(hits)
    
    stop(
      "Ambiguous source file. Keep only one final copy."
    )
  }
  
  
  hits[1]
}


d1_d6_file <-
  "results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv"


if (!file.exists(d1_d6_file)) {
  
  stop(
    paste0(
      "Missing frozen D1-D6 master:\n",
      d1_d6_file
    )
  )
}


d7a_file <- find_one_file(
  "D7A_M468_Hallmark_GSEA.csv"
)


d7b_file <- find_one_file(
  "D7B_SUM1315_Hallmark_GSEA.csv"
)


d8_file <- find_one_file(
  "D8_primary_Hallmark_GSEA.csv"
)


d5_file <- find_one_file(
  "D5_H3_P3_Hallmark_master_comparison.csv"
)


cat("D1-D6:", d1_d6_file, "\n")
cat("D7A   :", d7a_file, "\n")
cat("D7B   :", d7b_file, "\n")
cat("D8    :", d8_file, "\n")
cat("D5    :", d5_file, "\n\n")


# ============================================================
# 2. LOAD TABLES
# ============================================================

old <- read.csv(
  d1_d6_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


d7a <- read.csv(
  d7a_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


d7b <- read.csv(
  d7b_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


d8 <- read.csv(
  d8_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


d5 <- read.csv(
  d5_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 3. BASIC SOURCE QC
# ============================================================

stopifnot(
  nrow(old) == 50,
  nrow(d7a) == 50,
  nrow(d7b) == 50,
  nrow(d8) == 50,
  nrow(d5) == 50
)


stopifnot(
  !anyDuplicated(old$pathway),
  !anyDuplicated(d7a$pathway),
  !anyDuplicated(d7b$pathway),
  !anyDuplicated(d8$ID),
  !anyDuplicated(d5$pathway)
)


# Reorder every table to the frozen D1-D6 Hallmark order.

d7a <- d7a[
  match(
    old$pathway,
    d7a$pathway
  ),
  ,
  drop = FALSE
]


d7b <- d7b[
  match(
    old$pathway,
    d7b$pathway
  ),
  ,
  drop = FALSE
]


d8 <- d8[
  match(
    old$pathway,
    d8$ID
  ),
  ,
  drop = FALSE
]


d5 <- d5[
  match(
    old$pathway,
    d5$pathway
  ),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    old$pathway,
    d7a$pathway
  ),
  identical(
    old$pathway,
    d7b$pathway
  ),
  identical(
    old$pathway,
    d8$ID
  ),
  identical(
    old$pathway,
    d5$pathway
  )
)


cat("Source-table QC PASSED.\n\n")


# ============================================================
# 4. CLASSIFICATION FUNCTION
# ============================================================

classify_programs <- function(
    n_up,
    n_down,
    n_sig_up,
    n_sig_down,
    direction_threshold,
    strong_sig_threshold,
    specific_label
) {
  
  result <- rep(
    "Inconsistent",
    length(n_up)
  )
  
  
  # ----------------------------------------------------------
  # 1. CONTEXT-DEPENDENT
  # significant evidence exists in BOTH directions
  # ----------------------------------------------------------
  
  context <-
    n_sig_up >= 1 &
    n_sig_down >= 1
  
  
  result[context] <-
    "Context-dependent"
  
  
  # ----------------------------------------------------------
  # 2. STRONG RECURRENT
  # ----------------------------------------------------------
  
  unresolved <-
    result == "Inconsistent"
  
  
  strong_up <-
    unresolved &
    n_up >= direction_threshold &
    n_sig_up >= strong_sig_threshold &
    n_sig_down == 0
  
  
  strong_down <-
    unresolved &
    n_down >= direction_threshold &
    n_sig_down >= strong_sig_threshold &
    n_sig_up == 0
  
  
  result[
    strong_up |
      strong_down
  ] <- "Strong recurrent"
  
  
  # ----------------------------------------------------------
  # 3. DIRECTIONAL
  # ----------------------------------------------------------
  
  unresolved <-
    result == "Inconsistent"
  
  
  directional_up <-
    unresolved &
    n_up >= direction_threshold &
    n_sig_down == 0
  
  
  directional_down <-
    unresolved &
    n_down >= direction_threshold &
    n_sig_up == 0
  
  
  result[
    directional_up |
      directional_down
  ] <- "Directional"
  
  
  # ----------------------------------------------------------
  # 4. DATASET / UNIT SPECIFIC
  # ----------------------------------------------------------
  
  unresolved <-
    result == "Inconsistent"
  
  
  specific <-
    unresolved &
    (
      n_sig_up +
        n_sig_down ==
        1
    )
  
  
  result[specific] <-
    specific_label
  
  
  result
}


# ============================================================
# 5. D7 CONSERVATIVE ACCESSION VOTE
# ============================================================

d7_same_direction <-
  sign(
    d7a$NES
  ) ==
  sign(
    d7b$NES
  )


D7_NES <- ifelse(
  d7_same_direction,
  (
    d7a$NES +
      d7b$NES
  ) / 2,
  NA_real_
)


D7_direction <- ifelse(
  !d7_same_direction,
  0,
  sign(D7_NES)
)


D7_significant <-
  d7_same_direction &
  d7a$p.adjust < 0.05 &
  d7b$p.adjust < 0.05


cat("============================================================\n")
cat("D7 CHECKPOINT\n")
cat("============================================================\n")

cat(
  "Same direction:",
  sum(d7_same_direction),
  "/ 50\n"
)

cat(
  "Concordant significant:",
  sum(D7_significant),
  "/ 50\n\n"
)


stopifnot(
  sum(d7_same_direction) == 33,
  sum(D7_significant) == 6
)


# ============================================================
# 6. D8 PRIMARY ACCESSION VOTE
# ============================================================

D8_NES <-
  d8$NES


D8_direction <-
  sign(
    D8_NES
  )


D8_significant <-
  d8$p.adjust < 0.05


cat("============================================================\n")
cat("D8 CHECKPOINT\n")
cat("============================================================\n")

cat(
  "Primary significant Hallmarks:",
  sum(D8_significant),
  "/ 50\n\n"
)


stopifnot(
  sum(D8_significant) == 15
)


# ============================================================
# 7. REPRODUCE FROZEN D1-D6 ACCESSION CLASSIFICATION
# ============================================================

normalize_old_class <- function(x) {
  
  result <- rep(
    NA_character_,
    length(x)
  )
  
  
  result[
    grepl(
      "strong",
      x,
      ignore.case = TRUE
    )
  ] <- "Strong recurrent"
  
  
  result[
    grepl(
      "direction",
      x,
      ignore.case = TRUE
    )
  ] <- "Directional"
  
  
  result[
    grepl(
      "context",
      x,
      ignore.case = TRUE
    )
  ] <- "Context-dependent"
  
  
  result[
    grepl(
      "dataset",
      x,
      ignore.case = TRUE
    )
  ] <- "Dataset-specific"
  
  
  result[
    grepl(
      "inconsistent|weak",
      x,
      ignore.case = TRUE
    )
  ] <- "Inconsistent"
  
  
  result
}


old_accession_class <-
  normalize_old_class(
    old$classification
  )


accession_class_levels <- c(
  "Strong recurrent",
  "Directional",
  "Context-dependent",
  "Dataset-specific",
  "Inconsistent"
)


old_accession_counts <- as.integer(
  table(
    factor(
      old_accession_class,
      levels =
        accession_class_levels
    )
  )
)


expected_old_accession_counts <- c(
  1,
  12,
  15,
  12,
  10
)


cat("============================================================\n")
cat("D1-D6 ACCESSION CHECKPOINT\n")
cat("============================================================\n")


print(
  data.frame(
    class =
      accession_class_levels,
    observed =
      old_accession_counts,
    expected =
      expected_old_accession_counts
  ),
  row.names = FALSE
)


stopifnot(
  all(
    old_accession_counts ==
      expected_old_accession_counts
  )
)


cat("\nD1-D6 accession checkpoint PASSED.\n\n")


# ============================================================
# 8. FINAL D1-D8 ACCESSION COUNTS
# ============================================================

final_n_up <-
  old$n_UP +
  as.integer(
    D7_direction > 0
  ) +
  as.integer(
    D8_direction > 0
  )


final_n_down <-
  old$n_DOWN +
  as.integer(
    D7_direction < 0
  ) +
  as.integer(
    D8_direction < 0
  )


final_n_sig_up <-
  old$n_sig_UP +
  as.integer(
    D7_direction > 0 &
      D7_significant
  ) +
  as.integer(
    D8_direction > 0 &
      D8_significant
  )


final_n_sig_down <-
  old$n_sig_DOWN +
  as.integer(
    D7_direction < 0 &
      D7_significant
  ) +
  as.integer(
    D8_direction < 0 &
      D8_significant
  )


final_accession_class <-
  classify_programs(
    
    n_up =
      final_n_up,
    
    n_down =
      final_n_down,
    
    n_sig_up =
      final_n_sig_up,
    
    n_sig_down =
      final_n_sig_down,
    
    direction_threshold =
      6,
    
    strong_sig_threshold =
      4,
    
    specific_label =
      "Dataset-specific"
  )


final_accession_counts <- as.integer(
  table(
    factor(
      final_accession_class,
      levels =
        accession_class_levels
    )
  )
)


expected_final_accession_counts <- c(
  0,
  4,
  24,
  12,
  10
)


cat("============================================================\n")
cat("FINAL D1-D8 ACCESSION CHECKPOINT\n")
cat("============================================================\n")


print(
  data.frame(
    class =
      accession_class_levels,
    observed =
      final_accession_counts,
    expected =
      expected_final_accession_counts
  ),
  row.names = FALSE
)


stopifnot(
  all(
    final_accession_counts ==
      expected_final_accession_counts
  )
)


cat("\nFINAL D1-D8 accession checkpoint PASSED.\n\n")


# ============================================================
# 9. PANEL A DATA
# ============================================================

panelA_data <- rbind(
  
  data.frame(
    stage =
      "D1-D6\n6 accessions",
    class =
      accession_class_levels,
    n =
      old_accession_counts
  ),
  
  data.frame(
    stage =
      "D1-D8\n8 accessions",
    class =
      accession_class_levels,
    n =
      final_accession_counts
  )
)


panelA_data$stage <- factor(
  panelA_data$stage,
  levels = c(
    "D1-D6\n6 accessions",
    "D1-D8\n8 accessions"
  )
)


panelA_data$class <- factor(
  panelA_data$class,
  levels =
    accession_class_levels
)


# ============================================================
# 10. PANEL B — FROZEN D1 HYPOTHESIS HEATMAP
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


short_names <- c(
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" =
    "Oxidative phosphorylation",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" =
    "EMT",
  
  "HALLMARK_G2M_CHECKPOINT" =
    "G2M checkpoint",
  
  "HALLMARK_E2F_TARGETS" =
    "E2F targets",
  
  "HALLMARK_MTORC1_SIGNALING" =
    "mTORC1 signaling",
  
  "HALLMARK_MYC_TARGETS_V1" =
    "MYC targets V1",
  
  "HALLMARK_GLYCOLYSIS" =
    "Glycolysis",
  
  "HALLMARK_INFLAMMATORY_RESPONSE" =
    "Inflammatory response",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" =
    "TNFα / NF-κB",
  
  "HALLMARK_UV_RESPONSE_DN" =
    "UV response DN"
)


heat <- rbind(
  
  data.frame(
    pathway = old$pathway,
    dataset = "D1",
    NES = old$D1_NES,
    significant =
      old$D1_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D2",
    NES = old$D2_NES,
    significant =
      old$D2_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D3",
    NES = old$D3_NES,
    significant =
      old$D3_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D4",
    NES = old$D4_NES,
    significant =
      old$D4_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D5",
    NES = old$D5_NES,
    significant =
      old$D5_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D6",
    NES = old$D6_NES,
    significant =
      old$D6_sig_direction != 0,
    discordant = FALSE
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D7",
    NES = D7_NES,
    significant =
      D7_significant,
    discordant =
      !d7_same_direction
  ),
  
  data.frame(
    pathway = old$pathway,
    dataset = "D8",
    NES = D8_NES,
    significant =
      D8_significant,
    discordant = FALSE
  )
)


heat <- heat[
  heat$pathway %in%
    frozen_pathways,
  ,
  drop = FALSE
]


heat$pathway_label <-
  short_names[
    heat$pathway
  ]


heat$pathway_label <- factor(
  heat$pathway_label,
  levels = rev(
    unname(
      short_names[
        frozen_pathways
      ]
    )
  )
)


heat$dataset <- factor(
  heat$dataset,
  levels =
    paste0(
      "D",
      1:8
    )
)


heat$significance_symbol <- ifelse(
  heat$significant,
  "*",
  ""
)


# ============================================================
# 11. MODEL-LEVEL SENSITIVITY
#
# Frozen original model units:
# D1-def
# D1-prof
# D2
# D3
# D4
# HCC1806-H1
# HCC1806-H2
# PSN1-P1
# D6
#
# Expanded:
# + D7 M468
# + D7 SUM1315
# + D8 PH039
#
# IMPORTANT:
# P1 significance must survive P1 sensitivity analysis.
# This reproduces the locked old 9-unit checkpoint.
# ============================================================

P1_robust_significant <-
  old$P1_padj < 0.05 &
  old$P1_sens_padj < 0.05 &
  sign(
    old$P1_NES
  ) ==
  sign(
    old$P1_sens_NES
  )


unit_NES_9 <- cbind(
  
  D1_def =
    old$D1_def_NES,
  
  D1_prof =
    old$D1_prof_NES,
  
  D2 =
    old$D2_NES,
  
  D3 =
    old$D3_NES,
  
  D4 =
    old$D4_NES,
  
  HCC_H1 =
    old$H1_NES,
  
  HCC_H2 =
    old$H2_NES,
  
  PSN_P1 =
    old$P1_NES,
  
  D6 =
    old$D6_NES
)


unit_SIG_9 <- cbind(
  
  D1_def =
    old$D1_def_padj < 0.05,
  
  D1_prof =
    old$D1_prof_padj < 0.05,
  
  D2 =
    old$D2_padj < 0.05,
  
  D3 =
    old$D3_padj < 0.05,
  
  D4 =
    old$D4_padj < 0.05,
  
  HCC_H1 =
    old$H1_padj < 0.05,
  
  HCC_H2 =
    old$H2_padj < 0.05,
  
  PSN_P1 =
    P1_robust_significant,
  
  D6 =
    old$D6_padj < 0.05
)


model9_n_up <-
  rowSums(
    unit_NES_9 > 0
  )


model9_n_down <-
  rowSums(
    unit_NES_9 < 0
  )


model9_n_sig_up <-
  rowSums(
    unit_NES_9 > 0 &
      unit_SIG_9
  )


model9_n_sig_down <-
  rowSums(
    unit_NES_9 < 0 &
      unit_SIG_9
  )


model9_class <- classify_programs(
  
  n_up =
    model9_n_up,
  
  n_down =
    model9_n_down,
  
  n_sig_up =
    model9_n_sig_up,
  
  n_sig_down =
    model9_n_sig_down,
  
  direction_threshold =
    6,
  
  strong_sig_threshold =
    4,
  
  specific_label =
    "Unit-specific"
)


model_class_levels <- c(
  "Strong recurrent",
  "Directional",
  "Context-dependent",
  "Unit-specific",
  "Inconsistent"
)


model9_counts <- as.integer(
  table(
    factor(
      model9_class,
      levels =
        model_class_levels
    )
  )
)


expected_model9_counts <- c(
  4,
  14,
  24,
  5,
  3
)


cat("============================================================\n")
cat("OLD 9-UNIT MODEL CHECKPOINT\n")
cat("============================================================\n")


print(
  data.frame(
    class =
      model_class_levels,
    observed =
      model9_counts,
    expected =
      expected_model9_counts
  ),
  row.names = FALSE
)


stopifnot(
  all(
    model9_counts ==
      expected_model9_counts
  )
)


cat("\n9-unit model checkpoint PASSED.\n\n")


# ============================================================
# 12. EXPANDED 12-UNIT MODEL ANALYSIS
# ============================================================

unit_NES_12 <- cbind(
  
  unit_NES_9,
  
  D7_M468 =
    d7a$NES,
  
  D7_SUM1315 =
    d7b$NES,
  
  D8_PH039 =
    d8$NES
)


unit_SIG_12 <- cbind(
  
  unit_SIG_9,
  
  D7_M468 =
    d7a$p.adjust < 0.05,
  
  D7_SUM1315 =
    d7b$p.adjust < 0.05,
  
  D8_PH039 =
    d8$p.adjust < 0.05
)


model12_n_up <-
  rowSums(
    unit_NES_12 > 0
  )


model12_n_down <-
  rowSums(
    unit_NES_12 < 0
  )


model12_n_sig_up <-
  rowSums(
    unit_NES_12 > 0 &
      unit_SIG_12
  )


model12_n_sig_down <-
  rowSums(
    unit_NES_12 < 0 &
      unit_SIG_12
  )


model12_class <- classify_programs(
  
  n_up =
    model12_n_up,
  
  n_down =
    model12_n_down,
  
  n_sig_up =
    model12_n_sig_up,
  
  n_sig_down =
    model12_n_sig_down,
  
  direction_threshold =
    8,
  
  strong_sig_threshold =
    6,
  
  specific_label =
    "Unit-specific"
)


model12_counts <- as.integer(
  table(
    factor(
      model12_class,
      levels =
        model_class_levels
    )
  )
)


expected_model12_counts <- c(
  0,
  5,
  35,
  5,
  5
)


cat("============================================================\n")
cat("FINAL 12-UNIT MODEL CHECKPOINT\n")
cat("============================================================\n")


print(
  data.frame(
    class =
      model_class_levels,
    observed =
      model12_counts,
    expected =
      expected_model12_counts
  ),
  row.names = FALSE
)


stopifnot(
  all(
    model12_counts ==
      expected_model12_counts
  )
)


cat("\n12-unit model checkpoint PASSED.\n\n")


# ============================================================
# 13. PANEL C DATA
# ============================================================

panelC_data <- rbind(
  
  data.frame(
    stage =
      "D1-D6\n9 model units",
    class =
      model_class_levels,
    n =
      model9_counts
  ),
  
  data.frame(
    stage =
      "D1-D8\n12 model units",
    class =
      model_class_levels,
    n =
      model12_counts
  )
)


panelC_data$stage <- factor(
  panelC_data$stage,
  levels = c(
    "D1-D6\n9 model units",
    "D1-D8\n12 model units"
  )
)


panelC_data$class <- factor(
  panelC_data$class,
  levels =
    model_class_levels
)


# ============================================================
# 14. FINAL TREATMENT-STATE DECOMPOSITION
#
# OFF / DRUG-FREE:
#   D1, D4, D6
#
# ON:
#   D2, D3, D7
#
# D7 conservative handling:
#   if M468 and SUM1315 agree in direction:
#       include their mean as D7 accession NES
#   if discordant:
#       omit D7 for that pathway from ON mean
#
# D8:
#   excluded because treatment state at collection is uncertain
#
# D5:
#   H3 and P3 are DIRECT within-model treatment-withdrawal
#   contrasts.
# ============================================================

OFF_mean_NES <- rowMeans(
  cbind(
    old$D1_NES,
    old$D4_NES,
    old$D6_NES
  )
)


ON_mean_NES <- rowMeans(
  
  cbind(
    old$D2_NES,
    old$D3_NES,
    D7_NES
  ),
  
  na.rm = TRUE
)


observational_exposure_delta <-
  ON_mean_NES -
  OFF_mean_NES


D5_exposure_mean <- rowMeans(
  cbind(
    d5$H3_NES,
    d5$P3_NES
  )
)


exposure_rho <- cor(
  observational_exposure_delta,
  D5_exposure_mean,
  method = "spearman"
)


same_exposure_direction <-
  sign(
    observational_exposure_delta
  ) ==
  sign(
    D5_exposure_mean
  )


same_exposure_n <-
  sum(
    same_exposure_direction
  )


cat("============================================================\n")
cat("FINAL TREATMENT-STATE CHECKPOINT\n")
cat("============================================================\n")


cat(
  "Spearman rho:",
  round(
    exposure_rho,
    4
  ),
  "\n"
)


cat(
  "Same direction:",
  same_exposure_n,
  "/ 50\n\n"
)


stopifnot(
  abs(
    exposure_rho -
      (-0.2557023)
  ) < 0.005,
  same_exposure_n == 24
)


cat("Treatment-state checkpoint PASSED.\n\n")


treatment_data <- data.frame(
  
  pathway =
    old$pathway,
  
  OFF_mean_NES =
    OFF_mean_NES,
  
  ON_mean_NES =
    ON_mean_NES,
  
  observational_delta =
    observational_exposure_delta,
  
  D5_exposure_mean =
    D5_exposure_mean,
  
  robust_D5 =
    d5$robust_shared_exposure_candidate,
  
  same_exposure_direction =
    same_exposure_direction,
  
  stringsAsFactors = FALSE
)


# ============================================================
# 15. SELECT LABELS FOR PANEL D
# ============================================================

treatment_labels <- c(
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" =
    "TNFα/NF-κB",
  
  "HALLMARK_INFLAMMATORY_RESPONSE" =
    "Inflammatory",
  
  "HALLMARK_COMPLEMENT" =
    "Complement",
  
  "HALLMARK_APOPTOSIS" =
    "Apoptosis",
  
  "HALLMARK_MYC_TARGETS_V1" =
    "MYC V1",
  
  "HALLMARK_MTORC1_SIGNALING" =
    "mTORC1",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" =
    "OXPHOS",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" =
    "EMT"
)


treatment_data$label <- ""


label_match <- match(
  treatment_data$pathway,
  names(
    treatment_labels
  )
)


has_label <-
  !is.na(
    label_match
  )


treatment_data$label[
  has_label
] <- treatment_labels[
  label_match[
    has_label
  ]
]


# ============================================================
# 16. COLOUR PALETTE
# ============================================================

class_colours <- c(
  
  "Strong recurrent" =
    "#54278F",
  
  "Directional" =
    "#2B8CBE",
  
  "Context-dependent" =
    "#D7301F",
  
  "Dataset-specific" =
    "#FDAE6B",
  
  "Unit-specific" =
    "#FDAE6B",
  
  "Inconsistent" =
    "#BDBDBD"
)


# ============================================================
# 17. PANEL A — ACCESSION-LEVEL CLASSIFICATION
# ============================================================

pA <- ggplot(
  panelA_data,
  aes(
    x = stage,
    y = n,
    fill = class
  )
) +
  
  geom_col(
    width = 0.68,
    colour = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label =
        ifelse(
          n == 0,
          "",
          n
        )
    ),
    position =
      position_stack(
        vjust = 0.5
      ),
    colour = "white",
    fontface = "bold",
    size = 3.6
  ) +
  
  scale_fill_manual(
    values =
      class_colours[
        accession_class_levels
      ],
    drop = FALSE
  ) +
  
  scale_y_continuous(
    breaks =
      seq(
        0,
        50,
        10
      ),
    limits =
      c(
        0,
        50
      ),
    expand =
      c(
        0,
        0
      )
  ) +
  
  labs(
    title =
      "A  Accession-level reproducibility",
    subtitle =
      "Expansion eliminates the apparent strongly recurrent program",
    x = NULL,
    y =
      "Hallmark pathways",
    fill =
      "Classification"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    plot.subtitle =
      element_text(
        size = 9.5
      ),
    axis.text.x =
      element_text(
        face = "bold"
      ),
    legend.position =
      "bottom",
    legend.title =
      element_text(
        face = "bold"
      )
  )


# ============================================================
# 18. PANEL B — FROZEN-HYPOTHESIS HEATMAP
# ============================================================

pB <- ggplot(
  heat,
  aes(
    x = dataset,
    y = pathway_label,
    fill = NES
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.55,
    width = 0.96,
    height = 0.96
  ) +
  
  geom_text(
    aes(
      label =
        significance_symbol
    ),
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits =
      c(
        -3.5,
        3.5
      ),
    oob =
      scales::squish,
    na.value =
      "grey75"
  ) +
  
  labs(
    title =
      "B  Prospective validation of the frozen D1 hypothesis",
    subtitle =
      "* FDR < 0.05; grey D7 cells = discordant M468/SUM1315 direction",
    x =
      "Independent GEO accession",
    y = NULL,
    fill =
      "NES"
  ) +
  
  theme_minimal(
    base_size = 10.5
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    plot.subtitle =
      element_text(
        size = 9
      ),
    panel.grid =
      element_blank(),
    axis.text.x =
      element_text(
        face = "bold"
      ),
    axis.text.y =
      element_text(
        colour = "black"
      ),
    legend.position =
      "right"
  )


# ============================================================
# 19. PANEL C — MODEL-LEVEL SENSITIVITY
# ============================================================

pC <- ggplot(
  panelC_data,
  aes(
    x = stage,
    y = n,
    fill = class
  )
) +
  
  geom_col(
    width = 0.68,
    colour = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label =
        ifelse(
          n == 0,
          "",
          n
        )
    ),
    position =
      position_stack(
        vjust = 0.5
      ),
    colour = "white",
    fontface = "bold",
    size = 3.6
  ) +
  
  scale_fill_manual(
    values =
      class_colours[
        model_class_levels
      ],
    drop = FALSE
  ) +
  
  scale_y_continuous(
    breaks =
      seq(
        0,
        50,
        10
      ),
    limits =
      c(
        0,
        50
      ),
    expand =
      c(
        0,
        0
      )
  ) +
  
  labs(
    title =
      "C  Model-level sensitivity analysis",
    subtitle =
      "Secondary analysis; model units are not independent studies",
    x = NULL,
    y =
      "Hallmark pathways",
    fill =
      "Classification"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    plot.subtitle =
      element_text(
        size = 9.5
      ),
    axis.text.x =
      element_text(
        face = "bold"
      ),
    legend.position =
      "bottom",
    legend.title =
      element_text(
        face = "bold"
      )
  )


# ============================================================
# 20. PANEL D — TREATMENT-STATE DECOMPOSITION
# ============================================================

pD <- ggplot(
  treatment_data,
  aes(
    x =
      observational_delta,
    y =
      D5_exposure_mean
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      colour =
        robust_D5
    ),
    size = 2.8,
    alpha = 0.85
  ) +
  
  scale_colour_manual(
    values = c(
      "FALSE" =
        "grey60",
      "TRUE" =
        "#D73027"
    ),
    labels = c(
      "FALSE" =
        "Other Hallmarks",
      "TRUE" =
        "Robust D5 exposure response"
    ),
    name = NULL
  ) +
  
  labs(
    title =
      "D  Treatment-state decomposition",
    subtitle =
      paste0(
        "Observed ON−OFF vs direct D5 withdrawal response; ",
        "Spearman ρ = ",
        sprintf(
          "%.3f",
          exposure_rho
        ),
        "; ",
        same_exposure_n,
        "/50 direction agreement"
      ),
    x =
      "Observational ON − OFF mean NES",
    y =
      "D5 direct exposure-response mean NES"
  ) +
  
  theme_classic(
    base_size = 10.5
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    plot.subtitle =
      element_text(
        size = 9
      ),
    legend.position =
      "bottom"
  )


# Add pathway labels.
# Use ggrepel when available; otherwise standard text.

label_data <- treatment_data[
  treatment_data$label != "",
  ,
  drop = FALSE
]


if (
  requireNamespace(
    "ggrepel",
    quietly = TRUE
  )
) {
  
  pD <- pD +
    
    ggrepel::geom_text_repel(
      data =
        label_data,
      aes(
        label =
          label
      ),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    )
  
} else {
  
  pD <- pD +
    
    geom_text(
      data =
        label_data,
      aes(
        label =
          label
      ),
      size = 2.8,
      hjust = -0.08,
      vjust = -0.25,
      check_overlap = TRUE
    )
}


# ============================================================
# 21. FINAL FIGURE 8 ASSEMBLY
# ============================================================

figure8 <-
  
  (
    pA |
      pB
  ) /
  
  (
    pC |
      pD
  ) +
  
  plot_layout(
    widths =
      c(
        0.85,
        1.55
      ),
    heights =
      c(
        1,
        1
      ),
    guides =
      "collect"
  ) +
  
  plot_annotation(
    
    title =
      "Figure 8. Cross-dataset synthesis of acquired PARP-inhibitor resistance",
    
    subtitle =
      paste0(
        "Eight independent GEO accessions reveal predominantly context-dependent ",
        "transcriptional remodeling rather than a universal Hallmark resistance state."
      ),
    
    caption =
      paste0(
        "Primary inference is accession-level. ",
        "D7 contributes one conservative accession vote; ",
        "D7 models are separated only in the model-level sensitivity analysis. ",
        "D8 is excluded from ON/OFF treatment grouping because collection exposure is uncertain."
      ),
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 17
          ),
        
        plot.subtitle =
          element_text(
            size = 11
          ),
        
        plot.caption =
          element_text(
            size = 8.5,
            colour = "grey30",
            hjust = 0
          )
      )
  )


# ============================================================
# 22. SAVE FIGURE
# ============================================================

png_file <- file.path(
  out_dir,
  "Figure8_FINAL_D1_D8.png"
)


pdf_file <- file.path(
  out_dir,
  "Figure8_FINAL_D1_D8.pdf"
)


ggsave(
  filename =
    png_file,
  plot =
    figure8,
  width = 16,
  height = 11,
  units = "in",
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)


ggsave(
  filename =
    pdf_file,
  plot =
    figure8,
  width = 16,
  height = 11,
  units = "in",
  bg = "white",
  limitsize = FALSE
)


# ============================================================
# 23. SAVE FIGURE SOURCE TABLES
# ============================================================

accession_master <- data.frame(
  
  pathway =
    old$pathway,
  
  old_class =
    old_accession_class,
  
  final_class =
    final_accession_class,
  
  final_n_UP =
    final_n_up,
  
  final_n_DOWN =
    final_n_down,
  
  final_n_sig_UP =
    final_n_sig_up,
  
  final_n_sig_DOWN =
    final_n_sig_down,
  
  D7_NES =
    D7_NES,
  
  D7_significant =
    D7_significant,
  
  D8_NES =
    D8_NES,
  
  D8_significant =
    D8_significant,
  
  stringsAsFactors = FALSE
)


write.csv(
  accession_master,
  file.path(
    out_dir,
    "Figure8_accession_source_table.csv"
  ),
  row.names = FALSE
)


model_master <- data.frame(
  
  pathway =
    old$pathway,
  
  model9_class =
    model9_class,
  
  model12_class =
    model12_class,
  
  model12_n_UP =
    model12_n_up,
  
  model12_n_DOWN =
    model12_n_down,
  
  model12_n_sig_UP =
    model12_n_sig_up,
  
  model12_n_sig_DOWN =
    model12_n_sig_down,
  
  stringsAsFactors = FALSE
)


write.csv(
  model_master,
  file.path(
    out_dir,
    "Figure8_model_level_source_table.csv"
  ),
  row.names = FALSE
)


write.csv(
  treatment_data,
  file.path(
    out_dir,
    "Figure8_treatment_state_source_table.csv"
  ),
  row.names = FALSE
)


write.csv(
  heat,
  file.path(
    out_dir,
    "Figure8_frozen_hypothesis_heatmap_source.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. SAVE FINAL NUMERICAL CHECKPOINT
# ============================================================

checkpoint <- data.frame(
  
  metric = c(
    
    "Accession_strong",
    "Accession_directional",
    "Accession_context",
    "Accession_dataset_specific",
    "Accession_inconsistent",
    
    "Model_strong",
    "Model_directional",
    "Model_context",
    "Model_unit_specific",
    "Model_inconsistent",
    
    "Treatment_Spearman_rho",
    "Treatment_same_direction_n"
  ),
  
  value = c(
    
    final_accession_counts[1],
    final_accession_counts[2],
    final_accession_counts[3],
    final_accession_counts[4],
    final_accession_counts[5],
    
    model12_counts[1],
    model12_counts[2],
    model12_counts[3],
    model12_counts[4],
    model12_counts[5],
    
    exposure_rho,
    same_exposure_n
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  checkpoint,
  file.path(
    out_dir,
    "Figure8_FINAL_checkpoint.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 25. SESSION INFO
# ============================================================

sink(
  file.path(
    out_dir,
    "Figure8_sessionInfo.txt"
  )
)


sessionInfo()


sink()


# ============================================================
# 26. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 8 COMPLETE\n")
cat("============================================================\n\n")


cat("FINAL ACCESSION LEVEL\n")
cat("---------------------\n")
cat("Strong recurrent :", final_accession_counts[1], "\n")
cat("Directional      :", final_accession_counts[2], "\n")
cat("Context-dependent:", final_accession_counts[3], "\n")
cat("Dataset-specific :", final_accession_counts[4], "\n")
cat("Inconsistent     :", final_accession_counts[5], "\n\n")


cat("FINAL MODEL LEVEL\n")
cat("-----------------\n")
cat("Strong recurrent :", model12_counts[1], "\n")
cat("Directional      :", model12_counts[2], "\n")
cat("Context-dependent:", model12_counts[3], "\n")
cat("Unit-specific    :", model12_counts[4], "\n")
cat("Inconsistent     :", model12_counts[5], "\n\n")


cat("TREATMENT STATE\n")
cat("---------------\n")
cat(
  "Spearman rho:",
  round(
    exposure_rho,
    4
  ),
  "\n"
)
cat(
  "Same direction:",
  same_exposure_n,
  "/50\n\n"
)


cat("PNG:\n", png_file, "\n\n")
cat("PDF:\n", pdf_file, "\n\n")


cat(
  "All locked checkpoints reproduced successfully.\n"
)


cat("============================================================\n")