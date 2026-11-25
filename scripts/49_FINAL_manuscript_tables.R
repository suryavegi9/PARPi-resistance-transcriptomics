# ============================================================
# FINAL MANUSCRIPT TABLES
# PARP-INHIBITOR RESISTANCE TRANSCRIPTOMICS PROJECT
#
# Generates:
#
# Table 1 - Dataset/design overview
# Table 2 - Final accession-level Hallmark classification
# Table 3 - Robust context-dependent pathways
# Table 4 - Treatment-state evidence
#
# Supplementary:
# S1 - Full accession-level integration master
# S2 - Full model-level sensitivity master
# S3 - Full treatment-state master
#
# No new analysis is performed here.
# This script only packages LOCKED results.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. OUTPUT DIRECTORY
# ============================================================

outdir <- "manuscript/tables"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. INPUT FILES
# ============================================================

accession_master_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "FINAL_D1_D6_Hallmark_master.csv"
)


model_master_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "model_level_sensitivity/",
  "FINAL_model_level_Hallmark_master.csv"
)


robust_context_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "model_level_sensitivity/",
  "FINAL_robust_context_dependent_both_methods.csv"
)


treatment_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "treatment_state/",
  "FINAL_treatment_state_master.csv"
)


required_files <- c(
  accession_master_file,
  model_master_file,
  robust_context_file,
  treatment_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "One or more required integration outputs are missing."
  )
}


# ============================================================
# 2. HELPER
# ============================================================

clean_hallmark <- function(x) {
  
  x <- gsub(
    "^HALLMARK_",
    "",
    x
  )
  
  x <- gsub(
    "_",
    " ",
    x
  )
  
  x
}


# ============================================================
# ============================================================
# TABLE 1
# DATASET AND DESIGN OVERVIEW
# ============================================================
# ============================================================

table1 <- data.frame(
  
  Dataset = c(
    "D1A",
    "D1B",
    "D2",
    "D3",
    "D4",
    "D5A",
    "D5B",
    "D6"
  ),
  
  Accession = c(
    "GSE235980",
    "GSE235980",
    "GSE117765",
    "GSE153867",
    "GSE342642",
    "GSE202434",
    "GSE202434",
    "GSE191231"
  ),
  
  Model = c(
    "UWB1.289 BRCA1-deficient",
    "UWB1.289 BRCA1-proficient/restored",
    "PEO1 BRCA2-mutant",
    "A2780",
    "BRCA1-mutant HGSOC PDX",
    "HCC1806",
    "PSN1",
    "TP53-mutant, BRCA1/2-WT HGSOC PDX"
  ),
  
  System = c(
    "Ovarian cancer cell line",
    "Ovarian cancer cell line",
    "Ovarian cancer cell line",
    "Ovarian cancer cell line",
    "Patient-derived xenograft",
    "Breast cancer cell line",
    "Pancreatic cancer cell line",
    "Patient-derived xenograft"
  ),
  
  Biological_design = c(
    
    "Sensitive vs acquired-resistant; n=2 per group",
    
    "Sensitive vs acquired-resistant; n=2 per group",
    
    paste0(
      "2 parental samples vs 4 independently derived ",
      "olaparib-resistant clones"
    ),
    
    "8 parental vs 8 olaparib-resistant samples",
    
    paste0(
      "2 control tumors vs 2 resistant tumors; ",
      "3 spatial regions per tumor"
    ),
    
    paste0(
      "Parental, maintained talazoparib-resistant, ",
      "and drug-free persistent-resistant states; ",
      "3 replicates per state"
    ),
    
    paste0(
      "Parental, maintained talazoparib-resistant, ",
      "and withdrawn/re-sensitized states; ",
      "3 replicates per state"
    ),
    
    paste0(
      "2 vehicle-control mice vs 2 olaparib-selected ",
      "recurrent/resistant mice"
    )
  ),
  
  RNA_collection_state = c(
    
    "Drug-free",
    
    "Drug-free",
    
    "Ongoing olaparib exposure",
    
    "Ongoing olaparib exposure",
    
    "Untreated acquired-resistant tumor",
    
    paste0(
      "Maintained-on-drug and drug-free ",
      "persistent-resistant states"
    ),
    
    paste0(
      "Maintained-on-drug and withdrawn ",
      "re-sensitized states"
    ),
    
    paste0(
      "Off treatment; recurrent tumors collected ",
      "long after olaparib withdrawal"
    )
  ),
  
  Expression_input = c(
    "Raw counts",
    "Raw counts",
    "TPM",
    "FPKM",
    "Raw supplied counts",
    "log2(TPM+1)",
    "log2(expression+1)",
    "TPM"
  ),
  
  Primary_analysis = c(
    
    "DESeq2",
    
    "DESeq2",
    
    "limma-trend",
    
    "limma-trend",
    
    "voom + duplicateCorrelation",
    
    "limma-trend",
    
    "limma-trend",
    
    "limma-trend"
  ),
  
  Study_role = c(
    
    "Discovery",
    
    "Discovery",
    
    "External validation",
    
    "External validation",
    
    "In-vivo validation",
    
    "Treatment-state / persistent-resistance dissection",
    
    "Treatment-state / resistance-reversion dissection",
    
    "Independent off-treatment in-vivo validation"
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  table1,
  file.path(
    outdir,
    "Table1_dataset_design_overview.csv"
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================
# TABLE 2
# FINAL ACCESSION-LEVEL HALLMARK CLASSIFICATION
# ============================================================
# ============================================================

acc <- read.csv(
  accession_master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


table2 <- data.frame(
  
  Pathway =
    clean_hallmark(
      acc$pathway
    ),
  
  D1_NES =
    round(
      acc$D1_NES,
      3
    ),
  
  D2_NES =
    round(
      acc$D2_NES,
      3
    ),
  
  D3_NES =
    round(
      acc$D3_NES,
      3
    ),
  
  D4_NES =
    round(
      acc$D4_NES,
      3
    ),
  
  D5_NES =
    round(
      acc$D5_NES,
      3
    ),
  
  D6_NES =
    round(
      acc$D6_NES,
      3
    ),
  
  N_direction_UP =
    acc$n_UP,
  
  N_direction_DOWN =
    acc$n_DOWN,
  
  N_significant_UP =
    acc$n_sig_UP,
  
  N_significant_DOWN =
    acc$n_sig_DOWN,
  
  Majority_direction =
    acc$majority_direction,
  
  Final_classification =
    acc$classification,
  
  stringsAsFactors = FALSE
)


class_order <- c(
  
  "Strong recurrent",
  
  "Directionally recurrent",
  
  "Context-dependent",
  
  "Dataset-specific significant",
  
  "Inconsistent/weak"
)


table2$class_rank <- match(
  table2$Final_classification,
  class_order
)


table2 <- table2[
  order(
    table2$class_rank,
    -pmax(
      table2$N_significant_UP,
      table2$N_significant_DOWN
    ),
    table2$Pathway
  ),
]


table2$class_rank <- NULL


write.csv(
  table2,
  file.path(
    outdir,
    "Table2_final_accession_Hallmark_classification.csv"
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================
# TABLE 3
# ROBUST CONTEXT-DEPENDENT PATHWAYS
# ============================================================
# ============================================================

robust_context <- read.csv(
  robust_context_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


model_master <- read.csv(
  model_master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


table3 <- merge(
  
  robust_context[
    ,
    c(
      "pathway",
      "accession_level_classification",
      "model_level_classification"
    )
  ],
  
  model_master[
    ,
    c(
      "pathway",
      "n_UP",
      "n_DOWN",
      "n_sig_UP",
      "n_sig_DOWN"
    )
  ],
  
  by = "pathway",
  
  all.x = TRUE
)


table3 <- data.frame(
  
  Pathway =
    clean_hallmark(
      table3$pathway
    ),
  
  Model_level_UP =
    table3$n_UP,
  
  Model_level_DOWN =
    table3$n_DOWN,
  
  Significant_UP =
    table3$n_sig_UP,
  
  Significant_DOWN =
    table3$n_sig_DOWN,
  
  Accession_classification =
    table3$accession_level_classification,
  
  Model_level_classification =
    table3$model_level_classification,
  
  stringsAsFactors = FALSE
)


table3 <- table3[
  order(
    -(
      table3$Significant_UP +
        table3$Significant_DOWN
    ),
    table3$Pathway
  ),
]


write.csv(
  table3,
  file.path(
    outdir,
    "Table3_robust_context_dependent_pathways.csv"
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================
# TABLE 4
# TREATMENT-STATE EVIDENCE
# ============================================================
# ============================================================

treatment <- read.csv(
  treatment_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Main table:
# retain pathways with some classified exposure evidence.
# Weak/mixed pathways remain in Supplementary Table S3.

table4_raw <- treatment[
  treatment$exposure_evidence !=
    "Weak/mixed",
]


table4 <- data.frame(
  
  Pathway =
    clean_hallmark(
      table4_raw$pathway
    ),
  
  Observational_ON_minus_OFF =
    round(
      table4_raw$
        observational_exposure_delta,
      3
    ),
  
  D5_controlled_exposure_effect =
    round(
      table4_raw$
        D5_exposure_mean,
      3
    ),
  
  D5_robust_exposure_program =
    table4_raw$
    robust_exposure_program,
  
  Direction_agreement =
    table4_raw$
    exposure_direction_agreement,
  
  Clean_between_group_separation =
    table4_raw$
    clean_between_group_separation,
  
  Exposure_evidence =
    table4_raw$
    exposure_evidence,
  
  stringsAsFactors = FALSE
)


evidence_order <- c(
  
  "Strong convergent",
  
  "Convergent",
  
  "D5-specific exposure",
  
  "Observational-only"
)


table4$evidence_rank <- match(
  table4$Exposure_evidence,
  evidence_order
)


table4 <- table4[
  order(
    table4$evidence_rank,
    -abs(
      table4$
        D5_controlled_exposure_effect
    ),
    table4$Pathway
  ),
]


table4$evidence_rank <- NULL


write.csv(
  table4,
  file.path(
    outdir,
    "Table4_treatment_state_evidence.csv"
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================
# SUPPLEMENTARY TABLES
# ============================================================
# ============================================================


# ------------------------------------------------------------
# Supplementary Table S1
# Full accession-level master
# ------------------------------------------------------------

write.csv(
  acc,
  file.path(
    outdir,
    "TableS1_full_accession_level_master.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Supplementary Table S2
# Full model-level master
# ------------------------------------------------------------

write.csv(
  model_master,
  file.path(
    outdir,
    "TableS2_full_model_level_sensitivity_master.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Supplementary Table S3
# Full treatment-state master
# ------------------------------------------------------------

write.csv(
  treatment,
  file.path(
    outdir,
    "TableS3_full_treatment_state_master.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 3. OPTIONAL EXCEL WORKBOOK
# ============================================================

if (
  requireNamespace(
    "openxlsx",
    quietly = TRUE
  )
) {
  
  wb <- openxlsx::createWorkbook()
  
  
  openxlsx::addWorksheet(
    wb,
    "Table1_Design"
  )
  
  openxlsx::writeData(
    wb,
    "Table1_Design",
    table1
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Table2_Hallmark"
  )
  
  openxlsx::writeData(
    wb,
    "Table2_Hallmark",
    table2
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Table3_Context"
  )
  
  openxlsx::writeData(
    wb,
    "Table3_Context",
    table3
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "Table4_Treatment"
  )
  
  openxlsx::writeData(
    wb,
    "Table4_Treatment",
    table4
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "TableS1_Accessions"
  )
  
  openxlsx::writeData(
    wb,
    "TableS1_Accessions",
    acc
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "TableS2_Models"
  )
  
  openxlsx::writeData(
    wb,
    "TableS2_Models",
    model_master
  )
  
  
  openxlsx::addWorksheet(
    wb,
    "TableS3_Treatment"
  )
  
  openxlsx::writeData(
    wb,
    "TableS3_Treatment",
    treatment
  )
  
  
  openxlsx::saveWorkbook(
    
    wb,
    
    file.path(
      outdir,
      "FINAL_manuscript_tables.xlsx"
    ),
    
    overwrite = TRUE
  )
  
  
  cat(
    "\nExcel workbook generated.\n"
  )
  
} else {
  
  cat(
    paste0(
      "\nopenxlsx not installed. ",
      "CSV tables were generated normally.\n"
    )
  )
}


# ============================================================
# 4. QC SUMMARY
# ============================================================

cat("\n========================================\n")
cat("FINAL MANUSCRIPT TABLES COMPLETE\n")
cat("========================================\n\n")


cat(
  "Table 1 rows:",
  nrow(table1),
  "\n"
)


cat(
  "Table 2 Hallmarks:",
  nrow(table2),
  "\n"
)


cat(
  "Table 3 robust context-dependent pathways:",
  nrow(table3),
  "\n"
)


cat(
  "Table 4 treatment-state evidence pathways:",
  nrow(table4),
  "\n"
)


cat(
  "\nOutput directory:\n",
  outdir,
  "\n"
)


cat("\n========================================\n")