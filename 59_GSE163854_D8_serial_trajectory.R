# ============================================================
# GSE163854 - D8
# SECONDARY SERIAL-RESISTANCE TRAJECTORY ANALYSIS
# CLEAN VERSION
#
# IMPORTANT:
#   Exploratory / secondary analysis only.
#
#   It does NOT:
#     - replace the frozen S2 vs C3 primary comparison
#     - create extra accession-level votes
#     - modify Script 58 results
#
# Frozen gene universe:
#   13,410 genes
#
# Groups:
#
#   S2_control:
#     P2A1, P2A2, P2A3
#
#   C1_sensitive:
#     F1A3, F1A4
#
#   C1_mixed:
#     F1A5, F251
#
#   C2_resistant:
#     F351, F352
#
#   C3_resistant:
#     F451, F452, F453
#
# Contrasts:
#
#   C1_sensitive vs S2
#   C1_mixed vs S2
#   C2_resistant vs S2
#   C3_resistant vs S2
#   C1_mixed vs C1_sensitive
#
# Pipeline:
#   TMM -> voom -> limma -> moderated-t Hallmark GSEA
#
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

qc_dir <-
  "results/GSE163854/QC"

primary_dir <-
  "results/GSE163854/D8_analysis"

trajectory_dir <-
  "results/GSE163854/D8_trajectory"


dir.create(
  trajectory_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("GSE163854 - D8 SERIAL TRAJECTORY ANALYSIS\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

cran_packages <- c(
  "msigdbr",
  "ggplot2"
)

bioc_packages <- c(
  "edgeR",
  "limma",
  "clusterProfiler"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages("BiocManager")
}


for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


library(edgeR)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)


# ============================================================
# 2. LOAD FROZEN QC OBJECT
# ============================================================

qc_file <- file.path(
  qc_dir,
  "D8_QC_objects.rds"
)


if (!file.exists(qc_file)) {
  
  stop(
    "D8 QC object from Script 57 not found."
  )
}


qc <- readRDS(
  qc_file
)


# ============================================================
# 3. LOAD PRIMARY D8 HALLMARK RESULT
# ============================================================

primary_gsea_file <- file.path(
  primary_dir,
  "D8_primary_Hallmark_GSEA.csv"
)


if (!file.exists(primary_gsea_file)) {
  
  stop(
    "Primary D8 Hallmark result from Script 58 not found."
  )
}


primary_gsea <- read.csv(
  primary_gsea_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. VERIFY FROZEN GENE UNIVERSE
# ============================================================

primary_genes <- rownames(
  qc$primary_dge
)


cat(
  "Frozen primary gene universe:",
  length(primary_genes),
  "\n"
)


if (length(primary_genes) != 13410) {
  
  stop(
    paste0(
      "Expected 13,410 frozen genes, found ",
      length(primary_genes),
      "."
    )
  )
}


# ============================================================
# 5. DEFINE ALL 12 SAMPLES
# ============================================================

all_samples <- c(
  
  "P2A1",
  "P2A2",
  "P2A3",
  
  "F1A3",
  "F1A4",
  
  "F1A5",
  "F251",
  
  "F351",
  "F352",
  
  "F451",
  "F452",
  "F453"
)


if (!all(
  all_samples %in%
  colnames(qc$counts)
)) {
  
  stop(
    "Expected trajectory samples are missing from frozen count matrix."
  )
}


# ============================================================
# 6. RESTRICT TO FROZEN 13,410 GENES
# ============================================================

missing_genes <- setdiff(
  primary_genes,
  rownames(qc$counts)
)


if (length(missing_genes) > 0) {
  
  stop(
    "Some frozen primary genes are absent from full count matrix."
  )
}


counts <- qc$counts[
  primary_genes,
  all_samples,
  drop = FALSE
]


storage.mode(counts) <- "numeric"


cat(
  "Trajectory matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n"
)


# ============================================================
# 7. SAMPLE GROUPS
# ============================================================

group <- factor(
  
  c(
    
    "S2_control",
    "S2_control",
    "S2_control",
    
    "C1_sensitive",
    "C1_sensitive",
    
    "C1_mixed",
    "C1_mixed",
    
    "C2_resistant",
    "C2_resistant",
    
    "C3_resistant",
    "C3_resistant",
    "C3_resistant"
    
  ),
  
  levels = c(
    "S2_control",
    "C1_sensitive",
    "C1_mixed",
    "C2_resistant",
    "C3_resistant"
  )
)


sample_metadata <- data.frame(
  
  sample =
    all_samples,
  
  group =
    group,
  
  stringsAsFactors = FALSE
)


rownames(
  sample_metadata
) <- sample_metadata$sample


cat(
  "\nTrajectory sample design:\n\n"
)


print(
  sample_metadata,
  row.names = FALSE
)


write.csv(
  sample_metadata,
  file.path(
    trajectory_dir,
    "D8_trajectory_sample_design.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8. DGE OBJECT + TMM NORMALIZATION
# ============================================================

dge <- edgeR::DGEList(
  counts = counts
)


dge <- edgeR::normLibSizes(
  dge,
  method = "TMM"
)


norm_table <- data.frame(
  
  sample =
    colnames(dge),
  
  group =
    group,
  
  library_size =
    dge$samples$lib.size,
  
  norm_factor =
    dge$samples$norm.factors,
  
  effective_library_size =
    dge$samples$lib.size *
    dge$samples$norm.factors,
  
  stringsAsFactors = FALSE
)


cat(
  "\nTMM factors:\n\n"
)


print(
  norm_table,
  row.names = FALSE
)


write.csv(
  norm_table,
  file.path(
    trajectory_dir,
    "D8_trajectory_TMM_factors.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. DESIGN MATRIX
# ============================================================

design <- model.matrix(
  ~ 0 + group
)


colnames(design) <- levels(
  group
)


rownames(design) <- all_samples


cat(
  "\nDesign matrix:\n\n"
)


print(
  design
)


# ============================================================
# 10. CONTRAST MATRIX
# ============================================================

contrast_matrix <- limma::makeContrasts(
  
  C1_sensitive_vs_S2 =
    C1_sensitive - S2_control,
  
  C1_mixed_vs_S2 =
    C1_mixed - S2_control,
  
  C2_resistant_vs_S2 =
    C2_resistant - S2_control,
  
  C3_resistant_vs_S2 =
    C3_resistant - S2_control,
  
  C1_mixed_vs_C1_sensitive =
    C1_mixed - C1_sensitive,
  
  levels = design
)


cat(
  "\nContrast matrix:\n\n"
)


print(
  contrast_matrix
)


# ============================================================
# 11. VOOM
# ============================================================

pdf(
  file.path(
    trajectory_dir,
    "D8_trajectory_voom_mean_variance.pdf"
  ),
  width = 7,
  height = 6
)


v <- limma::voom(
  dge,
  design,
  plot = TRUE
)


dev.off()


# ============================================================
# 12. LIMMA MODEL
# ============================================================

fit <- limma::lmFit(
  v,
  design
)


fit <- limma::contrasts.fit(
  fit,
  contrast_matrix
)


fit <- limma::eBayes(
  fit
)


# ============================================================
# 13. DE EXTRACTION FUNCTION
# ============================================================

extract_de <- function(
    fitted_model,
    contrast_name
) {
  
  tab <- limma::topTable(
    
    fitted_model,
    
    coef =
      contrast_name,
    
    number =
      Inf,
    
    adjust.method =
      "BH",
    
    sort.by =
      "none"
  )
  
  
  tab$gene <- rownames(
    tab
  )
  
  
  tab <- tab[
    ,
    c(
      "gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    ),
    drop = FALSE
  ]
  
  
  rownames(tab) <- NULL
  
  
  return(tab)
}


# ============================================================
# 14. EXTRACT EACH CONTRAST EXPLICITLY
#
# No named-list indexing.
# ============================================================

de_C1_sensitive_vs_S2 <- extract_de(
  fit,
  "C1_sensitive_vs_S2"
)


de_C1_mixed_vs_S2 <- extract_de(
  fit,
  "C1_mixed_vs_S2"
)


de_C2_resistant_vs_S2 <- extract_de(
  fit,
  "C2_resistant_vs_S2"
)


de_C3_resistant_vs_S2 <- extract_de(
  fit,
  "C3_resistant_vs_S2"
)


de_C1_mixed_vs_C1_sensitive <- extract_de(
  fit,
  "C1_mixed_vs_C1_sensitive"
)


# ============================================================
# 15. SAVE DE TABLES
# ============================================================

write.csv(
  de_C1_sensitive_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_C1_sensitive_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  de_C1_mixed_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_C1_mixed_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  de_C2_resistant_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_C2_resistant_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  de_C3_resistant_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_C3_resistant_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  de_C1_mixed_vs_C1_sensitive,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_C1_mixed_vs_C1_sensitive.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. DE SUMMARY
# ============================================================

summarize_de <- function(
    tab,
    contrast_name
) {
  
  data.frame(
    
    contrast =
      contrast_name,
    
    genes_tested =
      nrow(tab),
    
    FDR_lt_0.05 =
      sum(
        tab$adj.P.Val < 0.05
      ),
    
    FDR_lt_0.05_abs_logFC_ge_1 =
      sum(
        tab$adj.P.Val < 0.05 &
          abs(tab$logFC) >= 1
      ),
    
    stringsAsFactors = FALSE
  )
}


de_summary <- rbind(
  
  summarize_de(
    de_C1_sensitive_vs_S2,
    "C1_sensitive_vs_S2"
  ),
  
  summarize_de(
    de_C1_mixed_vs_S2,
    "C1_mixed_vs_S2"
  ),
  
  summarize_de(
    de_C2_resistant_vs_S2,
    "C2_resistant_vs_S2"
  ),
  
  summarize_de(
    de_C3_resistant_vs_S2,
    "C3_resistant_vs_S2"
  ),
  
  summarize_de(
    de_C1_mixed_vs_C1_sensitive,
    "C1_mixed_vs_C1_sensitive"
  )
  
)


rownames(de_summary) <- NULL


cat("\n")
cat("============================================================\n")
cat("TRAJECTORY DE SUMMARY\n")
cat("============================================================\n\n")


print(
  de_summary,
  row.names = FALSE
)


write.csv(
  de_summary,
  file.path(
    trajectory_dir,
    "D8_trajectory_DE_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. HALLMARK GENE SETS
# ============================================================

hallmark_msig <- msigdbr::msigdbr(
  
  db_species =
    "HS",
  
  species =
    "Homo sapiens",
  
  collection =
    "H"
)


term2gene <- unique(
  hallmark_msig[
    ,
    c(
      "gs_name",
      "gene_symbol"
    )
  ]
)


term2gene <- term2gene[
  !is.na(term2gene$gene_symbol) &
    term2gene$gene_symbol != "",
  ,
  drop = FALSE
]


cat(
  "\nHallmark gene sets:",
  length(
    unique(
      term2gene$gs_name
    )
  ),
  "\n"
)


# ============================================================
# 18. GSEA FUNCTION
# ============================================================

run_gsea <- function(
    de_table,
    contrast_name
) {
  
  ranked <- de_table$t
  
  
  names(ranked) <- de_table$gene
  
  
  ranked <- ranked[
    is.finite(ranked)
  ]
  
  
  ranked <- sort(
    ranked,
    decreasing = TRUE
  )
  
  
  tie_fraction <-
    1 -
    length(
      unique(ranked)
    ) /
    length(ranked)
  
  
  cat(
    "\n",
    contrast_name,
    ": ranking genes = ",
    length(ranked),
    "; tie fraction = ",
    round(
      tie_fraction,
      6
    ),
    "\n",
    sep = ""
  )
  
  
  set.seed(
    20261110
  )
  
  
  gsea_object <- clusterProfiler::GSEA(
    
    geneList =
      ranked,
    
    TERM2GENE =
      term2gene,
    
    exponent =
      1,
    
    minGSSize =
      10,
    
    maxGSSize =
      500,
    
    eps =
      0,
    
    pvalueCutoff =
      1,
    
    pAdjustMethod =
      "BH",
    
    verbose =
      FALSE,
    
    by =
      "fgsea"
  )
  
  
  result <- as.data.frame(
    gsea_object@result
  )
  
  
  if (nrow(result) == 0) {
    
    stop(
      paste(
        "No Hallmark GSEA results for",
        contrast_name
      )
    )
  }
  
  
  result$contrast <-
    contrast_name
  
  
  return(result)
}


# ============================================================
# 19. RUN GSEA EXPLICITLY
# ============================================================

gsea_C1_sensitive_vs_S2 <- run_gsea(
  de_C1_sensitive_vs_S2,
  "C1_sensitive_vs_S2"
)


gsea_C1_mixed_vs_S2 <- run_gsea(
  de_C1_mixed_vs_S2,
  "C1_mixed_vs_S2"
)


gsea_C2_resistant_vs_S2 <- run_gsea(
  de_C2_resistant_vs_S2,
  "C2_resistant_vs_S2"
)


gsea_C3_resistant_vs_S2 <- run_gsea(
  de_C3_resistant_vs_S2,
  "C3_resistant_vs_S2"
)


gsea_C1_mixed_vs_C1_sensitive <- run_gsea(
  de_C1_mixed_vs_C1_sensitive,
  "C1_mixed_vs_C1_sensitive"
)


# ============================================================
# 20. SAVE GSEA TABLES
# ============================================================

write.csv(
  gsea_C1_sensitive_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_C1_sensitive_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_C1_mixed_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_C1_mixed_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_C2_resistant_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_C2_resistant_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_C3_resistant_vs_S2,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_C3_resistant_vs_S2.csv"
  ),
  row.names = FALSE
)


write.csv(
  gsea_C1_mixed_vs_C1_sensitive,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_C1_mixed_vs_C1_sensitive.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. GSEA SUMMARY
# ============================================================

summarize_gsea <- function(
    tab,
    contrast_name
) {
  
  data.frame(
    
    contrast =
      contrast_name,
    
    Hallmarks =
      nrow(tab),
    
    FDR_lt_0.05 =
      sum(
        tab$p.adjust < 0.05
      ),
    
    positive_sig =
      sum(
        tab$p.adjust < 0.05 &
          tab$NES > 0
      ),
    
    negative_sig =
      sum(
        tab$p.adjust < 0.05 &
          tab$NES < 0
      ),
    
    stringsAsFactors = FALSE
  )
}


gsea_summary <- rbind(
  
  summarize_gsea(
    gsea_C1_sensitive_vs_S2,
    "C1_sensitive_vs_S2"
  ),
  
  summarize_gsea(
    gsea_C1_mixed_vs_S2,
    "C1_mixed_vs_S2"
  ),
  
  summarize_gsea(
    gsea_C2_resistant_vs_S2,
    "C2_resistant_vs_S2"
  ),
  
  summarize_gsea(
    gsea_C3_resistant_vs_S2,
    "C3_resistant_vs_S2"
  ),
  
  summarize_gsea(
    gsea_C1_mixed_vs_C1_sensitive,
    "C1_mixed_vs_C1_sensitive"
  )
  
)


rownames(gsea_summary) <- NULL


cat("\n")
cat("============================================================\n")
cat("TRAJECTORY HALLMARK SUMMARY\n")
cat("============================================================\n\n")


print(
  gsea_summary,
  row.names = FALSE
)


write.csv(
  gsea_summary,
  file.path(
    trajectory_dir,
    "D8_trajectory_Hallmark_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. HELPER TO PREPARE STAGE TABLE
# ============================================================

prepare_stage <- function(
    gsea_table,
    stage_prefix
) {
  
  out <- gsea_table[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    ),
    drop = FALSE
  ]
  
  
  colnames(out) <- c(
    
    "pathway",
    
    paste0(
      stage_prefix,
      "_NES"
    ),
    
    paste0(
      stage_prefix,
      "_FDR"
    )
  )
  
  
  return(out)
}


# ============================================================
# 23. BUILD FOUR-STAGE TRAJECTORY TABLE
# ============================================================

stage_C1_sensitive <- prepare_stage(
  gsea_C1_sensitive_vs_S2,
  "C1_sensitive_vs_S2"
)


stage_C1_mixed <- prepare_stage(
  gsea_C1_mixed_vs_S2,
  "C1_mixed_vs_S2"
)


stage_C2_resistant <- prepare_stage(
  gsea_C2_resistant_vs_S2,
  "C2_resistant_vs_S2"
)


stage_C3_resistant <- prepare_stage(
  gsea_C3_resistant_vs_S2,
  "C3_resistant_vs_S2"
)


trajectory_table <- merge(
  stage_C1_sensitive,
  stage_C1_mixed,
  by = "pathway",
  all = TRUE
)


trajectory_table <- merge(
  trajectory_table,
  stage_C2_resistant,
  by = "pathway",
  all = TRUE
)


trajectory_table <- merge(
  trajectory_table,
  stage_C3_resistant,
  by = "pathway",
  all = TRUE
)


# ============================================================
# 24. ADD PRIMARY SCRIPT-58 C3 RESULT
# ============================================================

primary_small <- primary_gsea[
  ,
  c(
    "ID",
    "NES",
    "p.adjust"
  ),
  drop = FALSE
]


colnames(primary_small) <- c(
  "pathway",
  "Primary_C3_NES",
  "Primary_C3_FDR"
)


trajectory_table <- merge(
  trajectory_table,
  primary_small,
  by = "pathway",
  all = TRUE
)


# ============================================================
# 25. VALIDATE 12-SAMPLE C3 CONTRAST
# AGAINST FROZEN PRIMARY C3 RESULT
# ============================================================

c3_complete <- trajectory_table[
  is.finite(
    trajectory_table$C3_resistant_vs_S2_NES
  ) &
    is.finite(
      trajectory_table$Primary_C3_NES
    ),
  ,
  drop = FALSE
]


c3_nes_rho <- cor(
  c3_complete$C3_resistant_vs_S2_NES,
  c3_complete$Primary_C3_NES,
  method = "spearman"
)


c3_same_direction <- mean(
  sign(
    c3_complete$C3_resistant_vs_S2_NES
  ) ==
    sign(
      c3_complete$Primary_C3_NES
    )
)


cat("\n")
cat("============================================================\n")
cat("TRAJECTORY MODEL VALIDATION AGAINST PRIMARY D8\n")
cat("============================================================\n\n")


cat(
  "Hallmarks compared:",
  nrow(c3_complete),
  "\n"
)


cat(
  "C3 NES Spearman rho:",
  round(
    c3_nes_rho,
    5
  ),
  "\n"
)


cat(
  "C3 direction agreement:",
  round(
    100 *
      c3_same_direction,
    2
  ),
  "%\n"
)


# ============================================================
# 26. PRIMARY SIGNIFICANT HALLMARKS
# ============================================================

primary_significant <- primary_gsea[
  primary_gsea$p.adjust < 0.05,
  ,
  drop = FALSE
]


primary_trajectory <- trajectory_table[
  trajectory_table$pathway %in%
    primary_significant$ID,
  ,
  drop = FALSE
]


primary_trajectory$final_sign <- sign(
  primary_trajectory$Primary_C3_NES
)


# ============================================================
# 27. SAME-DIRECTION FLAGS AT EACH STAGE
# ============================================================

primary_trajectory$C1_sensitive_same_final_direction <-
  sign(
    primary_trajectory$C1_sensitive_vs_S2_NES
  ) ==
  primary_trajectory$final_sign


primary_trajectory$C1_mixed_same_final_direction <-
  sign(
    primary_trajectory$C1_mixed_vs_S2_NES
  ) ==
  primary_trajectory$final_sign


primary_trajectory$C2_resistant_same_final_direction <-
  sign(
    primary_trajectory$C2_resistant_vs_S2_NES
  ) ==
  primary_trajectory$final_sign


primary_trajectory$C3_resistant_same_final_direction <-
  sign(
    primary_trajectory$C3_resistant_vs_S2_NES
  ) ==
  primary_trajectory$final_sign


# ============================================================
# 28. SIGNIFICANT + SAME-DIRECTION FLAGS
# ============================================================

primary_trajectory$C1_sensitive_support <-
  primary_trajectory$C1_sensitive_same_final_direction &
  primary_trajectory$C1_sensitive_vs_S2_FDR < 0.05


primary_trajectory$C1_mixed_support <-
  primary_trajectory$C1_mixed_same_final_direction &
  primary_trajectory$C1_mixed_vs_S2_FDR < 0.05


primary_trajectory$C2_resistant_support <-
  primary_trajectory$C2_resistant_same_final_direction &
  primary_trajectory$C2_resistant_vs_S2_FDR < 0.05


primary_trajectory$C3_resistant_support <-
  primary_trajectory$C3_resistant_same_final_direction &
  primary_trajectory$C3_resistant_vs_S2_FDR < 0.05


# ============================================================
# 29. FIRST SIGNIFICANT SAME-DIRECTION STAGE
# ============================================================

first_stage <- rep(
  "No_stage_FDR_lt_0.05",
  nrow(primary_trajectory)
)


for (i in seq_len(
  nrow(primary_trajectory)
)) {
  
  if (
    !is.na(
      primary_trajectory$C1_sensitive_support[i]
    ) &&
    primary_trajectory$C1_sensitive_support[i]
  ) {
    
    first_stage[i] <- "C1_sensitive"
    
  } else if (
    
    !is.na(
      primary_trajectory$C1_mixed_support[i]
    ) &&
    primary_trajectory$C1_mixed_support[i]
    
  ) {
    
    first_stage[i] <- "C1_mixed"
    
  } else if (
    
    !is.na(
      primary_trajectory$C2_resistant_support[i]
    ) &&
    primary_trajectory$C2_resistant_support[i]
    
  ) {
    
    first_stage[i] <- "C2_resistant"
    
  } else if (
    
    !is.na(
      primary_trajectory$C3_resistant_support[i]
    ) &&
    primary_trajectory$C3_resistant_support[i]
    
  ) {
    
    first_stage[i] <- "C3_resistant"
  }
}


primary_trajectory$first_same_direction_sig_stage <-
  first_stage


# ============================================================
# 30. SAME-CYCLE C1 MIXED VS C1 SENSITIVE
# ============================================================

same_cycle_small <- gsea_C1_mixed_vs_C1_sensitive[
  ,
  c(
    "ID",
    "NES",
    "p.adjust"
  ),
  drop = FALSE
]


colnames(same_cycle_small) <- c(
  "pathway",
  "C1_mixed_vs_sensitive_NES",
  "C1_mixed_vs_sensitive_FDR"
)


primary_trajectory <- merge(
  primary_trajectory,
  same_cycle_small,
  by = "pathway",
  all.x = TRUE
)


primary_trajectory$C1_mixed_vs_sensitive_same_final_direction <-
  sign(
    primary_trajectory$C1_mixed_vs_sensitive_NES
  ) ==
  primary_trajectory$final_sign


primary_trajectory$C1_mixed_vs_sensitive_support <-
  primary_trajectory$C1_mixed_vs_sensitive_same_final_direction &
  primary_trajectory$C1_mixed_vs_sensitive_FDR < 0.05


# ============================================================
# 31. RESTORE PRIMARY SIGNIFICANCE ORDER
# ============================================================

primary_order_table <- primary_significant[
  order(
    primary_significant$p.adjust,
    -abs(
      primary_significant$NES
    )
  ),
  ,
  drop = FALSE
]


primary_order <- primary_order_table$ID


primary_trajectory$order_index <- match(
  primary_trajectory$pathway,
  primary_order
)


primary_trajectory <- primary_trajectory[
  order(
    primary_trajectory$order_index
  ),
  ,
  drop = FALSE
]


# ============================================================
# 32. SAVE PRIMARY HALLMARK TRAJECTORY
# ============================================================

write.csv(
  primary_trajectory,
  file.path(
    trajectory_dir,
    "D8_primary_significant_Hallmark_trajectory.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 33. PRINT PRIMARY HALLMARK TRAJECTORIES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PRIMARY D8 HALLMARK TRAJECTORIES\n")
cat("============================================================\n\n")


print_table <- primary_trajectory[
  ,
  c(
    
    "pathway",
    
    "C1_sensitive_vs_S2_NES",
    "C1_sensitive_vs_S2_FDR",
    
    "C1_mixed_vs_S2_NES",
    "C1_mixed_vs_S2_FDR",
    
    "C2_resistant_vs_S2_NES",
    "C2_resistant_vs_S2_FDR",
    
    "C3_resistant_vs_S2_NES",
    "C3_resistant_vs_S2_FDR",
    
    "Primary_C3_NES",
    "Primary_C3_FDR",
    
    "first_same_direction_sig_stage",
    
    "C1_mixed_vs_sensitive_NES",
    "C1_mixed_vs_sensitive_FDR",
    
    "C1_mixed_vs_sensitive_support"
    
  ),
  drop = FALSE
]


print(
  print_table,
  row.names = FALSE
)


# ============================================================
# 34. PRETTY PATHWAY LABELS
# ============================================================

pretty_pathway <- function(x) {
  
  x <- sub(
    "^HALLMARK_",
    "",
    x
  )
  
  
  x <- gsub(
    "_",
    " ",
    x
  )
  
  
  x <- tools::toTitleCase(
    tolower(x)
  )
  
  
  x <- gsub(
    "Interferon Alpha Response",
    "Interferon alpha response",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Interferon Gamma Response",
    "Interferon gamma response",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Tnfa Signaling Via Nfkb",
    "TNFa signaling via NF-kB",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Il6 Jak Stat3 Signaling",
    "IL6-JAK-STAT3 signaling",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "E2f Targets",
    "E2F targets",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Kras Signaling Up",
    "KRAS signaling up",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Il2 Stat5 Signaling",
    "IL2-STAT5 signaling",
    x,
    fixed = TRUE
  )
  
  
  x <- gsub(
    "Tgf Beta Signaling",
    "TGF-beta signaling",
    x,
    fixed = TRUE
  )
  
  
  return(x)
}


# ============================================================
# 35. HEATMAP DATA
# ============================================================

heatmap_C1_sensitive <- data.frame(
  
  pathway =
    primary_trajectory$pathway,
  
  stage =
    "C1 sensitive",
  
  NES =
    primary_trajectory$C1_sensitive_vs_S2_NES,
  
  stringsAsFactors = FALSE
)


heatmap_C1_mixed <- data.frame(
  
  pathway =
    primary_trajectory$pathway,
  
  stage =
    "C1 mixed",
  
  NES =
    primary_trajectory$C1_mixed_vs_S2_NES,
  
  stringsAsFactors = FALSE
)


heatmap_C2 <- data.frame(
  
  pathway =
    primary_trajectory$pathway,
  
  stage =
    "C2 resistant",
  
  NES =
    primary_trajectory$C2_resistant_vs_S2_NES,
  
  stringsAsFactors = FALSE
)


heatmap_C3 <- data.frame(
  
  pathway =
    primary_trajectory$pathway,
  
  stage =
    "C3 resistant",
  
  NES =
    primary_trajectory$C3_resistant_vs_S2_NES,
  
  stringsAsFactors = FALSE
)


heatmap_long <- rbind(
  
  heatmap_C1_sensitive,
  
  heatmap_C1_mixed,
  
  heatmap_C2,
  
  heatmap_C3
)


heatmap_long$pretty_pathway <- pretty_pathway(
  heatmap_long$pathway
)


path_order <- pretty_pathway(
  primary_trajectory$pathway
)


heatmap_long$pretty_pathway <- factor(
  
  heatmap_long$pretty_pathway,
  
  levels = rev(
    path_order
  )
)


heatmap_long$stage <- factor(
  
  heatmap_long$stage,
  
  levels = c(
    "C1 sensitive",
    "C1 mixed",
    "C2 resistant",
    "C3 resistant"
  )
)


# ============================================================
# 36. TRAJECTORY HEATMAP
# ============================================================

p_heatmap <- ggplot(
  
  heatmap_long,
  
  aes(
    x = stage,
    y = pretty_pathway,
    fill = NES
  )
  
) +
  
  geom_tile(
    linewidth = 0.4
  ) +
  
  geom_text(
    
    aes(
      label = sprintf(
        "%.2f",
        NES
      )
    ),
    
    size = 2.7
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    name = "NES"
  ) +
  
  labs(
    
    title =
      "D8 serial emergence of final resistance-associated programs",
    
    subtitle =
      "Each selected stage compared with passaged untreated S2 control",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8
    )
  )


ggsave(
  filename = file.path(
    trajectory_dir,
    "D8_primary_Hallmark_trajectory_heatmap.pdf"
  ),
  plot = p_heatmap,
  width = 9,
  height = 8,
  bg = "white"
)


ggsave(
  filename = file.path(
    trajectory_dir,
    "D8_primary_Hallmark_trajectory_heatmap.png"
  ),
  plot = p_heatmap,
  width = 9,
  height = 8,
  dpi = 600,
  bg = "white"
)


# ============================================================
# 37. SAVE FULL 50-HALLMARK TABLE
# ============================================================

write.csv(
  trajectory_table,
  file.path(
    trajectory_dir,
    "D8_all_Hallmark_stage_trajectory.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 38. SUMMARY COUNTS
# ============================================================

first_stage_counts <- table(
  primary_trajectory$first_same_direction_sig_stage
)


same_cycle_support_n <- sum(
  primary_trajectory$C1_mixed_vs_sensitive_support,
  na.rm = TRUE
)


trajectory_summary <- data.frame(
  
  metric = c(
    
    "Frozen_gene_universe",
    
    "Primary_significant_Hallmarks",
    
    "C3_vs_primary_NES_Spearman_rho",
    
    "C3_vs_primary_direction_agreement",
    
    "Primary_Hallmarks_supported_in_C1_mixed_vs_sensitive"
    
  ),
  
  value = c(
    
    nrow(counts),
    
    nrow(primary_trajectory),
    
    c3_nes_rho,
    
    c3_same_direction,
    
    same_cycle_support_n
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  trajectory_summary,
  file.path(
    trajectory_dir,
    "D8_trajectory_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("TRAJECTORY SUMMARY\n")
cat("============================================================\n\n")


print(
  trajectory_summary,
  row.names = FALSE
)


cat(
  "\nFirst same-direction significant stage counts:\n\n"
)


print(
  first_stage_counts
)


cat(
  "\nPrimary final-state Hallmarks also supported by ",
  "C1 mixed vs C1 sensitive: ",
  same_cycle_support_n,
  "\n",
  sep = ""
)


# ============================================================
# 39. SESSION INFO
# ============================================================

sink(
  file.path(
    trajectory_dir,
    "D8_trajectory_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 40. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D8 SERIAL TRAJECTORY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(
  "IMPORTANT:\n",
  "  This remains a secondary exploratory analysis.\n",
  "  No additional accession-level votes were created.\n",
  "  Frozen D8 primary result remains S2 vs C3.\n\n"
)


cat(
  "Interpretation questions:\n",
  "  1. Which C3 programs already exist in C1-sensitive tumors?\n",
  "  2. Which first emerge in C1-mixed tumors?\n",
  "  3. Which emerge only after resistance is established?\n",
  "  4. Which distinguish C1-mixed from C1-sensitive tumors?\n"
)


cat(
  "\n============================================================\n"
)