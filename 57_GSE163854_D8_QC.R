# ============================================================
# GSE163854 - D8
# QUALITY CONTROL BEFORE DIFFERENTIAL EXPRESSION
#
# Frozen primary comparison:
#
#   Sensitive/passaged control:
#     P2A1, P2A2, P2A3
#
#   Third-cycle resistant:
#     F451, F452, F453
#
# Primary expression filter frozen BEFORE DE:
#
#   CPM >= 1 in >= 3 of the 6 primary samples
#
# This script performs QC ONLY.
#
# NO differential expression
# NO GSEA
# NO biological result inspection
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


raw_dir <-
  "data_raw/GSE163854"


qc_dir <-
  "results/GSE163854/QC"


dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("GSE163854 - D8 QUALITY CONTROL\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

packages <- c(
  "edgeR",
  "ggplot2",
  "ggrepel"
)


for (pkg in packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    if (pkg == "edgeR") {
      
      if (!requireNamespace(
        "BiocManager",
        quietly = TRUE
      )) {
        install.packages("BiocManager")
      }
      
      BiocManager::install(
        "edgeR",
        ask = FALSE,
        update = FALSE
      )
      
    } else {
      
      install.packages(pkg)
      
    }
  }
}


library(edgeR)
library(ggplot2)
library(ggrepel)


# ============================================================
# 2. LOAD RAW COUNTS
# ============================================================

count_file <- file.path(
  raw_dir,
  "GSE163854_raw_counts.rds"
)


if (!file.exists(count_file)) {
  
  stop(
    "Raw-count RDS from Script 56 not found."
  )
}


counts_raw <- readRDS(
  count_file
)


cat(
  "Raw matrix:",
  nrow(counts_raw),
  "genes x",
  ncol(counts_raw) - 1,
  "samples\n"
)


# ============================================================
# 3. VERIFY EXPECTED STRUCTURE
# ============================================================

if (
  colnames(counts_raw)[1] !=
  "gene_symbol"
) {
  
  stop(
    "First column is not gene_symbol."
  )
}


expected_samples <- c(
  "F1A3",
  "F1A4",
  "F1A5",
  "F251",
  "F351",
  "F352",
  "F451",
  "F452",
  "F453",
  "P2A1",
  "P2A2",
  "P2A3"
)


observed_samples <- colnames(
  counts_raw
)[-1]


if (!setequal(
  expected_samples,
  observed_samples
)) {
  
  cat(
    "\nObserved sample columns:\n"
  )
  
  print(
    observed_samples
  )
  
  
  stop(
    "Sample columns differ from frozen D8 design."
  )
}


# ============================================================
# 4. BUILD COUNT MATRIX
# ============================================================

genes <- as.character(
  counts_raw$gene_symbol
)


counts <- as.matrix(
  counts_raw[
    ,
    expected_samples,
    drop = FALSE
  ]
)


storage.mode(counts) <- "numeric"


rownames(counts) <- genes


# ============================================================
# 5. FINAL COUNT VALIDATION
# ============================================================

if (anyNA(counts)) {
  
  stop(
    "Count matrix contains NA values."
  )
}


if (any(counts < 0)) {
  
  stop(
    "Count matrix contains negative values."
  )
}


integer_fraction <- mean(
  abs(
    counts -
      round(counts)
  ) < 1e-8
)


cat(
  "Integer-like count fraction:",
  integer_fraction,
  "\n"
)


if (integer_fraction != 1) {
  
  stop(
    "Counts are not completely integer-valued."
  )
}


# ============================================================
# 6. SAMPLE METADATA
# ============================================================

sample_metadata <- data.frame(
  
  sample = expected_samples,
  
  passage_group = c(
    "C1",
    "C1",
    "C1",
    "C1",
    "C2",
    "C2",
    "C3",
    "C3",
    "C3",
    "S2",
    "S2",
    "S2"
  ),
  
  phenotype = c(
    "Sensitive",
    "Sensitive",
    "Mixed",
    "Mixed",
    "Resistant",
    "Resistant",
    "Resistant",
    "Resistant",
    "Resistant",
    "Sensitive_control",
    "Sensitive_control",
    "Sensitive_control"
  ),
  
  treatment_history = c(
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Niraparib_selected",
    "Passaged_no_niraparib",
    "Passaged_no_niraparib",
    "Passaged_no_niraparib"
  ),
  
  primary_group = c(
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    "Resistant",
    "Resistant",
    "Resistant",
    "Sensitive",
    "Sensitive",
    "Sensitive"
  ),
  
  stringsAsFactors = FALSE
)


rownames(
  sample_metadata
) <- sample_metadata$sample


write.csv(
  sample_metadata,
  file.path(
    qc_dir,
    "D8_sample_metadata_frozen.csv"
  ),
  row.names = FALSE
)


cat(
  "\nFrozen D8 metadata:\n\n"
)


print(
  sample_metadata,
  row.names = FALSE
)


# ============================================================
# 7. LIBRARY SIZE QC
# ============================================================

library_sizes <- colSums(
  counts
)


library_table <- data.frame(
  
  sample =
    names(library_sizes),
  
  library_size =
    as.numeric(
      library_sizes
    ),
  
  million_reads =
    as.numeric(
      library_sizes
    ) / 1e6,
  
  fraction_of_median =
    as.numeric(
      library_sizes /
        median(library_sizes)
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  library_table,
  file.path(
    qc_dir,
    "D8_library_sizes.csv"
  ),
  row.names = FALSE
)


cat(
  "\nLibrary sizes:\n\n"
)


print(
  library_table,
  row.names = FALSE
)


# ============================================================
# 8. ALL-SAMPLE QC FILTER
#
# Visualization only.
#
# Because several biological stages contain n=2,
# retain genes with CPM >= 1 in >= 2 samples.
# ============================================================

dge_all <- DGEList(
  counts = counts
)


cpm_all <- cpm(
  dge_all
)


keep_all <- rowSums(
  cpm_all >= 1
) >= 2


cat(
  "\nAll-sample QC genes before filter:",
  nrow(counts),
  "\n"
)


cat(
  "All-sample QC genes retained:",
  sum(keep_all),
  "\n"
)


dge_all <- dge_all[
  keep_all,
  ,
  keep.lib.sizes = FALSE
]


dge_all <- calcNormFactors(
  dge_all,
  method = "TMM"
)


# ============================================================
# 9. TMM NORMALIZATION FACTORS
# ============================================================

norm_table <- data.frame(
  
  sample =
    colnames(dge_all),
  
  library_size =
    dge_all$samples$lib.size,
  
  norm_factor =
    dge_all$samples$norm.factors,
  
  effective_library_size =
    dge_all$samples$lib.size *
    dge_all$samples$norm.factors,
  
  stringsAsFactors = FALSE
)


cat(
  "\nTMM normalization factors:\n\n"
)


print(
  norm_table,
  row.names = FALSE
)


write.csv(
  norm_table,
  file.path(
    qc_dir,
    "D8_TMM_normalization_factors.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. LOG-CPM FOR QC
# ============================================================

logCPM_all <- cpm(
  dge_all,
  log = TRUE,
  prior.count = 2
)


# ============================================================
# 11. SAMPLE PEARSON CORRELATIONS
# ============================================================

sample_cor <- cor(
  logCPM_all,
  method = "pearson"
)


write.csv(
  sample_cor,
  file.path(
    qc_dir,
    "D8_sample_correlations.csv"
  )
)


cat(
  "\nSample Pearson correlations:\n\n"
)


print(
  round(
    sample_cor,
    3
  )
)


# ============================================================
# 12. CORRELATION HEATMAP DATA
# ============================================================

cor_long <- as.data.frame(
  as.table(
    sample_cor
  )
)


colnames(
  cor_long
) <- c(
  "Sample_1",
  "Sample_2",
  "Correlation"
)


p_cor <- ggplot(
  cor_long,
  aes(
    x = Sample_1,
    y = Sample_2,
    fill = Correlation
  )
) +
  
  geom_tile() +
  
  scale_fill_gradient(
    name = "Pearson r"
  ) +
  
  coord_equal() +
  
  labs(
    title =
      "GSE163854 sample correlation",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


ggsave(
  file.path(
    qc_dir,
    "D8_sample_correlation_heatmap.pdf"
  ),
  p_cor,
  width = 8,
  height = 7,
  bg = "white"
)


# ============================================================
# 13. ALL-SAMPLE PCA
#
# Top 5000 most variable expressed genes.
# ============================================================

gene_variance_all <- apply(
  logCPM_all,
  1,
  var
)


n_top_all <- min(
  5000,
  length(
    gene_variance_all
  )
)


top_all <- order(
  gene_variance_all,
  decreasing = TRUE
)[
  seq_len(
    n_top_all
  )
]


pca_all <- prcomp(
  t(
    logCPM_all[
      top_all,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)


variance_all <- 100 *
  pca_all$sdev^2 /
  sum(
    pca_all$sdev^2
  )


pca_all_df <- data.frame(
  
  sample =
    rownames(
      pca_all$x
    ),
  
  PC1 =
    pca_all$x[
      ,
      1
    ],
  
  PC2 =
    pca_all$x[
      ,
      2
    ],
  
  stringsAsFactors = FALSE
)


pca_all_df <- merge(
  pca_all_df,
  sample_metadata,
  by = "sample",
  sort = FALSE
)


cat(
  "\nALL-SAMPLE PCA:\n"
)


cat(
  "PC1:",
  round(
    variance_all[1],
    2
  ),
  "%\n"
)


cat(
  "PC2:",
  round(
    variance_all[2],
    2
  ),
  "%\n"
)


print(
  pca_all_df
)


write.csv(
  pca_all_df,
  file.path(
    qc_dir,
    "D8_all_sample_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


p_all <- ggplot(
  pca_all_df,
  aes(
    x = PC1,
    y = PC2,
    shape = phenotype,
    label = sample
  )
) +
  
  geom_point(
    size = 3.5
  ) +
  
  ggrepel::geom_text_repel(
    size = 3,
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "GSE163854 - all-sample PCA",
    
    subtitle =
      "Serial PH039 niraparib-resistance trajectory",
    
    x = paste0(
      "PC1 (",
      round(
        variance_all[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        variance_all[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "Phenotype"
  ) +
  
  theme_classic(
    base_size = 11
  )


ggsave(
  file.path(
    qc_dir,
    "D8_all_sample_PCA.pdf"
  ),
  p_all,
  width = 7,
  height = 6,
  bg = "white"
)


# ============================================================
# 14. PRIMARY CONTRAST
# ============================================================

sensitive_primary <- c(
  "P2A1",
  "P2A2",
  "P2A3"
)


resistant_primary <- c(
  "F451",
  "F452",
  "F453"
)


primary_samples <- c(
  sensitive_primary,
  resistant_primary
)


primary_counts <- counts[
  ,
  primary_samples,
  drop = FALSE
]


primary_group <- factor(
  c(
    rep(
      "Sensitive",
      3
    ),
    rep(
      "Resistant",
      3
    )
  ),
  levels = c(
    "Sensitive",
    "Resistant"
  )
)


# ============================================================
# 15. FROZEN PRIMARY EXPRESSION FILTER
#
# CPM >= 1 in >= 3 of the 6 primary samples.
# ============================================================

dge_primary_unfiltered <- DGEList(
  counts = primary_counts,
  group = primary_group
)


primary_cpm <- cpm(
  dge_primary_unfiltered
)


keep_primary <- rowSums(
  primary_cpm >= 1
) >= 3


cat(
  "\n============================================================\n"
)


cat(
  "FROZEN PRIMARY FILTER\n"
)


cat(
  "============================================================\n\n"
)


cat(
  "Rule: CPM >= 1 in >= 3 of 6 primary samples\n"
)


cat(
  "Genes before filtering:",
  nrow(primary_counts),
  "\n"
)


cat(
  "Genes retained:",
  sum(keep_primary),
  "\n"
)


cat(
  "Genes removed:",
  sum(!keep_primary),
  "\n"
)


filter_table <- data.frame(
  
  gene =
    rownames(primary_counts),
  
  retained =
    keep_primary,
  
  stringsAsFactors = FALSE
)


write.csv(
  filter_table,
  file.path(
    qc_dir,
    "D8_primary_filter_status.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. PRIMARY TMM NORMALIZATION
# ============================================================

dge_primary <- dge_primary_unfiltered[
  keep_primary,
  ,
  keep.lib.sizes = FALSE
]


dge_primary <- calcNormFactors(
  dge_primary,
  method = "TMM"
)


primary_norm_table <- data.frame(
  
  sample =
    colnames(dge_primary),
  
  group =
    primary_group,
  
  library_size =
    dge_primary$samples$lib.size,
  
  norm_factor =
    dge_primary$samples$norm.factors,
  
  effective_library_size =
    dge_primary$samples$lib.size *
    dge_primary$samples$norm.factors,
  
  stringsAsFactors = FALSE
)


cat(
  "\nPrimary TMM factors:\n\n"
)


print(
  primary_norm_table,
  row.names = FALSE
)


write.csv(
  primary_norm_table,
  file.path(
    qc_dir,
    "D8_primary_TMM_factors.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. PRIMARY LOG-CPM
# ============================================================

primary_logCPM <- cpm(
  dge_primary,
  log = TRUE,
  prior.count = 2
)


# ============================================================
# 18. PRIMARY PCA
# ============================================================

primary_variance <- apply(
  primary_logCPM,
  1,
  var
)


n_top_primary <- min(
  5000,
  length(
    primary_variance
  )
)


top_primary <- order(
  primary_variance,
  decreasing = TRUE
)[
  seq_len(
    n_top_primary
  )
]


pca_primary <- prcomp(
  t(
    primary_logCPM[
      top_primary,
      ,
      drop = FALSE
    ]
  ),
  center = TRUE,
  scale. = FALSE
)


variance_primary <- 100 *
  pca_primary$sdev^2 /
  sum(
    pca_primary$sdev^2
  )


primary_pca_df <- data.frame(
  
  sample =
    rownames(
      pca_primary$x
    ),
  
  PC1 =
    pca_primary$x[
      ,
      1
    ],
  
  PC2 =
    pca_primary$x[
      ,
      2
    ],
  
  group =
    primary_group,
  
  stringsAsFactors = FALSE
)


cat(
  "\nPRIMARY PCA:\n"
)


cat(
  "PC1:",
  round(
    variance_primary[1],
    2
  ),
  "%\n"
)


cat(
  "PC2:",
  round(
    variance_primary[2],
    2
  ),
  "%\n\n"
)


print(
  primary_pca_df,
  row.names = FALSE
)


write.csv(
  primary_pca_df,
  file.path(
    qc_dir,
    "D8_primary_PCA_coordinates.csv"
  ),
  row.names = FALSE
)


p_primary <- ggplot(
  primary_pca_df,
  aes(
    x = PC1,
    y = PC2,
    shape = group,
    label = sample
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  ggrepel::geom_text_repel(
    size = 3.2,
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "GSE163854 - frozen primary contrast",
    
    subtitle =
      "Passaged sensitive control vs third-cycle resistant PH039",
    
    x = paste0(
      "PC1 (",
      round(
        variance_primary[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        variance_primary[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "State"
  ) +
  
  theme_classic(
    base_size = 11
  )


ggsave(
  file.path(
    qc_dir,
    "D8_primary_PCA.pdf"
  ),
  p_primary,
  width = 7,
  height = 6,
  bg = "white"
)


# ============================================================
# 19. PRIMARY WITHIN-GROUP CORRELATIONS
# ============================================================

primary_cor <- cor(
  primary_logCPM,
  method = "pearson"
)


sensitive_cor <- primary_cor[
  sensitive_primary,
  sensitive_primary,
  drop = FALSE
]


resistant_cor <- primary_cor[
  resistant_primary,
  resistant_primary,
  drop = FALSE
]


upper_values <- function(x) {
  
  x[
    upper.tri(x)
  ]
}


sensitive_cor_values <- upper_values(
  sensitive_cor
)


resistant_cor_values <- upper_values(
  resistant_cor
)


cat(
  "\n============================================================\n"
)


cat(
  "PRIMARY WITHIN-GROUP CORRELATIONS\n"
)


cat(
  "============================================================\n\n"
)


cat(
  "Sensitive controls:\n"
)


cat(
  "Mean =",
  round(
    mean(
      sensitive_cor_values
    ),
    5
  ),
  "\n"
)


cat(
  "Min  =",
  round(
    min(
      sensitive_cor_values
    ),
    5
  ),
  "\n"
)


cat(
  "Max  =",
  round(
    max(
      sensitive_cor_values
    ),
    5
  ),
  "\n\n"
)


cat(
  "C3 resistant:\n"
)


cat(
  "Mean =",
  round(
    mean(
      resistant_cor_values
    ),
    5
  ),
  "\n"
)


cat(
  "Min  =",
  round(
    min(
      resistant_cor_values
    ),
    5
  ),
  "\n"
)


cat(
  "Max  =",
  round(
    max(
      resistant_cor_values
    ),
    5
  ),
  "\n"
)


# ============================================================
# 20. PRIMARY MDS
# ============================================================

pdf(
  file.path(
    qc_dir,
    "D8_primary_MDS.pdf"
  ),
  width = 7,
  height = 6
)


plotMDS(
  dge_primary,
  labels = colnames(
    dge_primary
  )
)


dev.off()


# ============================================================
# 21. SAVE QC OBJECTS FOR LATER ANALYSIS / FIGURES
# ============================================================

qc_objects <- list(
  
  counts =
    counts,
  
  sample_metadata =
    sample_metadata,
  
  all_keep =
    keep_all,
  
  all_dge =
    dge_all,
  
  all_logCPM =
    logCPM_all,
  
  sample_correlations =
    sample_cor,
  
  all_pca = list(
    coordinates =
      pca_all_df,
    variance =
      variance_all
  ),
  
  primary_samples =
    primary_samples,
  
  primary_group =
    primary_group,
  
  primary_keep =
    keep_primary,
  
  primary_dge =
    dge_primary,
  
  primary_logCPM =
    primary_logCPM,
  
  primary_correlations =
    primary_cor,
  
  primary_pca = list(
    coordinates =
      primary_pca_df,
    variance =
      variance_primary
  )
)


saveRDS(
  qc_objects,
  file.path(
    qc_dir,
    "D8_QC_objects.rds"
  )
)


# ============================================================
# 22. QC SUMMARY
# ============================================================

qc_summary <- data.frame(
  
  metric = c(
    
    "Raw_gene_rows",
    
    "All_QC_genes_retained",
    
    "Primary_genes_retained",
    
    "Primary_PC1_percent",
    
    "Primary_PC2_percent",
    
    "Sensitive_mean_correlation",
    
    "Sensitive_min_correlation",
    
    "Resistant_mean_correlation",
    
    "Resistant_min_correlation"
    
  ),
  
  value = c(
    
    nrow(counts),
    
    sum(keep_all),
    
    sum(keep_primary),
    
    variance_primary[1],
    
    variance_primary[2],
    
    mean(
      sensitive_cor_values
    ),
    
    min(
      sensitive_cor_values
    ),
    
    mean(
      resistant_cor_values
    ),
    
    min(
      resistant_cor_values
    )
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  qc_summary,
  file.path(
    qc_dir,
    "D8_QC_summary.csv"
  ),
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)


cat(
  "D8 QC SUMMARY\n"
)


cat(
  "============================================================\n\n"
)


print(
  qc_summary,
  row.names = FALSE
)


# ============================================================
# 23. SESSION INFO
# ============================================================

sink(
  file.path(
    qc_dir,
    "D8_QC_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. END
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D8 QC COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Primary comparison remains frozen as:\n"
)


cat(
  "Sensitive: P2A1, P2A2, P2A3\n"
)


cat(
  "Resistant: F451, F452, F453\n\n"
)


cat(
  "NO differential expression or GSEA was performed.\n"
)


cat(
  "If QC is acceptable, the next script will perform:\n",
  "  TMM -> voom -> limma -> Hallmark GSEA\n",
  "  plus independent edgeR-QL sensitivity analysis.\n"
)


cat(
  "\n============================================================\n"
)