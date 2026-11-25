# ============================================================
# GSE163854 - D8
# FINAL PUBLICATION FIGURES
#
# Uses ONLY frozen outputs from:
#   Script 57 - QC
#   Script 58 - primary DE/GSEA + edgeR sensitivity
#   Script 59 - serial trajectory
#
# NO DE re-analysis
# NO GSEA re-analysis
# NO threshold changes
#
# MAIN FIGURE
#   A - Primary S2 vs C3 PCA
#   B - Robust final-state Hallmark programs
#   C - Serial trajectory of final-state Hallmarks
#   D - D1 frozen-hypothesis comparison
#
# SUPPLEMENTARY FIGURE
#   S1 - Full 12-sample serial PCA
#   S2 - Sample correlation heatmap
#   S3 - limma vs edgeR gene-level validation
#   S4 - limma vs edgeR Hallmark validation
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

qc_dir <-
  "results/GSE163854/QC"

analysis_dir <-
  "results/GSE163854/D8_analysis"

trajectory_dir <-
  "results/GSE163854/D8_trajectory"

figure_dir <-
  "results/GSE163854/figures"


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n")
cat("============================================================\n")
cat("GSE163854 - D8 FINAL FIGURES\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "ggplot2",
  "ggrepel",
  "patchwork"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


library(ggplot2)
library(ggrepel)
library(patchwork)


# ============================================================
# 2. REQUIRED FILES
# ============================================================

qc_file <- file.path(
  qc_dir,
  "D8_QC_objects.rds"
)


hallmark_validation_file <- file.path(
  analysis_dir,
  "D8_Hallmark_limma_vs_edgeR.csv"
)


gene_validation_file <- file.path(
  analysis_dir,
  "D8_limma_vs_edgeR_gene_comparison.csv"
)


d1_validation_file <- file.path(
  analysis_dir,
  "D8_vs_D1_frozen_validation.csv"
)


trajectory_file <- file.path(
  trajectory_dir,
  "D8_primary_significant_Hallmark_trajectory.csv"
)


required_files <- c(
  qc_file,
  hallmark_validation_file,
  gene_validation_file,
  d1_validation_file,
  trajectory_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  cat(
    "Missing files:\n"
  )
  
  print(
    missing_files
  )
  
  stop(
    "Required frozen D8 outputs are missing."
  )
}


cat(
  "All required D8 frozen outputs found.\n"
)


# ============================================================
# 3. LOAD FROZEN OBJECTS
# ============================================================

qc <- readRDS(
  qc_file
)


hallmark_validation <- read.csv(
  hallmark_validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


gene_validation <- read.csv(
  gene_validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


d1_validation <- read.csv(
  d1_validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


trajectory <- read.csv(
  trajectory_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. PATHWAY LABEL FUNCTION
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
    "IFN-alpha response",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Interferon Gamma Response",
    "IFN-gamma response",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Tnfa Signaling Via Nfkb",
    "TNFa/NF-kB signaling",
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
  
  x <- gsub(
    "E2f Targets",
    "E2F targets",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "G2m Checkpoint",
    "G2M checkpoint",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Mtorc1 Signaling",
    "mTORC1 signaling",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Myc Targets V1",
    "MYC targets V1",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Myc Targets V2",
    "MYC targets V2",
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
    "Uv Response Up",
    "UV response up",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Uv Response Dn",
    "UV response down",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Oxidative Phosphorylation",
    "OXPHOS",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Epithelial Mesenchymal Transition",
    "EMT",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "Inflammatory Response",
    "Inflammatory response",
    x,
    fixed = TRUE
  )
  
  return(x)
}


# ============================================================
# 5. GLOBAL THEME
# ============================================================

publication_theme <- theme_classic(
  base_size = 11
) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.subtitle = element_text(
      size = 9
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 6. MAIN PANEL A
# PRIMARY S2 VS C3 PCA
# ============================================================

primary_pca <- qc$primary_pca$coordinates

primary_variance <- qc$primary_pca$variance


primary_pca$group <- factor(
  primary_pca$group,
  levels = c(
    "Sensitive",
    "Resistant"
  )
)


p_A <- ggplot(
  primary_pca,
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
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  
  labs(
    
    title =
      "A. Primary D8 resistance contrast",
    
    subtitle =
      "Passaged sensitive control vs third-cycle resistant PH039 PDX",
    
    x = paste0(
      "PC1 (",
      round(
        primary_variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        primary_variance[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "State"
  ) +
  
  publication_theme


# ============================================================
# 7. MAIN PANEL B
# ROBUST FINAL-STATE HALLMARKS
#
# Significant in both limma and edgeR
# with the same direction.
# ============================================================

robust <- hallmark_validation[
  !is.na(
    hallmark_validation$robust_significant
  ) &
    hallmark_validation$robust_significant,
  ,
  drop = FALSE
]


if (nrow(robust) == 0) {
  
  stop(
    "No robust Hallmark programs detected."
  )
}


robust$pretty <- pretty_pathway(
  robust$ID
)


robust <- robust[
  order(
    robust$NES_limma,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


robust$pretty <- factor(
  robust$pretty,
  levels = rev(
    robust$pretty
  )
)


robust_long <- rbind(
  
  data.frame(
    
    pathway =
      robust$pretty,
    
    method =
      "limma-voom",
    
    NES =
      robust$NES_limma,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    
    pathway =
      robust$pretty,
    
    method =
      "edgeR-QL",
    
    NES =
      robust$NES_edgeR,
    
    stringsAsFactors = FALSE
  )
)


robust_long$pathway <- factor(
  robust_long$pathway,
  levels = levels(
    robust$pretty
  )
)


p_B <- ggplot() +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  
  geom_segment(
    
    data = robust,
    
    aes(
      x = NES_limma,
      xend = NES_edgeR,
      y = pretty,
      yend = pretty
    ),
    
    linewidth = 0.5
  ) +
  
  geom_point(
    
    data = robust_long,
    
    aes(
      x = NES,
      y = pathway,
      shape = method
    ),
    
    size = 2.8
  ) +
  
  labs(
    
    title =
      "B. Robust final-state Hallmark programs",
    
    subtitle =
      "Significant in both limma-voom and edgeR-QL with concordant direction",
    
    x =
      "Normalized enrichment score",
    
    y =
      NULL,
    
    shape =
      "Analysis"
  ) +
  
  publication_theme +
  
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      size = 8
    )
  )


# ============================================================
# 8. MAIN PANEL C
# SERIAL TRAJECTORY OF FINAL-STATE PROGRAMS
# ============================================================

trajectory_C1_sensitive <- data.frame(
  
  pathway =
    trajectory$pathway,
  
  stage =
    "C1 sensitive",
  
  NES =
    trajectory$C1_sensitive_vs_S2_NES,
  
  FDR =
    trajectory$C1_sensitive_vs_S2_FDR,
  
  stringsAsFactors = FALSE
)


trajectory_C1_mixed <- data.frame(
  
  pathway =
    trajectory$pathway,
  
  stage =
    "C1 mixed",
  
  NES =
    trajectory$C1_mixed_vs_S2_NES,
  
  FDR =
    trajectory$C1_mixed_vs_S2_FDR,
  
  stringsAsFactors = FALSE
)


trajectory_C2 <- data.frame(
  
  pathway =
    trajectory$pathway,
  
  stage =
    "C2 resistant",
  
  NES =
    trajectory$C2_resistant_vs_S2_NES,
  
  FDR =
    trajectory$C2_resistant_vs_S2_FDR,
  
  stringsAsFactors = FALSE
)


trajectory_C3 <- data.frame(
  
  pathway =
    trajectory$pathway,
  
  stage =
    "C3 resistant",
  
  NES =
    trajectory$C3_resistant_vs_S2_NES,
  
  FDR =
    trajectory$C3_resistant_vs_S2_FDR,
  
  stringsAsFactors = FALSE
)


trajectory_long <- rbind(
  trajectory_C1_sensitive,
  trajectory_C1_mixed,
  trajectory_C2,
  trajectory_C3
)


trajectory_long$pretty <- pretty_pathway(
  trajectory_long$pathway
)


pathway_order <- pretty_pathway(
  trajectory$pathway
)


trajectory_long$pretty <- factor(
  trajectory_long$pretty,
  levels = rev(
    pathway_order
  )
)


trajectory_long$stage <- factor(
  trajectory_long$stage,
  levels = c(
    "C1 sensitive",
    "C1 mixed",
    "C2 resistant",
    "C3 resistant"
  )
)


trajectory_long$label <- paste0(
  
  sprintf(
    "%.2f",
    trajectory_long$NES
  ),
  
  ifelse(
    trajectory_long$FDR < 0.05,
    "*",
    ""
  )
)


p_C <- ggplot(
  trajectory_long,
  aes(
    x = stage,
    y = pretty,
    fill = NES
  )
) +
  
  geom_tile(
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    size = 2.7
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    name = "NES"
  ) +
  
  labs(
    
    title =
      "C. Serial emergence of final resistance-associated programs",
    
    subtitle =
      "Each stage compared with passaged untreated S2 controls; * FDR < 0.05",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  publication_theme +
  
  theme(
    
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8
    )
  )


# ============================================================
# 9. MAIN PANEL D
# D1 FROZEN HYPOTHESIS COMPARISON
# ============================================================

d1_validation$pretty <- pretty_pathway(
  d1_validation$pathway
)


d1_validation$display_label <- paste0(
  d1_validation$pretty,
  " [D1 expected ",
  d1_validation$expected_direction,
  "]"
)


d1_validation$display_label <- factor(
  d1_validation$display_label,
  levels = rev(
    d1_validation$display_label
  )
)


d1_validation$significance <- ifelse(
  d1_validation$p.adjust_limma < 0.05,
  "FDR < 0.05",
  "Not significant"
)


d1_validation$direction_status <- ifelse(
  d1_validation$primary_direction_match,
  "D1 direction matched",
  "Opposite to D1"
)


p_D <- ggplot(
  d1_validation,
  aes(
    x = NES_limma,
    y = display_label,
    fill = direction_status,
    alpha = significance
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  
  geom_col(
    width = 0.68
  ) +
  
  scale_alpha_manual(
    values = c(
      "FDR < 0.05" = 1,
      "Not significant" = 0.45
    )
  ) +
  
  labs(
    
    title =
      "D. D1 frozen hypotheses tested in D8",
    
    subtitle =
      "Observed D8 NES; bracketed text shows the direction predicted before D8 analysis",
    
    x =
      "D8 Hallmark NES",
    
    y =
      NULL,
    
    fill =
      "Direction",
    
    alpha =
      "Significance"
  ) +
  
  publication_theme +
  
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(
      size = 8
    )
  )


# ============================================================
# 10. COMBINE MAIN FIGURE
# ============================================================

main_figure <- (
  
  p_A |
    p_B
  
) /
  
  p_C /
  
  p_D +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.45,
      1.15
    )
  ) &
  
  theme(
    plot.background =
      element_rect(
        fill = "white",
        color = NA
      )
  )


# ============================================================
# 11. SAVE MAIN FIGURE
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D8_GSE163854_main.pdf"
  ),
  plot = main_figure,
  width = 14,
  height = 17,
  units = "in",
  bg = "white"
)


ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D8_GSE163854_main.png"
  ),
  plot = main_figure,
  width = 14,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# 12. SAVE INDIVIDUAL MAIN PANELS
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "Figure_D8_A_primary_PCA.pdf"
  ),
  p_A,
  width = 6,
  height = 5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_B_robust_Hallmarks.pdf"
  ),
  p_B,
  width = 7,
  height = 7,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_C_serial_trajectory.pdf"
  ),
  p_C,
  width = 9,
  height = 8,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_D_D1_validation.pdf"
  ),
  p_D,
  width = 9,
  height = 6.5,
  bg = "white"
)


# ============================================================
# 13. SUPPLEMENTARY PANEL S1
# FULL 12-SAMPLE PCA
# ============================================================

all_pca <- qc$all_pca$coordinates

all_variance <- qc$all_pca$variance


p_S1 <- ggplot(
  all_pca,
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
    size = 2.8,
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "S1. Full serial-selection PCA",
    
    subtitle =
      "Sensitive, mixed and resistant stages across the PH039 selection experiment",
    
    x = paste0(
      "PC1 (",
      round(
        all_variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        all_variance[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "Phenotype"
  ) +
  
  publication_theme


# ============================================================
# 14. SUPPLEMENTARY PANEL S2
# SAMPLE CORRELATION HEATMAP
# ============================================================

sample_cor <- qc$sample_correlations


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


p_S2 <- ggplot(
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
      "S2. Sample correlation matrix",
    x = NULL,
    y = NULL
  ) +
  
  publication_theme +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    
    axis.text.y = element_text(
      size = 7
    )
  )


# ============================================================
# 15. SUPPLEMENTARY PANEL S3
# GENE-LEVEL METHOD VALIDATION
# ============================================================

gene_rho <- cor(
  gene_validation$logFC_limma,
  gene_validation$logFC_edgeR,
  method = "spearman"
)


gene_nonzero <- gene_validation[
  gene_validation$logFC_limma != 0 &
    gene_validation$logFC_edgeR != 0,
  ,
  drop = FALSE
]


gene_sign_agreement <- mean(
  sign(
    gene_nonzero$logFC_limma
  ) ==
    sign(
      gene_nonzero$logFC_edgeR
    )
)


p_S3 <- ggplot(
  gene_validation,
  aes(
    x = logFC_edgeR,
    y = logFC_limma
  )
) +
  
  geom_point(
    alpha = 0.14,
    size = 0.65
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  
  coord_equal() +
  
  labs(
    
    title =
      "S3. Gene-level method validation",
    
    subtitle = paste0(
      "Spearman rho = ",
      round(
        gene_rho,
        3
      ),
      "; direction agreement = ",
      round(
        100 *
          gene_sign_agreement,
        1
      ),
      "%"
    ),
    
    x =
      "edgeR-QL log2FC",
    
    y =
      "limma-voom log2FC"
  ) +
  
  publication_theme


# ============================================================
# 16. SUPPLEMENTARY PANEL S4
# HALLMARK METHOD VALIDATION
# ============================================================

hallmark_complete <- hallmark_validation[
  is.finite(
    hallmark_validation$NES_limma
  ) &
    is.finite(
      hallmark_validation$NES_edgeR
    ),
  ,
  drop = FALSE
]


hallmark_rho <- cor(
  hallmark_complete$NES_limma,
  hallmark_complete$NES_edgeR,
  method = "spearman"
)


hallmark_direction_agreement <- mean(
  sign(
    hallmark_complete$NES_limma
  ) ==
    sign(
      hallmark_complete$NES_edgeR
    )
)


hallmark_complete$pretty <- pretty_pathway(
  hallmark_complete$ID
)


hallmark_labels <- hallmark_complete[
  !is.na(
    hallmark_complete$robust_significant
  ) &
    hallmark_complete$robust_significant,
  ,
  drop = FALSE
]


p_S4 <- ggplot(
  hallmark_complete,
  aes(
    x = NES_edgeR,
    y = NES_limma
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.3
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 2.3,
    alpha = 0.75
  ) +
  
  ggrepel::geom_text_repel(
    
    data =
      hallmark_labels,
    
    aes(
      label = pretty
    ),
    
    size = 2.3,
    
    max.overlaps = Inf,
    
    box.padding = 0.3,
    
    point.padding = 0.25
  ) +
  
  coord_equal() +
  
  labs(
    
    title =
      "S4. Hallmark method validation",
    
    subtitle = paste0(
      "NES Spearman rho = ",
      round(
        hallmark_rho,
        3
      ),
      "; ",
      round(
        100 *
          hallmark_direction_agreement,
        1
      ),
      "% direction agreement"
    ),
    
    x =
      "edgeR-QL NES",
    
    y =
      "limma-voom NES"
  ) +
  
  publication_theme


# ============================================================
# 17. COMBINE SUPPLEMENTARY FIGURE
# ============================================================

supplementary_figure <- (
  
  p_S1 |
    p_S2
  
) /
  
  (
    
    p_S3 |
      p_S4
    
  ) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1
    )
  ) &
  
  theme(
    plot.background =
      element_rect(
        fill = "white",
        color = NA
      )
  )


# ============================================================
# 18. SAVE SUPPLEMENTARY FIGURE
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D8_GSE163854_supplementary.pdf"
  ),
  plot = supplementary_figure,
  width = 14,
  height = 12,
  units = "in",
  bg = "white"
)


ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D8_GSE163854_supplementary.png"
  ),
  plot = supplementary_figure,
  width = 14,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# 19. SAVE INDIVIDUAL SUPPLEMENTARY PANELS
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "Figure_D8_S1_full_serial_PCA.pdf"
  ),
  p_S1,
  width = 7,
  height = 6,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_S2_sample_correlations.pdf"
  ),
  p_S2,
  width = 7,
  height = 6,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_S3_gene_method_validation.pdf"
  ),
  p_S3,
  width = 6,
  height = 5.5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D8_S4_Hallmark_method_validation.pdf"
  ),
  p_S4,
  width = 7,
  height = 6,
  bg = "white"
)


# ============================================================
# 20. SAVE FIGURE SOURCE TABLES
# ============================================================

write.csv(
  robust_long,
  file.path(
    figure_dir,
    "D8_Figure_B_robust_Hallmark_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  trajectory_long,
  file.path(
    figure_dir,
    "D8_Figure_C_trajectory_data.csv"
  ),
  row.names = FALSE
)


write.csv(
  d1_validation,
  file.path(
    figure_dir,
    "D8_Figure_D_D1_validation_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. FIGURE SUMMARY
# ============================================================

figure_summary <- data.frame(
  
  metric = c(
    
    "Primary_PC1_percent",
    
    "Primary_PC2_percent",
    
    "Robust_final_Hallmarks",
    
    "Trajectory_final_Hallmarks",
    
    "Gene_method_Spearman_rho",
    
    "Gene_direction_agreement",
    
    "Hallmark_method_Spearman_rho",
    
    "Hallmark_direction_agreement",
    
    "D1_directional_matches",
    
    "D1_robust_significant_matches"
  ),
  
  value = c(
    
    primary_variance[1],
    
    primary_variance[2],
    
    nrow(
      robust
    ),
    
    nrow(
      trajectory
    ),
    
    gene_rho,
    
    gene_sign_agreement,
    
    hallmark_rho,
    
    hallmark_direction_agreement,
    
    sum(
      d1_validation$primary_direction_match,
      na.rm = TRUE
    ),
    
    sum(
      d1_validation$robust_significant_match,
      na.rm = TRUE
    )
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  figure_summary,
  file.path(
    figure_dir,
    "D8_figure_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. SESSION INFO
# ============================================================

sink(
  file.path(
    figure_dir,
    "D8_figures_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 23. FINAL OUTPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D8 FIGURE GENERATION COMPLETE\n")
cat("============================================================\n\n")


cat(
  "MAIN FIGURE:\n"
)

cat(
  file.path(
    figure_dir,
    "Figure_D8_GSE163854_main.pdf"
  ),
  "\n\n"
)


cat(
  "SUPPLEMENTARY FIGURE:\n"
)

cat(
  file.path(
    figure_dir,
    "Figure_D8_GSE163854_supplementary.pdf"
  ),
  "\n\n"
)


cat(
  "FIGURE SUMMARY:\n\n"
)


print(
  figure_summary,
  row.names = FALSE
)


cat("\n")
cat(
  "No DE, GSEA, filtering or pathway classification was rerun.\n"
)

cat(
  "Frozen D8 accession-level vote remains S2 vs C3 only.\n"
)

cat("============================================================\n")