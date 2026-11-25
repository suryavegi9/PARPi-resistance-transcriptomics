# ============================================================
# GSE165914 - D7
# PUBLICATION-QUALITY FIGURES
#
# Uses ONLY frozen outputs from:
#   Script 52 - QC / PCA
#   Script 53 - DE + Hallmark GSEA
#   Script 54 - independent author-DE validation
#
# NO filtering
# NO DE re-analysis
# NO GSEA re-analysis
# NO threshold changes
#
# Outputs:
#
# MAIN D7 FIGURE
#   A1 - MDA-MB-468 PCA
#   A2 - SUM1315 PCA
#   B  - Hallmark NES agreement scatter
#   C  - Concordant significant Hallmark heatmap
#   D  - Original D1 frozen pathways in D7
#
# SUPPLEMENTARY D7 FIGURE
#   S1 - sample correlation heatmap
#   S2 - MDA-MB-468 author-vs-our DE validation
#   S3 - SUM1315 author-vs-our DE validation
#
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

qc_dir <- "results/GSE165914/QC"

analysis_dir <-
  "results/GSE165914/D7_analysis"

validation_dir <-
  "results/GSE165914/D7_validation"

figure_dir <-
  "results/GSE165914/figures"


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


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


cat("\n")
cat("============================================================\n")
cat("GSE165914 - D7 FIGURE GENERATION\n")
cat("============================================================\n\n")


# ============================================================
# 2. REQUIRED FILES
# ============================================================

qc_file <- file.path(
  qc_dir,
  "D7_QC_objects.rds"
)


integration_file <- file.path(
  analysis_dir,
  "D7_accession_Hallmark_integration.csv"
)


frozen_file <- file.path(
  analysis_dir,
  "D7_vs_D1_frozen_validation.csv"
)


m468_validation_file <- file.path(
  validation_dir,
  "D7A_M468_author_vs_our_DE.csv"
)


sum_validation_file <- file.path(
  validation_dir,
  "D7B_SUM1315_author_vs_our_DE.csv"
)


required_files <- c(
  qc_file,
  integration_file,
  frozen_file,
  m468_validation_file,
  sum_validation_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  cat(
    "Missing required files:\n"
  )
  
  print(
    missing_files
  )
  
  stop(
    "Cannot generate D7 figures because required files are missing."
  )
}


cat(
  "All required D7 result files found.\n"
)


# ============================================================
# 3. LOAD FROZEN RESULTS
# ============================================================

qc <- readRDS(
  qc_file
)


hallmark <- read.csv(
  integration_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


frozen <- read.csv(
  frozen_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


m468_validation <- read.csv(
  m468_validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


sum_validation <- read.csv(
  sum_validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. PRETTY PATHWAY LABEL FUNCTION
#
# ASCII-only to avoid Windows encoding problems.
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
  
  
  replacements_from <- c(
    "Myc",
    "E2f",
    "Mtorc1",
    "G2m",
    "Tnfa",
    "Nfkb",
    "Ros",
    "Uv",
    "Il2",
    "Il6",
    "Jak",
    "Stat3",
    "Kras",
    "P53"
  )
  
  
  replacements_to <- c(
    "MYC",
    "E2F",
    "mTORC1",
    "G2M",
    "TNFa",
    "NF-kB",
    "ROS",
    "UV",
    "IL2",
    "IL6",
    "JAK",
    "STAT3",
    "KRAS",
    "P53"
  )
  
  
  for (i in seq_along(replacements_from)) {
    
    x <- gsub(
      replacements_from[i],
      replacements_to[i],
      x,
      fixed = TRUE
    )
    
  }
  
  
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
# 5. GLOBAL FIGURE THEME
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
    
    legend.position = "right",
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 6. PANEL A1 - MDA-MB-468 PCA
# ============================================================

m468_pca <- qc$m468_pca$coordinates

m468_variance <- qc$m468_pca$variance


m468_pca$state <- factor(
  m468_pca$state,
  levels = c(
    "Parental",
    "OlaR"
  )
)


p_A1 <- ggplot(
  m468_pca,
  aes(
    x = PC1,
    y = PC2,
    shape = state,
    label = sample
  )
) +
  
  geom_point(
    size = 3.5
  ) +
  
  ggrepel::geom_text_repel(
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  
  labs(
    
    title =
      "A1. MDA-MB-468",
    
    subtitle =
      "Parental vs acquired olaparib-resistant",
    
    x = paste0(
      "PC1 (",
      round(
        m468_variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        m468_variance[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "State"
  ) +
  
  publication_theme


# ============================================================
# 7. PANEL A2 - SUM1315 PCA
# ============================================================

sum_pca <- qc$sum_pca$coordinates

sum_variance <- qc$sum_pca$variance


sum_pca$state <- factor(
  sum_pca$state,
  levels = c(
    "Parental",
    "OlaR"
  )
)


p_A2 <- ggplot(
  sum_pca,
  aes(
    x = PC1,
    y = PC2,
    shape = state,
    label = sample
  )
) +
  
  geom_point(
    size = 3.5
  ) +
  
  ggrepel::geom_text_repel(
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  
  labs(
    
    title =
      "A2. SUM1315",
    
    subtitle =
      "Parental vs acquired olaparib-resistant",
    
    x = paste0(
      "PC1 (",
      round(
        sum_variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        sum_variance[2],
        1
      ),
      "%)"
    ),
    
    shape =
      "State"
  ) +
  
  publication_theme


# ============================================================
# 8. PANEL B - HALLMARK NES AGREEMENT
# ============================================================

complete_hallmark <- hallmark[
  is.finite(
    hallmark$M468_NES
  ) &
    is.finite(
      hallmark$SUM1315_NES
    ),
  ,
  drop = FALSE
]


hallmark_rho <- cor(
  complete_hallmark$M468_NES,
  complete_hallmark$SUM1315_NES,
  method = "spearman"
)


same_direction_n <- sum(
  sign(
    complete_hallmark$M468_NES
  ) ==
    sign(
      complete_hallmark$SUM1315_NES
    )
)


complete_hallmark$category <- "Other"


concordant_sig <- (
  !is.na(
    complete_hallmark$concordant_significant
  ) &
    complete_hallmark$concordant_significant
)


opposite_both_sig <- (
  !is.na(
    complete_hallmark$same_direction
  ) &
    !complete_hallmark$same_direction &
    !is.na(
      complete_hallmark$M468_FDR
    ) &
    !is.na(
      complete_hallmark$SUM1315_FDR
    ) &
    complete_hallmark$M468_FDR < 0.05 &
    complete_hallmark$SUM1315_FDR < 0.05
)


complete_hallmark$category[
  concordant_sig
] <- "Concordant significant"


complete_hallmark$category[
  opposite_both_sig
] <- "Opposite significant"


complete_hallmark$category <- factor(
  complete_hallmark$category,
  levels = c(
    "Other",
    "Concordant significant",
    "Opposite significant"
  )
)


label_data <- complete_hallmark[
  concordant_sig |
    opposite_both_sig,
  ,
  drop = FALSE
]


label_data$label <- pretty_pathway(
  label_data$pathway
)


p_B <- ggplot(
  complete_hallmark,
  aes(
    x = M468_NES,
    y = SUM1315_NES,
    shape = category
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  
  geom_point(
    size = 2.8,
    alpha = 0.80
  ) +
  
  ggrepel::geom_text_repel(
    
    data = label_data,
    
    aes(
      label = label
    ),
    
    size = 2.6,
    
    max.overlaps = Inf,
    
    box.padding = 0.35,
    
    point.padding = 0.25,
    
    show.legend = FALSE
  ) +
  
  labs(
    
    title =
      "B. Hallmark pathway agreement",
    
    subtitle = paste0(
      "Spearman rho = ",
      round(
        hallmark_rho,
        3
      ),
      "; ",
      same_direction_n,
      "/",
      nrow(
        complete_hallmark
      ),
      " pathways have the same direction"
    ),
    
    x =
      "MDA-MB-468 Hallmark NES",
    
    y =
      "SUM1315 Hallmark NES",
    
    shape =
      "Pathway class"
  ) +
  
  publication_theme


# ============================================================
# 9. SAVE PANEL-B SOURCE DATA
# ============================================================

write.csv(
  complete_hallmark,
  file.path(
    figure_dir,
    "D7_Figure_B_Hallmark_scatter_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. PANEL C - CONCORDANT SIGNIFICANT PATHWAYS
# ============================================================

concordant <- hallmark[
  !is.na(
    hallmark$concordant_significant
  ) &
    hallmark$concordant_significant,
  ,
  drop = FALSE
]


if (nrow(concordant) == 0) {
  
  stop(
    "No concordant significant Hallmark pathways found."
  )
}


concordant$pretty <- pretty_pathway(
  concordant$pathway
)


# Order by average D7 NES from highest to lowest.

concordant <- concordant[
  order(
    concordant$D7_NES,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


path_order <- concordant$pretty


concordant_long <- rbind(
  
  data.frame(
    
    pathway =
      concordant$pretty,
    
    model =
      "MDA-MB-468",
    
    NES =
      concordant$M468_NES,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    
    pathway =
      concordant$pretty,
    
    model =
      "SUM1315",
    
    NES =
      concordant$SUM1315_NES,
    
    stringsAsFactors = FALSE
  )
  
)


concordant_long$pathway <- factor(
  concordant_long$pathway,
  levels = rev(
    path_order
  )
)


concordant_long$model <- factor(
  concordant_long$model,
  levels = c(
    "MDA-MB-468",
    "SUM1315"
  )
)


p_C <- ggplot(
  concordant_long,
  aes(
    x = model,
    y = pathway,
    fill = NES
  )
) +
  
  geom_tile(
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        NES
      )
    ),
    size = 3.2
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    name = "NES"
  ) +
  
  labs(
    
    title =
      "C. Concordant significant programs",
    
    subtitle =
      "FDR < 0.05 in both acquired-resistant models",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  publication_theme +
  
  theme(
    axis.text.x = element_text(
      face = "bold"
    )
  )


write.csv(
  concordant_long,
  file.path(
    figure_dir,
    "D7_Figure_C_concordant_Hallmark_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. PANEL D - ORIGINAL D1 FROZEN PATHWAYS
# ============================================================

frozen$pretty <- pretty_pathway(
  frozen$pathway
)


frozen$pathway_label <- paste0(
  frozen$pretty,
  " [D1 expected ",
  frozen$D1_expected_direction,
  "]"
)


# Preserve the original frozen hypothesis order.

frozen$pathway_label <- factor(
  frozen$pathway_label,
  levels = rev(
    frozen$pathway_label
  )
)


frozen_long <- rbind(
  
  data.frame(
    
    pathway =
      frozen$pathway,
    
    pathway_label =
      frozen$pathway_label,
    
    model =
      "MDA-MB-468",
    
    NES =
      frozen$M468_NES,
    
    FDR =
      frozen$M468_FDR,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    
    pathway =
      frozen$pathway,
    
    pathway_label =
      frozen$pathway_label,
    
    model =
      "SUM1315",
    
    NES =
      frozen$SUM1315_NES,
    
    FDR =
      frozen$SUM1315_FDR,
    
    stringsAsFactors = FALSE
  )
  
)


frozen_long$pathway_label <- factor(
  frozen_long$pathway_label,
  levels = levels(
    frozen$pathway_label
  )
)


frozen_long$model <- factor(
  frozen_long$model,
  levels = c(
    "MDA-MB-468",
    "SUM1315"
  )
)


frozen_long$significance <- ifelse(
  frozen_long$FDR < 0.05,
  "FDR < 0.05",
  "Not significant"
)


p_D <- ggplot(
  frozen_long,
  aes(
    x = pathway_label,
    y = NES,
    fill = model,
    alpha = significance
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  
  geom_col(
    position = position_dodge(
      width = 0.78
    ),
    width = 0.68
  ) +
  
  coord_flip() +
  
  scale_alpha_manual(
    values = c(
      "FDR < 0.05" = 1,
      "Not significant" = 0.45
    )
  ) +
  
  labs(
    
    title =
      "D. Original D1 frozen pathways in D7",
    
    subtitle =
      paste0(
        "Bar direction = D7 NES; ",
        "bracketed text = direction predicted from D1"
      ),
    
    x =
      NULL,
    
    y =
      "Normalized enrichment score",
    
    fill =
      "D7 model",
    
    alpha =
      "Significance"
  ) +
  
  publication_theme +
  
  theme(
    legend.position = "bottom"
  )


write.csv(
  frozen_long,
  file.path(
    figure_dir,
    "D7_Figure_D_frozen_pathway_data.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. MAIN D7 FIGURE
# ============================================================

main_figure <- (
  
  p_A1 |
    p_A2
  
) /
  
  (
    
    p_B |
      p_C
    
  ) /
  
  p_D +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1.20,
      1.45
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
# 13. SAVE MAIN FIGURE
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D7_GSE165914_main.pdf"
  ),
  plot = main_figure,
  width = 14,
  height = 15,
  units = "in",
  bg = "white"
)


ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D7_GSE165914_main.png"
  ),
  plot = main_figure,
  width = 14,
  height = 15,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# 14. SAVE INDIVIDUAL MAIN PANELS
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "Figure_D7_A1_M468_PCA.pdf"
  ),
  p_A1,
  width = 6,
  height = 5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_A2_SUM1315_PCA.pdf"
  ),
  p_A2,
  width = 6,
  height = 5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_B_Hallmark_agreement.pdf"
  ),
  p_B,
  width = 7,
  height = 6,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_C_concordant_Hallmarks.pdf"
  ),
  p_C,
  width = 6,
  height = 5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_D_D1_frozen_pathways.pdf"
  ),
  p_D,
  width = 10,
  height = 7,
  bg = "white"
)


# ============================================================
# 15. SUPPLEMENTARY PANEL S1
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


p_S1 <- ggplot(
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
      "S1. D7 sample correlation matrix",
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
# 16. SUPPLEMENTARY VALIDATION FUNCTION
# ============================================================

make_validation_plot <- function(
    dat,
    panel_title
) {
  
  
  rho <- cor(
    dat$our_logFC,
    dat$author_logFC_aligned,
    method = "spearman"
  )
  
  
  nonzero <- dat[
    dat$our_logFC != 0 &
      dat$author_logFC_aligned != 0,
    ,
    drop = FALSE
  ]
  
  
  sign_agreement <- mean(
    sign(
      nonzero$our_logFC
    ) ==
      sign(
        nonzero$author_logFC_aligned
      )
  )
  
  
  p <- ggplot(
    dat,
    aes(
      x = author_logFC_aligned,
      y = our_logFC
    )
  ) +
    
    geom_point(
      alpha = 0.18,
      size = 0.65
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    
    labs(
      
      title =
        panel_title,
      
      subtitle = paste0(
        "Spearman rho = ",
        round(
          rho,
          3
        ),
        "; direction agreement = ",
        round(
          100 *
            sign_agreement,
          1
        ),
        "%"
      ),
      
      x =
        "Author edgeR log2FC",
      
      y =
        "Independent limma logFC"
    ) +
    
    publication_theme
  
  
  return(p)
}


# ============================================================
# 17. SUPPLEMENTARY S2 / S3
# ============================================================

p_S2 <- make_validation_plot(
  m468_validation,
  "S2. MDA-MB-468 independent DE validation"
)


p_S3 <- make_validation_plot(
  sum_validation,
  "S3. SUM1315 independent DE validation"
)


# ============================================================
# 18. COMBINE SUPPLEMENTARY FIGURE
# ============================================================

supplementary_figure <- (
  
  p_S1 |
    p_S2 |
    p_S3
  
) &
  
  theme(
    plot.background =
      element_rect(
        fill = "white",
        color = NA
      )
  )


# ============================================================
# 19. SAVE SUPPLEMENTARY FIGURE
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D7_supplementary_QC_validation.pdf"
  ),
  plot = supplementary_figure,
  width = 18,
  height = 6,
  units = "in",
  bg = "white"
)


ggsave(
  filename = file.path(
    figure_dir,
    "Figure_D7_supplementary_QC_validation.png"
  ),
  plot = supplementary_figure,
  width = 18,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# 20. SAVE INDIVIDUAL SUPPLEMENTARY PANELS
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "Figure_D7_S1_sample_correlations.pdf"
  ),
  p_S1,
  width = 8,
  height = 7,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_S2_M468_author_validation.pdf"
  ),
  p_S2,
  width = 6,
  height = 5,
  bg = "white"
)


ggsave(
  file.path(
    figure_dir,
    "Figure_D7_S3_SUM1315_author_validation.pdf"
  ),
  p_S3,
  width = 6,
  height = 5,
  bg = "white"
)


# ============================================================
# 21. FIGURE SUMMARY
# ============================================================

figure_summary <- data.frame(
  
  metric = c(
    "M468_PC1_percent",
    "M468_PC2_percent",
    "SUM1315_PC1_percent",
    "SUM1315_PC2_percent",
    "Hallmark_NES_Spearman_rho",
    "Hallmark_same_direction_n",
    "Hallmark_total_n",
    "Concordant_significant_Hallmarks",
    "Opposite_significant_Hallmarks"
  ),
  
  value = c(
    m468_variance[1],
    m468_variance[2],
    sum_variance[1],
    sum_variance[2],
    hallmark_rho,
    same_direction_n,
    nrow(
      complete_hallmark
    ),
    sum(
      concordant_sig
    ),
    sum(
      opposite_both_sig
    )
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  figure_summary,
  file.path(
    figure_dir,
    "D7_figure_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. SESSION INFO
# ============================================================

sink(
  file.path(
    figure_dir,
    "D7_figures_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 23. FINAL OUTPUT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("D7 FIGURE GENERATION COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Main figure:\n",
  file.path(
    figure_dir,
    "Figure_D7_GSE165914_main.pdf"
  ),
  "\n\n"
)


cat(
  "Supplementary figure:\n",
  file.path(
    figure_dir,
    "Figure_D7_supplementary_QC_validation.pdf"
  ),
  "\n\n"
)


cat(
  "Main-panel statistics:\n\n"
)


print(
  figure_summary,
  row.names = FALSE
)


cat("\n")
cat(
  "No biological analysis was rerun or modified.\n"
)

cat("============================================================\n")