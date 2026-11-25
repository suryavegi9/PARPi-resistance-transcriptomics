# ============================================================
# GSE191231 - DATASET 6
# PUBLICATION FIGURES
#
# Figure 7A - PCA
# Figure 7B - Primary vs sensitivity Hallmark concordance
# Figure 7C - Robust Hallmark programs
# Figure 7D - Frozen D1 hypothesis validation
#
# Main message:
# D6 shows a robust but mixed off-treatment resistance state,
# supporting selective pathway recurrence plus substantial
# cross-model heterogeneity.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)


# ============================================================
# 0. DIRECTORIES
# ============================================================

figdir <- "figures/GSE191231_D6"

datadir <- file.path(
  figdir,
  "figure_data"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  datadir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. COMMON HELPERS
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
# FIGURE 7A
# PCA OF D6
# ============================================================
# ============================================================


tpm_file <- paste0(
  "results/GSE191231/",
  "GSE191231_filtered_TPM_matrix.csv"
)


tpm_dat <- read.csv(
  tpm_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


ensembl <- tpm_dat$ensembl_id

tpm <- as.matrix(
  tpm_dat[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(tpm) <- "numeric"


log_tpm <- log2(
  tpm + 1
)


pca <- prcomp(
  t(log_tpm),
  center = TRUE,
  scale. = FALSE
)


variance <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )


pca_df <- data.frame(
  
  sample = rownames(
    pca$x
  ),
  
  group = c(
    "Vehicle",
    "Vehicle",
    "Olaparib-selected",
    "Olaparib-selected"
  ),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  stringsAsFactors = FALSE
)


pca_df$label <- c(
  "834",
  "846",
  "827",
  "878"
)


write.csv(
  pca_df,
  file.path(
    datadir,
    "Figure7A_PCA_data.csv"
  ),
  row.names = FALSE
)


p7a <- ggplot(
  
  pca_df,
  
  aes(
    x = PC1,
    y = PC2,
    colour = group
  )
) +
  
  geom_point(
    size = 5
  ) +
  
  geom_text_repel(
    
    aes(
      label = label
    ),
    
    size = 4,
    
    box.padding = 0.5,
    
    point.padding = 0.3,
    
    min.segment.length = 0,
    
    max.overlaps = Inf,
    
    seed = 123,
    
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    
    values = c(
      
      "Vehicle" =
        "#D55E00",
      
      "Olaparib-selected" =
        "#0072B2"
    )
  ) +
  
  labs(
    
    title =
      "A  Transcriptomic separation of olaparib-selected PDX tumors",
    
    subtitle =
      paste0(
        "PC1 = ",
        round(
          variance[1],
          1
        ),
        "%; PC2 = ",
        round(
          variance[2],
          1
        ),
        "%"
      ),
    
    x = paste0(
      "PC1 (",
      round(
        variance[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        variance[2],
        1
      ),
      "%)"
    ),
    
    colour = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),
    
    legend.position =
      "bottom"
  )


ggsave(
  file.path(
    figdir,
    "Figure7A_D6_PCA.png"
  ),
  p7a,
  width = 7.5,
  height = 6,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure7A_D6_PCA.pdf"
  ),
  p7a,
  width = 7.5,
  height = 6
)


# ============================================================
# ============================================================
# FIGURE 7B
# PRIMARY vs SENSITIVITY HALLMARK CONCORDANCE
# ============================================================
# ============================================================


gsea_compare <- read.csv(
  
  paste0(
    "results/GSE191231/Hallmark_GSEA/",
    "GSE191231_Hallmark_GSEA_primary_vs_sensitivity.csv"
  ),
  
  stringsAsFactors = FALSE
)


gsea_compare$label <- clean_hallmark(
  gsea_compare$pathway
)


gsea_compare$category <- ifelse(
  
  gsea_compare$robust_significant,
  
  "Robust significant",
  
  "Other"
)


label_pathways <- c(
  
  "HALLMARK_HYPOXIA",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_GLYCOLYSIS",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)


label_df <- gsea_compare[
  gsea_compare$pathway %in%
    label_pathways,
]


write.csv(
  gsea_compare,
  file.path(
    datadir,
    "Figure7B_GSEA_concordance_data.csv"
  ),
  row.names = FALSE
)


p7b <- ggplot(
  
  gsea_compare,
  
  aes(
    x = primary_NES,
    y = sensitivity_NES
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  
  geom_point(
    
    aes(
      colour = category
    ),
    
    size = 3.4
  ) +
  
  geom_text_repel(
    
    data = label_df,
    
    aes(
      label = label
    ),
    
    size = 3,
    
    box.padding = 0.45,
    
    point.padding = 0.3,
    
    min.segment.length = 0,
    
    max.overlaps = Inf,
    
    seed = 123,
    
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    
    values = c(
      
      "Other" =
        "grey70",
      
      "Robust significant" =
        "#009E73"
    )
  ) +
  
  labs(
    
    title =
      "B  Hallmark enrichment is robust to analysis strategy",
    
    subtitle =
      paste0(
        "Spearman rho = 0.879; ",
        "21/50 Hallmarks significant and directionally consistent"
      ),
    
    x =
      "Primary TPM / limma-trend NES",
    
    y =
      "RSEM / voom sensitivity NES",
    
    colour = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),
    
    legend.position =
      "bottom",
    
    plot.margin =
      margin(
        15,
        25,
        15,
        25
      )
  )


ggsave(
  file.path(
    figdir,
    "Figure7B_D6_Hallmark_robustness.png"
  ),
  p7b,
  width = 8,
  height = 7,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure7B_D6_Hallmark_robustness.pdf"
  ),
  p7b,
  width = 8,
  height = 7
)


# ============================================================
# ============================================================
# FIGURE 7C
# ROBUST SIGNIFICANT HALLMARKS
# ============================================================
# ============================================================


robust <- read.csv(
  
  paste0(
    "results/GSE191231/Hallmark_GSEA/",
    "GSE191231_Hallmark_GSEA_robust_significant.csv"
  ),
  
  stringsAsFactors = FALSE
)


robust$pathway_label <- clean_hallmark(
  robust$pathway
)


robust$mean_NES <- rowMeans(
  robust[
    ,
    c(
      "primary_NES",
      "sensitivity_NES"
    )
  ]
)


robust <- robust[
  order(
    robust$mean_NES
  ),
]


robust$pathway_label <- factor(
  
  robust$pathway_label,
  
  levels =
    robust$pathway_label
)


write.csv(
  robust,
  file.path(
    datadir,
    "Figure7C_robust_Hallmarks.csv"
  ),
  row.names = FALSE
)


p7c <- ggplot(
  robust
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_segment(
    
    aes(
      x = primary_NES,
      xend = sensitivity_NES,
      y = pathway_label,
      yend = pathway_label
    ),
    
    linewidth = 0.7,
    colour = "grey60"
  ) +
  
  geom_point(
    
    aes(
      x = primary_NES,
      y = pathway_label
    ),
    
    size = 3.3,
    shape = 16,
    colour = "#0072B2"
  ) +
  
  geom_point(
    
    aes(
      x = sensitivity_NES,
      y = pathway_label
    ),
    
    size = 3.5,
    shape = 17,
    colour = "#D55E00"
  ) +
  
  labs(
    
    title =
      "C  Robust transcriptional programs in persistent resistance",
    
    subtitle =
      paste0(
        "Positive NES = higher in olaparib-selected tumors; ",
        "negative NES = higher in vehicle controls"
      ),
    
    x =
      "Normalized enrichment score (NES)",
    
    y = NULL,
    
    caption =
      "Circle = TPM/limma primary; triangle = RSEM/voom sensitivity"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    
    axis.text.y =
      element_text(
        size = 8.5
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9
      )
  )


ggsave(
  file.path(
    figdir,
    "Figure7C_D6_robust_Hallmarks.png"
  ),
  p7c,
  width = 9.5,
  height = 8.5,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure7C_D6_robust_Hallmarks.pdf"
  ),
  p7c,
  width = 9.5,
  height = 8.5
)


# ============================================================
# ============================================================
# FIGURE 7D
# FROZEN D1 HYPOTHESIS VALIDATION
# ============================================================
# ============================================================


frozen <- read.csv(
  
  paste0(
    "results/GSE191231/Hallmark_GSEA/",
    "GSE191231_D1_frozen_validation.csv"
  ),
  
  stringsAsFactors = FALSE
)


frozen$pathway_label <- clean_hallmark(
  frozen$pathway
)


# Add expected direction to label

frozen$expected_symbol <- ifelse(
  frozen$expected_direction == "UP",
  "\u2191",
  "\u2193"
)


frozen$display_label <- paste0(
  frozen$pathway_label,
  "  [expected ",
  frozen$expected_symbol,
  "]"
)


# ============================================================
# Long-format heatmap
# ============================================================

frozen_long <- bind_rows(
  
  data.frame(
    
    pathway =
      frozen$pathway,
    
    display_label =
      frozen$display_label,
    
    analysis =
      "Primary",
    
    NES =
      frozen$primary_NES,
    
    padj =
      frozen$primary_padj,
    
    expected_direction =
      frozen$expected_direction,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    
    pathway =
      frozen$pathway,
    
    display_label =
      frozen$display_label,
    
    analysis =
      "Sensitivity",
    
    NES =
      frozen$sensitivity_NES,
    
    padj =
      frozen$sensitivity_padj,
    
    expected_direction =
      frozen$expected_direction,
    
    stringsAsFactors = FALSE
  )
)


frozen_long$significant <-
  frozen_long$padj < 0.05


order_paths <- c(
  
  "OXIDATIVE PHOSPHORYLATION",
  
  "EPITHELIAL MESENCHYMAL TRANSITION",
  
  "G2M CHECKPOINT",
  
  "E2F TARGETS",
  
  "MTORC1 SIGNALING",
  
  "MYC TARGETS V1",
  
  "GLYCOLYSIS",
  
  "INFLAMMATORY RESPONSE",
  
  "TNFA SIGNALING VIA NFKB",
  
  "UV RESPONSE DN"
)


label_order <- frozen$display_label[
  match(
    order_paths,
    frozen$pathway_label
  )
]


frozen_long$display_label <- factor(
  
  frozen_long$display_label,
  
  levels = rev(
    label_order
  )
)


frozen_long$analysis <- factor(
  
  frozen_long$analysis,
  
  levels = c(
    "Primary",
    "Sensitivity"
  )
)


max_abs <- max(
  abs(
    frozen_long$NES
  ),
  na.rm = TRUE
)


write.csv(
  frozen_long,
  file.path(
    datadir,
    "Figure7D_frozen_D1_validation.csv"
  ),
  row.names = FALSE
)


p7d <- ggplot(
  
  frozen_long,
  
  aes(
    x = analysis,
    y = display_label,
    fill = NES
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  
  geom_text(
    
    aes(
      label = ifelse(
        significant,
        "*",
        ""
      )
    ),
    
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    
    low = "#2166AC",
    
    mid = "white",
    
    high = "#B2182B",
    
    midpoint = 0,
    
    limits = c(
      -max_abs,
      max_abs
    ),
    
    name = "NES"
  ) +
  
  labs(
    
    title =
      "D  D6 provides selective rather than complete validation of D1",
    
    subtitle =
      paste0(
        "3/10 frozen pathways match the expected direction ",
        "and remain significant in both analyses"
      ),
    
    x = NULL,
    
    y = NULL,
    
    caption =
      paste0(
        "* BH-adjusted P < 0.05. ",
        "Arrows indicate the pre-specified direction expected from D1."
      )
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 10
      ),
    
    axis.text.y =
      element_text(
        size = 8.5
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9
      )
  )


ggsave(
  file.path(
    figdir,
    "Figure7D_D6_frozen_D1_validation.png"
  ),
  p7d,
  width = 8,
  height = 6.8,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure7D_D6_frozen_D1_validation.pdf"
  ),
  p7d,
  width = 8,
  height = 6.8
)


# ============================================================
# ============================================================
# COMBINED FIGURE 7
# ============================================================
# ============================================================

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {
  
  library(patchwork)
  
  
  combined <- (
    
    p7a |
      p7b
    
  ) /
    
    (
      
      p7c |
        p7d
    ) +
    
    plot_layout(
      heights = c(
        0.9,
        1.1
      )
    ) +
    
    plot_annotation(
      
      title =
        paste0(
          "Figure 7. GSE191231 reveals a robust but ",
          "context-specific off-treatment resistance state"
        ),
      
      theme = theme(
        
        plot.title =
          element_text(
            face = "bold",
            size = 18,
            hjust = 0
          )
      )
    )
  
  
  ggsave(
    
    file.path(
      figdir,
      "Figure7_D6_combined.png"
    ),
    
    combined,
    
    width = 18,
    
    height = 14,
    
    dpi = 400
  )
  
  
  ggsave(
    
    file.path(
      figdir,
      "Figure7_D6_combined.pdf"
    ),
    
    combined,
    
    width = 18,
    
    height = 14
  )
}


# ============================================================
# FINAL
# ============================================================

cat("\n========================================\n")
cat("D6 FIGURES COMPLETE\n")
cat("========================================\n\n")

cat("Figure 7A - PCA\n")
cat("Figure 7B - Hallmark robustness\n")
cat("Figure 7C - robust Hallmark programs\n")
cat("Figure 7D - frozen D1 validation\n")

cat(
  "\nOutput directory:",
  figdir,
  "\n"
)