# ============================================================
# GSE153867 - PUBLICATION FIGURES
# Dataset 3: A2780 parental vs olaparib-resistant
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(ggplot2)
library(ggrepel)

dir.create(
  "figures/GSE153867",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# FIGURE 3A - PCA
# ============================================================

log_expr <- read.csv(
  "results/GSE153867/GSE153867_log2FPKMplus1.csv",
  row.names = 1,
  check.names = FALSE
)

log_expr <- as.matrix(log_expr)


metadata <- data.frame(
  
  sample = colnames(log_expr),
  
  phenotype = ifelse(
    grepl("^C-", colnames(log_expr)),
    "Resistant",
    "Parental"
  ),
  
  row.names = colnames(log_expr),
  
  stringsAsFactors = FALSE
)


pca <- prcomp(
  t(log_expr),
  center = TRUE,
  scale. = FALSE
)


variance <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100


pca_df <- data.frame(
  
  sample = rownames(pca$x),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  phenotype = metadata[
    rownames(pca$x),
    "phenotype"
  ]
)


p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = phenotype
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text_repel(
    aes(label = sample),
    size = 3.5,
    show.legend = FALSE,
    box.padding = 0.35,
    point.padding = 0.3
  ) +
  labs(
    title = "Transcriptomic separation of parental and olaparib-resistant A2780 cells",
    subtitle = "Eight parental and eight resistant samples from GSE153867",
    x = paste0(
      "PC1 (",
      round(variance[1], 1),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(variance[2], 1),
      "% variance)"
    ),
    color = "Phenotype"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    plot.margin = margin(
      15, 20, 15, 15
    )
  )


print(p_pca)


ggsave(
  "figures/GSE153867/Figure_3A_PCA.png",
  p_pca,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE153867/Figure_3A_PCA.pdf",
  p_pca,
  width = 9,
  height = 6
)


# ============================================================
# FIGURE 3B - PRIMARY SIGNIFICANT HALLMARKS
# ============================================================

gsea <- read.csv(
  "results/GSE153867/GSE153867_Hallmark_GSEA_all.csv"
)


gsea_sig <- subset(
  gsea,
  p.adjust < 0.05
)


gsea_sig$Pathway <- gsub(
  "^HALLMARK_",
  "",
  gsea_sig$ID
)

gsea_sig$Pathway <- gsub(
  "_",
  " ",
  gsea_sig$Pathway
)


gsea_sig$Direction <- ifelse(
  gsea_sig$NES > 0,
  "Enriched in resistant",
  "Depleted in resistant"
)


gsea_sig <- gsea_sig[
  order(gsea_sig$NES),
]


gsea_sig$Pathway <- factor(
  gsea_sig$Pathway,
  levels = gsea_sig$Pathway
)


p_hallmark <- ggplot(
  gsea_sig,
  aes(
    x = Pathway,
    y = NES,
    fill = Direction
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "Significant Hallmark programs in GSE153867",
    subtitle = "Primary analysis; BH-adjusted P < 0.05",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      15, 25, 15, 15
    )
  )


print(p_hallmark)


ggsave(
  "figures/GSE153867/Figure_3B_significant_Hallmarks.png",
  p_hallmark,
  width = 10,
  height = 7.5,
  dpi = 300
)

ggsave(
  "figures/GSE153867/Figure_3B_significant_Hallmarks.pdf",
  p_hallmark,
  width = 10,
  height = 7.5
)


# ============================================================
# FIGURE 3C - FROZEN DATASET-1 HYPOTHESIS VALIDATION
# ============================================================

validation <- read.csv(
  "results/GSE153867/GSE153867_frozen_hypothesis_validation.csv"
)


validation$Pathway <- gsub(
  "^HALLMARK_",
  "",
  validation$Hallmark
)

validation$Pathway <- gsub(
  "_",
  " ",
  validation$Pathway
)


validation$Result <- ifelse(
  
  validation$direction_concordant &
    validation$significant,
  
  "Concordant + significant",
  
  ifelse(
    
    validation$direction_concordant &
      !validation$significant,
    
    "Concordant direction, not significant",
    
    ifelse(
      
      !validation$direction_concordant &
        validation$significant,
      
      "Opposite + significant",
      
      "Opposite direction, not significant"
    )
  )
)


validation$Expected <- paste0(
  validation$Pathway,
  "  [expected ",
  validation$expected_direction,
  "]"
)


validation <- validation[
  order(validation$NES),
]


validation$Expected <- factor(
  validation$Expected,
  levels = validation$Expected
)


p_validation <- ggplot(
  validation,
  aes(
    x = Expected,
    y = NES,
    fill = Result
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "External validation of the GSE235980 resistance signature",
    subtitle = "1/10 pathways retained the predicted direction; none were concordant and significant",
    x = NULL,
    y = "GSE153867 normalized enrichment score (NES)",
    fill = "Validation result"
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 11,
      margin = margin(b = 10)
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(
      15, 25, 15, 15
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )


print(p_validation)


ggsave(
  "figures/GSE153867/Figure_3C_frozen_hypothesis_validation.png",
  p_validation,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/GSE153867/Figure_3C_frozen_hypothesis_validation.pdf",
  p_validation,
  width = 11,
  height = 7
)


# ============================================================
# FIGURE 3D - FILTER ROBUSTNESS
# ============================================================

sensitivity <- read.csv(
  "results/GSE153867/sensitivity/GSE153867_primary_vs_strict_filter.csv"
)


frozen_pathways <- validation$Hallmark


sensitivity$Frozen <- ifelse(
  sensitivity$ID %in% frozen_pathways,
  "Frozen pathway",
  "Other Hallmark"
)


rho <- cor(
  sensitivity$NES_primary,
  sensitivity$NES_strict,
  method = "spearman",
  use = "complete.obs"
)


same_direction <- sum(
  sign(sensitivity$NES_primary) ==
    sign(sensitivity$NES_strict),
  na.rm = TRUE
)


p_sensitivity <- ggplot(
  sensitivity,
  aes(
    x = NES_primary,
    y = NES_strict,
    color = Frozen
  )
) +
  geom_point(
    size = 2.6,
    alpha = 0.85
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.8,
    label = paste0(
      "Spearman rho = ",
      round(rho, 3),
      "\nSame NES direction = ",
      same_direction,
      "/",
      nrow(sensitivity)
    ),
    size = 4
  ) +
  labs(
    title = "Hallmark enrichment is robust to expression-filter choice",
    subtitle = "Primary filter compared with the pre-specified stricter sensitivity filter",
    x = "NES: primary analysis",
    y = "NES: stricter-filter analysis",
    color = NULL
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      15, 20, 15, 15
    )
  )


print(p_sensitivity)


ggsave(
  "figures/GSE153867/Figure_3D_filter_robustness.png",
  p_sensitivity,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  "figures/GSE153867/Figure_3D_filter_robustness.pdf",
  p_sensitivity,
  width = 7.5,
  height = 6.5
)


# ============================================================
# COMPLETE
# ============================================================

cat("\n========================================\n")
cat("GSE153867 FIGURE PACKAGE COMPLETE\n")
cat("========================================\n")

cat("\nCreated:\n")
cat("Figure 3A - PCA\n")
cat("Figure 3B - Significant Hallmarks\n")
cat("Figure 3C - Frozen hypothesis validation\n")
cat("Figure 3D - Filter robustness\n")

cat(
  "\nSaved under:\n",
  "figures/GSE153867/\n"
)