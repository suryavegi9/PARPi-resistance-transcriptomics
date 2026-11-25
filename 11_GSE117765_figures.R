# ============================================================
# GSE117765 - PUBLICATION-QUALITY FIGURES
#
# Dataset 2:
# PEO1 parental vs independently derived olaparib-resistant clones
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(ggplot2)
library(ggrepel)

dir.create(
  "figures/GSE117765",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# FIGURE 2A - PCA
# ============================================================

log_expr <- read.csv(
  "results/GSE117765/GSE117765_log2TPMplus1.csv",
  row.names = 1,
  check.names = FALSE
)

log_expr <- as.matrix(log_expr)

metadata <- data.frame(
  
  sample = colnames(log_expr),
  
  phenotype = ifelse(
    grepl("Adherent", colnames(log_expr)),
    "Parental",
    "Resistant"
  ),
  
  label = c(
    "Parent-1",
    "Parent-2",
    "Clone-10",
    "Clone-18",
    "Clone-4",
    "Clone-5"
  ),
  
  row.names = colnames(log_expr)
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
  ],
  
  label = metadata[
    rownames(pca$x),
    "label"
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
    aes(label = label),
    size = 4,
    show.legend = FALSE,
    box.padding = 0.4,
    point.padding = 0.4
  ) +
  labs(
    title = "Transcriptomic separation of parental and olaparib-resistant PEO1 cells",
    subtitle = "Four independently derived resistant clones show a strong resistance-associated shift",
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
  "figures/GSE117765/Figure_2A_PCA_refined.png",
  p_pca,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE117765/Figure_2A_PCA_refined.pdf",
  p_pca,
  width = 9,
  height = 6
)


# ============================================================
# FIGURE 2B - SIGNIFICANT HALLMARK GSEA
# ============================================================

gsea <- read.csv(
  "results/GSE117765/GSE117765_Hallmark_GSEA_all.csv"
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
    width = 0.65
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "Significant Hallmark programs in GSE117765",
    subtitle = "Hallmark GSEA at BH-adjusted P < 0.05",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = NULL
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


print(p_hallmark)


ggsave(
  "figures/GSE117765/Figure_2B_significant_Hallmarks.png",
  p_hallmark,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  "figures/GSE117765/Figure_2B_significant_Hallmarks.pdf",
  p_hallmark,
  width = 9,
  height = 5.5
)


# ============================================================
# FIGURE 2C - FROZEN HYPOTHESIS VALIDATION
# ============================================================

validation <- read.csv(
  "results/GSE117765/GSE117765_frozen_hypothesis_validation.csv"
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
    subtitle = paste0(
      "5/10 pathways retained the predicted direction; ",
      "none were concordant and significant"
    ),
    x = NULL,
    y = "GSE117765 normalized enrichment score (NES)",
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
  "figures/GSE117765/Figure_2C_frozen_hypothesis_validation.png",
  p_validation,
  width = 11,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/GSE117765/Figure_2C_frozen_hypothesis_validation.pdf",
  p_validation,
  width = 11,
  height = 7
)

 # ============================================================
# FIGURE 2D - SENSITIVITY ANALYSIS CONCORDANCE
# ============================================================

sensitivity <- read.csv(
  "results/GSE117765/sensitivity/Hallmark_primary_vs_complete_case.csv"
)


sensitivity$Frozen <- ifelse(
  
  sensitivity$ID %in%
    validation$Hallmark,
  
  "Frozen pathway",
  "Other Hallmark"
)


rho_sensitivity <- cor(
  sensitivity$NES_primary,
  sensitivity$NES_complete_case,
  method = "spearman",
  use = "complete.obs"
)


same_direction <- sum(
  sign(sensitivity$NES_primary) ==
    sign(sensitivity$NES_complete_case)
)


p_sensitivity <- ggplot(
  sensitivity,
  aes(
    x = NES_primary,
    y = NES_complete_case,
    color = Frozen
  )
) +
  geom_point(
    alpha = 0.8,
    size = 2.5
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
      round(rho_sensitivity, 3),
      "\nSame NES direction = ",
      same_direction,
      "/",
      nrow(sensitivity)
    ),
    size = 4
  ) +
  labs(
    title = "Hallmark enrichment is robust to missing-value handling",
    subtitle = "Primary NA→0 analysis compared with complete-case sensitivity analysis",
    x = "NES: primary analysis",
    y = "NES: complete-case analysis",
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
  "figures/GSE117765/Figure_2D_sensitivity_concordance.png",
  p_sensitivity,
  width = 7,
  height = 6.5,
  dpi = 300
)

ggsave(
  "figures/GSE117765/Figure_2D_sensitivity_concordance.pdf",
  p_sensitivity,
  width = 7,
  height = 6.5
)


# ============================================================
# COMPLETE
# ============================================================

cat("\n========================================\n")
cat("GSE117765 FIGURE PACKAGE COMPLETE\n")
cat("========================================\n")

cat("\nCreated:\n")

cat(
  "Figure 2A - PCA\n",
  "Figure 2B - Significant Hallmarks\n",
  "Figure 2C - Frozen hypothesis validation\n",
  "Figure 2D - Sensitivity concordance\n"
)

cat(
  "\nSaved under:\n",
  "figures/GSE117765/\n"
)