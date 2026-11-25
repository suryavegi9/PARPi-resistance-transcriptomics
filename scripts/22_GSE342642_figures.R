# ============================================================
# GSE342642 - DATASET 4
# FINAL FIGURES 4A-4D
#
# 4A: PCA / sample structure
# 4B: Significant Hallmark GSEA
# 4C: Frozen Dataset-1 validation
# 4D: Hallmark robustness across sensitivity analyses
#
# ASCII-safe version for Windows PDF/PNG rendering
# ============================================================

library(ggplot2)

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

figdir <- "figures/GSE342642"
datadir <- file.path(figdir, "figure_data")

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
# HELPER: CLEAN HALLMARK NAMES
# ============================================================

clean_hallmark <- function(x) {
  
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
  
  # Restore biological capitalization
  x <- gsub("E2f", "E2F", x)
  x <- gsub("G2m", "G2M", x)
  x <- gsub("Myc", "MYC", x)
  x <- gsub("Mtorc1", "mTORC1", x)
  x <- gsub("Dna", "DNA", x)
  
  # ASCII-safe versions
  x <- gsub("Tnfa", "TNFa", x)
  x <- gsub("Nfkb", "NF-kB", x)
  
  x <- gsub("Uv", "UV", x)
  x <- gsub("Il2", "IL2", x)
  x <- gsub(
    "Pi3k Akt Mtor",
    "PI3K-AKT-mTOR",
    x
  )
  
  x
}


# ============================================================
# FIGURE 4A
# PCA OF 12 REGIONAL SAMPLES
# ============================================================

pca <- read.csv(
  "results/GSE342642/QC/PCA_coordinates.csv",
  stringsAsFactors = FALSE
)


# Short sample labels
pca$label <- paste0(
  pca$tumor,
  "-R",
  pca$region
)


pca$condition <- factor(
  pca$condition,
  levels = c(
    "Control",
    "Resistant"
  )
)


pca$tumor <- factor(
  pca$tumor,
  levels = c(
    "437",
    "440",
    "433",
    "438"
  )
)


write.csv(
  pca,
  file.path(
    datadir,
    "Figure4A_PCA_data.csv"
  ),
  row.names = FALSE
)


p4A <- ggplot(
  pca,
  aes(
    x = PC1,
    y = PC2,
    color = tumor,
    shape = condition
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    vjust = -1,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  labs(
    title = "GSE342642 PDX transcriptomic structure",
    subtitle = "Three spatial regions sampled from each of four tumors",
    x = "PC1 (49.7%)",
    y = "PC2 (21.0%)",
    color = "Tumor",
    shape = "Condition"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )


ggsave(
  file.path(
    figdir,
    "Figure4A_PCA.pdf"
  ),
  p4A,
  width = 8,
  height = 6.5
)


ggsave(
  file.path(
    figdir,
    "Figure4A_PCA.png"
  ),
  p4A,
  width = 8,
  height = 6.5,
  dpi = 400
)


# ============================================================
# FIGURE 4B
# SIGNIFICANT HALLMARK GSEA
# ============================================================

sig <- read.csv(
  paste0(
    "results/GSE342642/Hallmark_GSEA/",
    "GSE342642_Hallmark_GSEA_significant.csv"
  ),
  stringsAsFactors = FALSE
)


sig$Pathway <- clean_hallmark(
  sig$ID
)


sig$Direction <- ifelse(
  sig$NES > 0,
  "Higher in resistant",
  "Lower in resistant"
)


# Order from most negative to most positive
sig$Pathway <- factor(
  sig$Pathway,
  levels = sig$Pathway[
    order(sig$NES)
  ]
)


write.csv(
  sig,
  file.path(
    datadir,
    "Figure4B_significant_Hallmarks.csv"
  ),
  row.names = FALSE
)


p4B <- ggplot(
  sig,
  aes(
    x = NES,
    y = Pathway,
    fill = Direction
  )
) +
  
  geom_col(
    width = 0.72
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Hallmark programs associated with acquired olaparib resistance",
    subtitle = "Significant pathways: BH-adjusted P < 0.05",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top",
    axis.text.y = element_text(
      size = 9
    )
  )


ggsave(
  file.path(
    figdir,
    "Figure4B_significant_Hallmarks.pdf"
  ),
  p4B,
  width = 9,
  height = 8
)


ggsave(
  file.path(
    figdir,
    "Figure4B_significant_Hallmarks.png"
  ),
  p4B,
  width = 9,
  height = 8,
  dpi = 400
)


# ============================================================
# FIGURE 4C
# FROZEN DATASET-1 VALIDATION
# ============================================================

frozen <- read.csv(
  paste0(
    "results/GSE342642/Hallmark_GSEA/",
    "GSE342642_frozen_D1_validation.csv"
  ),
  stringsAsFactors = FALSE
)


frozen$Pathway <- clean_hallmark(
  frozen$pathway
)


# ASCII-safe expected directions
frozen$Expected <- ifelse(
  frozen$expected_direction == "UP",
  "UP",
  "DOWN"
)


frozen$Pathway_label <- paste0(
  frozen$Pathway,
  "   [D1 ",
  frozen$Expected,
  "]"
)


frozen$Validation <- ifelse(
  frozen$concordant_and_significant,
  "Concordant + significant",
  ifelse(
    frozen$direction_concordant,
    "Direction concordant",
    "Discordant"
  )
)


# Preserve original frozen-hypothesis order
frozen$Pathway_label <- factor(
  frozen$Pathway_label,
  levels = rev(
    frozen$Pathway_label
  )
)


write.csv(
  frozen,
  file.path(
    datadir,
    "Figure4C_frozen_D1_validation.csv"
  ),
  row.names = FALSE
)


p4C <- ggplot(
  frozen,
  aes(
    x = NES,
    y = Pathway_label,
    color = Validation
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = Pathway_label
    ),
    linewidth = 0.7,
    color = "grey70"
  ) +
  
  geom_point(
    size = 4
  ) +
  
  labs(
    title = "Pre-specified Dataset-1 pathway validation in the PDX model",
    subtitle = "7/10 direction concordant; 5/10 concordant and significant",
    x = "Dataset-4 NES",
    y = NULL,
    color = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top",
    axis.text.y = element_text(
      size = 9
    )
  )


ggsave(
  file.path(
    figdir,
    "Figure4C_frozen_D1_validation.pdf"
  ),
  p4C,
  width = 9,
  height = 6.5
)


ggsave(
  file.path(
    figdir,
    "Figure4C_frozen_D1_validation.png"
  ),
  p4C,
  width = 9,
  height = 6.5,
  dpi = 400
)


# ============================================================
# FIGURE 4D
# ROBUSTNESS OF CORE PATHWAYS
# ============================================================

robust <- read.csv(
  paste0(
    "results/GSE342642/Hallmark_sensitivity/",
    "Hallmark_three_analysis_comparison.csv"
  ),
  stringsAsFactors = FALSE
)


focus_pathways <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)


robust <- robust[
  match(
    focus_pathways,
    robust$ID
  ),
]


# Convert wide table to long format
robust_long <- rbind(
  
  data.frame(
    ID = robust$ID,
    Analysis = "Primary blocked",
    NES = robust$NES_primary
  ),
  
  data.frame(
    ID = robust$ID,
    Analysis = "Without 440-R3",
    NES = robust$NES_without440R3
  ),
  
  data.frame(
    ID = robust$ID,
    Analysis = "Tumor pseudobulk",
    NES = robust$NES_pseudobulk
  )
)


robust_long$Pathway <- clean_hallmark(
  robust_long$ID
)


robust_long$Analysis <- factor(
  robust_long$Analysis,
  levels = c(
    "Primary blocked",
    "Without 440-R3",
    "Tumor pseudobulk"
  )
)


robust_long$Pathway <- factor(
  robust_long$Pathway,
  levels = rev(
    clean_hallmark(
      focus_pathways
    )
  )
)


write.csv(
  robust_long,
  file.path(
    datadir,
    "Figure4D_Hallmark_robustness.csv"
  ),
  row.names = FALSE
)


p4D <- ggplot(
  robust_long,
  aes(
    x = NES,
    y = Pathway,
    shape = Analysis
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_point(
    position = position_dodge(
      width = 0.55
    ),
    size = 3.8
  ) +
  
  labs(
    title = "Key Dataset-4 pathway signals are robust",
    subtitle = paste0(
      "Primary vs -440-R3: rho = 0.996; ",
      "primary vs pseudobulk: rho = 0.983"
    ),
    x = "Normalized enrichment score (NES)",
    y = NULL,
    shape = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )


ggsave(
  file.path(
    figdir,
    "Figure4D_Hallmark_robustness.pdf"
  ),
  p4D,
  width = 8.5,
  height = 5.5
)


ggsave(
  file.path(
    figdir,
    "Figure4D_Hallmark_robustness.png"
  ),
  p4D,
  width = 8.5,
  height = 5.5,
  dpi = 400
)


# ============================================================
# FINISHED
# ============================================================

cat("\n========================================\n")
cat("DATASET 4 FIGURES COMPLETE\n")
cat("========================================\n\n")

cat("Created:\n")
cat("Figure 4A - PCA\n")
cat("Figure 4B - Significant Hallmark GSEA\n")
cat("Figure 4C - Frozen D1 validation\n")
cat("Figure 4D - Hallmark robustness\n\n")

cat(
  "Output directory:\n",
  normalizePath(figdir),
  "\n"
)