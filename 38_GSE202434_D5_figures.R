# ============================================================
# GSE202434 - DATASET 5
# FINAL D5 FIGURES
#
# Figure 6A - PCA of HCC1806 and PSN1
# Figure 6B - H3 vs P3 Hallmark concordance
# Figure 6C - Robust exposure-associated programs
# Figure 6D - Six-contrast mechanistic heatmap
#
# D5 main message:
# Resistant states differ strongly between models,
# but transcriptional response to drug exposure/withdrawal
# is substantially conserved.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(ggplot2)
library(dplyr)
library(tidyr)


# ============================================================
# 0. OUTPUT DIRECTORIES
# ============================================================

figdir <- "figures/GSE202434_D5"

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
# FIGURE 6A
# HCC1806 + PSN1 PCA
# ============================================================

hcc_pca <- read.csv(
  "results/GSE202434/HCC1806/HCC1806_PCA_coordinates.csv",
  stringsAsFactors = FALSE
)

psn_pca <- read.csv(
  "results/GSE202434/PSN1/PSN1_PCA_coordinates.csv",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Standardize condition labels
# ------------------------------------------------------------

hcc_pca$model <- "HCC1806"

hcc_pca$state <- ifelse(
  hcc_pca$group == "Parental",
  "Parental",
  ifelse(
    hcc_pca$group == "TalaR_M",
    "Maintained on drug",
    "Drug-free / withdrawn"
  )
)


psn_pca$model <- "PSN1"

psn_pca$state <- ifelse(
  psn_pca$group == "Parental",
  "Parental",
  ifelse(
    psn_pca$group == "TalaR_M",
    "Maintained on drug",
    "Drug-free / withdrawn"
  )
)


pca_data <- bind_rows(
  hcc_pca,
  psn_pca
)


pca_data$state <- factor(
  pca_data$state,
  levels = c(
    "Parental",
    "Drug-free / withdrawn",
    "Maintained on drug"
  )
)


# Short sample labels

pca_data$sample_short <- gsub(
  "^HCC1806_",
  "",
  pca_data$sample
)

pca_data$sample_short <- gsub(
  "^PSN1-",
  "",
  pca_data$sample_short
)


write.csv(
  pca_data,
  file.path(
    datadir,
    "Figure6A_PCA_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Percent variance is different between models.
# Use facet-specific axis labels in subtitle instead.
# ------------------------------------------------------------

p6a <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = state
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    aes(
      label = sample_short
    ),
    nudge_y = 2,
    size = 2.7,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ model,
    scales = "free"
  ) +
  
  labs(
    title = "A  Transcriptomic structure of the two D5 resistance models",
    
    subtitle = paste0(
      "HCC1806: PC1 = 88.9%, PC2 = 5.7%; ",
      "PSN1: PC1 = 53.8%, PC2 = 29.6%"
    ),
    
    x = "PC1",
    y = "PC2",
    colour = "State"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.position = "bottom"
  )


ggsave(
  file.path(
    figdir,
    "Figure6A_D5_PCA.png"
  ),
  p6a,
  width = 11,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure6A_D5_PCA.pdf"
  ),
  p6a,
  width = 11,
  height = 5.5
)


# ============================================================
# FIGURE 6B
# HCC1806 H3 vs PSN1 P3 HALLMARK CONCORDANCE
# ============================================================

d5_compare <- read.csv(
  "results/GSE202434/D5_integration/D5_H3_P3_Hallmark_master_comparison.csv",
  stringsAsFactors = FALSE
)


clean_name <- function(x) {
  
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


d5_compare$label <- clean_name(
  d5_compare$pathway
)


d5_compare$category <- ifelse(
  d5_compare$robust_shared_exposure_candidate,
  "Robust shared exposure program",
  ifelse(
    d5_compare$strong_divergent_candidate,
    "Strong divergent",
    "Other"
  )
)


# ------------------------------------------------------------
# Label selected strongest programs only
# ------------------------------------------------------------

label_pathways <- c(
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)


d5_compare$plot_label <- ifelse(
  d5_compare$pathway %in% label_pathways,
  d5_compare$label,
  ""
)


write.csv(
  d5_compare,
  file.path(
    datadir,
    "Figure6B_H3_P3_concordance_data.csv"
  ),
  row.names = FALSE
)


p6b <- ggplot(
  d5_compare,
  aes(
    x = H3_NES,
    y = P3_NES
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
    size = 3
  ) +
  
  geom_text(
    aes(
      label = plot_label
    ),
    size = 2.8,
    nudge_y = 0.12,
    check_overlap = TRUE
  ) +
  
  labs(
    title = "B  Conserved pathway response to talazoparib exposure state",
    
    subtitle = paste0(
      "HCC1806 H3 vs PSN1 P3: Spearman rho = 0.663; ",
      "39/50 Hallmarks share direction"
    ),
    
    x = "HCC1806 H3 NES\nMaintained vs drug-free resistant",
    
    y = "PSN1 P3 NES\nMaintained vs withdrawn",
    
    colour = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.position = "bottom"
  )


ggsave(
  file.path(
    figdir,
    "Figure6B_D5_H3_P3_concordance.png"
  ),
  p6b,
  width = 8.5,
  height = 7,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure6B_D5_H3_P3_concordance.pdf"
  ),
  p6b,
  width = 8.5,
  height = 7
)


# ============================================================
# FIGURE 6C
# ROBUST SHARED EXPOSURE-ASSOCIATED PROGRAMS
# ============================================================

exposure <- read.csv(
  "results/GSE202434/D5_integration/D5_robust_exposure_programs.csv",
  stringsAsFactors = FALSE
)


exposure$pathway_label <- clean_name(
  exposure$pathway
)


exposure_long <- exposure %>%
  
  select(
    pathway,
    pathway_label,
    exposure_state,
    H3_NES,
    P3_NES
  ) %>%
  
  pivot_longer(
    cols = c(
      H3_NES,
      P3_NES
    ),
    names_to = "model",
    values_to = "NES"
  )


exposure_long$model <- factor(
  exposure_long$model,
  levels = c(
    "H3_NES",
    "P3_NES"
  ),
  labels = c(
    "HCC1806",
    "PSN1"
  )
)


# Order by average NES

exposure_order <- exposure_long %>%
  
  group_by(
    pathway_label
  ) %>%
  
  summarise(
    mean_NES = mean(NES),
    .groups = "drop"
  ) %>%
  
  arrange(
    mean_NES
  ) %>%
  
  pull(
    pathway_label
  )


exposure_long$pathway_label <- factor(
  exposure_long$pathway_label,
  levels = exposure_order
)


write.csv(
  exposure_long,
  file.path(
    datadir,
    "Figure6C_exposure_programs_data.csv"
  ),
  row.names = FALSE
)


p6c <- ggplot(
  exposure_long,
  aes(
    x = NES,
    y = pathway_label,
    shape = model
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_line(
    aes(
      group = pathway_label
    ),
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 3.3
  ) +
  
  labs(
    title = "C  Robust exposure-associated Hallmark programs",
    
    subtitle = paste0(
      "Positive NES = higher during ongoing talazoparib exposure; ",
      "negative NES = higher after withdrawal"
    ),
    
    x = "Normalized enrichment score (NES)",
    y = NULL,
    shape = "Model"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 8.5
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.position = "bottom"
  )


ggsave(
  file.path(
    figdir,
    "Figure6C_D5_exposure_programs.png"
  ),
  p6c,
  width = 9,
  height = 8,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure6C_D5_exposure_programs.pdf"
  ),
  p6c,
  width = 9,
  height = 8
)


# ============================================================
# FIGURE 6D
# SIX-CONTRAST MECHANISTIC HEATMAP
# ============================================================

mech <- read.csv(
  "results/GSE202434/D5_integration/D5_six_contrast_mechanistic_master.csv",
  stringsAsFactors = FALSE
)


focus <- c(
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  
  "HALLMARK_GLYCOLYSIS",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)


heat <- mech %>%
  
  filter(
    pathway %in% focus
  ) %>%
  
  select(
    pathway,
    
    H1_NES,
    H1_padj,
    
    H2_NES,
    H2_padj,
    
    H3_NES,
    H3_padj,
    
    P1_NES,
    P1_padj,
    
    P2_NES,
    P2_padj,
    
    P3_NES,
    P3_padj
  )


heat$pathway_label <- clean_name(
  heat$pathway
)


# ------------------------------------------------------------
# Convert six contrasts to long format
# ------------------------------------------------------------

heat_long <- bind_rows(
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "H1",
    NES = heat$H1_NES,
    padj = heat$H1_padj
  ),
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "H2",
    NES = heat$H2_NES,
    padj = heat$H2_padj
  ),
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "H3",
    NES = heat$H3_NES,
    padj = heat$H3_padj
  ),
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "P1",
    NES = heat$P1_NES,
    padj = heat$P1_padj
  ),
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "P2",
    NES = heat$P2_NES,
    padj = heat$P2_padj
  ),
  
  data.frame(
    pathway = heat$pathway,
    pathway_label = heat$pathway_label,
    contrast = "P3",
    NES = heat$P3_NES,
    padj = heat$P3_padj
  )
)


heat_long$significant <-
  heat_long$padj < 0.05


heat_long$contrast <- factor(
  
  heat_long$contrast,
  
  levels = c(
    "H1",
    "H2",
    "H3",
    "P1",
    "P2",
    "P3"
  ),
  
  labels = c(
    "HCC\nM vs P",
    "HCC\nDF vs P",
    "HCC\nM vs DF",
    "PSN1\nM vs P",
    "PSN1\nWD vs P",
    "PSN1\nM vs WD"
  )
)


heat_long$pathway_label <- factor(
  
  heat_long$pathway_label,
  
  levels = rev(
    c(
      "E2F TARGETS",
      "G2M CHECKPOINT",
      "MYC TARGETS V1",
      "MTORC1 SIGNALING",
      "OXIDATIVE PHOSPHORYLATION",
      "GLYCOLYSIS",
      "TNFA SIGNALING VIA NFKB",
      "EPITHELIAL MESENCHYMAL TRANSITION"
    )
  )
)


write.csv(
  heat_long,
  file.path(
    datadir,
    "Figure6D_six_contrast_heatmap_data.csv"
  ),
  row.names = FALSE
)


p6d <- ggplot(
  heat_long,
  aes(
    x = contrast,
    y = pathway_label,
    fill = NES
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        significant,
        "*",
        ""
      )
    ),
    size = 5
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "NES"
  ) +
  
  labs(
    title = "D  Resistance states diverge while exposure responses converge",
    
    subtitle = paste0(
      "HCC1806 and PSN1 show opposing resistant-state pathway profiles, ",
      "but similar maintained-versus-withdrawn responses. ",
      "* BH-adjusted P < 0.05."
    ),
    
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      face = "bold",
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


ggsave(
  file.path(
    figdir,
    "Figure6D_D5_mechanistic_heatmap.png"
  ),
  p6d,
  width = 10,
  height = 6.5,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure6D_D5_mechanistic_heatmap.pdf"
  ),
  p6d,
  width = 10,
  height = 6.5
)


# ============================================================
# OPTIONAL COMBINED FIGURE
# ============================================================

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {
  
  library(patchwork)
  
  combined <- (
    p6a |
      p6b
  ) /
    (
      p6c |
        p6d
    ) +
    
    plot_annotation(
      title = paste0(
        "Figure 6. Treatment-state effects dominate the reproducible ",
        "transcriptional response in GSE202434"
      )
    )
  
  
  ggsave(
    file.path(
      figdir,
      "Figure6_D5_combined.png"
    ),
    combined,
    width = 18,
    height = 14,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      figdir,
      "Figure6_D5_combined.pdf"
    ),
    combined,
    width = 18,
    height = 14
  )
}


# ============================================================
# FINAL QC
# ============================================================

cat("\n========================================\n")
cat("D5 FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Figure 6A - PCA\n"
)

cat(
  "Figure 6B - H3 vs P3 Hallmark concordance\n"
)

cat(
  "Figure 6C - Robust exposure-associated programs\n"
)

cat(
  "Figure 6D - Six-contrast mechanistic heatmap\n"
)

cat(
  "\nOutput:",
  figdir,
  "\n"
)