# ============================================================
# GSE202434 - DATASET 5
# ENHANCED PUBLICATION-READY FIGURES
#
# Figure 6A - PCA: HCC1806 + PSN1
# Figure 6B - H3 vs P3 Hallmark concordance
# Figure 6C - Robust exposure-associated programs
# Figure 6D - Six-contrast mechanistic heatmap
#
# Main D5 message:
# Resistant-state transcriptomes diverge between models,
# whereas pathway responses to talazoparib exposure/
# withdrawal show substantial convergence.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. PACKAGES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# patchwork is used only for the combined figure
# install.packages("patchwork") if needed


# ============================================================
# 1. OUTPUT DIRECTORIES
# ============================================================

figdir <- "figures/GSE202434_D5/enhanced"

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
# 2. COMMON SETTINGS
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


# Accessible / publication-friendly palette

state_colors <- c(
  
  "Parental" =
    "#D55E00",
  
  "Drug-free / withdrawn" =
    "#009E73",
  
  "Maintained on drug" =
    "#0072B2"
)


category_colors <- c(
  
  "Other" =
    "grey70",
  
  "Robust shared exposure program" =
    "#009E73",
  
  "Strong divergent" =
    "#0072B2"
)


# ============================================================
# ============================================================
# FIGURE 6A
# PCA OF HCC1806 + PSN1
# ============================================================
# ============================================================


hcc_pca <- read.csv(
  
  "results/GSE202434/HCC1806/HCC1806_PCA_coordinates.csv",
  
  stringsAsFactors = FALSE
)


psn_pca <- read.csv(
  
  "results/GSE202434/PSN1/PSN1_PCA_coordinates.csv",
  
  stringsAsFactors = FALSE
)


# ============================================================
# 6A.1 STANDARDIZE BIOLOGICAL STATE
# ============================================================

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


# ============================================================
# 6A.2 CREATE SHORT, UNAMBIGUOUS LABELS
# ============================================================

get_rep <- function(x) {
  
  sub(
    ".*(?:rep|-)([123])$",
    "\\1",
    x,
    perl = TRUE
  )
}


pca_data$replicate <- get_rep(
  pca_data$sample
)


pca_data$sample_short <- case_when(
  
  pca_data$group == "Parental" ~
    paste0(
      "P",
      pca_data$replicate
    ),
  
  pca_data$model == "HCC1806" &
    pca_data$group == "TalaR_DF" ~
    paste0(
      "DF",
      pca_data$replicate
    ),
  
  pca_data$model == "PSN1" &
    pca_data$group == "Withdraw" ~
    paste0(
      "WD",
      pca_data$replicate
    ),
  
  pca_data$group == "TalaR_M" ~
    paste0(
      "M",
      pca_data$replicate
    ),
  
  TRUE ~
    pca_data$sample
)


write.csv(
  
  pca_data,
  
  file.path(
    datadir,
    "Figure6A_PCA_data_enhanced.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 6A.3 PCA PLOT
# ============================================================

p6a <- ggplot(
  
  pca_data,
  
  aes(
    x = PC1,
    y = PC2,
    colour = state
  )
) +
  
  geom_point(
    size = 4.3,
    alpha = 0.95
  ) +
  
  geom_text_repel(
    
    aes(
      label = sample_short
    ),
    
    size = 3.5,
    
    box.padding = 0.35,
    
    point.padding = 0.25,
    
    min.segment.length = 0,
    
    segment.colour = "grey60",
    
    max.overlaps = Inf,
    
    show.legend = FALSE,
    
    seed = 123
  ) +
  
  facet_wrap(
    ~ model,
    scales = "free"
  ) +
  
  scale_colour_manual(
    values = state_colors
  ) +
  
  labs(
    
    title =
      "A  Transcriptomic structure of the two D5 resistance models",
    
    subtitle =
      paste0(
        "HCC1806: PC1 = 88.9%, PC2 = 5.7%; ",
        "PSN1: PC1 = 53.8%, PC2 = 29.6%"
      ),
    
    x = "PC1",
    
    y = "PC2",
    
    colour = "State"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    strip.background =
      element_rect(
        fill = "grey95",
        colour = "grey25"
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 11
      ),
    
    legend.position =
      "bottom",
    
    legend.title =
      element_text(
        face = "bold"
      ),
    
    plot.margin =
      margin(
        12,
        20,
        12,
        12
      )
  )


ggsave(
  
  file.path(
    figdir,
    "Figure6A_D5_PCA_enhanced.png"
  ),
  
  p6a,
  
  width = 11,
  
  height = 6,
  
  dpi = 400
)


ggsave(
  
  file.path(
    figdir,
    "Figure6A_D5_PCA_enhanced.pdf"
  ),
  
  p6a,
  
  width = 11,
  
  height = 6
)


# ============================================================
# ============================================================
# FIGURE 6B
# H3 vs P3 HALLMARK CONCORDANCE
# ============================================================
# ============================================================


d5_compare <- read.csv(
  
  paste0(
    "results/GSE202434/D5_integration/",
    "D5_H3_P3_Hallmark_master_comparison.csv"
  ),
  
  stringsAsFactors = FALSE
)


d5_compare$label <- clean_hallmark(
  d5_compare$pathway
)


# ============================================================
# 6B.1 CLASSIFICATION
# ============================================================

d5_compare$category <- ifelse(
  
  d5_compare$robust_shared_exposure_candidate,
  
  "Robust shared exposure program",
  
  ifelse(
    
    d5_compare$strong_divergent_candidate,
    
    "Strong divergent",
    
    "Other"
  )
)


d5_compare$category <- factor(
  
  d5_compare$category,
  
  levels = c(
    "Other",
    "Robust shared exposure program",
    "Strong divergent"
  )
)


# ============================================================
# 6B.2 PATHWAYS TO LABEL
# ============================================================

label_pathways <- c(
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_MTORC1_SIGNALING",
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)


label_data <- d5_compare[
  d5_compare$pathway %in%
    label_pathways,
]


write.csv(
  
  d5_compare,
  
  file.path(
    datadir,
    "Figure6B_H3_P3_concordance_enhanced.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 6B.3 SCATTER
# ============================================================

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
    
    linewidth = 0.45,
    
    colour = "grey45"
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.45,
    
    colour = "grey45"
  ) +
  
  geom_abline(
    
    slope = 1,
    
    intercept = 0,
    
    linetype = "dotted",
    
    linewidth = 0.55,
    
    colour = "grey30"
  ) +
  
  geom_point(
    
    aes(
      colour = category
    ),
    
    size = 3.4,
    
    alpha = 0.95
  ) +
  
  geom_text_repel(
    
    data = label_data,
    
    aes(
      label = label
    ),
    
    size = 3.2,
    
    box.padding = 0.45,
    
    point.padding = 0.35,
    
    min.segment.length = 0,
    
    segment.colour = "grey50",
    
    seed = 123,
    
    max.overlaps = Inf,
    
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = category_colors
  ) +
  
  scale_x_continuous(
    
    expand = expansion(
      mult = c(
        0.15,
        0.12
      )
    )
  ) +
  
  scale_y_continuous(
    
    expand = expansion(
      mult = c(
        0.12,
        0.15
      )
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    
    title =
      "B  Conserved pathway response to talazoparib exposure state",
    
    subtitle =
      paste0(
        "HCC1806 H3 vs PSN1 P3: Spearman rho = 0.663; ",
        "39/50 Hallmarks share direction"
      ),
    
    x =
      paste0(
        "HCC1806 H3 NES\n",
        "Maintained vs drug-free resistant"
      ),
    
    y =
      paste0(
        "PSN1 P3 NES\n",
        "Maintained vs withdrawn"
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
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 11
      ),
    
    legend.position =
      "bottom",
    
    plot.margin =
      margin(
        15,
        35,
        15,
        35
      )
  )


ggsave(
  
  file.path(
    figdir,
    "Figure6B_D5_H3_P3_concordance_enhanced.png"
  ),
  
  p6b,
  
  width = 9,
  
  height = 7.5,
  
  dpi = 400
)


ggsave(
  
  file.path(
    figdir,
    "Figure6B_D5_H3_P3_concordance_enhanced.pdf"
  ),
  
  p6b,
  
  width = 9,
  
  height = 7.5
)


# ============================================================
# ============================================================
# FIGURE 6C
# ROBUST SHARED EXPOSURE PROGRAMS
# ============================================================
# ============================================================


exposure <- read.csv(
  
  paste0(
    "results/GSE202434/D5_integration/",
    "D5_robust_exposure_programs.csv"
  ),
  
  stringsAsFactors = FALSE
)


exposure$pathway_label <- clean_hallmark(
  exposure$pathway
)


# ============================================================
# 6C.1 ORDER PATHWAYS
# ============================================================

exposure$mean_NES <- rowMeans(
  
  exposure[
    ,
    c(
      "H3_NES",
      "P3_NES"
    )
  ]
)


exposure <- exposure[
  order(
    exposure$mean_NES
  ),
]


exposure$pathway_label <- factor(
  
  exposure$pathway_label,
  
  levels =
    exposure$pathway_label
)


exposure$exposure_state <- factor(
  
  exposure$exposure_state,
  
  levels = c(
    "Higher_after_withdrawal",
    "Higher_on_drug"
  ),
  
  labels = c(
    "Higher after withdrawal",
    "Higher during ongoing drug exposure"
  )
)


write.csv(
  
  exposure,
  
  file.path(
    datadir,
    "Figure6C_exposure_programs_enhanced.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 6C.2 DUMBBELL PLOT
# ============================================================

p6c <- ggplot(
  exposure
) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.45,
    
    colour = "grey45"
  ) +
  
  geom_segment(
    
    aes(
      x = H3_NES,
      xend = P3_NES,
      y = pathway_label,
      yend = pathway_label
    ),
    
    colour = "grey55",
    
    linewidth = 0.75
  ) +
  
  geom_point(
    
    aes(
      x = H3_NES,
      y = pathway_label
    ),
    
    shape = 16,
    
    size = 3.3,
    
    colour = "#0072B2"
  ) +
  
  geom_point(
    
    aes(
      x = P3_NES,
      y = pathway_label
    ),
    
    shape = 17,
    
    size = 3.6,
    
    colour = "#D55E00"
  ) +
  
  facet_grid(
    
    exposure_state ~ .,
    
    scales = "free_y",
    
    space = "free_y",
    
    switch = "y"
  ) +
  
  labs(
    
    title =
      "C  Robust exposure-associated Hallmark programs",
    
    subtitle =
      paste0(
        "Concordant and significant in HCC1806 and PSN1, ",
        "with PSN1 outlier sensitivity support"
      ),
    
    x =
      "Normalized enrichment score (NES)",
    
    y = NULL,
    
    caption =
      "Circle = HCC1806 H3; triangle = PSN1 P3"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    strip.background =
      element_blank(),
    
    strip.placement =
      "outside",
    
    strip.text.y.left =
      element_text(
        angle = 0,
        face = "bold",
        size = 9
      ),
    
    axis.text.y =
      element_text(
        size = 9
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 10.5
      ),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9
      ),
    
    plot.margin =
      margin(
        12,
        15,
        12,
        20
      )
  )


ggsave(
  
  file.path(
    figdir,
    "Figure6C_D5_exposure_programs_enhanced.png"
  ),
  
  p6c,
  
  width = 10,
  
  height = 8.5,
  
  dpi = 400
)


ggsave(
  
  file.path(
    figdir,
    "Figure6C_D5_exposure_programs_enhanced.pdf"
  ),
  
  p6c,
  
  width = 10,
  
  height = 8.5
)


# ============================================================
# ============================================================
# FIGURE 6D
# SIX-CONTRAST MECHANISTIC HEATMAP
# ============================================================
# ============================================================


mech <- read.csv(
  
  paste0(
    "results/GSE202434/D5_integration/",
    "D5_six_contrast_mechanistic_master.csv"
  ),
  
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
  )


# ============================================================
# 6D.1 BUILD LONG TABLE
# ============================================================

heat_long <- bind_rows(
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "H1",
    
    NES = heat$H1_NES,
    
    padj = heat$H1_padj
  ),
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "H2",
    
    NES = heat$H2_NES,
    
    padj = heat$H2_padj
  ),
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "H3",
    
    NES = heat$H3_NES,
    
    padj = heat$H3_padj
  ),
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "P1",
    
    NES = heat$P1_NES,
    
    padj = heat$P1_padj
  ),
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "P2",
    
    NES = heat$P2_NES,
    
    padj = heat$P2_padj
  ),
  
  data.frame(
    
    pathway = heat$pathway,
    
    contrast = "P3",
    
    NES = heat$P3_NES,
    
    padj = heat$P3_padj
  )
)


heat_long$pathway_label <- clean_hallmark(
  heat_long$pathway
)


heat_long$significant <-
  heat_long$padj < 0.05


# ============================================================
# 6D.2 FIX ORDER
# ============================================================

path_order <- c(
  
  "E2F TARGETS",
  
  "G2M CHECKPOINT",
  
  "MYC TARGETS V1",
  
  "MTORC1 SIGNALING",
  
  "OXIDATIVE PHOSPHORYLATION",
  
  "GLYCOLYSIS",
  
  "TNFA SIGNALING VIA NFKB",
  
  "EPITHELIAL MESENCHYMAL TRANSITION"
)


heat_long$pathway_label <- factor(
  
  heat_long$pathway_label,
  
  levels = rev(
    path_order
  )
)


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


# Symmetric NES scale

max_abs_nes <- max(
  abs(
    heat_long$NES
  )
)


write.csv(
  
  heat_long,
  
  file.path(
    datadir,
    "Figure6D_six_contrast_heatmap_enhanced.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 6D.3 HEATMAP
# ============================================================

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
  
  # Visual separator between HCC1806 and PSN1
  
  geom_vline(
    
    xintercept = 3.5,
    
    linewidth = 0.8,
    
    colour = "grey25"
  ) +
  
  scale_fill_gradient2(
    
    low = "#2166AC",
    
    mid = "white",
    
    high = "#B2182B",
    
    midpoint = 0,
    
    limits = c(
      -max_abs_nes,
      max_abs_nes
    ),
    
    name = "NES"
  ) +
  
  labs(
    
    title =
      "D  Resistance states diverge while exposure responses converge",
    
    subtitle =
      paste0(
        "Opposing resistant-state profiles, but concordant ",
        "maintained-versus-withdrawn responses"
      ),
    
    x = NULL,
    
    y = NULL,
    
    caption =
      paste0(
        "* BH-adjusted P < 0.05. ",
        "M = maintained on drug; DF = drug-free resistant; ",
        "WD = withdrawn; P = parental."
      )
  ) +
  
  theme_minimal(
    base_size = 12
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
        size = 10
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 10.5
      ),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9
      ),
    
    legend.title =
      element_text(
        face = "bold"
      ),
    
    plot.margin =
      margin(
        15,
        25,
        15,
        15
      )
  )


ggsave(
  
  file.path(
    figdir,
    "Figure6D_D5_mechanistic_heatmap_enhanced.png"
  ),
  
  p6d,
  
  width = 10.5,
  
  height = 7,
  
  dpi = 400
)


ggsave(
  
  file.path(
    figdir,
    "Figure6D_D5_mechanistic_heatmap_enhanced.pdf"
  ),
  
  p6d,
  
  width = 10.5,
  
  height = 7
)


# ============================================================
# ============================================================
# COMBINED FIGURE 6
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
    
    p6a + p6b +
      
      plot_layout(
        widths = c(
          1.15,
          1
        )
      )
    
  ) /
    
    (
      
      p6c + p6d +
        
        plot_layout(
          widths = c(
            0.9,
            1.1
          )
        )
    ) +
    
    plot_annotation(
      
      title =
        paste0(
          "Figure 6. Treatment-state effects dominate the ",
          "reproducible transcriptional response in GSE202434"
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
      "Figure6_D5_combined_enhanced.png"
    ),
    
    combined,
    
    width = 19,
    
    height = 14,
    
    dpi = 400
  )
  
  
  ggsave(
    
    file.path(
      figdir,
      "Figure6_D5_combined_enhanced.pdf"
    ),
    
    combined,
    
    width = 19,
    
    height = 14
  )
  
  
  cat(
    "\nEnhanced combined Figure 6 generated.\n"
  )
  
} else {
  
  cat(
    paste0(
      "\npatchwork is not installed.\n",
      "Panels 6A-6D were generated successfully.\n",
      "Run install.packages('patchwork') if you want the combined figure.\n"
    )
  )
}


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n========================================\n")
cat("ENHANCED D5 FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Figure 6A - enhanced PCA\n"
)

cat(
  "Figure 6B - enhanced H3/P3 concordance\n"
)

cat(
  "Figure 6C - enhanced exposure-program dumbbell plot\n"
)

cat(
  "Figure 6D - enhanced mechanistic heatmap\n"
)

cat(
  "\nOutput directory:\n",
  figdir,
  "\n"
)

cat("\n========================================\n")