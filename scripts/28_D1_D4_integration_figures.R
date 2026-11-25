# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 6: INTEGRATION FIGURES
#
# Figure 5A - 50-Hallmark NES heatmap
# Figure 5B - Pairwise cross-dataset concordance
# Figure 5C - Pathway classification summary
# Figure 5D - Context / treatment-state examples
#
# IMPORTANT:
# These are INTERIM D1-D4 integration figures.
# Final manuscript figures will be regenerated after D5/D6.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. PACKAGES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)


# ============================================================
# 1. OUTPUT DIRECTORIES
# ============================================================

figdir <- "figures/integration_D1_D4"

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
# 2. LOAD CLASSIFIED MASTER TABLE
# ============================================================

master <- read.csv(
  "results/integration_D1_D4/D1_D4_Hallmark_classification.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


stopifnot(
  nrow(master) == 50
)


# ============================================================
# 3. CLEAN PATHWAY LABELS
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


master$pathway_label <- clean_hallmark(
  master$pathway
)


# ============================================================
# 4. LOCK CLASS ORDER
# ============================================================

class_order <- c(
  "Strong recurrent",
  "Context-dependent",
  "Directionally recurrent",
  "Dataset-specific significant",
  "Inconsistent / weak"
)


master$class <- factor(
  master$class,
  levels = class_order
)


# ============================================================
# FIGURE 5A
# 50-HALLMARK CROSS-DATASET NES HEATMAP
#
# D1 NES = descriptive mean of two BRCA-context contrasts.
#
# D1 gets a significance marker ONLY when both D1 contrasts
# are significant AND concordant.
# ============================================================


# ------------------------------------------------------------
# 5A.1 Make long NES table
# ------------------------------------------------------------

heatmap_data <- master %>%
  
  select(
    pathway,
    pathway_label,
    class,
    
    D1_NES,
    D2_NES,
    D3_NES,
    D4_NES,
    
    D1_formal_significant,
    D2_significant,
    D3_significant,
    D4_significant
  ) %>%
  
  pivot_longer(
    cols = c(
      D1_NES,
      D2_NES,
      D3_NES,
      D4_NES
    ),
    names_to = "dataset",
    values_to = "NES"
  )


# ------------------------------------------------------------
# 5A.2 Add significance flag
# ------------------------------------------------------------

heatmap_data$significant <- FALSE


heatmap_data$significant[
  heatmap_data$dataset == "D1_NES"
] <-
  master$D1_formal_significant[
    match(
      heatmap_data$pathway[
        heatmap_data$dataset == "D1_NES"
      ],
      master$pathway
    )
  ]


heatmap_data$significant[
  heatmap_data$dataset == "D2_NES"
] <-
  master$D2_significant[
    match(
      heatmap_data$pathway[
        heatmap_data$dataset == "D2_NES"
      ],
      master$pathway
    )
  ]


heatmap_data$significant[
  heatmap_data$dataset == "D3_NES"
] <-
  master$D3_significant[
    match(
      heatmap_data$pathway[
        heatmap_data$dataset == "D3_NES"
      ],
      master$pathway
    )
  ]


heatmap_data$significant[
  heatmap_data$dataset == "D4_NES"
] <-
  master$D4_significant[
    match(
      heatmap_data$pathway[
        heatmap_data$dataset == "D4_NES"
      ],
      master$pathway
    )
  ]


# ------------------------------------------------------------
# 5A.3 Friendly dataset names
# ------------------------------------------------------------

heatmap_data$dataset <- factor(
  heatmap_data$dataset,
  levels = c(
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES"
  ),
  labels = c(
    "D1\nGSE235980",
    "D2\nGSE117765",
    "D3\nGSE153867",
    "D4\nGSE342642"
  )
)


# ------------------------------------------------------------
# 5A.4 Order pathways within classes
#
# Stronger signals appear higher within each category.
# ------------------------------------------------------------

master$mean_abs_NES <- rowMeans(
  abs(
    as.matrix(
      master[
        ,
        c(
          "D1_NES",
          "D2_NES",
          "D3_NES",
          "D4_NES"
        )
      ]
    )
  )
)


pathway_order <- master %>%
  
  arrange(
    class,
    desc(mean_abs_NES)
  ) %>%
  
  pull(
    pathway_label
  )


heatmap_data$pathway_label <- factor(
  heatmap_data$pathway_label,
  levels = rev(
    pathway_order
  )
)


# ------------------------------------------------------------
# 5A.5 Save source data
# ------------------------------------------------------------

write.csv(
  heatmap_data,
  file.path(
    datadir,
    "Figure5A_Hallmark_heatmap_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5A.6 Plot
# ------------------------------------------------------------

p5a <- ggplot(
  heatmap_data,
  aes(
    x = dataset,
    y = pathway_label,
    fill = NES
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        significant,
        "*",
        ""
      )
    ),
    size = 4
  ) +
  
  facet_grid(
    class ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "NES"
  ) +
  
  labs(
    title = "A  Cross-dataset Hallmark enrichment landscape",
    subtitle = paste0(
      "Positive NES = enriched in resistant state; ",
      "negative NES = depleted in resistant state. ",
      "* adjusted P < 0.05."
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 7.5
    ),
    
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 8
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 9
    ),
    
    legend.position = "right"
  )


ggsave(
  file.path(
    figdir,
    "Figure5A_D1_D4_Hallmark_heatmap.png"
  ),
  p5a,
  width = 9,
  height = 15,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure5A_D1_D4_Hallmark_heatmap.pdf"
  ),
  p5a,
  width = 9,
  height = 15
)


# ============================================================
# FIGURE 5B
# PAIRWISE DATASET CONCORDANCE
#
# Cell fill = Spearman correlation of 50 Hallmark NES values
# Text = Spearman rho + directional agreement percentage
# ============================================================


nes_matrix <- master %>%
  
  select(
    D1_NES,
    D2_NES,
    D3_NES,
    D4_NES
  )


dataset_names <- c(
  "D1",
  "D2",
  "D3",
  "D4"
)


rho_matrix <- cor(
  nes_matrix,
  method = "spearman"
)


direction_matrix <- matrix(
  NA_real_,
  nrow = 4,
  ncol = 4
)


for (i in 1:4) {
  
  for (j in 1:4) {
    
    direction_matrix[i, j] <-
      
      100 * mean(
        sign(nes_matrix[[i]]) ==
          sign(nes_matrix[[j]])
      )
  }
}


concordance_data <- expand.grid(
  dataset_x = dataset_names,
  dataset_y = dataset_names,
  stringsAsFactors = FALSE
)


concordance_data$rho <- NA_real_
concordance_data$direction_percent <- NA_real_


for (k in seq_len(nrow(concordance_data))) {
  
  i <- match(
    concordance_data$dataset_x[k],
    dataset_names
  )
  
  j <- match(
    concordance_data$dataset_y[k],
    dataset_names
  )
  
  concordance_data$rho[k] <-
    rho_matrix[i, j]
  
  concordance_data$direction_percent[k] <-
    direction_matrix[i, j]
}


concordance_data$label <- paste0(
  "rho = ",
  sprintf(
    "%.2f",
    concordance_data$rho
  ),
  "\n",
  round(
    concordance_data$direction_percent
  ),
  "% direction"
)


concordance_data$dataset_x <- factor(
  concordance_data$dataset_x,
  levels = dataset_names
)

concordance_data$dataset_y <- factor(
  concordance_data$dataset_y,
  levels = rev(
    dataset_names
  )
)


write.csv(
  concordance_data,
  file.path(
    datadir,
    "Figure5B_pairwise_concordance_data.csv"
  ),
  row.names = FALSE
)


p5b <- ggplot(
  concordance_data,
  aes(
    x = dataset_x,
    y = dataset_y,
    fill = rho
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 1
  ) +
  
  geom_text(
    aes(
      label = label
    ),
    size = 3.7
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Spearman\nrho"
  ) +
  
  coord_equal() +
  
  labs(
    title = "B  Pairwise pathway-profile concordance",
    subtitle = "50 Hallmark gene sets per dataset",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text = element_text(
      face = "bold"
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


ggsave(
  file.path(
    figdir,
    "Figure5B_D1_D4_pairwise_concordance.png"
  ),
  p5b,
  width = 7,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure5B_D1_D4_pairwise_concordance.pdf"
  ),
  p5b,
  width = 7,
  height = 6
)


# ============================================================
# FIGURE 5C
# PATHWAY CLASSIFICATION SUMMARY
# ============================================================


class_data <- master %>%
  
  count(
    class,
    name = "n_pathways"
  )


class_data$class <- factor(
  class_data$class,
  levels = class_order
)


write.csv(
  class_data,
  file.path(
    datadir,
    "Figure5C_pathway_classification_counts.csv"
  ),
  row.names = FALSE
)


p5c <- ggplot(
  class_data,
  aes(
    x = class,
    y = n_pathways,
    fill = class
  )
) +
  
  geom_col(
    width = 0.72,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(
      label = n_pathways
    ),
    vjust = -0.4,
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  
  labs(
    title = "C  Cross-dataset pathway classification",
    subtitle = "Classification integrates direction and statistical support",
    x = NULL,
    y = "Number of Hallmark pathways"
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 14
    )
  )


ggsave(
  file.path(
    figdir,
    "Figure5C_D1_D4_classification_summary.png"
  ),
  p5c,
  width = 8,
  height = 5.5,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure5C_D1_D4_classification_summary.pdf"
  ),
  p5c,
  width = 8,
  height = 5.5
)
# ============================================================
# FIGURE 5D - REVISED
# ============================================================


# ------------------------------------------------------------
# Fixed x positions for individual datasets
#
# Prevents random jitter and ensures significance markers
# remain aligned with their corresponding points.
# ------------------------------------------------------------

context_data$x_position <- NA_real_

context_data$x_position[
  context_data$dataset == "D1"
] <- 0.90

context_data$x_position[
  context_data$dataset == "D4"
] <- 1.10

context_data$x_position[
  context_data$dataset == "D2"
] <- 1.90

context_data$x_position[
  context_data$dataset == "D3"
] <- 2.10


# ------------------------------------------------------------
# Exposure-group means
#
# Descriptive only - no connecting line because these are
# independent datasets, not paired longitudinal measurements.
# ------------------------------------------------------------

mean_data <- context_data %>%
  
  group_by(
    pathway_label,
    exposure_group
  ) %>%
  
  summarise(
    mean_NES = mean(NES),
    .groups = "drop"
  )


mean_data$x_position <- ifelse(
  mean_data$exposure_group ==
    "Drug-free / untreated",
  1,
  2
)


# ------------------------------------------------------------
# Save revised plotting data
# ------------------------------------------------------------

write.csv(
  context_data,
  file.path(
    datadir,
    "Figure5D_context_pathways_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Revised plot
# ------------------------------------------------------------

p5d <- ggplot() +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  # Individual datasets
  geom_point(
    data = context_data,
    aes(
      x = x_position,
      y = NES,
      shape = dataset,
      colour = dataset
    ),
    size = 3.3
  ) +
  
  # Significance markers aligned with points
  geom_text(
    data = context_data,
    aes(
      x = x_position,
      y = NES,
      label = ifelse(
        significant,
        "*",
        ""
      )
    ),
    nudge_x = 0.055,
    nudge_y = 0.10,
    size = 4,
    fontface = "bold"
  ) +
  
  # Descriptive exposure-group mean
  geom_point(
    data = mean_data,
    aes(
      x = x_position,
      y = mean_NES
    ),
    shape = 18,
    size = 4.5,
    colour = "black"
  ) +
  
  facet_wrap(
    ~ pathway_label,
    ncol = 2
  ) +
  
  scale_x_continuous(
    breaks = c(
      1,
      2
    ),
    labels = c(
      "Drug-free / untreated",
      "Ongoing olaparib"
    ),
    limits = c(
      0.72,
      2.28
    )
  ) +
  
  labs(
    title =
      "D  Context dependence of selected resistance-associated programs",
    
    subtitle = paste0(
      "Points = individual datasets; black diamonds = descriptive ",
      "exposure-group means. * adjusted P < 0.05."
    ),
    
    x = NULL,
    
    y = "Normalized enrichment score (NES)",
    
    shape = "Dataset",
    
    colour = "Dataset"
  ) +
  
  theme_minimal(
    base_size = 10.5
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    
    axis.text.x = element_text(
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
    "Figure5D_D1_D4_context_dependence.png"
  ),
  p5d,
  width = 10,
  height = 9,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure5D_D1_D4_context_dependence.pdf"
  ),
  p5d,
  width = 10,
  height = 9
)

 
# ------------------------------------------------------------
# Add significance
# ------------------------------------------------------------

context_data$significant <- FALSE


context_data$significant[
  context_data$dataset == "D1"
] <-
  master$D1_formal_significant[
    match(
      context_data$pathway[
        context_data$dataset == "D1"
      ],
      master$pathway
    )
  ]


context_data$significant[
  context_data$dataset == "D2"
] <-
  master$D2_significant[
    match(
      context_data$pathway[
        context_data$dataset == "D2"
      ],
      master$pathway
    )
  ]


context_data$significant[
  context_data$dataset == "D3"
] <-
  master$D3_significant[
    match(
      context_data$pathway[
        context_data$dataset == "D3"
      ],
      master$pathway
    )
  ]


context_data$significant[
  context_data$dataset == "D4"
] <-
  master$D4_significant[
    match(
      context_data$pathway[
        context_data$dataset == "D4"
      ],
      master$pathway
    )
  ]


# ------------------------------------------------------------
# Clean facet labels
# ------------------------------------------------------------

context_labels <- c(
  "TNFA SIGNALING VIA NFKB",
  "E2F TARGETS",
  "G2M CHECKPOINT",
  "MYC TARGETS V1",
  "OXIDATIVE PHOSPHORYLATION",
  "EPITHELIAL MESENCHYMAL TRANSITION"
)


context_data$pathway_label <- factor(
  context_data$pathway_label,
  levels = context_labels
)


write.csv(
  context_data,
  file.path(
    datadir,
    "Figure5D_context_pathways_data.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Plot individual datasets + exposure-group mean
# ------------------------------------------------------------

p5d <- ggplot(
  context_data,
  aes(
    x = exposure_group,
    y = NES
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  stat_summary(
    aes(
      group = 1
    ),
    fun = mean,
    geom = "line",
    linewidth = 0.7,
    colour = "black"
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4,
    colour = "black"
  ) +
  
  geom_point(
    aes(
      shape = dataset,
      colour = dataset
    ),
    size = 3.2,
    position = position_jitter(
      width = 0.06,
      height = 0
    )
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        significant,
        "*",
        ""
      )
    ),
    position = position_jitter(
      width = 0.06,
      height = 0
    ),
    nudge_y = 0.15,
    size = 4
  ) +
  
  facet_wrap(
    ~ pathway_label,
    ncol = 2
  ) +
  
  labs(
    title = "D  Context dependence of selected resistance-associated programs",
    subtitle = paste0(
      "Points = individual datasets; diamonds/lines = exposure-group means. ",
      "* adjusted P < 0.05."
    ),
    x = NULL,
    y = "Normalized enrichment score (NES)",
    shape = "Dataset",
    colour = "Dataset"
  ) +
  
  theme_minimal(
    base_size = 10.5
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    
    axis.text.x = element_text(
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
    "Figure5D_D1_D4_context_dependence.png"
  ),
  p5d,
  width = 10,
  height = 9,
  dpi = 300
)


ggsave(
  file.path(
    figdir,
    "Figure5D_D1_D4_context_dependence.pdf"
  ),
  p5d,
  width = 10,
  height = 9
)


# ============================================================
# OPTIONAL COMPOSITE FIGURE
#
# Only generated if patchwork is already installed.
# No automatic package installation.
# ============================================================

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {
  
  library(patchwork)
  
  combined <- (
    p5a |
      (
        p5b /
          p5c /
          p5d
      )
  ) +
    
    plot_layout(
      widths = c(
        1.25,
        1
      )
    ) +
    
    plot_annotation(
      title = paste0(
        "Figure 5. Cross-dataset transcriptomic integration ",
        "of acquired PARP-inhibitor resistance"
      )
    )
  
  
  ggsave(
    file.path(
      figdir,
      "Figure5_D1_D4_integration_combined.png"
    ),
    combined,
    width = 18,
    height = 16,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      figdir,
      "Figure5_D1_D4_integration_combined.pdf"
    ),
    combined,
    width = 18,
    height = 16
  )
  
  
  cat(
    "\nCombined Figure 5 generated using patchwork.\n"
  )
  
} else {
  
  cat(
    "\npatchwork not installed - individual panels were still generated successfully.\n"
  )
  
}


# ============================================================
# FINAL QC SUMMARY
# ============================================================

cat("\n========================================\n")
cat("INTEGRATION FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Figure 5A: 50-Hallmark NES heatmap\n"
)

cat(
  "Figure 5B: Pairwise concordance matrix\n"
)

cat(
  "Figure 5C: Pathway classification summary\n"
)

cat(
  "Figure 5D: Selected context-dependent programs\n"
)

cat(
  "\nOutput directory:\n",
  figdir,
  "\n"
)

cat(
  "\nFigure-data directory:\n",
  datadir,
  "\n"
)

cat("\n========================================\n")