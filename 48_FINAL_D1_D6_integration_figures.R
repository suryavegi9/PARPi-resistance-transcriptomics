# ============================================================
# FINAL D1-D6 CROSS-DATASET SYNTHESIS FIGURE
#
# Figure 8A
#   Six-accession Hallmark NES landscape
#
# Figure 8B
#   Accession-level vs model-level classification counts
#
# Figure 8C
#   Observational treatment-state signal vs controlled D5
#   exposure signal
#
# Figure 8D
#   Significant UP/DOWN support for the original frozen
#   D1 pathways under model-level sensitivity integration
#
# Central manuscript message:
#
# Acquired PARPi resistance is dominated by context-dependent
# transcriptional remodeling, with limited recurrent programs.
# Glycolysis suppression is the strongest recurrent signal.
# Treatment state explains a subset, but not the majority,
# of cross-dataset heterogeneity.
# ============================================================


setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ============================================================
# 0. PACKAGES
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)


# ============================================================
# 1. FILES
# ============================================================

master_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "FINAL_D1_D6_Hallmark_master.csv"
)


accession_counts_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "FINAL_D1_D6_classification_counts.csv"
)


model_counts_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "model_level_sensitivity/",
  "FINAL_model_level_classification_counts.csv"
)


treatment_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "treatment_state/",
  "FINAL_treatment_state_master.csv"
)


frozen_file <- paste0(
  "results/FINAL_D1_D6_integration/",
  "model_level_sensitivity/",
  "FINAL_frozen_pathways_accession_vs_model.csv"
)


required_files <- c(
  master_file,
  accession_counts_file,
  model_counts_file,
  treatment_file,
  frozen_file
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  
  print(missing_files)
  
  stop(
    "Required integration files are missing."
  )
}


# ============================================================
# 2. OUTPUT
# ============================================================

figdir <- "figures/FINAL_D1_D6_integration"

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
# 3. HELPERS
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
# FIGURE 8A
# SIX-ACCESSION HALLMARK LANDSCAPE
# ============================================================
# ============================================================


master <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_master <- c(
  
  "pathway",
  
  "D1_NES",
  "D2_NES",
  "D3_NES",
  "D4_NES",
  "D5_NES",
  "D6_NES",
  
  "D1_sig_direction",
  "D2_sig_direction",
  "D3_sig_direction",
  "D4_sig_direction",
  "D5_sig_direction",
  "D6_sig_direction",
  
  "classification",
  
  "n_UP",
  "n_DOWN",
  "n_sig_UP",
  "n_sig_DOWN"
)


stopifnot(
  all(
    required_master %in%
      colnames(master)
  )
)


# ============================================================
# 8A.1 ORDER PATHWAYS
#
# Put the most scientifically important categories first.
# ============================================================

class_order <- c(
  
  "Strong recurrent",
  
  "Context-dependent",
  
  "Directionally recurrent",
  
  "Dataset-specific significant",
  
  "Inconsistent/weak"
)


master$class_rank <- match(
  master$classification,
  class_order
)


master$max_direction_support <- pmax(
  master$n_UP,
  master$n_DOWN
)


master$max_sig_support <- pmax(
  master$n_sig_UP,
  master$n_sig_DOWN
)


master <- master %>%
  
  arrange(
    class_rank,
    desc(max_sig_support),
    desc(max_direction_support),
    pathway
  )


master$pathway_label <- clean_hallmark(
  master$pathway
)


pathway_order <- rev(
  master$pathway_label
)


# ============================================================
# 8A.2 BUILD LONG NES TABLE
# ============================================================

nes_long <- bind_rows(
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D1",
    NES = master$D1_NES,
    significant = master$D1_sig_direction != 0,
    classification = master$classification
  ),
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D2",
    NES = master$D2_NES,
    significant = master$D2_sig_direction != 0,
    classification = master$classification
  ),
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D3",
    NES = master$D3_NES,
    significant = master$D3_sig_direction != 0,
    classification = master$classification
  ),
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D4",
    NES = master$D4_NES,
    significant = master$D4_sig_direction != 0,
    classification = master$classification
  ),
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D5",
    NES = master$D5_NES,
    significant = master$D5_sig_direction != 0,
    classification = master$classification
  ),
  
  data.frame(
    pathway = master$pathway,
    pathway_label = master$pathway_label,
    dataset = "D6",
    NES = master$D6_NES,
    significant = master$D6_sig_direction != 0,
    classification = master$classification
  )
)


nes_long$dataset <- factor(
  nes_long$dataset,
  levels = paste0(
    "D",
    1:6
  )
)


nes_long$pathway_label <- factor(
  nes_long$pathway_label,
  levels = pathway_order
)


max_abs_nes <- max(
  abs(
    nes_long$NES
  ),
  na.rm = TRUE
)


write.csv(
  nes_long,
  file.path(
    datadir,
    "Figure8A_accession_Hallmark_landscape.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 8A.3 HEATMAP
# ============================================================

p8a <- ggplot(
  
  nes_long,
  
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
    
    size = 2.5,
    fontface = "bold"
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
      "A  Cross-dataset Hallmark landscape of acquired PARPi resistance",
    
    subtitle =
      paste0(
        "One formal vote per GEO accession; ",
        "* indicates accession-level significance"
      ),
    
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10
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
        size = 6.6
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5
      ),
    
    legend.title =
      element_text(
        face = "bold"
      ),
    
    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


ggsave(
  file.path(
    figdir,
    "Figure8A_final_Hallmark_landscape.png"
  ),
  p8a,
  width = 9,
  height = 13,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure8A_final_Hallmark_landscape.pdf"
  ),
  p8a,
  width = 9,
  height = 13
)


# ============================================================
# ============================================================
# FIGURE 8B
# ACCESSION vs MODEL-LEVEL CLASSIFICATION
# ============================================================
# ============================================================


acc_counts <- read.csv(
  accession_counts_file,
  stringsAsFactors = FALSE
)


model_counts <- read.csv(
  model_counts_file,
  stringsAsFactors = FALSE
)


acc_counts$analysis <-
  "Accession-level primary"


model_counts$analysis <-
  "Model-level sensitivity"


# Normalize category naming

acc_counts$display_class <- acc_counts$classification

model_counts$display_class <- model_counts$classification


acc_counts$display_class[
  acc_counts$display_class ==
    "Dataset-specific significant"
] <- "Specific significant"


model_counts$display_class[
  model_counts$display_class ==
    "Unit-specific significant"
] <- "Specific significant"


classification_plot <- bind_rows(
  
  acc_counts[
    ,
    c(
      "display_class",
      "n",
      "analysis"
    )
  ],
  
  model_counts[
    ,
    c(
      "display_class",
      "n",
      "analysis"
    )
  ]
)


class_levels <- c(
  
  "Strong recurrent",
  
  "Directionally recurrent",
  
  "Context-dependent",
  
  "Specific significant",
  
  "Inconsistent/weak"
)


classification_plot$display_class <- factor(
  
  classification_plot$display_class,
  
  levels = class_levels
)


classification_plot$analysis <- factor(
  
  classification_plot$analysis,
  
  levels = c(
    "Accession-level primary",
    "Model-level sensitivity"
  )
)


write.csv(
  classification_plot,
  file.path(
    datadir,
    "Figure8B_classification_comparison.csv"
  ),
  row.names = FALSE
)


p8b <- ggplot(
  
  classification_plot,
  
  aes(
    x = display_class,
    y = n,
    fill = analysis
  )
) +
  
  geom_col(
    
    position = position_dodge(
      width = 0.75
    ),
    
    width = 0.68
  ) +
  
  geom_text(
    
    aes(
      label = n
    ),
    
    position = position_dodge(
      width = 0.75
    ),
    
    vjust = -0.4,
    
    size = 3.5
  ) +
  
  scale_fill_manual(
    
    values = c(
      "Accession-level primary" =
        "#0072B2",
      
      "Model-level sensitivity" =
        "#D55E00"
    )
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
    
    title =
      "B  Context dependence persists across integration strategies",
    
    subtitle =
      paste0(
        "All 15 accession-level context-dependent Hallmarks ",
        "remain context-dependent at model level"
      ),
    
    x = NULL,
    
    y =
      "Number of Hallmark pathways",
    
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      ),
    
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
    "Figure8B_integration_strategy_comparison.png"
  ),
  p8b,
  width = 8.5,
  height = 6.5,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure8B_integration_strategy_comparison.pdf"
  ),
  p8b,
  width = 8.5,
  height = 6.5
)


# ============================================================
# ============================================================
# FIGURE 8C
# TREATMENT-STATE EVIDENCE
# ============================================================
# ============================================================


treatment <- read.csv(
  treatment_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_treatment <- c(
  
  "pathway",
  
  "observational_exposure_delta",
  
  "D5_exposure_mean",
  
  "exposure_evidence"
)


stopifnot(
  all(
    required_treatment %in%
      colnames(treatment)
  )
)


treatment$label <- clean_hallmark(
  treatment$pathway
)


# ============================================================
# Selected labels:
#
# Strong convergent pathways +
# major biologically important discordant programs.
# ============================================================

label_pathways <- c(
  
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  
  "HALLMARK_INFLAMMATORY_RESPONSE",
  
  "HALLMARK_COMPLEMENT",
  
  "HALLMARK_APOPTOSIS",
  
  "HALLMARK_E2F_TARGETS",
  
  "HALLMARK_G2M_CHECKPOINT",
  
  "HALLMARK_MYC_TARGETS_V1",
  
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  
  "HALLMARK_MTORC1_SIGNALING"
)


label_treatment <- treatment[
  treatment$pathway %in%
    label_pathways,
]


evidence_levels <- c(
  
  "Strong convergent",
  
  "Convergent",
  
  "D5-specific exposure",
  
  "Observational-only",
  
  "Weak/mixed"
)


treatment$exposure_evidence <- factor(
  
  treatment$exposure_evidence,
  
  levels = evidence_levels
)


write.csv(
  treatment,
  file.path(
    datadir,
    "Figure8C_treatment_state_evidence.csv"
  ),
  row.names = FALSE
)


p8c <- ggplot(
  
  treatment,
  
  aes(
    x = observational_exposure_delta,
    y = D5_exposure_mean
  )
) +
  
  geom_hline(
    
    yintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.45,
    
    colour = "grey50"
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.45,
    
    colour = "grey50"
  ) +
  
  geom_point(
    
    aes(
      colour = exposure_evidence
    ),
    
    size = 3.5,
    alpha = 0.95
  ) +
  
  geom_text_repel(
    
    data = label_treatment,
    
    aes(
      label = label
    ),
    
    size = 2.8,
    
    box.padding = 0.45,
    
    point.padding = 0.3,
    
    min.segment.length = 0,
    
    max.overlaps = Inf,
    
    seed = 123,
    
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    
    values = c(
      
      "Strong convergent" =
        "#009E73",
      
      "Convergent" =
        "#56B4E9",
      
      "D5-specific exposure" =
        "#E69F00",
      
      "Observational-only" =
        "#CC79A7",
      
      "Weak/mixed" =
        "grey70"
    )
  ) +
  
  labs(
    
    title =
      "C  Treatment state explains only a subset of pathway heterogeneity",
    
    subtitle =
      paste0(
        "Across all 50 Hallmarks: rho = -0.064; ",
        "52% directional agreement"
      ),
    
    x =
      paste0(
        "Observational ON vs OFF exposure difference\n",
        "(ON mean NES - OFF mean NES)"
      ),
    
    y =
      paste0(
        "Controlled D5 exposure effect\n",
        "(mean H3/P3 NES)"
      ),
    
    colour =
      "Exposure evidence"
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
    "Figure8C_treatment_state_integration.png"
  ),
  p8c,
  width = 9,
  height = 7.5,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure8C_treatment_state_integration.pdf"
  ),
  p8c,
  width = 9,
  height = 7.5
)


# ============================================================
# ============================================================
# FIGURE 8D
# FROZEN D1 PATHWAYS:
# MODEL-LEVEL SIGNIFICANT UP/DOWN SUPPORT
# ============================================================
# ============================================================


frozen <- read.csv(
  frozen_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_frozen <- c(
  
  "pathway",
  
  "n_sig_UP",
  
  "n_sig_DOWN",
  
  "model_level_classification",
  
  "accession_level_classification"
)


stopifnot(
  all(
    required_frozen %in%
      colnames(frozen)
  )
)


frozen$pathway_label <- clean_hallmark(
  frozen$pathway
)


# ============================================================
# Display order
# ============================================================

frozen_order <- c(
  
  "GLYCOLYSIS",
  
  "MTORC1 SIGNALING",
  
  "INFLAMMATORY RESPONSE",
  
  "TNFA SIGNALING VIA NFKB",
  
  "E2F TARGETS",
  
  "G2M CHECKPOINT",
  
  "MYC TARGETS V1",
  
  "EPITHELIAL MESENCHYMAL TRANSITION",
  
  "OXIDATIVE PHOSPHORYLATION",
  
  "UV RESPONSE DN"
)


frozen$pathway_label <- factor(
  
  frozen$pathway_label,
  
  levels = rev(
    frozen_order
  )
)


# ============================================================
# Diverging bar table
#
# UP support = positive
# DOWN support = negative
# ============================================================

frozen_long <- bind_rows(
  
  data.frame(
    
    pathway_label =
      frozen$pathway_label,
    
    direction =
      "Significant UP",
    
    support =
      frozen$n_sig_UP,
    
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    
    pathway_label =
      frozen$pathway_label,
    
    direction =
      "Significant DOWN",
    
    support =
      -frozen$n_sig_DOWN,
    
    stringsAsFactors = FALSE
  )
)


write.csv(
  frozen_long,
  file.path(
    datadir,
    "Figure8D_frozen_pathway_vote_balance.csv"
  ),
  row.names = FALSE
)


p8d <- ggplot(
  
  frozen_long,
  
  aes(
    x = support,
    y = pathway_label,
    fill = direction
  )
) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linewidth = 0.6,
    
    colour = "grey30"
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "Significant UP" =
        "#B2182B",
      
      "Significant DOWN" =
        "#2166AC"
    )
  ) +
  
  scale_x_continuous(
    
    breaks = -5:5,
    
    labels = abs(
      -5:5
    )
  ) +
  
  labs(
    
    title =
      "D  Original candidate pathways resolve into recurrence or context dependence",
    
    subtitle =
      paste0(
        "Model-level sensitivity analysis; ",
        "bar length = number of significant supporting models"
      ),
    
    x =
      "Number of significant model-level observations",
    
    y = NULL,
    
    fill = NULL
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
    
    axis.text.y =
      element_text(
        size = 9
      ),
    
    legend.position =
      "bottom"
  )


ggsave(
  file.path(
    figdir,
    "Figure8D_frozen_pathway_vote_balance.png"
  ),
  p8d,
  width = 9,
  height = 7,
  dpi = 400
)


ggsave(
  file.path(
    figdir,
    "Figure8D_frozen_pathway_vote_balance.pdf"
  ),
  p8d,
  width = 9,
  height = 7
)


# ============================================================
# ============================================================
# COMBINED FINAL FIGURE 8
# ============================================================
# ============================================================

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {
  
  library(patchwork)
  
  
  # Heatmap gets more vertical space.
  
  right_column <- (
    
    p8b /
      p8c /
      p8d
  ) +
    
    plot_layout(
      heights = c(
        0.85,
        1.05,
        1
      )
    )
  
  
  combined <- (
    
    p8a |
      right_column
    
  ) +
    
    plot_layout(
      widths = c(
        1.05,
        1
      )
    ) +
    
    plot_annotation(
      
      title =
        paste0(
          "Figure 8. Acquired PARP-inhibitor resistance is dominated ",
          "by context-dependent transcriptional remodeling"
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
      "Figure8_FINAL_D1_D6_synthesis.png"
    ),
    
    combined,
    
    width = 20,
    
    height = 16,
    
    dpi = 400
  )
  
  
  ggsave(
    
    file.path(
      figdir,
      "Figure8_FINAL_D1_D6_synthesis.pdf"
    ),
    
    combined,
    
    width = 20,
    
    height = 16
  )
}


# ============================================================
# FINAL QC
# ============================================================

cat("\n========================================\n")
cat("FINAL D1-D6 SYNTHESIS FIGURES COMPLETE\n")
cat("========================================\n\n")

cat(
  "Figure 8A - accession-level 50-Hallmark landscape\n"
)

cat(
  "Figure 8B - accession vs model-level classifications\n"
)

cat(
  "Figure 8C - treatment-state integration\n"
)

cat(
  "Figure 8D - frozen pathway significant-vote balance\n"
)

cat(
  "\nOutput directory:\n",
  figdir,
  "\n"
)

cat("\n========================================\n")