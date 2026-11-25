# ============================================================
# GSE235980 - REFINED FIGURES
# Cleaner final version before moving to Dataset 2
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(DESeq2)
library(ggplot2)
library(ggrepel)

dir.create(
  "figures/GSE235980_refined",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# COMMON INPUTS
# ============================================================

counts_raw <- read.delim(
  "data_raw/GSE235980_CountReads.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- counts_raw$Geneid
storage.mode(counts) <- "integer"

keep <- rowSums(counts >= 10) >= 2
counts_filtered <- counts[keep, ]

metadata <- data.frame(
  BRCA1_status = c(
    "Deficient", "Deficient", "Deficient", "Deficient",
    "Proficient", "Proficient", "Proficient", "Proficient"
  ),
  resistance = c(
    "Resistant", "Resistant", "Sensitive", "Sensitive",
    "Resistant", "Resistant", "Sensitive", "Sensitive"
  ),
  replicate = c(1, 2, 1, 2, 1, 2, 1, 2),
  row.names = colnames(counts_filtered)
)

metadata$group <- factor(
  paste(metadata$BRCA1_status, metadata$resistance, sep = "_")
)

metadata$sample_label <- c(
  "Def-R1", "Def-R2", "Def-S1", "Def-S2",
  "Prof-R1", "Prof-R2", "Prof-S1", "Prof-S2"
)

stopifnot(
  identical(rownames(metadata), colnames(counts_filtered))
)

res_def <- read.csv(
  "results/differential_expression/BRCA1_deficient_all_DESeq2_results.csv",
  row.names = 1,
  check.names = FALSE
)

res_prof <- read.csv(
  "results/differential_expression/BRCA1_proficient_all_DESeq2_results.csv",
  row.names = 1,
  check.names = FALSE
)

shared_results <- read.csv(
  "results/differential_expression/shared_DEGs_comparison.csv"
)

hallmark <- read.csv(
  "results/pathway_analysis/Hallmark_GSEA_shared_same_direction.csv"
)


# ============================================================
# FIGURE 1A - PCA (cleaner labels)
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ group
)

dds <- estimateSizeFactors(dds)

vsd <- vst(
  dds,
  blind = TRUE
)

pca_data <- plotPCA(
  vsd,
  intgroup = c("BRCA1_status", "resistance"),
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar")
)

pca_data$sample <- rownames(pca_data)
pca_data$sample_label <- metadata[pca_data$sample, "sample_label"]

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = BRCA1_status,
    shape = resistance
  )
) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = sample_label),
    size = 3.8,
    show.legend = FALSE,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  labs(
    title = "Transcriptomic structure of GSE235980",
    subtitle = "Olaparib-sensitive and resistant cells across BRCA1 backgrounds",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "BRCA1 status",
    shape = "Olaparib phenotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 10)
  )

print(p_pca)

ggsave(
  "figures/GSE235980_refined/Figure_1A_PCA_refined.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE235980_refined/Figure_1A_PCA_refined.pdf",
  p_pca,
  width = 8,
  height = 6
)


# ============================================================
# FIGURE 1B - CONCORDANCE (grey background + up/down highlights)
# ============================================================

all_changes <- data.frame(
  gene = rownames(res_def),
  log2FC_deficient = res_def$log2FoldChange,
  log2FC_proficient = res_prof$log2FoldChange
)

all_changes <- all_changes[
  complete.cases(
    all_changes$log2FC_deficient,
    all_changes$log2FC_proficient
  ),
]

shared_consistent <- subset(
  shared_results,
  same_direction == TRUE
)

up_consistent <- shared_consistent$gene[
  shared_consistent$log2FC_deficient > 0
]

down_consistent <- shared_consistent$gene[
  shared_consistent$log2FC_deficient < 0
]

all_changes$category <- "Other gene"
all_changes$category[all_changes$gene %in% up_consistent] <-
  "Shared consistent up"
all_changes$category[all_changes$gene %in% down_consistent] <-
  "Shared consistent down"

all_changes$category <- factor(
  all_changes$category,
  levels = c(
    "Other gene",
    "Shared consistent up",
    "Shared consistent down"
  )
)

rho <- cor(
  all_changes$log2FC_deficient,
  all_changes$log2FC_proficient,
  method = "spearman"
)

p_concordance <- ggplot(
  all_changes,
  aes(
    x = log2FC_deficient,
    y = log2FC_proficient,
    color = category
  )
) +
  geom_point(
    data = subset(all_changes, category == "Other gene"),
    alpha = 0.25,
    size = 0.7
  ) +
  geom_point(
    data = subset(all_changes, category != "Other gene"),
    alpha = 0.7,
    size = 1.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.7
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
      "\nShared consistent DEGs = 692"
    ),
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Other gene" = "grey75",
      "Shared consistent up" = "#D55E00",
      "Shared consistent down" = "#0072B2"
    )
  ) +
  labs(
    title = "Concordance of olaparib resistance-associated expression changes",
    subtitle = "Shared directionally consistent DEGs highlighted by regulation direction",
    x = "log2 fold change: BRCA1-deficient",
    y = "log2 fold change: BRCA1-proficient",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_concordance)

ggsave(
  "figures/GSE235980_refined/Figure_1B_log2FC_concordance_refined.png",
  p_concordance,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE235980_refined/Figure_1B_log2FC_concordance_refined.pdf",
  p_concordance,
  width = 8,
  height = 6
)


# ============================================================
# FIGURE 1C - SHARED DEG SUMMARY (better order + percentages)
# ============================================================

n_up <- sum(shared_consistent$log2FC_deficient > 0)
n_down <- sum(shared_consistent$log2FC_deficient < 0)
n_opposite <- sum(shared_results$same_direction == FALSE)

deg_summary <- data.frame(
  category = factor(
    c(
      "Consistently upregulated",
      "Consistently downregulated",
      "Opposite direction"
    ),
    levels = c(
      "Consistently upregulated",
      "Consistently downregulated",
      "Opposite direction"
    )
  ),
  genes = c(n_up, n_down, n_opposite)
)

total_shared <- nrow(shared_results)
consistent_total <- n_up + n_down
consistent_pct <- round(100 * consistent_total / total_shared, 1)

p_deg <- ggplot(
  deg_summary,
  aes(
    x = category,
    y = genes,
    fill = category
  )
) +
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = genes),
    vjust = -0.4,
    size = 4.5
  ) +
  labs(
    title = "Shared resistance-associated DEGs",
    subtitle = paste0(
      consistent_total, "/", total_shared,
      " shared DEGs (", consistent_pct,
      "%) changed in the same direction"
    ),
    x = NULL,
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 18, hjust = 1),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_deg)

ggsave(
  "figures/GSE235980_refined/Figure_1C_shared_DEGs_refined.png",
  p_deg,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

ggsave(
  "figures/GSE235980_refined/Figure_1C_shared_DEGs_refined.pdf",
  p_deg,
  width = 7.5,
  height = 5.5
)


# ============================================================
# FIGURE 2 - SHARED HALLMARK GSEA (cleaner margins/text)
# ============================================================

hallmark$label <- gsub("^HALLMARK_", "", hallmark$Hallmark)
hallmark$label <- gsub("_", " ", hallmark$label)

hallmark <- hallmark[
  order(hallmark$mean_abs_NES, decreasing = FALSE),
]

hallmark$label <- factor(
  hallmark$label,
  levels = hallmark$label
)

hallmark_long <- rbind(
  data.frame(
    pathway = hallmark$label,
    background = "BRCA1-deficient",
    NES = hallmark$NES_deficient
  ),
  data.frame(
    pathway = hallmark$label,
    background = "BRCA1-proficient",
    NES = hallmark$NES_proficient
  )
)

p_hallmark <- ggplot(
  hallmark_long,
  aes(
    x = pathway,
    y = NES,
    fill = background
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "Shared Hallmark programs associated with olaparib resistance",
    subtitle = "Positive NES = enriched in resistant cells",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = "BRCA1 background"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(18, 25, 12, 12)
  )

print(p_hallmark)

ggsave(
  "figures/GSE235980_refined/Figure_2_Hallmark_GSEA_refined.png",
  p_hallmark,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/GSE235980_refined/Figure_2_Hallmark_GSEA_refined.pdf",
  p_hallmark,
  width = 10,
  height = 7
)


# ============================================================
# DONE
# ============================================================

cat("\n===========================================\n")
cat("REFINED GSE235980 FIGURES GENERATED\n")
cat("===========================================\n")
cat("\nSaved in: figures/GSE235980_refined/\n")