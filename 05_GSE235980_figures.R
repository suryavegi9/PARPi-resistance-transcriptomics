# ============================================================
# GSE235980 - Publication-Quality Figures
# Independent PARPi Resistance Transcriptomics Project
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

library(DESeq2)
library(ggplot2)

# Create dedicated output folder
dir.create(
  "figures/GSE235980",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# FIGURE 1A - PCA
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

# Low-count filtering
keep <- rowSums(counts >= 10) >= 2

counts_filtered <- counts[keep, ]


# Metadata
metadata <- data.frame(
  
  BRCA1_status = c(
    "Deficient",
    "Deficient",
    "Deficient",
    "Deficient",
    "Proficient",
    "Proficient",
    "Proficient",
    "Proficient"
  ),
  
  resistance = c(
    "Resistant",
    "Resistant",
    "Sensitive",
    "Sensitive",
    "Resistant",
    "Resistant",
    "Sensitive",
    "Sensitive"
  ),
  
  replicate = c(
    1, 2, 1, 2,
    1, 2, 1, 2
  ),
  
  row.names = colnames(counts_filtered)
)

metadata$group <- factor(
  paste(
    metadata$BRCA1_status,
    metadata$resistance,
    sep = "_"
  )
)

stopifnot(
  identical(
    rownames(metadata),
    colnames(counts_filtered)
  )
)

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
  intgroup = c(
    "BRCA1_status",
    "resistance"
  ),
  returnData = TRUE
)

percent_var <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)

pca_data$sample <- rownames(pca_data)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = BRCA1_status,
    shape = resistance
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text(
    aes(label = sample),
    vjust = -0.9,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "Transcriptomic structure of GSE235980",
    subtitle = "Olaparib-sensitive and resistant cells across BRCA1 backgrounds",
    x = paste0(
      "PC1 (",
      percent_var[1],
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      "% variance)"
    ),
    color = "BRCA1 status",
    shape = "Olaparib phenotype"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

print(p_pca)

ggsave(
  "figures/GSE235980/Figure_1A_PCA.png",
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE235980/Figure_1A_PCA.pdf",
  p_pca,
  width = 8,
  height = 6
)


# ============================================================
# FIGURE 1B - GENOME-WIDE LOG2FC CONCORDANCE
# ============================================================

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

shared_results <- read.csv(
  "results/differential_expression/shared_DEGs_comparison.csv"
)

consistent_genes <- shared_results$gene[
  shared_results$same_direction == TRUE
]

all_changes$category <-
  ifelse(
    all_changes$gene %in% consistent_genes,
    "Shared directionally consistent DEG",
    "Other gene"
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
    alpha = 0.35,
    size = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -1,
    label = paste0(
      "Spearman rho = ",
      round(rho, 3)
    ),
    size = 4
  ) +
  labs(
    title = "Concordance of olaparib resistance-associated expression changes",
    subtitle = "692 shared directionally consistent DEGs highlighted",
    x = "log2 fold change: BRCA1-deficient",
    y = "log2 fold change: BRCA1-proficient",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

print(p_concordance)

ggsave(
  "figures/GSE235980/Figure_1B_log2FC_concordance.png",
  p_concordance,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/GSE235980/Figure_1B_log2FC_concordance.pdf",
  p_concordance,
  width = 7,
  height = 6
)


# ============================================================
# FIGURE 1C - SHARED DEG SUMMARY
# ============================================================

consistent <- subset(
  shared_results,
  same_direction == TRUE
)

n_up <- sum(
  consistent$log2FC_deficient > 0
)

n_down <- sum(
  consistent$log2FC_deficient < 0
)

n_opposite <- sum(
  shared_results$same_direction == FALSE
)

deg_summary <- data.frame(
  category = c(
    "Consistently upregulated",
    "Consistently downregulated",
    "Opposite direction"
  ),
  genes = c(
    n_up,
    n_down,
    n_opposite
  )
)

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
    size = 4
  ) +
  labs(
    title = "Shared resistance-associated DEGs",
    subtitle = "Genes significant in both BRCA1 backgrounds",
    x = NULL,
    y = "Number of genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

print(p_deg)

ggsave(
  "figures/GSE235980/Figure_1C_shared_DEGs.png",
  p_deg,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/GSE235980/Figure_1C_shared_DEGs.pdf",
  p_deg,
  width = 7,
  height = 5
)


# ============================================================
# FIGURE 2 - SHARED HALLMARK GSEA
# ============================================================

hallmark <- read.csv(
  "results/pathway_analysis/Hallmark_GSEA_shared_same_direction.csv"
)

# Clean Hallmark names
hallmark$label <- gsub(
  "^HALLMARK_",
  "",
  hallmark$Hallmark
)

hallmark$label <- gsub(
  "_",
  " ",
  hallmark$label
)

# Order by mean absolute NES
hallmark <- hallmark[
  order(
    hallmark$mean_abs_NES,
    decreasing = FALSE
  ),
]

hallmark$label <- factor(
  hallmark$label,
  levels = hallmark$label
)


# Convert wide table to long format
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
    position = position_dodge(
      width = 0.75
    ),
    width = 0.65
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  labs(
    title = "Shared Hallmark programs associated with olaparib resistance",
    subtitle = "Positive NES = enriched in resistant cells; negative NES = depleted in resistant cells",
    x = NULL,
    y = "Normalized enrichment score (NES)",
    fill = "BRCA1 background"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )

print(p_hallmark)

ggsave(
  "figures/GSE235980/Figure_2_Hallmark_GSEA.png",
  p_hallmark,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/GSE235980/Figure_2_Hallmark_GSEA.pdf",
  p_hallmark,
  width = 9,
  height = 7
)


# ============================================================
# FINAL CHECK
# ============================================================

cat("\n======================================\n")
cat("GSE235980 FIGURE GENERATION COMPLETE\n")
cat("======================================\n")

cat("\nFigure 1A: PCA\n")
cat("Figure 1B: log2FC concordance\n")
cat("Figure 1C: shared DEG summary\n")
cat("Figure 2: Hallmark GSEA\n")

cat(
  "\nAll figures saved under:\n",
  "figures/GSE235980/\n"
)