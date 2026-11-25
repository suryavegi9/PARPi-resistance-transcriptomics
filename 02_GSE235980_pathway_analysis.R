# ============================================================
# GSE235980 - Pathway Enrichment Analysis
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

# ------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 2. LOAD SAVED SHARED DEG RESULTS
# ------------------------------------------------------------

shared_up <- read.csv(
  "results/differential_expression/shared_consistent_upregulated.csv"
)

shared_down <- read.csv(
  "results/differential_expression/shared_consistent_downregulated.csv"
)

up_genes <- shared_up$gene
down_genes <- shared_down$gene

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")


# ------------------------------------------------------------
# 3. RECREATE BACKGROUND GENE SET
# ------------------------------------------------------------

counts_raw <- read.delim(
  "data_raw/GSE235980_CountReads.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts_raw[, -1])
rownames(counts) <- counts_raw$Geneid

keep <- rowSums(counts >= 10) >= 2

counts_filtered <- counts[keep, ]

background_genes <- rownames(counts_filtered)

cat("Background genes:", length(background_genes), "\n")


# ------------------------------------------------------------
# 4. GO ENRICHMENT - UPREGULATED GENES
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_up_df <- as.data.frame(ego_up)

cat("Significant GO terms - UP:", nrow(ego_up_df), "\n")

print(
  head(ego_up_df, 10)
)


# ------------------------------------------------------------
# 5. GO ENRICHMENT - DOWNREGULATED GENES
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_down_df <- as.data.frame(ego_down)

cat("Significant GO terms - DOWN:", nrow(ego_down_df), "\n")


# ------------------------------------------------------------
# 6. CHECK DOWNREGULATED GENE MAPPING
# ------------------------------------------------------------

down_mapped <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat(
  "Mapped downregulated genes:",
  length(unique(down_mapped$SYMBOL)),
  "out of",
  length(down_genes),
  "\n"
)


# ------------------------------------------------------------
# 7. EXPLORATORY DOWNREGULATED GO TERMS
#    Not considered statistically significant
# ------------------------------------------------------------

ego_down_all <- enrichGO(
  gene = down_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

down_all_df <- as.data.frame(ego_down_all)

down_all_df <- down_all_df[
  order(down_all_df$p.adjust),
]

print(
  head(down_all_df, 10)
)


# ------------------------------------------------------------
# 8. SAVE RESULTS
# ------------------------------------------------------------

write.csv(
  ego_up_df,
  "results/pathway_analysis/GO_upregulated.csv",
  row.names = FALSE
)

write.csv(
  ego_down_df,
  "results/pathway_analysis/GO_downregulated.csv",
  row.names = FALSE
)

write.csv(
  down_all_df,
  "results/pathway_analysis/GO_downregulated_exploratory.csv",
  row.names = FALSE
)