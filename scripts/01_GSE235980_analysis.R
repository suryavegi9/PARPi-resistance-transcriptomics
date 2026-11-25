# ============================================================
# GSE235980 - Olaparib Resistance RNA-seq Analysis
# Project: PARP Inhibitor Resistance
# ============================================================


# ------------------------------------------------------------
# 1. SET WORKING DIRECTORY
# ------------------------------------------------------------

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))


# ------------------------------------------------------------
# 2. LOAD RAW COUNT DATA
# ------------------------------------------------------------

counts_raw <- read.delim(
  "data_raw/GSE235980_CountReads.txt",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(counts_raw)
head(counts_raw)
colnames(counts_raw)


# ------------------------------------------------------------
# 3. CREATE COUNT MATRIX
# ------------------------------------------------------------

counts <- as.matrix(counts_raw[, -1])

rownames(counts) <- counts_raw$Geneid

storage.mode(counts) <- "integer"

dim(counts)
head(counts)


# ------------------------------------------------------------
# 4. BASIC DATA CHECKS
# ------------------------------------------------------------

sum(duplicated(rownames(counts)))

sum(is.na(rownames(counts)))

all(
  sapply(
    as.data.frame(counts),
    function(x) all(x >= 0 & x %% 1 == 0)
  )
)


# ------------------------------------------------------------
# 5. CREATE SAMPLE METADATA
# ------------------------------------------------------------

metadata <- data.frame(
  
  BRCA1_status = c(
    "deficient",
    "deficient",
    "deficient",
    "deficient",
    "proficient",
    "proficient",
    "proficient",
    "proficient"
  ),
  
  resistance = c(
    "resistant",
    "resistant",
    "sensitive",
    "sensitive",
    "resistant",
    "resistant",
    "sensitive",
    "sensitive"
  ),
  
  replicate = c(
    1, 2,
    1, 2,
    1, 2,
    1, 2
  ),
  
  row.names = colnames(counts)
)

metadata


# Confirm sample order

stopifnot(
  identical(
    rownames(metadata),
    colnames(counts)
  )
)


# ------------------------------------------------------------
# 6. CREATE COMBINED EXPERIMENTAL GROUP
# ------------------------------------------------------------

metadata$group <- paste(
  metadata$BRCA1_status,
  metadata$resistance,
  sep = "_"
)

metadata$group <- factor(metadata$group)

table(metadata$group)


# Save metadata

write.csv(
  metadata,
  "metadata/GSE235980_metadata.csv"
)


# ------------------------------------------------------------
# 7. LIBRARY SIZE QC
# ------------------------------------------------------------

library_sizes <- colSums(counts)

library_sizes

round(
  library_sizes / 1e6,
  2
)

barplot(
  library_sizes / 1e6,
  las = 2,
  ylab = "Library size (millions of reads)",
  main = "GSE235980 Library Sizes"
)


# ------------------------------------------------------------
# 8. LOW-COUNT FILTERING
# ------------------------------------------------------------

nrow(counts)

sum(
  rowSums(counts) == 0
)

# Keep genes with >=10 reads in at least 2 samples

keep <- rowSums(counts >= 10) >= 2

table(keep)

counts_filtered <- counts[keep, ]

dim(counts_filtered)


# ------------------------------------------------------------
# 9. LOAD DESEQ2
# ------------------------------------------------------------

library(DESeq2)


# ------------------------------------------------------------
# 10. CREATE DESEQ2 DATASET
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  
  countData = counts_filtered,
  
  colData = metadata,
  
  design = ~ group
)


# ------------------------------------------------------------
# 11. RUN DESEQ2
# ------------------------------------------------------------

dds <- DESeq(dds)

sizeFactors(dds)


# ------------------------------------------------------------
# 12. PCA
# ------------------------------------------------------------

vsd <- vst(
  dds,
  blind = TRUE
)

plotPCA(
  vsd,
  intgroup = c(
    "BRCA1_status",
    "resistance"
  )
)


# ------------------------------------------------------------
# 13. BRCA1-DEFICIENT:
#     RESISTANT VS SENSITIVE
# ------------------------------------------------------------

res_def <- results(
  
  dds,
  
  contrast = c(
    "group",
    "deficient_resistant",
    "deficient_sensitive"
  )
)

summary(res_def)


# ------------------------------------------------------------
# 14. BRCA1-PROFICIENT:
#     RESISTANT VS SENSITIVE
# ------------------------------------------------------------

res_prof <- results(
  
  dds,
  
  contrast = c(
    "group",
    "proficient_resistant",
    "proficient_sensitive"
  )
)

summary(res_prof)


# ------------------------------------------------------------
# 15. SIGNIFICANT DEGs
# Criteria:
# adjusted p-value < 0.05
# |log2 fold change| >= 1
# ------------------------------------------------------------

sig_def <- subset(
  as.data.frame(res_def),
  
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)


sig_prof <- subset(
  as.data.frame(res_prof),
  
  !is.na(padj) &
    padj < 0.05 &
    abs(log2FoldChange) >= 1
)


# Number of significant DEGs

nrow(sig_def)

nrow(sig_prof)


# Direction of change

table(
  sig_def$log2FoldChange > 0
)

table(
  sig_prof$log2FoldChange > 0
)


# ------------------------------------------------------------
# 16. FIND SHARED SIGNIFICANT GENES
# ------------------------------------------------------------

shared_genes <- intersect(
  rownames(sig_def),
  rownames(sig_prof)
)

length(shared_genes)


# ------------------------------------------------------------
# 17. COMPARE DIRECTION OF SHARED GENES
# ------------------------------------------------------------

shared_results <- data.frame(
  
  gene = shared_genes,
  
  log2FC_deficient =
    sig_def[
      shared_genes,
      "log2FoldChange"
    ],
  
  log2FC_proficient =
    sig_prof[
      shared_genes,
      "log2FoldChange"
    ]
)


shared_results$same_direction <-
  
  sign(
    shared_results$log2FC_deficient
  ) ==
  
  sign(
    shared_results$log2FC_proficient
  )


table(
  shared_results$same_direction
)


# ------------------------------------------------------------
# 18. CONSISTENT SHARED GENES
# ------------------------------------------------------------

shared_consistent <- subset(
  shared_results,
  same_direction == TRUE
)

nrow(shared_consistent)


# TRUE = increased in resistant cells
# FALSE = decreased in resistant cells

table(
  shared_consistent$log2FC_deficient > 0
)


# ------------------------------------------------------------
# 19. SAVE RESULTS
# ------------------------------------------------------------

write.csv(
  as.data.frame(res_def),
  "results/differential_expression/BRCA1_deficient_all_DESeq2_results.csv"
)

write.csv(
  as.data.frame(res_prof),
  "results/differential_expression/BRCA1_proficient_all_DESeq2_results.csv"
)

write.csv(
  sig_def,
  "results/differential_expression/BRCA1_deficient_significant_DEGs.csv"
)

write.csv(
  sig_prof,
  "results/differential_expression/BRCA1_proficient_significant_DEGs.csv"
)

write.csv(
  shared_results,
  "results/differential_expression/shared_DEGs_comparison.csv",
  row.names = FALSE
)

write.csv(
  shared_consistent,
  "results/differential_expression/shared_consistent_DEGs.csv",
  row.names = FALSE
)
nrow(sig_def)

nrow(sig_prof)

length(shared_genes)

table(shared_results$same_direction)
shared_consistent <- subset(
  shared_results,
  same_direction == TRUE
)

table(shared_consistent$log2FC_deficient > 0)
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

cor.test(
  all_changes$log2FC_deficient,
  all_changes$log2FC_proficient,
  method = "spearman"
)
exists("all_changes")
dim(all_changes)
cor.test(
  all_changes$log2FC_deficient,
  all_changes$log2FC_proficient,
  method = "spearman"
)
table(shared_consistent$log2FC_deficient > 0)
# ------------------------------------------------------------
# 20. GLOBAL CONCORDANCE BETWEEN BRCA1 BACKGROUNDS
# ------------------------------------------------------------

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

cor.test(
  all_changes$log2FC_deficient,
  all_changes$log2FC_proficient,
  method = "spearman"
)

table(shared_consistent$log2FC_deficient > 0)
nrow(shared_consistent)
table(shared_consistent$log2FC_deficient > 0)
library(ggplot2)

p_concordance <- ggplot(
  all_changes,
  aes(
    x = log2FC_deficient,
    y = log2FC_proficient
  )
) +
  geom_point(
    alpha = 0.25,
    size = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Concordance of Olaparib Resistance-Associated Expression Changes",
    subtitle = "BRCA1-deficient vs BRCA1-proficient backgrounds",
    x = "log2 fold change: BRCA1-deficient",
    y = "log2 fold change: BRCA1-proficient"
  ) +
  theme_minimal()

p_concordance
exists("p_concordance")
p_concordance
shared_up <- subset(
  shared_consistent,
  log2FC_deficient > 0
)

shared_down <- subset(
  shared_consistent,
  log2FC_deficient < 0
)

nrow(shared_up)
nrow(shared_down)
shared_up <- subset(
  shared_consistent,
  log2FC_deficient > 0
)

shared_down <- subset(
  shared_consistent,
  log2FC_deficient < 0
)
exists("shared_up")
exists("shared_down")
exists("shared_up")
nrow(shared_up)
nrow(shared_down)
write.csv(
  shared_up,
  "results/differential_expression/shared_consistent_upregulated.csv",
  row.names = FALSE
)

write.csv(
  shared_down,
  "results/differential_expression/shared_consistent_downregulated.csv",
  row.names = FALSE
)
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db"
))
n
library(clusterProfiler)
library(org.Hs.eg.db)
BiocManager::install("GO.db")
n
library(clusterProfiler)
library(org.Hs.eg.db)
BiocManager::install("HDO.db")
n
library(clusterProfiler)
install.packages(
  "shadowtext",
  dependencies = TRUE
)
BiocManager::install(
  c(
    "clusterProfiler",
    "enrichplot",
    "DOSE",
    "GO.db",
    "HDO.db",
    "org.Hs.eg.db"
  ),
  dependencies = TRUE,
  ask = FALSE,
  update = TRUE,
  force = TRUE
)
library(clusterProfiler)
library(org.Hs.eg.db)
"shadowtext" %in% rownames(installed.packages())
install.packages(
  "shadowtext",
  repos = "https://cloud.r-project.org",
  dependencies = TRUE
)
packageVersion("shadowtext")
R.version.string
R.version.string
library(BiocManager)
BiocManager::install(
  c(
    "DESeq2",
    "clusterProfiler",
    "org.Hs.eg.db"
  ),
  ask = FALSE,
  update = FALSE
)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
BiocManager::install(
  "clusterProfiler",
  ask = FALSE,
  update = FALSE,
  force = TRUE
)

library(clusterProfiler)
n
n
n
library(clusterProfiler)
library(org.Hs.eg.db)
nrow(shared_up)
nrow(shared_down)
up_genes <- shared_up$gene
down_genes <- shared_down$gene
background_genes <- rownames(counts_filtered)
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
head(as.data.frame(ego_up), 10)
up_results <- as.data.frame(ego_up)

up_results[
  1:10,
  c(
    "Description",
    "GeneRatio",
    "p.adjust",
    "geneID"
  )
]
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
head(as.data.frame(ego_down), 10)
down_mapped <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

nrow(down_mapped)
length(unique(down_mapped$SYMBOL))
length(down_genes)
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

head(
  down_all_df[
    order(down_all_df$p.adjust),
  ],
  10
)
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
exists("ego_down_all")
down_all_df <- as.data.frame(ego_down_all)
head(
  down_all_df[
    order(down_all_df$p.adjust),
  ],
  10
)
down_mapped <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

nrow(down_mapped)
length(unique(down_mapped$SYMBOL))
length(down_genes)
nrow(down_mapped)

length(unique(down_mapped$SYMBOL))

length(down_genes)

head(
  down_all_df[
    order(down_all_df$p.adjust),
  ],
  10
)
