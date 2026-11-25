# ============================================================
# GSE342642 - DATASET 4
# HALLMARK GSEA SENSITIVITY ANALYSES
#
# A. Primary regional model WITHOUT low-depth sample 440-R3
# B. Tumor-level pseudobulk Hallmark GSEA
#
# These analyses test robustness of the PRIMARY result.
# ============================================================

library(edgeR)
library(limma)
library(clusterProfiler)
library(msigdbr)

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/GSE342642/Hallmark_sensitivity"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. LOAD ORIGINAL GEO DATA
# ============================================================

geo_file <- file.choose()

dat <- read.delim(
  geo_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_cols <- grep(
  "_Read_Count$",
  colnames(dat),
  value = TRUE
)

sample_names <- sub(
  "_Read_Count$",
  "",
  count_cols
)

counts <- as.matrix(
  dat[, count_cols]
)

storage.mode(counts) <- "integer"

colnames(counts) <- sample_names
rownames(counts) <- as.character(dat$Gene_ID)


condition <- ifelse(
  grepl("^Cont_", sample_names),
  "Control",
  "Resistant"
)

tumor <- sub(
  "^(Cont|Resi)_([0-9]+).*",
  "\\2",
  sample_names
)

metadata <- data.frame(
  sample = sample_names,
  condition = factor(
    condition,
    levels = c("Control", "Resistant")
  ),
  tumor = factor(tumor)
)

rownames(metadata) <- sample_names


# ============================================================
# 2. HALLMARK DATABASE
# ============================================================

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

term2gene <- unique(
  data.frame(
    term = hallmark$gs_name,
    gene = as.character(hallmark$ncbi_gene),
    stringsAsFactors = FALSE
  )
)

term2gene <- term2gene[
  !is.na(term2gene$gene) &
    term2gene$gene != "",
]


run_hallmark_gsea <- function(t_stats) {
  
  t_stats <- t_stats[
    !is.na(t_stats)
  ]
  
  t_stats <- sort(
    t_stats,
    decreasing = TRUE
  )
  
  set.seed(20261023)
  
  x <- GSEA(
    geneList = t_stats,
    TERM2GENE = term2gene,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 1e-10,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  )
  
  as.data.frame(x)
}


# ============================================================
# PART A
# REGION-LEVEL MODEL WITHOUT 440-R3
# ============================================================

cat("\n========================================\n")
cat("A. REMOVE 440-R3 SENSITIVITY\n")
cat("========================================\n\n")

remove_sample <- "Cont_440_P3_2-3"

keep_samples <- colnames(counts) != remove_sample

counts_A <- counts[, keep_samples]

meta_A <- metadata[
  colnames(counts_A),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Tumor-aware filtering
# Same biological rule as primary analysis
# ------------------------------------------------------------

dge_A <- DGEList(
  counts = counts_A
)

dge_A <- normLibSizes(
  dge_A,
  method = "TMM"
)

cpm_A <- cpm(
  dge_A,
  log = FALSE
)

tumor_levels_A <- levels(
  droplevels(meta_A$tumor)
)

tumor_median_cpm_A <- sapply(
  tumor_levels_A,
  function(t) {
    
    idx <- meta_A$tumor == t
    
    apply(
      cpm_A[, idx, drop = FALSE],
      1,
      median
    )
  }
)

keep_A <- rowSums(
  tumor_median_cpm_A >= 1
) >= 2

dge_A <- dge_A[
  keep_A,
  ,
  keep.lib.sizes = FALSE
]

dge_A <- normLibSizes(
  dge_A,
  method = "TMM"
)


# ------------------------------------------------------------
# Blocked limma-voom
# ------------------------------------------------------------

design_A <- model.matrix(
  ~ condition,
  data = meta_A
)

v1_A <- voom(
  dge_A,
  design_A,
  plot = FALSE
)

cor1_A <- duplicateCorrelation(
  v1_A,
  design_A,
  block = meta_A$tumor
)

v2_A <- voom(
  dge_A,
  design_A,
  plot = FALSE,
  block = meta_A$tumor,
  correlation = cor1_A$consensus.correlation
)

cor2_A <- duplicateCorrelation(
  v2_A,
  design_A,
  block = meta_A$tumor
)

fit_A <- lmFit(
  v2_A,
  design_A,
  block = meta_A$tumor,
  correlation = cor2_A$consensus.correlation
)

fit_A <- eBayes(
  fit_A
)

tt_A <- topTable(
  fit_A,
  coef = "conditionResistant",
  number = Inf,
  sort.by = "none"
)

stats_A <- tt_A$t
names(stats_A) <- rownames(tt_A)

gsea_A <- run_hallmark_gsea(
  stats_A
)

write.csv(
  gsea_A,
  file.path(
    outdir,
    "GSE342642_GSEA_without_440R3.csv"
  ),
  row.names = FALSE
)


# ============================================================
# PART B
# TUMOR-LEVEL PSEUDOBULK GSEA
# ============================================================

cat("\n========================================\n")
cat("B. PSEUDOBULK GSEA\n")
cat("========================================\n\n")

tumor_ids <- c(
  "437",
  "440",
  "433",
  "438"
)

pb_counts <- sapply(
  tumor_ids,
  function(t) {
    
    idx <- metadata$tumor == t
    
    rowSums(
      counts[, idx, drop = FALSE]
    )
  }
)

colnames(pb_counts) <- paste0(
  "Tumor_",
  tumor_ids
)

pb_condition <- factor(
  c(
    "Control",
    "Control",
    "Resistant",
    "Resistant"
  ),
  levels = c(
    "Control",
    "Resistant"
  )
)

pb_meta <- data.frame(
  condition = pb_condition
)

rownames(pb_meta) <- colnames(pb_counts)


# Use genes tested in PRIMARY model
primary <- read.csv(
  "results/GSE342642/primary_DEG/GSE342642_primary_DEG_results.csv",
  stringsAsFactors = FALSE
)

tested_genes <- as.character(
  primary$Gene_ID
)

pb_counts <- pb_counts[
  intersect(
    tested_genes,
    rownames(pb_counts)
  ),
  ,
  drop = FALSE
]


dge_B <- DGEList(
  counts = pb_counts
)

dge_B <- normLibSizes(
  dge_B,
  method = "TMM"
)

design_B <- model.matrix(
  ~ condition,
  data = pb_meta
)

v_B <- voom(
  dge_B,
  design_B,
  plot = FALSE
)

fit_B <- lmFit(
  v_B,
  design_B
)

fit_B <- eBayes(
  fit_B
)

tt_B <- topTable(
  fit_B,
  coef = "conditionResistant",
  number = Inf,
  sort.by = "none"
)

stats_B <- tt_B$t
names(stats_B) <- rownames(tt_B)

gsea_B <- run_hallmark_gsea(
  stats_B
)

write.csv(
  gsea_B,
  file.path(
    outdir,
    "GSE342642_GSEA_pseudobulk.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 3. LOAD PRIMARY GSEA
# ============================================================

primary_gsea <- read.csv(
  "results/GSE342642/Hallmark_GSEA/GSE342642_Hallmark_GSEA_all.csv",
  stringsAsFactors = FALSE
)


# ============================================================
# 4. MERGE ALL THREE ANALYSES
# ============================================================

comparison <- merge(
  primary_gsea[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  gsea_A[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  by = "ID",
  suffixes = c(
    "_primary",
    "_without440R3"
  )
)

comparison <- merge(
  comparison,
  gsea_B[
    ,
    c(
      "ID",
      "NES",
      "p.adjust"
    )
  ],
  by = "ID"
)

colnames(comparison)[
  colnames(comparison) == "NES"
] <- "NES_pseudobulk"

colnames(comparison)[
  colnames(comparison) == "p.adjust"
] <- "padj_pseudobulk"


comparison$same_direction_without440R3 <-
  sign(comparison$NES_primary) ==
  sign(comparison$NES_without440R3)

comparison$same_direction_pseudobulk <-
  sign(comparison$NES_primary) ==
  sign(comparison$NES_pseudobulk)


# ============================================================
# 5. FROZEN D1 PATHWAYS
# ============================================================

frozen_pathways <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_UV_RESPONSE_DN"
)

frozen_compare <- comparison[
  match(
    frozen_pathways,
    comparison$ID
  ),
]


# ============================================================
# 6. FIVE REPLICATED D1/D4 PATHWAYS
# ============================================================

core5 <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_GLYCOLYSIS"
)

core5_compare <- comparison[
  match(
    core5,
    comparison$ID
  ),
]


# ============================================================
# 7. SUMMARY
# ============================================================

cat("\n========================================\n")
cat("GLOBAL HALLMARK ROBUSTNESS\n")
cat("========================================\n\n")

cat(
  "Same direction after removing 440-R3:",
  sum(
    comparison$same_direction_without440R3
  ),
  "/",
  nrow(comparison),
  "\n"
)

cat(
  "Same direction in pseudobulk:",
  sum(
    comparison$same_direction_pseudobulk
  ),
  "/",
  nrow(comparison),
  "\n"
)


rho_A <- cor(
  comparison$NES_primary,
  comparison$NES_without440R3,
  method = "spearman"
)

rho_B <- cor(
  comparison$NES_primary,
  comparison$NES_pseudobulk,
  method = "spearman"
)


cat(
  "Spearman rho primary vs without 440-R3:",
  round(rho_A, 4),
  "\n"
)

cat(
  "Spearman rho primary vs pseudobulk:",
  round(rho_B, 4),
  "\n"
)


cat("\n========================================\n")
cat("FROZEN D1 PATHWAYS\n")
cat("========================================\n\n")

print(
  frozen_compare[
    ,
    c(
      "ID",
      "NES_primary",
      "NES_without440R3",
      "NES_pseudobulk",
      "same_direction_without440R3",
      "same_direction_pseudobulk"
    )
  ]
)


cat("\n========================================\n")
cat("CORE FIVE D1/D4 PATHWAYS\n")
cat("========================================\n\n")

print(
  core5_compare[
    ,
    c(
      "ID",
      "NES_primary",
      "p.adjust_primary",
      "NES_without440R3",
      "p.adjust_without440R3",
      "NES_pseudobulk",
      "padj_pseudobulk"
    )
  ]
)


write.csv(
  comparison,
  file.path(
    outdir,
    "Hallmark_three_analysis_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  frozen_compare,
  file.path(
    outdir,
    "Frozen_D1_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  core5_compare,
  file.path(
    outdir,
    "Core5_sensitivity_comparison.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("HALLMARK SENSITIVITY COMPLETE\n")
cat("========================================\n")
