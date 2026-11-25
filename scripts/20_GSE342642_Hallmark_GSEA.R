# ============================================================
# GSE342642 - DATASET 4
# HALLMARK GSEA
#
# Input:
# Full moderated-t ranking from PRIMARY blocked model
#
# Goal:
# 1. Identify Hallmark programs associated with resistance
# 2. Test the FROZEN Dataset-1 validation hypothesis
# ============================================================

library(clusterProfiler)
library(msigdbr)

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/GSE342642/Hallmark_GSEA"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load ranked primary result
# ------------------------------------------------------------

rank_file <- paste0(
  "results/GSE342642/primary_DEG/",
  "GSE342642_moderated_t_ranking.csv"
)

rank_df <- read.csv(
  rank_file,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2. Build ranked gene vector
#
# Gene_ID = NCBI / Entrez-like ID
# Statistic = moderated t
#
# Positive t:
# higher in resistant
#
# Negative t:
# lower in resistant
# ------------------------------------------------------------

rank_df <- rank_df[
  !is.na(rank_df$t) &
    !is.na(rank_df$Gene_ID),
]

rank_df$Gene_ID <- as.character(
  rank_df$Gene_ID
)


# Check duplicate gene IDs
cat("\n========================================\n")
cat("RANKING QC\n")
cat("========================================\n\n")

cat(
  "Ranked rows:",
  nrow(rank_df),
  "\n"
)

cat(
  "Duplicated Gene_IDs:",
  sum(duplicated(rank_df$Gene_ID)),
  "\n"
)


# If duplicates unexpectedly exist:
# retain the row with greatest absolute t-statistic.

rank_df <- rank_df[
  order(
    abs(rank_df$t),
    decreasing = TRUE
  ),
]

rank_df <- rank_df[
  !duplicated(rank_df$Gene_ID),
]

gene_list <- rank_df$t

names(gene_list) <- rank_df$Gene_ID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


cat(
  "Unique ranked genes:",
  length(gene_list),
  "\n"
)


# ------------------------------------------------------------
# 3. Retrieve current human Hallmark gene sets
# ------------------------------------------------------------

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

cat(
  "\nMSigDB version:",
  unique(hallmark$db_version),
  "\n"
)

cat(
  "Hallmark gene sets:",
  length(unique(hallmark$gs_name)),
  "\n"
)


# ------------------------------------------------------------
# 4. TERM2GENE table
#
# Match using NCBI Gene IDs
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 5. Check overlap before GSEA
# ------------------------------------------------------------

mapped_ranked_genes <- intersect(
  names(gene_list),
  term2gene$gene
)

cat("\n========================================\n")
cat("HALLMARK MAPPING\n")
cat("========================================\n\n")

cat(
  "Ranked genes represented in Hallmark database:",
  length(mapped_ranked_genes),
  "\n"
)

cat(
  "Mapping percentage:",
  round(
    100 *
      length(mapped_ranked_genes) /
      length(gene_list),
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 6. Run ranked Hallmark GSEA
# ------------------------------------------------------------

set.seed(20261023)

gsea <- GSEA(
  geneList = gene_list,
  
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


gsea_res <- as.data.frame(
  gsea
)


# ------------------------------------------------------------
# 7. Sort by adjusted p-value
# ------------------------------------------------------------

gsea_res <- gsea_res[
  order(
    gsea_res$p.adjust,
    -abs(gsea_res$NES)
  ),
]


write.csv(
  gsea_res,
  file.path(
    outdir,
    "GSE342642_Hallmark_GSEA_all.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Significant Hallmarks
#
# Locked criterion:
# BH adjusted p < 0.05
# ------------------------------------------------------------

sig <- gsea_res[
  gsea_res$p.adjust < 0.05,
]

cat("\n========================================\n")
cat("SIGNIFICANT HALLMARKS\n")
cat("========================================\n\n")

cat(
  "Significant Hallmarks:",
  nrow(sig),
  "\n\n"
)

if (nrow(sig) > 0) {
  
  print(
    sig[
      ,
      c(
        "ID",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
}


write.csv(
  sig,
  file.path(
    outdir,
    "GSE342642_Hallmark_GSEA_significant.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Frozen Dataset-1 validation hypothesis
#
# IMPORTANT:
# These predictions were fixed BEFORE Dataset 4 analysis.
# Do not change them based on Dataset 4 results.
# ------------------------------------------------------------

frozen <- data.frame(
  
  pathway = c(
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
  ),
  
  expected_direction = c(
    "UP",
    "UP",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN",
    "DOWN"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. Extract frozen pathways
# ------------------------------------------------------------

frozen_result <- merge(
  frozen,
  gsea_res[
    ,
    c(
      "ID",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ],
  by.x = "pathway",
  by.y = "ID",
  all.x = TRUE
)


# Observed direction from NES
frozen_result$observed_direction <- ifelse(
  frozen_result$NES > 0,
  "UP",
  "DOWN"
)


frozen_result$direction_concordant <-
  frozen_result$expected_direction ==
  frozen_result$observed_direction


frozen_result$significant <-
  frozen_result$p.adjust < 0.05


frozen_result$concordant_and_significant <-
  frozen_result$direction_concordant &
  frozen_result$significant


# Restore original frozen pathway order
frozen_result <- frozen_result[
  match(
    frozen$pathway,
    frozen_result$pathway
  ),
]


cat("\n========================================\n")
cat("FROZEN DATASET-1 VALIDATION TEST\n")
cat("========================================\n\n")

print(
  frozen_result
)


cat(
  "\nDirection concordant:",
  sum(
    frozen_result$direction_concordant,
    na.rm = TRUE
  ),
  "/",
  nrow(frozen_result),
  "\n"
)

cat(
  "Concordant AND significant:",
  sum(
    frozen_result$concordant_and_significant,
    na.rm = TRUE
  ),
  "/",
  nrow(frozen_result),
  "\n"
)


write.csv(
  frozen_result,
  file.path(
    outdir,
    "GSE342642_frozen_D1_validation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Show strongest positive / negative pathways
# ------------------------------------------------------------

positive <- gsea_res[
  gsea_res$NES > 0,
]

negative <- gsea_res[
  gsea_res$NES < 0,
]


positive <- positive[
  order(
    -positive$NES
  ),
]

negative <- negative[
  order(
    negative$NES
  ),
]


cat("\n========================================\n")
cat("TOP POSITIVE HALLMARKS\n")
cat("========================================\n\n")

print(
  head(
    positive[
      ,
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    10
  )
)


cat("\n========================================\n")
cat("TOP NEGATIVE HALLMARKS\n")
cat("========================================\n\n")

print(
  head(
    negative[
      ,
      c(
        "ID",
        "NES",
        "p.adjust"
      )
    ],
    10
  )
)


cat("\n========================================\n")
cat("HALLMARK GSEA COMPLETE\n")
cat("========================================\n")