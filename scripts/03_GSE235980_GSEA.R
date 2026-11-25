# ============================================================
# GSE235980 - Ranked GSEA
# Shared olaparib resistance programs across BRCA1 contexts
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

# ------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)


# ------------------------------------------------------------
# 2. LOAD SAVED DESEQ2 RESULTS
# ------------------------------------------------------------

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

cat("BRCA1-deficient genes:", nrow(res_def), "\n")
cat("BRCA1-proficient genes:", nrow(res_prof), "\n")


# ------------------------------------------------------------
# 3. BUILD RANKED GENE LISTS
#
# We use the DESeq2 Wald statistic rather than log2FC alone.
# Positive stat = higher in resistant cells.
# Negative stat = lower in resistant cells.
# ------------------------------------------------------------

rank_def <- res_def$stat
names(rank_def) <- rownames(res_def)

rank_prof <- res_prof$stat
names(rank_prof) <- rownames(res_prof)


# Remove missing values

rank_def <- rank_def[!is.na(rank_def)]
rank_prof <- rank_prof[!is.na(rank_prof)]


# Sort from most positive to most negative

rank_def <- sort(
  rank_def,
  decreasing = TRUE
)

rank_prof <- sort(
  rank_prof,
  decreasing = TRUE
)


cat("Ranked deficient genes:", length(rank_def), "\n")
cat("Ranked proficient genes:", length(rank_prof), "\n")


# ------------------------------------------------------------
# 4. GSEA - BRCA1 DEFICIENT
# ------------------------------------------------------------

gsea_def <- gseGO(
  geneList = rank_def,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE
)

gsea_def_df <- as.data.frame(gsea_def)

cat(
  "Significant GSEA terms - BRCA1 deficient:",
  nrow(gsea_def_df),
  "\n"
)


# ------------------------------------------------------------
# 5. GSEA - BRCA1 PROFICIENT
# ------------------------------------------------------------

gsea_prof <- gseGO(
  geneList = rank_prof,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE
)

gsea_prof_df <- as.data.frame(gsea_prof)

cat(
  "Significant GSEA terms - BRCA1 proficient:",
  nrow(gsea_prof_df),
  "\n"
)


# ------------------------------------------------------------
# 6. SAVE FULL GSEA RESULTS
# ------------------------------------------------------------

write.csv(
  gsea_def_df,
  "results/pathway_analysis/GSEA_BRCA1_deficient_GO_BP.csv",
  row.names = FALSE
)

write.csv(
  gsea_prof_df,
  "results/pathway_analysis/GSEA_BRCA1_proficient_GO_BP.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. FIND GO TERMS SIGNIFICANT IN BOTH BACKGROUNDS
# ------------------------------------------------------------

shared_ids <- intersect(
  gsea_def_df$ID,
  gsea_prof_df$ID
)

cat(
  "Significant GO terms shared between backgrounds:",
  length(shared_ids),
  "\n"
)


shared_gsea <- data.frame(
  ID = shared_ids,
  
  Description =
    gsea_def_df[
      match(shared_ids, gsea_def_df$ID),
      "Description"
    ],
  
  NES_deficient =
    gsea_def_df[
      match(shared_ids, gsea_def_df$ID),
      "NES"
    ],
  
  padj_deficient =
    gsea_def_df[
      match(shared_ids, gsea_def_df$ID),
      "p.adjust"
    ],
  
  NES_proficient =
    gsea_prof_df[
      match(shared_ids, gsea_prof_df$ID),
      "NES"
    ],
  
  padj_proficient =
    gsea_prof_df[
      match(shared_ids, gsea_prof_df$ID),
      "p.adjust"
    ]
)


# ------------------------------------------------------------
# 8. CHECK WHETHER ENRICHMENT DIRECTION AGREES
# ------------------------------------------------------------

shared_gsea$same_direction <-
  sign(shared_gsea$NES_deficient) ==
  sign(shared_gsea$NES_proficient)


cat(
  "Shared GO terms with same NES direction:",
  sum(shared_gsea$same_direction),
  "\n"
)

cat(
  "Shared GO terms with opposite NES direction:",
  sum(!shared_gsea$same_direction),
  "\n"
)


# ------------------------------------------------------------
# 9. KEEP DIRECTIONALLY CONSISTENT SHARED PATHWAYS
# ------------------------------------------------------------

shared_gsea_consistent <- subset(
  shared_gsea,
  same_direction == TRUE
)

# Strongest first
shared_gsea_consistent$mean_abs_NES <-
  (
    abs(shared_gsea_consistent$NES_deficient) +
      abs(shared_gsea_consistent$NES_proficient)
  ) / 2

shared_gsea_consistent <-
  shared_gsea_consistent[
    order(
      shared_gsea_consistent$mean_abs_NES,
      decreasing = TRUE
    ),
  ]


# ------------------------------------------------------------
# 10. DISPLAY TOP SHARED PATHWAYS
# ------------------------------------------------------------

print(
  head(
    shared_gsea_consistent[
      ,
      c(
        "Description",
        "NES_deficient",
        "NES_proficient",
        "padj_deficient",
        "padj_proficient"
      )
    ],
    20
  )
)


# ------------------------------------------------------------
# 11. SAVE SHARED PATHWAYS
# ------------------------------------------------------------

write.csv(
  shared_gsea,
  "results/pathway_analysis/GSEA_shared_GO_BP_all.csv",
  row.names = FALSE
)

write.csv(
  shared_gsea_consistent,
  "results/pathway_analysis/GSEA_shared_GO_BP_same_direction.csv",
  row.names = FALSE
)
# ------------------------------------------------------------
# 12. INSPECT LEADING-EDGE GENES OF SHARED GSEA TERMS
# ------------------------------------------------------------

top_shared_ids <- head(
  shared_gsea_consistent$ID,
  10
)

leading_edge <- data.frame()

for (go_id in top_shared_ids) {
  
  def_row <- gsea_def_df[gsea_def_df$ID == go_id, ]
  prof_row <- gsea_prof_df[gsea_prof_df$ID == go_id, ]
  
  def_genes <- unlist(
    strsplit(def_row$core_enrichment, "/")
  )
  
  prof_genes <- unlist(
    strsplit(prof_row$core_enrichment, "/")
  )
  
  common_genes <- intersect(
    def_genes,
    prof_genes
  )
  
  leading_edge <- rbind(
    leading_edge,
    data.frame(
      ID = go_id,
      Description = def_row$Description,
      deficient_core_n = length(def_genes),
      proficient_core_n = length(prof_genes),
      shared_core_n = length(common_genes),
      shared_core_genes = paste(
        common_genes,
        collapse = "/"
      )
    )
  )
}

print(leading_edge)