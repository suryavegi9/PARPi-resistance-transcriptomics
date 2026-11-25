# ============================================================
# GSE235980 - Hallmark GSEA
# Shared olaparib-resistance programs across BRCA1 contexts
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

library(clusterProfiler)
library(msigdbr)


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


# ------------------------------------------------------------
# 3. CREATE RANKED GENE LISTS
# Wald statistic:
# positive = enriched toward resistant state
# negative = enriched toward sensitive state
# ------------------------------------------------------------

rank_def <- res_def$stat
names(rank_def) <- rownames(res_def)
rank_def <- rank_def[!is.na(rank_def)]
rank_def <- sort(rank_def, decreasing = TRUE)

rank_prof <- res_prof$stat
names(rank_prof) <- rownames(res_prof)
rank_prof <- rank_prof[!is.na(rank_prof)]
rank_prof <- sort(rank_prof, decreasing = TRUE)

cat("Ranked deficient genes:", length(rank_def), "\n")
cat("Ranked proficient genes:", length(rank_prof), "\n")


# ------------------------------------------------------------
# 4. GET HUMAN MSigDB HALLMARK GENE SETS
# ------------------------------------------------------------

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_t2g <- unique(
  hallmark[, c("gs_name", "gene_symbol")]
)

cat(
  "Number of Hallmark gene sets:",
  length(unique(hallmark_t2g$gs_name)),
  "\n"
)


# ------------------------------------------------------------
# 5. HALLMARK GSEA - BRCA1 DEFICIENT
# ------------------------------------------------------------

hall_def <- GSEA(
  geneList = rank_def,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

hall_def_df <- as.data.frame(hall_def)


# ------------------------------------------------------------
# 6. HALLMARK GSEA - BRCA1 PROFICIENT
# ------------------------------------------------------------

hall_prof <- GSEA(
  geneList = rank_prof,
  TERM2GENE = hallmark_t2g,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

hall_prof_df <- as.data.frame(hall_prof)


# ------------------------------------------------------------
# 7. DEFINE SIGNIFICANT HALLMARKS
# ------------------------------------------------------------

sig_hall_def <- subset(
  hall_def_df,
  !is.na(p.adjust) & p.adjust < 0.05
)

sig_hall_prof <- subset(
  hall_prof_df,
  !is.na(p.adjust) & p.adjust < 0.05
)

cat(
  "Significant Hallmarks - BRCA1 deficient:",
  nrow(sig_hall_def),
  "\n"
)

cat(
  "Significant Hallmarks - BRCA1 proficient:",
  nrow(sig_hall_prof),
  "\n"
)


# ------------------------------------------------------------
# 8. SHARED SIGNIFICANT HALLMARKS
# ------------------------------------------------------------

shared_hall_ids <- intersect(
  sig_hall_def$ID,
  sig_hall_prof$ID
)

cat(
  "Shared significant Hallmarks:",
  length(shared_hall_ids),
  "\n"
)


shared_hall <- data.frame(
  Hallmark = shared_hall_ids,
  
  NES_deficient =
    sig_hall_def[
      match(shared_hall_ids, sig_hall_def$ID),
      "NES"
    ],
  
  padj_deficient =
    sig_hall_def[
      match(shared_hall_ids, sig_hall_def$ID),
      "p.adjust"
    ],
  
  NES_proficient =
    sig_hall_prof[
      match(shared_hall_ids, sig_hall_prof$ID),
      "NES"
    ],
  
  padj_proficient =
    sig_hall_prof[
      match(shared_hall_ids, sig_hall_prof$ID),
      "p.adjust"
    ]
)


# ------------------------------------------------------------
# 9. DIRECTIONAL CONSISTENCY
# ------------------------------------------------------------

shared_hall$same_direction <-
  sign(shared_hall$NES_deficient) ==
  sign(shared_hall$NES_proficient)

cat(
  "Shared Hallmarks with same direction:",
  sum(shared_hall$same_direction),
  "\n"
)

cat(
  "Shared Hallmarks with opposite direction:",
  sum(!shared_hall$same_direction),
  "\n"
)


shared_hall_consistent <- subset(
  shared_hall,
  same_direction == TRUE
)

shared_hall_consistent$mean_abs_NES <-
  (
    abs(shared_hall_consistent$NES_deficient) +
      abs(shared_hall_consistent$NES_proficient)
  ) / 2

shared_hall_consistent <-
  shared_hall_consistent[
    order(
      shared_hall_consistent$mean_abs_NES,
      decreasing = TRUE
    ),
  ]


# ------------------------------------------------------------
# 10. PRINT SHARED HALLMARK RESULTS
# ------------------------------------------------------------

print(shared_hall_consistent)


# ------------------------------------------------------------
# 11. CHECK SPECIFIC BIOLOGICAL PROGRAMS
# ------------------------------------------------------------

programs_of_interest <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_HYPOXIA"
)

check_def <- hall_def_df[
  hall_def_df$ID %in% programs_of_interest,
  c("ID", "NES", "p.adjust")
]

check_prof <- hall_prof_df[
  hall_prof_df$ID %in% programs_of_interest,
  c("ID", "NES", "p.adjust")
]

names(check_def)[2:3] <-
  c("NES_deficient", "padj_deficient")

names(check_prof)[2:3] <-
  c("NES_proficient", "padj_proficient")

program_check <- merge(
  check_def,
  check_prof,
  by = "ID",
  all = TRUE
)

print(program_check)


# ------------------------------------------------------------
# 12. SAVE EVERYTHING
# ------------------------------------------------------------

write.csv(
  hall_def_df,
  "results/pathway_analysis/Hallmark_GSEA_BRCA1_deficient_all.csv",
  row.names = FALSE
)

write.csv(
  hall_prof_df,
  "results/pathway_analysis/Hallmark_GSEA_BRCA1_proficient_all.csv",
  row.names = FALSE
)

write.csv(
  shared_hall,
  "results/pathway_analysis/Hallmark_GSEA_shared.csv",
  row.names = FALSE
)

write.csv(
  shared_hall_consistent,
  "results/pathway_analysis/Hallmark_GSEA_shared_same_direction.csv",
  row.names = FALSE
)

write.csv(
  program_check,
  "results/pathway_analysis/Hallmark_GSEA_programs_of_interest.csv",
  row.names = FALSE
)