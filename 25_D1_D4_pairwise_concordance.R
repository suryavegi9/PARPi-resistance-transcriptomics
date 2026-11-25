# ============================================================
# D1-D4 CROSS-DATASET INTEGRATION
# STEP 3: PAIRWISE HALLMARK CONCORDANCE
#
# Questions:
# 1. Which datasets have similar global Hallmark profiles?
# 2. Which pairs show high/low directional agreement?
# 3. Is D4 globally more similar to D1, D2, or D3?
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

outdir <- "results/integration_D1_D4"

master <- read.csv(
  file.path(
    outdir,
    "D1_D4_Hallmark_master_matrix.csv"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 1. NES MATRIX
# ============================================================

nes <- master[
  ,
  c(
    "D1_NES",
    "D2_NES",
    "D3_NES",
    "D4_NES"
  )
]

rownames(nes) <- master$pathway


# ============================================================
# 2. SPEARMAN CORRELATION
#
# Spearman asks:
# Do pathways that rank high/low in one dataset tend to
# rank similarly in another dataset?
# ============================================================

rho <- cor(
  nes,
  method = "spearman"
)

cat("\n========================================\n")
cat("PAIRWISE SPEARMAN NES CORRELATIONS\n")
cat("========================================\n\n")

print(
  round(
    rho,
    3
  )
)

write.csv(
  rho,
  file.path(
    outdir,
    "D1_D4_pairwise_Spearman_NES.csv"
  )
)


# ============================================================
# 3. DIRECTIONAL AGREEMENT
# ============================================================

datasets <- colnames(nes)

direction_matrix <- matrix(
  NA_real_,
  nrow = length(datasets),
  ncol = length(datasets),
  dimnames = list(
    datasets,
    datasets
  )
)

for (i in seq_along(datasets)) {
  
  for (j in seq_along(datasets)) {
    
    direction_matrix[i, j] <- mean(
      sign(nes[[i]]) ==
        sign(nes[[j]])
    )
    
  }
}


cat("\n========================================\n")
cat("PAIRWISE DIRECTION AGREEMENT\n")
cat("========================================\n\n")

print(
  round(
    100 * direction_matrix,
    1
  )
)

write.csv(
  100 * direction_matrix,
  file.path(
    outdir,
    "D1_D4_pairwise_direction_agreement_percent.csv"
  )
)


# ============================================================
# 4. LONG-FORM PAIRWISE SUMMARY
# ============================================================

pairs <- combn(
  datasets,
  2
)

pair_summary <- data.frame()

for (k in seq_len(ncol(pairs))) {
  
  a <- pairs[1, k]
  b <- pairs[2, k]
  
  x <- nes[[a]]
  y <- nes[[b]]
  
  pair_summary <- rbind(
    
    pair_summary,
    
    data.frame(
      dataset_1 = a,
      dataset_2 = b,
      
      spearman_rho = cor(
        x,
        y,
        method = "spearman"
      ),
      
      same_direction_n = sum(
        sign(x) == sign(y)
      ),
      
      opposite_direction_n = sum(
        sign(x) != sign(y)
      ),
      
      direction_agreement_percent =
        100 * mean(
          sign(x) == sign(y)
        ),
      
      stringsAsFactors = FALSE
    )
  )
}


pair_summary <- pair_summary[
  order(
    -pair_summary$spearman_rho
  ),
]


cat("\n========================================\n")
cat("PAIRWISE SUMMARY RANKED BY SIMILARITY\n")
cat("========================================\n\n")

print(
  pair_summary,
  row.names = FALSE
)


write.csv(
  pair_summary,
  file.path(
    outdir,
    "D1_D4_pairwise_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. PATHWAYS DRIVING D1-D4 AGREEMENT / DISAGREEMENT
#
# D1 vs D4 is especially interesting because D4 appeared
# partially similar to D1 during validation.
# ============================================================

d1d4 <- data.frame(
  
  pathway = master$pathway,
  
  D1_NES = master$D1_NES,
  
  D4_NES = master$D4_NES,
  
  same_direction =
    sign(master$D1_NES) ==
    sign(master$D4_NES),
  
  abs_difference =
    abs(
      master$D1_NES -
        master$D4_NES
    ),
  
  stringsAsFactors = FALSE
)


d1d4 <- d1d4[
  order(
    d1d4$same_direction,
    d1d4$abs_difference
  ),
]


write.csv(
  d1d4,
  file.path(
    outdir,
    "D1_vs_D4_pathway_concordance.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("D1 vs D4 SUMMARY\n")
cat("========================================\n\n")

cat(
  "Same direction:",
  sum(d1d4$same_direction),
  "/ 50\n"
)

cat(
  "Opposite direction:",
  sum(!d1d4$same_direction),
  "/ 50\n"
)

cat(
  "Direction agreement:",
  round(
    100 * mean(d1d4$same_direction),
    1
  ),
  "%\n"
)


cat("\n========================================\n")
cat("PAIRWISE CONCORDANCE COMPLETE\n")
cat("========================================\n")