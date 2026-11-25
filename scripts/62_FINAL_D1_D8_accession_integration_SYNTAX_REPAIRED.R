# ============================================================
# 62_FINAL_D1_D8_accession_integration.R
# FINAL CLEAN VERSION
#
# D1-D6: frozen; never re-analysed here
# D7: GSE165914, one conservative accession vote
# D8: GSE163854, frozen S2 vs C3 primary vote
# GSE287274: excluded
#
# FINAL FROZEN N = 8 RULES
# Directional support: ceil(2N/3) = 6/8
# Strong significant support: ceil(N/2) = 4/8
#
# Classification priority:
# 1. Context-dependent
# 2. Strong recurrent
# 3. Directional
# 4. Dataset-specific
# 5. Inconsistent
# ============================================================


# ============================================================
# 0. PROJECT SETUP
# ============================================================

setwd(Sys.getenv("PARPI_PROJECT_ROOT", unset = getwd()))

old_master_file <-
  "results/FINAL_D1_D6_integration/FINAL_D1_D6_Hallmark_master.csv"

d7_root <-
  "results/GSE165914"

d8_root <-
  "results/GSE163854"

out_dir <-
  "results/FINAL_D1_D8_integration/accession_level"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("============================================================\n")
cat("FINAL D1-D8 ACCESSION-LEVEL INTEGRATION\n")
cat("============================================================\n\n")


# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

norm_name <- function(x) {
  
  tolower(
    gsub(
      "[^a-zA-Z0-9]",
      "",
      x
    )
  )
}


find_col <- function(
    column_names,
    aliases
) {
  
  hit <- which(
    norm_name(column_names) %in%
      norm_name(aliases)
  )
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  column_names[hit[1]]
}


classify_pathways <- function(
    nUP,
    nDOWN,
    nSigUP,
    nSigDown,
    direction_threshold,
    strong_sig_threshold
) {
  
  result <- rep(
    "Inconsistent",
    length(nUP)
  )
  
  
  # ----------------------------------------------------------
  # CONTEXT-DEPENDENT
  # ----------------------------------------------------------
  
  context <-
    nSigUP >= 1 &
    nSigDown >= 1
  
  result[context] <-
    "Context-dependent"
  
  
  # ----------------------------------------------------------
  # STRONG RECURRENT
  # ----------------------------------------------------------
  
  strong_up <-
    !context &
    nUP >= direction_threshold &
    nSigUP >= strong_sig_threshold &
    nSigDown == 0
  
  strong_down <-
    !context &
    nDOWN >= direction_threshold &
    nSigDown >= strong_sig_threshold &
    nSigUP == 0
  
  result[
    strong_up |
      strong_down
  ] <- "Strong recurrent"
  
  
  # ----------------------------------------------------------
  # DIRECTIONAL
  # ----------------------------------------------------------
  
  unresolved <-
    result == "Inconsistent"
  
  directional_up <-
    unresolved &
    nUP >= direction_threshold &
    nSigDown == 0
  
  directional_down <-
    unresolved &
    nDOWN >= direction_threshold &
    nSigUP == 0
  
  result[
    directional_up |
      directional_down
  ] <- "Directional"
  
  
  # ----------------------------------------------------------
  # DATASET-SPECIFIC
  # ----------------------------------------------------------
  
  unresolved <-
    result == "Inconsistent"
  
  dataset_specific <-
    unresolved &
    (
      nSigUP +
        nSigDown ==
        1
    )
  
  result[dataset_specific] <-
    "Dataset-specific"
  
  
  result
}


read_matching_csv <- function(
    root,
    schema,
    expected_rows = 50
) {
  
  files <- list.files(
    root,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop(
      paste(
        "No CSV files found under:",
        root
      )
    )
  }
  
  
  for (f in files) {
    
    x <- tryCatch(
      
      read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    if (is.null(x)) {
      next
    }
    
    
    if (
      !is.null(expected_rows) &&
      nrow(x) != expected_rows
    ) {
      next
    }
    
    
    found <- vapply(
      
      schema,
      
      function(aliases) {
        
        find_col(
          colnames(x),
          aliases
        )
      },
      
      character(1)
    )
    
    
    if (all(!is.na(found))) {
      
      return(
        list(
          data = x,
          file = f,
          cols = found
        )
      )
    }
  }
  
  
  stop(
    paste(
      "No compatible CSV found under:",
      root
    )
  )
}


# ============================================================
# 2. FROZEN D1-D6 MASTER
# ============================================================

if (!file.exists(old_master_file)) {
  
  stop(
    paste(
      "Frozen D1-D6 master not found:",
      old_master_file
    )
  )
}


old_raw <- read.csv(
  old_master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


old_cols <- c(
  
  pathway =
    find_col(
      colnames(old_raw),
      c(
        "pathway",
        "hallmark",
        "ID"
      )
    ),
  
  nUP =
    find_col(
      colnames(old_raw),
      c(
        "nUP",
        "n_up",
        "up_count",
        "num_up",
        "n_up_accessions",
        "n_up_datasets"
      )
    ),
  
  nDOWN =
    find_col(
      colnames(old_raw),
      c(
        "nDOWN",
        "n_down",
        "down_count",
        "num_down",
        "n_down_accessions",
        "n_down_datasets"
      )
    ),
  
  nSigUP =
    find_col(
      colnames(old_raw),
      c(
        "nSigUP",
        "n_sig_up",
        "sig_up_count",
        "num_sig_up",
        "significant_up_n",
        "n_significant_up"
      )
    ),
  
  nSigDown =
    find_col(
      colnames(old_raw),
      c(
        "nSigDown",
        "n_sig_down",
        "sig_down_count",
        "num_sig_down",
        "significant_down_n",
        "n_significant_down"
      )
    )
)


if (any(is.na(old_cols))) {
  
  cat(
    "Columns present in frozen D1-D6 master:\n"
  )
  
  print(
    colnames(old_raw)
  )
  
  stop(
    "Could not identify required D1-D6 accession-count columns."
  )
}


old <- data.frame(
  
  pathway =
    as.character(
      old_raw[[old_cols["pathway"]]]
    ),
  
  old_nUP =
    as.integer(
      old_raw[[old_cols["nUP"]]]
    ),
  
  old_nDOWN =
    as.integer(
      old_raw[[old_cols["nDOWN"]]]
    ),
  
  old_nSigUP =
    as.integer(
      old_raw[[old_cols["nSigUP"]]]
    ),
  
  old_nSigDown =
    as.integer(
      old_raw[[old_cols["nSigDown"]]]
    ),
  
  stringsAsFactors = FALSE
)


if (nrow(old) != 50) {
  
  stop(
    "Frozen D1-D6 master does not contain 50 Hallmark pathways."
  )
}


if (
  anyDuplicated(
    old$pathway
  ) > 0
) {
  
  stop(
    "Duplicated pathways detected in frozen D1-D6 master."
  )
}


if (
  anyNA(old)
) {
  
  stop(
    "Missing values detected in required frozen D1-D6 fields."
  )
}


if (
  any(
    old$old_nUP +
    old$old_nDOWN >
    6
  )
) {
  
  stop(
    "Frozen D1-D6 direction counts exceed six accessions."
  )
}


if (
  any(
    old$old_nSigUP >
    old$old_nUP
  ) ||
  any(
    old$old_nSigDown >
    old$old_nDOWN
  )
) {
  
  stop(
    "Frozen D1-D6 significant counts exceed direction counts."
  )
}


# ============================================================
# 3. REPRODUCE LOCKED D1-D6 CLASSIFICATION
# ============================================================

old$old_class <- classify_pathways(
  
  nUP =
    old$old_nUP,
  
  nDOWN =
    old$old_nDOWN,
  
  nSigUP =
    old$old_nSigUP,
  
  nSigDown =
    old$old_nSigDown,
  
  direction_threshold =
    4,
  
  strong_sig_threshold =
    3
)


class_levels <- c(
  "Strong recurrent",
  "Directional",
  "Context-dependent",
  "Dataset-specific",
  "Inconsistent"
)


expected_old <- c(
  1,
  12,
  15,
  12,
  10
)


observed_old <- as.integer(
  table(
    factor(
      old$old_class,
      levels = class_levels
    )
  )
)


old_checkpoint <- data.frame(
  
  class =
    class_levels,
  
  expected =
    expected_old,
  
  observed =
    observed_old,
  
  stringsAsFactors = FALSE
)


cat(
  "============================================================\n"
)

cat(
  "D1-D6 FROZEN CHECKPOINT\n"
)

cat(
  "============================================================\n"
)


print(
  old_checkpoint,
  row.names = FALSE
)


if (
  !all(
    observed_old ==
    expected_old
  )
) {
  
  stop(
    "D1-D6 frozen checkpoint failed."
  )
}


cat(
  "\nD1-D6 checkpoint PASSED.\n"
)

cat(
  "Master:",
  old_master_file,
  "\n\n"
)


# ============================================================
# 4. LOAD D7
# ============================================================

d7_schema <- list(
  
  pathway =
    c(
      "pathway"
    ),
  
  M468_NES =
    c(
      "M468_NES"
    ),
  
  M468_FDR =
    c(
      "M468_FDR"
    ),
  
  SUM1315_NES =
    c(
      "SUM1315_NES"
    ),
  
  SUM1315_FDR =
    c(
      "SUM1315_FDR"
    )
)


d7_match <- read_matching_csv(
  root = d7_root,
  schema = d7_schema,
  expected_rows = 50
)


d7_raw <- d7_match$data

d7_cols <- d7_match$cols


d7 <- data.frame(
  
  pathway =
    as.character(
      d7_raw[[d7_cols["pathway"]]]
    ),
  
  M468_NES =
    as.numeric(
      d7_raw[[d7_cols["M468_NES"]]]
    ),
  
  M468_FDR =
    as.numeric(
      d7_raw[[d7_cols["M468_FDR"]]]
    ),
  
  SUM1315_NES =
    as.numeric(
      d7_raw[[d7_cols["SUM1315_NES"]]]
    ),
  
  SUM1315_FDR =
    as.numeric(
      d7_raw[[d7_cols["SUM1315_FDR"]]]
    ),
  
  stringsAsFactors = FALSE
)


if (
  anyNA(d7) ||
  anyDuplicated(
    d7$pathway
  ) > 0
) {
  
  stop(
    "D7 input failed integrity checks."
  )
}


# ============================================================
# 5. CONSERVATIVE D7 ACCESSION VOTE
# ============================================================

d7_same_direction <-
  sign(
    d7$M468_NES
  ) ==
  sign(
    d7$SUM1315_NES
  )


d7$D7_direction <-
  "Discordant"


d7$D7_direction[
  d7_same_direction &
    d7$M468_NES > 0
] <- "UP"


d7$D7_direction[
  d7_same_direction &
    d7$M468_NES < 0
] <- "DOWN"


d7$D7_significant <-
  d7_same_direction &
  d7$M468_FDR < 0.05 &
  d7$SUM1315_FDR < 0.05


d7$D7_NES <-
  NA_real_


d7$D7_NES[
  d7_same_direction
] <-
  (
    d7$M468_NES[
      d7_same_direction
    ] +
      d7$SUM1315_NES[
        d7_same_direction
      ]
  ) / 2


cat(
  "============================================================\n"
)

cat(
  "D7 CHECKPOINT\n"
)

cat(
  "============================================================\n"
)

cat(
  "D7 file:",
  d7_match$file,
  "\n"
)

cat(
  "Same-direction pathways:",
  sum(
    d7_same_direction
  ),
  "/ 50\n"
)

cat(
  "Concordant significant pathways:",
  sum(
    d7$D7_significant
  ),
  "/ 50\n\n"
)


if (
  sum(
    d7_same_direction
  ) != 33
) {
  
  stop(
    "D7 checkpoint failed: expected 33 same-direction Hallmarks."
  )
}


if (
  sum(
    d7$D7_significant
  ) != 6
) {
  
  stop(
    "D7 checkpoint failed: expected 6 concordant-significant Hallmarks."
  )
}


# ============================================================
# 6. LOAD D8
#
# First try dedicated accession-vote schema.
# If not available, use frozen primary Hallmark GSEA table.
# ============================================================

d8_vote_schema <- list(
  
  pathway =
    c(
      "pathway"
    ),
  
  NES =
    c(
      "D8_NES"
    ),
  
  FDR =
    c(
      "D8_FDR"
    )
)


d8_gsea_schema <- list(
  
  pathway =
    c(
      "ID",
      "pathway"
    ),
  
  NES =
    c(
      "NES_limma",
      "NES"
    ),
  
  FDR =
    c(
      "p.adjust_limma",
      "p.adjust",
      "FDR"
    )
)


d8_match <- tryCatch(
  
  read_matching_csv(
    root = d8_root,
    schema = d8_vote_schema,
    expected_rows = 50
  ),
  
  error = function(e) {
    NULL
  }
)


if (
  is.null(
    d8_match
  )
) {
  
  d8_match <- read_matching_csv(
    root = d8_root,
    schema = d8_gsea_schema,
    expected_rows = 50
  )
}


d8_raw <- d8_match$data

d8_cols <- d8_match$cols


d8 <- data.frame(
  
  pathway =
    as.character(
      d8_raw[[d8_cols["pathway"]]]
    ),
  
  D8_NES =
    as.numeric(
      d8_raw[[d8_cols["NES"]]]
    ),
  
  D8_FDR =
    as.numeric(
      d8_raw[[d8_cols["FDR"]]]
    ),
  
  stringsAsFactors = FALSE
)


if (
  anyNA(d8) ||
  anyDuplicated(
    d8$pathway
  ) > 0
) {
  
  stop(
    "D8 input failed integrity checks."
  )
}


d8$D8_direction <- ifelse(
  
  d8$D8_NES > 0,
  
  "UP",
  
  ifelse(
    d8$D8_NES < 0,
    "DOWN",
    "NONE"
  )
)


d8$D8_significant <-
  d8$D8_FDR < 0.05


cat(
  "============================================================\n"
)

cat(
  "D8 CHECKPOINT\n"
)

cat(
  "============================================================\n"
)

cat(
  "D8 file:",
  d8_match$file,
  "\n"
)

cat(
  "Significant primary Hallmarks:",
  sum(
    d8$D8_significant
  ),
  "/ 50\n\n"
)


if (
  sum(
    d8$D8_significant
  ) != 15
) {
  
  stop(
    "D8 checkpoint failed: expected 15 significant primary Hallmarks."
  )
}


# ============================================================
# 7. VERIFY PATHWAY SETS
# ============================================================

if (
  !setequal(
    old$pathway,
    d7$pathway
  )
) {
  
  stop(
    "D7 pathway set does not match D1-D6."
  )
}


if (
  !setequal(
    old$pathway,
    d8$pathway
  )
) {
  
  stop(
    "D8 pathway set does not match D1-D6."
  )
}


# ============================================================
# 8. MERGE IN FROZEN D1-D6 ORDER
# ============================================================

d7_keep <- d7[
  ,
  c(
    "pathway",
    "D7_NES",
    "D7_direction",
    "D7_significant"
  )
]


d8_keep <- d8[
  ,
  c(
    "pathway",
    "D8_NES",
    "D8_FDR",
    "D8_direction",
    "D8_significant"
  )
]


master <- merge(
  old,
  d7_keep,
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)


master <- merge(
  master,
  d8_keep,
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)


master$order_index <- match(
  master$pathway,
  old$pathway
)


master <- master[
  order(
    master$order_index
  ),
  ,
  drop = FALSE
]


if (
  nrow(master) != 50 ||
  anyNA(
    master$order_index
  )
) {
  
  stop(
    "Final D1-D8 pathway merge failed."
  )
}


# ============================================================
# 9. ADD D7 VOTES
# ============================================================

d7_up <-
  master$D7_direction ==
  "UP"


d7_down <-
  master$D7_direction ==
  "DOWN"


d7_sig_up <-
  d7_up &
  master$D7_significant


d7_sig_down <-
  d7_down &
  master$D7_significant


# ============================================================
# 10. ADD D8 VOTES
# ============================================================

d8_up <-
  master$D8_direction ==
  "UP"


d8_down <-
  master$D8_direction ==
  "DOWN"


d8_sig_up <-
  d8_up &
  master$D8_significant


d8_sig_down <-
  d8_down &
  master$D8_significant


# ============================================================
# 11. FINAL D1-D8 ACCESSION COUNTS
# ============================================================

master$nUP <-
  master$old_nUP +
  as.integer(
    d7_up
  ) +
  as.integer(
    d8_up
  )


master$nDOWN <-
  master$old_nDOWN +
  as.integer(
    d7_down
  ) +
  as.integer(
    d8_down
  )


master$nSigUP <-
  master$old_nSigUP +
  as.integer(
    d7_sig_up
  ) +
  as.integer(
    d8_sig_up
  )


master$nSigDown <-
  master$old_nSigDown +
  as.integer(
    d7_sig_down
  ) +
  as.integer(
    d8_sig_down
  )


if (
  any(
    master$nUP +
    master$nDOWN >
    8
  )
) {
  
  stop(
    "Final direction counts exceed eight accessions."
  )
}


if (
  any(
    master$nSigUP >
    master$nUP
  ) ||
  any(
    master$nSigDown >
    master$nDOWN
  )
) {
  
  stop(
    "Final significant counts exceed direction counts."
  )
}


# ============================================================
# 12. FROZEN FINAL N=8 THRESHOLDS
# ============================================================

N_accessions <-
  8


direction_threshold <-
  ceiling(
    2 *
      N_accessions /
      3
  )


strong_sig_threshold <-
  ceiling(
    N_accessions /
      2
  )


if (
  direction_threshold != 6 ||
  strong_sig_threshold != 4
) {
  
  stop(
    "Frozen N=8 thresholds were not reproduced."
  )
}


cat(
  "============================================================\n"
)

cat(
  "FINAL FROZEN THRESHOLDS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Eligible accessions:",
  N_accessions,
  "\n"
)

cat(
  "Directional threshold:",
  direction_threshold,
  "/ 8\n"
)

cat(
  "Strong significant threshold:",
  strong_sig_threshold,
  "/ 8\n\n"
)


# ============================================================
# 13. FINAL CLASSIFICATION
# ============================================================

master$final_class <- classify_pathways(
  
  nUP =
    master$nUP,
  
  nDOWN =
    master$nDOWN,
  
  nSigUP =
    master$nSigUP,
  
  nSigDown =
    master$nSigDown,
  
  direction_threshold =
    direction_threshold,
  
  strong_sig_threshold =
    strong_sig_threshold
)


master$dominant_direction <- ifelse(
  
  master$nUP >
    master$nDOWN,
  
  "UP",
  
  ifelse(
    master$nDOWN >
      master$nUP,
    "DOWN",
    "Tie"
  )
)


master$max_direction_votes <- pmax(
  master$nUP,
  master$nDOWN
)


master$direction_support_fraction <-
  master$max_direction_votes /
  N_accessions


master$total_significant_accessions <-
  master$nSigUP +
  master$nSigDown


# ============================================================
# 14. FINAL CLASS SUMMARY
# ============================================================

final_counts <- as.integer(
  table(
    factor(
      master$final_class,
      levels = class_levels
    )
  )
)


final_summary <- data.frame(
  
  class =
    class_levels,
  
  n =
    final_counts,
  
  stringsAsFactors = FALSE
)


cat(
  "============================================================\n"
)

cat(
  "FINAL D1-D8 CLASSIFICATION SUMMARY\n"
)

cat(
  "============================================================\n"
)


print(
  final_summary,
  row.names = FALSE
)


# ============================================================
# 15. STRONG RECURRENT PATHWAYS
# ============================================================

strong <- master[
  master$final_class ==
    "Strong recurrent",
  ,
  drop = FALSE
]


if (
  nrow(strong) > 0
) {
  
  strong <- strong[
    order(
      -strong$max_direction_votes,
      -strong$total_significant_accessions
    ),
    ,
    drop = FALSE
  ]
}


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "STRONG RECURRENT PATHWAYS\n"
)

cat(
  "============================================================\n"
)


if (
  nrow(strong) == 0
) {
  
  cat(
    "NONE\n"
  )
  
} else {
  
  print(
    
    strong[
      ,
      c(
        "pathway",
        "dominant_direction",
        "nUP",
        "nDOWN",
        "nSigUP",
        "nSigDown",
        "old_class",
        "final_class"
      )
    ],
    
    row.names = FALSE
  )
}


# ============================================================
# 16. DIRECTIONAL PATHWAYS
# ============================================================

directional <- master[
  master$final_class ==
    "Directional",
  ,
  drop = FALSE
]


if (
  nrow(directional) > 0
) {
  
  directional <- directional[
    order(
      -directional$max_direction_votes,
      -directional$total_significant_accessions
    ),
    ,
    drop = FALSE
  ]
}


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "DIRECTIONAL PATHWAYS\n"
)

cat(
  "============================================================\n"
)


if (
  nrow(directional) == 0
) {
  
  cat(
    "NONE\n"
  )
  
} else {
  
  print(
    
    directional[
      ,
      c(
        "pathway",
        "dominant_direction",
        "nUP",
        "nDOWN",
        "nSigUP",
        "nSigDown",
        "old_class",
        "final_class"
      )
    ],
    
    row.names = FALSE
  )
}


# ============================================================
# 17. CONTEXT-DEPENDENT PATHWAYS
# ============================================================

context <- master[
  master$final_class ==
    "Context-dependent",
  ,
  drop = FALSE
]


if (
  nrow(context) > 0
) {
  
  context <- context[
    order(
      -context$total_significant_accessions,
      -context$max_direction_votes
    ),
    ,
    drop = FALSE
  ]
}


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "CONTEXT-DEPENDENT PATHWAYS\n"
)

cat(
  "============================================================\n"
)


if (
  nrow(context) == 0
) {
  
  cat(
    "NONE\n"
  )
  
} else {
  
  print(
    
    context[
      ,
      c(
        "pathway",
        "nUP",
        "nDOWN",
        "nSigUP",
        "nSigDown",
        "old_class",
        "final_class"
      )
    ],
    
    row.names = FALSE
  )
}


# ============================================================
# 18. OLD -> NEW CLASS TRANSITIONS
# ============================================================

transition_table <- as.data.frame(
  table(
    
    Old_class =
      factor(
        master$old_class,
        levels = class_levels
      ),
    
    Final_class =
      factor(
        master$final_class,
        levels = class_levels
      )
  )
)


transition_table <- transition_table[
  transition_table$Freq >
    0,
  ,
  drop = FALSE
]


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "D1-D6 -> D1-D8 CLASS TRANSITIONS\n"
)

cat(
  "============================================================\n"
)


print(
  transition_table,
  row.names = FALSE
)


# ============================================================
# 19. PATHWAYS WHOSE CLASS CHANGED
# ============================================================

changed <- master[
  master$old_class !=
    master$final_class,
  ,
  drop = FALSE
]


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "PATHWAYS WHOSE CLASS CHANGED AFTER D7 + D8\n"
)

cat(
  "============================================================\n"
)

cat(
  "Changed:",
  nrow(changed),
  "/ 50\n\n"
)


if (
  nrow(changed) > 0
) {
  
  print(
    
    changed[
      ,
      c(
        "pathway",
        "old_class",
        "final_class",
        "old_nUP",
        "old_nDOWN",
        "old_nSigUP",
        "old_nSigDown",
        "D7_direction",
        "D7_significant",
        "D8_direction",
        "D8_significant",
        "nUP",
        "nDOWN",
        "nSigUP",
        "nSigDown"
      )
    ],
    
    row.names = FALSE
  )
}


# ============================================================
# 20. ORIGINAL D1 FROZEN PATHWAYS
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


frozen_index <- match(
  frozen_pathways,
  master$pathway
)


if (
  anyNA(
    frozen_index
  )
) {
  
  stop(
    "One or more frozen D1 pathways are missing."
  )
}


frozen <- master[
  frozen_index,
  ,
  drop = FALSE
]


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "D1 FROZEN PATHWAYS - FINAL D1-D8 STATUS\n"
)

cat(
  "============================================================\n"
)


print(
  
  frozen[
    ,
    c(
      "pathway",
      "old_nUP",
      "old_nDOWN",
      "old_nSigUP",
      "old_nSigDown",
      "D7_NES",
      "D7_direction",
      "D7_significant",
      "D8_NES",
      "D8_direction",
      "D8_significant",
      "nUP",
      "nDOWN",
      "nSigUP",
      "nSigDown",
      "dominant_direction",
      "old_class",
      "final_class"
    )
  ],
  
  row.names = FALSE
)


# ============================================================
# 21. GLYCOLYSIS SPECIFIC CHECK
# ============================================================

glycolysis <- master[
  master$pathway ==
    "HALLMARK_GLYCOLYSIS",
  ,
  drop = FALSE
]


cat("\n")
cat(
  "============================================================\n"
)

cat(
  "GLYCOLYSIS CHECK\n"
)

cat(
  "============================================================\n"
)


print(
  
  glycolysis[
    ,
    c(
      "pathway",
      "old_nUP",
      "old_nDOWN",
      "old_nSigUP",
      "old_nSigDown",
      "D7_NES",
      "D7_direction",
      "D7_significant",
      "D8_NES",
      "D8_direction",
      "D8_significant",
      "nUP",
      "nDOWN",
      "nSigUP",
      "nSigDown",
      "old_class",
      "final_class"
    )
  ],
  
  row.names = FALSE
)


# ============================================================
# 22. SAVE OUTPUTS
# ============================================================

write.csv(
  
  master,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_accession_Hallmark_master.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  final_summary,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_class_summary.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  transition_table,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_class_transitions.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  changed,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_changed_pathways.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  
  frozen,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_D1_frozen_pathways.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 23. DATASET MANIFEST
# ============================================================

dataset_manifest <- data.frame(
  
  dataset =
    c(
      "D1",
      "D2",
      "D3",
      "D4",
      "D5",
      "D6",
      "D7",
      "D8",
      "D9_candidate"
    ),
  
  accession =
    c(
      "GSE235980",
      "GSE117765",
      "GSE153867",
      "GSE342642",
      "GSE202434",
      "GSE191231",
      "GSE165914",
      "GSE163854",
      "GSE287274"
    ),
  
  included =
    c(
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      FALSE
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  dataset_manifest,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_dataset_manifest.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 24. SAVE PRE-SPECIFIED THRESHOLDS
# ============================================================

threshold_table <- data.frame(
  
  parameter =
    c(
      "Eligible_accessions",
      "Directional_support_threshold",
      "Strong_significant_support_threshold"
    ),
  
  value =
    c(
      8,
      6,
      4
    ),
  
  rule =
    c(
      "N",
      "ceil(2N/3)",
      "ceil(N/2)"
    ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  threshold_table,
  
  file.path(
    out_dir,
    "FINAL_D1_D8_frozen_thresholds.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 25. SESSION INFO
# ============================================================

sink(
  file.path(
    out_dir,
    "FINAL_D1_D8_sessionInfo.txt"
  )
)


sessionInfo()


sink()


# ============================================================
# 26. FINAL MESSAGE
# ============================================================

cat("\n")
cat(
  "============================================================\n"
)

cat(
  "FINAL D1-D8 ACCESSION INTEGRATION COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "D1-D6 master:",
  old_master_file,
  "\n"
)

cat(
  "D7 source:",
  d7_match$file,
  "\n"
)

cat(
  "D8 source:",
  d8_match$file,
  "\n"
)

cat(
  "Changed pathway classes:",
  nrow(changed),
  "/ 50\n"
)

cat(
  "Output directory:",
  out_dir,
  "\n"
)

cat(
  "============================================================\n"
)
