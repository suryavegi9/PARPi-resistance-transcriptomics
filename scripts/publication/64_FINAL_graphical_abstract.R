
# ============================================================
# 64_FINAL_graphical_abstract.R
# PARPi resistance — Graphical abstract
# Generates artwork from frozen publication outputs
# No differential expression or GSEA is performed
# ============================================================

library(grid)

# ------------------------------------------------------------
# 1. PROJECT DIRECTORY
# ------------------------------------------------------------

# Run from the original project root
project_root <- Sys.getenv(
  "PARPI_PROJECT_ROOT",
  unset = getwd()
)

source_dir <- file.path(
  project_root,
  "results",
  "FINAL_D1_D8_integration",
  "Figure8_FINAL"
)

output_dir <- file.path(
  project_root,
  "publication",
  "graphical_abstract"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. LOAD FROZEN SOURCE TABLES
# ------------------------------------------------------------

accession <- read.csv(
  file.path(
    source_dir,
    "Figure8_accession_source_table.csv"
  )
)

model <- read.csv(
  file.path(
    source_dir,
    "Figure8_model_level_source_table.csv"
  )
)

treatment <- read.csv(
  file.path(
    source_dir,
    "Figure8_treatment_state_source_table.csv"
  )
)

checkpoint <- read.csv(
  file.path(
    source_dir,
    "Figure8_FINAL_checkpoint.csv"
  )
)

# ------------------------------------------------------------
# 3. SOURCE INTEGRITY CHECKS
# ------------------------------------------------------------

stopifnot(
  nrow(accession) == 50,
  nrow(model) == 50,
  nrow(treatment) == 50,
  !anyDuplicated(accession$pathway),
  !anyDuplicated(model$pathway),
  !anyDuplicated(treatment$pathway),
  setequal(accession$pathway, model$pathway),
  setequal(accession$pathway, treatment$pathway)
)

# ------------------------------------------------------------
# 4. EXTRACT FROZEN RESULTS
# ------------------------------------------------------------

n_accessions <- 8
n_units <- 12
n_hypotheses <- 10
n_hallmarks <- nrow(accession)

n_context_accession <- sum(
  accession$final_class == "Context-dependent"
)

n_context_model <- sum(
  model$model12_class == "Context-dependent"
)

n_strong_accession <- checkpoint$value[
  checkpoint$metric == "Accession_strong"
]

n_strong_model <- checkpoint$value[
  checkpoint$metric == "Model_strong"
]

n_robust_D5 <- sum(treatment$robust_D5)

n_agree_D5 <- sum(
  treatment$robust_D5 &
    treatment$same_exposure_direction
)

rho <- checkpoint$value[
  checkpoint$metric == "Treatment_Spearman_rho"
]

# Confirm agreement with frozen publication results

stopifnot(
  n_context_accession == 24,
  n_context_model == 35,
  n_strong_accession == 0,
  n_strong_model == 0,
  n_robust_D5 == 18,
  n_agree_D5 == 9,
  abs(rho - (-0.255702280912365)) < 1e-10
)

# ------------------------------------------------------------
# 5. COLOUR PALETTE
# ------------------------------------------------------------

col <- list(
  navy = "#17384A",
  blue = "#226A91",
  teal = "#267C83",
  text = "#17283B",
  muted = "#526678",
  panel = "#F3F6F9",
  border = "#CBD5E1",
  track = "#DCE5EB",
  white = "#FFFFFF"
)

# ------------------------------------------------------------
# 6. GRAPHICAL FUNCTIONS
# ------------------------------------------------------------

txt <- function(
    label, x, y,
    size = 11,
    colour = col$text,
    bold = FALSE,
    just = "centre"
) {
  grid.text(
    label,
    x = x,
    y = y,
    just = just,
    gp = gpar(
      fontsize = size,
      col = colour,
      fontfamily = "sans",
      fontface = if (bold) "bold" else "plain"
    )
  )
}

panel <- function(x, heading) {
  
  grid.roundrect(
    x = x,
    y = 0.485,
    width = 0.29,
    height = 0.65,
    r = unit(0.025, "snpc"),
    gp = gpar(
      fill = col$panel,
      col = col$border,
      lwd = 1
    )
  )
  
  txt(
    heading,
    x = x,
    y = 0.765,
    size = 13,
    bold = TRUE
  )
}

progress_bar <- function(
    x, y, count, total,
    width = 0.24
) {
  grid.rect(
    x = x,
    y = y,
    width = width,
    height = 0.032,
    just = c("left", "centre"),
    gp = gpar(
      fill = col$track,
      col = NA
    )
  )
  
  grid.rect(
    x = x,
    y = y,
    width = width * count / total,
    height = 0.032,
    just = c("left", "centre"),
    gp = gpar(
      fill = col$teal,
      col = NA
    )
  )
}

connector <- function(x1, x2) {
  grid.lines(
    x = c(x1, x2),
    y = c(0.485, 0.485),
    arrow = arrow(
      length = unit(2, "mm"),
      type = "closed"
    ),
    gp = gpar(
      col = col$muted,
      fill = col$muted,
      lwd = 1.3
    )
  )
}

# ------------------------------------------------------------
# 7. DRAW GRAPHICAL ABSTRACT
# ------------------------------------------------------------

draw_abstract <- function() {
  
  grid.newpage()
  
  grid.rect(
    gp = gpar(
      fill = col$white,
      col = NA
    )
  )
  
  # Main heading
  
  txt(
    "Acquired PARP inhibitor resistance",
    x = 0.5,
    y = 0.945,
    size = 19,
    bold = TRUE,
    colour = col$navy
  )
  
  txt(
    "Cross-dataset transcriptomic analysis",
    x = 0.5,
    y = 0.885,
    size = 11,
    colour = col$muted
  )
  
  # PANEL 1 — STUDY DESIGN
  
  panel(0.165, "Study design")
  
  txt(n_accessions, 0.165, 0.675,
      size = 27, bold = TRUE,
      colour = col$blue)
  
  txt("GEO accessions", 0.165, 0.615)
  
  txt(n_units, 0.165, 0.505,
      size = 24, bold = TRUE)
  
  txt("Resistance-state units",
      0.165, 0.445, size = 10)
  
  txt(n_hypotheses, 0.165, 0.335,
      size = 24, bold = TRUE)
  
  txt("Frozen discovery hypotheses",
      0.165, 0.275, size = 10)
  
  # PANEL 2 — RECURRENCE
  
  panel(0.5, "Hallmark recurrence")
  
  txt(
    paste0(n_strong_accession, " / ", n_hallmarks),
    0.5, 0.675,
    size = 27,
    bold = TRUE
  )
  
  txt("Strongly recurrent pathways",
      0.5, 0.615, size = 10)
  
  txt("Context-dependent pathways",
      0.5, 0.535,
      size = 11,
      bold = TRUE)
  
  txt("8 accessions",
      0.38, 0.465,
      size = 10, just = "left")
  
  txt(
    paste0(n_context_accession, " / 50"),
    0.62, 0.465,
    size = 11,
    bold = TRUE,
    just = "right"
  )
  
  progress_bar(
    0.38, 0.420,
    n_context_accession, n_hallmarks
  )
  
  txt("12 units",
      0.38, 0.335,
      size = 10, just = "left")
  
  txt(
    paste0(n_context_model, " / 50"),
    0.62, 0.335,
    size = 11,
    bold = TRUE,
    just = "right"
  )
  
  progress_bar(
    0.38, 0.290,
    n_context_model, n_hallmarks
  )
  
  # PANEL 3 — TREATMENT STATE
  
  panel(0.835, "Treatment state")
  
  txt(
    n_robust_D5, 0.835, 0.675,
    size = 27,
    bold = TRUE,
    colour = col$blue
  )
  
  txt("Robust exposure programs",
      0.835, 0.615, size = 10)
  
  txt(
    paste0(n_agree_D5, " / ", n_robust_D5),
    0.835, 0.505,
    size = 24, bold = TRUE
  )
  
  txt(
    "Agree with observational direction",
    0.835, 0.445,
    size = 9
  )
  
  txt(
    paste0("Spearman rho = ",
           sprintf("%.3f", rho)),
    0.835, 0.335,
    size = 15,
    bold = TRUE
  )
  
  txt("Across all 50 Hallmarks",
      0.835, 0.275, size = 10)
  
  # Directional connectors
  
  connector(0.313, 0.352)
  connector(0.648, 0.687)
  
  # MAIN CONCLUSION
  
  grid.roundrect(
    x = 0.5,
    y = 0.085,
    width = 0.96,
    height = 0.105,
    r = unit(0.018, "snpc"),
    gp = gpar(
      fill = col$navy,
      col = NA
    )
  )
  
  txt(
    "Context-dependent transcriptional remodeling predominates",
    0.5, 0.085,
    size = 13,
    colour = col$white,
    bold = TRUE
  )
}

# ------------------------------------------------------------
# 8. EXPORT VECTOR PDF
# ------------------------------------------------------------

if (!capabilities("cairo")) {
  stop("Cairo graphics support is required for PDF export.")
}

cairo_pdf(
  filename = file.path(
    output_dir,
    "Graphical_abstract.pdf"
  ),
  width = 10,
  height = 4.7,
  family = "sans",
  bg = "white"
)

draw_abstract()
dev.off()

# ------------------------------------------------------------
# 9. EXPORT HIGH-RESOLUTION PNG
# ------------------------------------------------------------

png(
  filename = file.path(
    output_dir,
    "Graphical_abstract.png"
  ),
  width = 10,
  height = 4.7,
  units = "in",
  res = 600,
  type = "cairo",
  bg = "white"
)

draw_abstract()
dev.off()

# ------------------------------------------------------------
# 10. EXPORT TIFF
# ------------------------------------------------------------

tiff(
  filename = file.path(
    output_dir,
    "Graphical_abstract.tiff"
  ),
  width = 10,
  height = 4.7,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

draw_abstract()
dev.off()

# Record the actual generation environment

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_graphical_abstract.txt"
  )
)

cat("\nGraphical abstract generated successfully.\n")
cat("Output directory:", output_dir, "\n")

print(list.files(output_dir))
