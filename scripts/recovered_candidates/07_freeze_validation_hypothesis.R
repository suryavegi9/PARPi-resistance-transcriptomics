dir.create(
  "hypotheses",
  showWarnings = FALSE
)

hypothesis <- c(
  
  "FROZEN BEFORE EXTERNAL VALIDATION",
  
  "Discovery dataset: GSE235980",
  
  "",
  
  "Primary resistance-enriched predictions:",
  "1. HALLMARK_OXIDATIVE_PHOSPHORYLATION: positive NES",
  "2. HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION: positive NES",
  
  "",
  
  "Primary resistance-depleted predictions:",
  "3. HALLMARK_G2M_CHECKPOINT: negative NES",
  "4. HALLMARK_E2F_TARGETS: negative NES",
  "5. HALLMARK_MTORC1_SIGNALING: negative NES",
  "6. HALLMARK_MYC_TARGETS_V1: negative NES",
  "7. HALLMARK_GLYCOLYSIS: negative NES",
  
  "",
  
  "Additional discovery signals:",
  "8. HALLMARK_INFLAMMATORY_RESPONSE: negative NES",
  "9. HALLMARK_TNFA_SIGNALING_VIA_NFKB: negative NES",
  "10. HALLMARK_UV_RESPONSE_DN: negative NES",
  
  "",
  
  "Negative/context-specific observations:",
  "- HALLMARK_DNA_REPAIR was not significantly shared.",
  "- HALLMARK_HYPOXIA was context-specific.",
  
  "",
  
  "Validation rule:",
  "Do not alter these predictions after viewing external datasets.",
  "Report concordant and discordant pathways."
)

writeLines(
  hypothesis,
  "hypotheses/GSE235980_frozen_validation_hypothesis.txt"
)
