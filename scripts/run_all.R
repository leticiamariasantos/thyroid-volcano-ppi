# ═══════════════════════════════════════════════════════════════════════════════
# run_all.R — Orquestrador principal (ETAPA 25). Idempotente; para em erro crítico.
#
# Fluxo (discovery já congelado): 01 acquire -> 02 QC -> 03 DEG -> 04 GSEA ->
# 05 redundancy -> 06 PPI -> 07 prioritization -> 09 composition -> 10 comp sensitivity
# -> 11 counts -> 12 counts DE -> 13 comparison -> 19 composition formal ->
# 20 ppi audit -> 21 prioritization audit -> 25 hallmark/reactome ->
# 26 composition covariate -> 15 validate outputs -> (validation pipeline separado)
#
# NOTA: scripts 16-24 (validação GEO) são orquestrados por run_validation_pipeline.R.
# ═══════════════════════════════════════════════════════════════════════════════════

RS <- file.path(R.home("bin"), "Rscript.exe")
root <- getwd()

scripts <- c(
  "01_acquire_global_data.R", "02_qc.R", "03_deg.R", "04_gsea.R",
  "05_pathway_redundancy.R", "06_ppi.R", "07_prioritization.R",
  "09_composition.R", "10_composition_sensitivity.R",
  "11_count_recount3.R", "12_counts_deg.R", "13_comparison.R",
  "19_composition_formal.R", "20_ppi_audit.R", "21_prioritization_audit.R",
  "25_gsea_hallmark_reactome.R", "26_composition_covariate.R",
  "15_validate_outputs.R"
)

run <- function(s) {
  p <- file.path("scripts", s)
  if (!file.exists(p)) { cat(sprintf("SKIP (ausente): %s\n", s)); return(invisible()) }
  cat(sprintf("\n════ RUN %s ════\n", s))
  rc <- system2(RS, p, wait = TRUE)
  if (rc != 0) { cat(sprintf("ERRO CRÍTICO em %s (exit %d)\n", s, rc)); quit(status = rc) }
}

for (s in scripts) run(s)
cat("\n══ run_all.R CONCLUÍDO ══\n")
