# ═══════════════════════════════════════════════════════════════════════════════
# run_validation_pipeline.R — PARTE 27: orquestrador do pipeline de validação
#
# Fluxo (idempotente; para em erro crítico):
#   16 download/import -> 17 validação DE -> 22 metadata/composição/integração ->
#   23 vias -> 24 integridade
# (19–21 são auditorias discovery, não dependem de GEO)
# ═══════════════════════════════════════════════════════════════════════════════

RS <- file.path(R.home("bin"), "Rscript.exe")
scripts <- c(
  "16_geo_download.R",
  "17_geo_validate.R",
  "22_validation_enhance.R",
  "23_validation_pathway.R",
  "24_integrity_tests.R"
)
root <- getwd()
for (s in scripts) {
  cat(sprintf("\n════ RUN %s ════\n", s))
  rc <- system2(RS, file.path("scripts", s), wait = TRUE)
  if (rc != 0) {
    cat(sprintf("ERRO CRÍTICO em %s (exit %d). Pipeline interrompido.\n", s, rc))
    quit(status = rc)
  }
}
cat("\n══ VALIDATION PIPELINE COMPLETO ══\n")
