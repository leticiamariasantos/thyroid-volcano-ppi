# ═══════════════════════════════════════════════════════════════════════════════
# run_phase2.R — Orquestrador da FASE 2 (reboot) — painel de 30 vias
# Executa os scripts scripts/phase2/ em ordem. Para em erro crítico.
# Uso: Rscript scripts/run_phase2.R
# ═══════════════════════════════════════════════════════════════════════════════
RS <- file.path(R.home("bin"), "Rscript.exe")
steps <- c(
  "00_config.R", "01_audit_data.R", "02_qc.R", "03_de.R", "04_panel.R",
  "05_gsea_global.R", "06_gsea_panel.R", "07_redundancy.R", "08_composition.R",
  "09_ppi.R", "10_validation.R", "11_singlecell.R", "12_protein_mutation.R",
  "13_robustness.R", "14_convergence.R", "15_figures.R", "16_validate.R"
)
for (s in steps) {
  p <- file.path("scripts", "phase2", s)
  cat("══ Executando", p, "══\n")
  st <- system2(RS, args = p, stdout = "", stderr = "")
  if (st != 0) stop("Falha em ", p, " (exit ", st, ")")
}
cat("══ FASE 2 concluída (0 falhas) ══\n")
