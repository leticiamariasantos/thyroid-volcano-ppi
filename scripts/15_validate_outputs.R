# ═══════════════════════════════════════════════════════════════════════════════
# 15_validate_outputs.R — Validação técnica (FASE 2 — reboot)
#
# ADAPTADO para a nova estrutura da fase 2. Este script agora delega à validação
# canônica em scripts/phase2/17_validate.R. A versão anterior validava a estrutura
# antiga (painel de 10 vias), removida nesta reconstrução.
#
# Execução: Rscript scripts/15_validate_outputs.R
# Objetivo: 0 failures.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))

cat("══ Validação de outputs (FASE 2 — reboot) ══\n")
cat("Delegando para scripts/phase2/17_validate.R ...\n\n")

source(here("scripts", "phase2", "17_validate.R"))
