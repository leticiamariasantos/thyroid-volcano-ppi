# ═══════════════════════════════════════════════════════════════════════════════
# setup_renv.R — Restaura o ambiente R a partir do renv.lock (NÃO-destrutivo)
#
# ATENÇÃO: este script NÃO re-inicializa o renv nem sobrescreve o renv.lock.
# Ele apenas restaura os pacotes nas versões fixadas no lockfile.
# ═══════════════════════════════════════════════════════════════════════════════

if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv", repos = "https://cloud.r-project.org")

renv::restore()

lock <- renv::lockfile_read("renv.lock")
cat("Ambiente restaurado a partir de renv.lock (", length(lock$Packages), " pacotes).\n", sep = "")
