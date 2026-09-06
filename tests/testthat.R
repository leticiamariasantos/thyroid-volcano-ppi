# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat.R — Test runner da Fase 2 (reboot)
# Uso: testthat::test_dir("tests/testthat")
# ═══════════════════════════════════════════════════════════════════════════════
library(testthat)
library(here)
test_dir(here::here("tests", "testthat"))
