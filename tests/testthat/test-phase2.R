# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat/test-phase2.R — Sanity checks da Fase 2 (reboot)
# ═══════════════════════════════════════════════════════════════════════════════
context("Fase 2 — sanity checks")

test_that("outputs-chave da Fase 2 existem", {
  expect_true(file.exists(here::here("results/phase2/reports/PHASE2_FINAL_REPORT.md")))
  expect_true(file.exists(here::here("results/phase2/pathways/PANEL_30_PATHWAYS.tsv")))
  expect_true(file.exists(here::here("results/phase2/pathways/PANEL_SELECTION_RATIONALE.md")))
  expect_true(file.exists(here::here("results/phase2/ppi/PPI_summary.tsv")))
  expect_true(file.exists(here::here("results/phase2/reports/PHASE2_LANGUAGE_LOCK_AUDIT.md")))
})

test_that("painel tem exatamente 30 vias (10 + 20)", {
  p <- utils::read.delim(here::here("results/phase2/pathways/PANEL_30_PATHWAYS.tsv"))
  expect_equal(nrow(p), 30)
  expect_equal(sum(p$panel_group == "ORIGINAL_10"), 10)
  expect_equal(sum(p$panel_group == "ADDITIONAL_20"), 20)
})

test_that("PPI tem FN1 como hub principal", {
  s <- utils::read.delim(here::here("results/phase2/ppi/PPI_summary.tsv"))
  expect_equal(s$value[s$metric == "top_hub_gene"], "FN1")
})
