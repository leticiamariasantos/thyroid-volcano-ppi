# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat/test-phase2.R — Sanity + regressão da Fase 2 (reboot)
# Testes de regressão dos números-chave (sem internet, sem dados grandes).
# ═══════════════════════════════════════════════════════════════════════════════
context("Fase 2 — sanity + regressão")

test_that("outputs-chave da Fase 2 existem", {
  expect_true(file.exists(here::here("results/phase2/reports/PHASE2_FINAL_REPORT.md")))
  expect_true(file.exists(here::here("results/phase2/pathways/PANEL_30_PATHWAYS.tsv")))
  expect_true(file.exists(here::here("results/phase2/pathways/PANEL_SELECTION_RATIONALE.md")))
  expect_true(file.exists(here::here("results/phase2/ppi/PPI_summary.tsv")))
  expect_true(file.exists(here::here("results/phase2/reports/PHASE2_LANGUAGE_LOCK_AUDIT.md")))
})

test_that("painel tem exatamente 30 vias (10 + 20)", {
  p <- utils::read.delim(here::here("results/phase2/pathways/PANEL_30_PATHWAYS.tsv"),
                         check.names = FALSE)
  expect_equal(nrow(p), 30)
  expect_equal(sum(p$panel_group == "ORIGINAL_10"), 10)
  expect_equal(sum(p$panel_group == "ADDITIONAL_20"), 20)
})

test_that("regressão: limma reporta 12.200 DEGs", {
  d <- utils::read.delim(here::here("results/phase2/differential_expression/DEG_summary.tsv"),
                         check.names = FALSE)
  expect_equal(d$total_deg[d$method == "limma"], 12200)
  expect_equal(d$genes_tested[d$method == "limma"], 26011)
})

test_that("regressão: 6 vias robustas (escore 4/4)", {
  r <- utils::read.delim(here::here("results/phase2/pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv"),
                         check.names = FALSE)
  expect_equal(sum(r$robustness_class == "ROBUSTA"), 6)
  robustas <- sort(r$database_id[r$robustness_class == "ROBUSTA"])
  expect_equal(robustas, sort(c("hsa00190", "hsa03030", "hsa03050",
                                "hsa04110", "hsa04115", "hsa04612")))
})

test_that("regressão: PPI tem FN1 como hub principal", {
  s <- utils::read.delim(here::here("results/phase2/ppi/PPI_summary.tsv"),
                         check.names = FALSE)
  expect_equal(s$value[s$metric == "top_hub_gene"], "FN1")
  expect_equal(as.integer(s$value[s$metric == "nodes"]), 70)
})

test_that("candidatos ITGA2/FN1/CCND1 presentes na convergência", {
  c <- utils::read.delim(here::here("results/phase2/validation/convergence_candidates.tsv"),
                         check.names = FALSE)
  expect_true(all(c("ITGA2", "FN1", "CCND1") %in% c$gene))
  expect_true(all(sign(c$logFC_limma) > 0))  # os três são Up no tumor
})
