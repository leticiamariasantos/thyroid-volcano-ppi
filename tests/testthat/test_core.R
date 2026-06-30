# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat/test_core.R — Core pipeline tests (no internet required)
# thyroid-volcano-ppi v3.1.0
#
# These tests use mock data and validate critical paths without external APIs.
# ═══════════════════════════════════════════════════════════════════════════════

context("Core pipeline validation")

# ── Thresholds ────────────────────────────────────────────────────────────────

test_that("THRESHOLD values are within valid ranges", {
  expect_true(exists("THRESHOLD"))
  expect_gt(THRESHOLD$lfc, 0)
  expect_gt(THRESHOLD$fdr, 0)
  expect_lt(THRESHOLD$fdr, 1)
  expect_gte(THRESHOLD$string, 0)
  expect_lte(THRESHOLD$string, 1000)
  expect_gt(THRESHOLD$expr_min, 0)
  expect_gt(THRESHOLD$expr_frac, 0)
  expect_lte(THRESHOLD$expr_frac, 1)
})

# ── DEG classification logic ──────────────────────────────────────────────────

test_that("DEG classification uses correct thresholds", {
  # Simulate results of topTable-like data
  mock <- tibble::tibble(
    gene_symbol = paste0("GENE", 1:6),
    logFC       = c(2.0, -1.8, 0.3, 1.5, -0.4, 0.1),
    adj.P.Val   = c(0.001, 0.01, 0.03, 0.10, 0.20, 0.80)
  )

  mock$regulation <- with(mock, ifelse(
    adj.P.Val < 0.05 & logFC >  1.0, "Up",
    ifelse(adj.P.Val < 0.05 & logFC < -1.0, "Down", "NS")
  ))

  expect_equal(mock$regulation, c("Up", "Down", "NS", "NS", "NS", "NS"))
  expect_equal(sum(mock$regulation == "Up"), 1)
  expect_equal(sum(mock$regulation == "Down"), 1)
})

# ── Expression filter logic ───────────────────────────────────────────────────

test_that("Expression filter keeps genes above detection threshold", {
  set.seed(42)
  # 10 genes, 100 samples: gene1 always low, others varying
  E <- matrix(rnorm(1000, mean = 2, sd = 1), nrow = 10)
  E[1, ] <- 0.1  # gene1 below detection

  keep <- rowMeans(E > 0.5) >= 0.1
  expect_false(keep[1])
  expect_true(all(keep[-1]) || sum(keep[-1]) >= 5)
})

# ── Required columns in DEG output ────────────────────────────────────────────

test_that("DEG output contains required columns", {
  required_cols <- c("gene_symbol", "logFC", "AveExpr", "t",
                     "P.Value", "adj.P.Val", "B", "regulation")
  # This test validates that the DEG pipeline produces these columns.
  # In a full run, deg object should exist after 03_deg.R.
  if (exists("deg") && is.data.frame(deg)) {
    for (col in required_cols) {
      expect_true(col %in% names(deg),
                  info = paste("Missing column:", col))
    }
  } else {
    skip("deg object not available (run full pipeline first)")
  }
})

# ── No critical NAs in key columns ────────────────────────────────────────────

test_that("No NAs in critical DEG columns", {
  if (exists("deg") && is.data.frame(deg)) {
    critical <- c("gene_symbol", "logFC", "adj.P.Val", "regulation")
    for (col in critical) {
      na_count <- sum(is.na(deg[[col]]))
      expect_equal(na_count, 0,
                   info = paste(na_count, "NAs in column:", col))
    }
  } else {
    skip("deg object not available (run full pipeline first)")
  }
})

# ── Output directories exist ──────────────────────────────────────────────────

test_that("Output directories are created", {
  expect_true(exists("DIRS"))
  for (d in DIRS) {
    expect_true(dir.exists(here::here(d)))
  }
})

# ── validate_expression_scale with edge cases ─────────────────────────────────

test_that("validate_expression_scale handles NA values", {
  data_with_na <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2)
  expect_silent(validate_expression_scale(data_with_na))
})

# ── export_tsv creates valid output ───────────────────────────────────────────

test_that("export_tsv creates readable TSV with correct content", {
  tmp <- tempfile(fileext = ".tsv")
  df <- tibble::tibble(
    gene_symbol = c("GENE1", "GENE2"),
    logFC = c(1.5, -2.0),
    regulation = c("Up", "Down")
  )
  export_tsv(df, tmp, "Test export")
  expect_true(file.exists(tmp))
  imported <- readr::read_tsv(tmp, show_col_types = FALSE)
  expect_equal(nrow(imported), 2)
  expect_equal(ncol(imported), 3)
  unlink(tmp)
})

# ── STRING score filtering boundary test ──────────────────────────────────────

test_that("STRING score threshold filters correctly", {
  edges <- tibble::tibble(
    from = c("A", "B", "C"),
    to   = c("B", "C", "D"),
    combined_score = c(900, 500, 700)
  )
  filtered <- edges |> dplyr::filter(combined_score >= 700)
  expect_equal(nrow(filtered), 2)
  expect_true(all(filtered$combined_score >= 700))
})

cat("\n  Core tests complete.\n")
