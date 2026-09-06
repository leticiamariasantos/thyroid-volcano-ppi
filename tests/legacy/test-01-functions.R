# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat/test-01-functions.R — Unit tests for core functions
# thyroid-volcano-ppi
# ═══════════════════════════════════════════════════════════════════════════════

context("Core utility functions")

# ── validate_expression_scale ─────────────────────────────────────────────────

test_that("validate_expression_scale detects linear-scale data", {
  linear_data <- matrix(rnorm(100, mean = 500, sd = 100), nrow = 10)
  expect_warning(
    validate_expression_scale(linear_data),
    "appear integer"
  )
})

test_that("validate_expression_scale rejects extreme values", {
  extreme_data <- matrix(c(0, 5, 100), nrow = 1)
  expect_error(
    validate_expression_scale(extreme_data),
    "not in log2"
  )
})

test_that("validate_expression_scale accepts valid log2 data", {
  valid_log2 <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 10)
  expect_silent(validate_expression_scale(valid_log2))
})

# ── fetch_kegg_genes ──────────────────────────────────────────────────────────

test_that("fetch_kegg_genes returns tibble with correct columns", {
  skip_if_offline()
  result <- fetch_kegg_genes("hsa04919")
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("pathway_id", "gene_symbol"))
  expect_gt(nrow(result), 50)
})

test_that("fetch_kegg_genes filters fusion artifacts", {
  skip_if_offline()
  result <- fetch_kegg_genes("hsa04919")
  has_fusion <- any(grepl("-.+-", result$gene_symbol))
  expect_false(has_fusion)
})

# ── export_tsv ────────────────────────────────────────────────────────────────

test_that("export_tsv creates TSV file correctly", {
  tmp <- tempfile(fileext = ".tsv")
  df <- tibble::tibble(a = 1:3, b = letters[1:3])
  export_tsv(df, tmp)
  expect_true(file.exists(tmp))
  imported <- readr::read_tsv(tmp, show_col_types = FALSE)
  expect_equal(nrow(imported), 3)
  expect_equal(ncol(imported), 2)
  unlink(tmp)
})

# ── build_igraph ──────────────────────────────────────────────────────────────

test_that("build_igraph creates undirected graph correctly", {
  edges <- tibble::tibble(
    from = c("A", "B"),
    to   = c("B", "C"),
    combined_score = c(900, 800)
  )
  nodes <- tibble::tibble(
    STRING_id   = c("A", "B", "C"),
    gene_symbol = c("GENE1", "GENE2", "GENE3"),
    regulation  = c("Up", "Down", "NS")
  )
  g <- build_igraph(edges, nodes)
  expect_s3_class(g, "igraph")
  expect_false(igraph::is_directed(g))
  expect_equal(igraph::vcount(g), 3)
  expect_equal(igraph::ecount(g), 2)
})

# ── compute_centrality ────────────────────────────────────────────────────────

test_that("compute_centrality returns metrics for all nodes", {
  edges <- tibble::tibble(
    from = c("A", "B", "C"),
    to   = c("B", "C", "A"),
    combined_score = c(900, 800, 950)
  )
  nodes <- tibble::tibble(
    STRING_id   = c("A", "B", "C"),
    gene_symbol = c("G1", "G2", "G3"),
    regulation  = c("Up", "Down", "Up")
  )
  g <- build_igraph(edges, nodes)
  cent <- compute_centrality(g)
  expect_equal(nrow(cent), 3)
  expect_named(cent, c("gene_symbol", "STRING_id", "regulation",
                       "degree", "betweenness", "closeness", "hub_score"))
  expect_equal(cent$degree, c(2, 2, 2))
})
