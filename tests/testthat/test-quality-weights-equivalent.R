source(here::here("R", "quality_weights_equivalent.R"))

test_that("sample-specific updates agree with limma for paired and composition designs", {
  skip_if_not_installed("limma")
  set.seed(42)
  n <- 118L
  group <- factor(rep(c("Normal", "Tumor"), n / 2L))
  paired <- model.matrix(~factor(rep(seq_len(n / 2L), each = 2L)) + group)
  composition <- model.matrix(~group + matrix(rnorm(n * 8L), n, 8L))
  E <- matrix(rnorm(24L * n), 24L, n)
  W <- matrix(runif(length(E), 0.1, 3), nrow(E), n)
  for (design in list(paired, composition)) {
    original <- limma::arrayWeights(E, design, weights = W, method = "genebygene")
    actual <- array_weights_samplewise_fast(E, design, W)
    expect_lt(max(abs(log(as.numeric(actual) / original))), 1e-7)
  }
})

test_that("both voom passes and the fitted expression model remain equivalent", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  set.seed(42)
  counts <- matrix(rnbinom(300L * 24L, mu = 40, size = 3), 300L, 24L)
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  dge <- edgeR::normLibSizes(edgeR::DGEList(counts))
  design <- model.matrix(~factor(rep(c("Normal", "Tumor"), each = 12L)))
  original <- limma::voomWithQualityWeights(dge, design, plot = FALSE)
  actual <- voom_quality_weights_equivalent(dge, design)
  expect_equal(actual$E, original$E, tolerance = 1e-7)
  expect_equal(actual$weights, original$weights, tolerance = 1e-7)
  expect_equal(actual$targets$sample.weights, original$targets$sample.weights, tolerance = 1e-7)
  ref_fit <- limma::eBayes(limma::lmFit(original, design), robust = TRUE)
  fast_fit <- limma::eBayes(limma::lmFit(actual, design), robust = TRUE)
  expect_equal(fast_fit$coefficients, ref_fit$coefficients, tolerance = 1e-6)
  expect_equal(fast_fit$t, ref_fit$t, tolerance = 1e-6)
})

test_that("zero-residual genes follow the reference skip rule", {
  set.seed(42)
  design <- matrix(1, 16L, 1L)
  E <- matrix(rnorm(20L * 16L), 20L, 16L)
  E[1, ] <- 4
  W <- matrix(1, nrow(E), ncol(E))
  original <- limma::arrayWeights(E, design, weights = W, method = "genebygene")
  actual <- array_weights_samplewise_fast(E, design, W)
  expect_equal(as.numeric(actual), as.numeric(original), tolerance = 1e-7)
})
