# Equivalent sample-specific gene-by-gene updates for limma 3.68.0.
# See docs/STAGE24_ACCELERATION_2026.md for the algebra and acceptance gates.
array_weights_samplewise_fast <- function(E, design, weights, prior.n = 10,
                                         solve_tol = 1e-12, progress = NULL) {
  n <- ncol(E)
  k <- n - 1L
  info <- prior.n * (diag(k) + 1)
  diagonal <- rep(prior.n, k)
  common <- prior.n
  gam <- numeric(k)
  aw <- rep(1, n)
  iterations <- integer(nrow(E))
  for (i in seq_len(nrow(E))) {
    w <- aw * weights[i, ]
    fit <- lm.wfit(design, E[i, ], w)
    h <- 1 - hat(fit$qr)
    s2 <- mean(fit$effects[-seq_len(fit$rank)]^2)
    if (s2 < 1e-15) next
    b <- h[seq_len(k)] - h[n]
    diagonal <- diagonal + h[seq_len(k)]
    common <- common + h[n]
    info <- info + h[n] - tcrossprod(b) / sum(h)
    diag(info) <- diag(info) + h[seq_len(k)]
    z <- w * fit$residuals^2 / s2 - h
    rhs <- z[seq_len(k)] - z[n]
    # Woodbury preconditioner for diag(diagonal) + common * 11'.
    precondition <- function(r) {
      scaled <- r / diagonal
      scaled - (sum(scaled) / (1 / common + sum(1 / diagonal))) / diagonal
    }
    step <- numeric(k)
    residual <- rhs
    target <- solve_tol * max(1, sqrt(sum(rhs^2)))
    zz <- precondition(residual)
    direction <- zz
    rz <- sum(residual * zz)
    for (it in seq_len(50L)) {
      if (sqrt(sum(residual^2)) <= target) break
      product <- drop(info %*% direction)
      alpha <- rz / sum(direction * product)
      step <- step + alpha * direction
      residual <- residual - alpha * product
      zz <- precondition(residual)
      rz_new <- sum(residual * zz)
      direction <- zz + (rz_new / rz) * direction
      rz <- rz_new
    }
    # A true-residual check avoids accepting recursively accumulated error.
    if (any(!is.finite(step)) || sqrt(sum((rhs - drop(info %*% step))^2)) > target * 5) {
      step <- solve(info, rhs)
    }
    gam <- gam + step
    aw <- exp(c(-gam, sum(gam)))
    iterations[i] <- it
    if (is.function(progress) && (i %% 1000L == 0L || i == nrow(E))) progress(i, nrow(E))
  }
  attr(aw, "iterations") <- iterations
  aw
}


# Change only the linear algebra used by the two arrayWeights calls; the voom
# transformations, trends, empirical Bayes and DE models remain those of limma.
voom_quality_weights_equivalent <- function(counts, design, progress = NULL) {
  reference <- limma::voomWithQualityWeights
  if (as.character(utils::packageVersion("limma")) != "3.68.0")
    return(reference(counts, design, plot = FALSE))
  scope <- new.env(parent = environment(reference))
  scope$arrayWeights <- function(object, design, method = "genebygene", maxiter = 50L,
                                 tol = 1e-5, var.design = NULL, var.group = NULL,
                                 trace = FALSE, ...) {
    E <- object$E
    weights <- object$weights
    supported <- identical(method, "genebygene") && is.null(var.design) &&
      is.null(var.group) && is.matrix(E) && is.matrix(weights) &&
      identical(dim(E), dim(weights)) && all(is.finite(E)) &&
      all(is.finite(weights)) && all(weights > 0) && nrow(E) >= 2L &&
      ncol(E) - qr(design)$rank >= 2L && qr(design)$rank == ncol(design)
    if (!supported) {
      return(limma::arrayWeights(object, design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design, var.group = var.group, trace = trace, ...))
    }
    setNames(as.numeric(array_weights_samplewise_fast(E, design, weights,
      progress = progress)), colnames(E))
  }
  environment(reference) <- scope
  reference(counts, design, plot = FALSE)
}
