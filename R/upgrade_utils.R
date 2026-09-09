# Shared helpers for the 2026 methodological upgrade.

assert_full_rank <- function(design, label = "design") {
  rank <- qr(design)$rank
  if (rank < ncol(design)) {
    stop(label, " is not identifiable (rank ", rank, "/", ncol(design),
         "). Do not label this analysis as batch-corrected.")
  }
  invisible(TRUE)
}

sample_code <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- grepl("^TCGA-[^-]+-[^-]+-[0-9]{2}", x)
  out[hit] <- substr(x[hit], 14L, 15L)
  out
}

tcga_key <- function(x) {
  ifelse(grepl("^TCGA-", x), paste0(substr(x, 1L, 12L), "-", substr(x, 14L, 15L)), x)
}

collapse_counts <- function(counts, symbols) {
  symbols <- as.character(symbols)
  symbols[is.na(symbols) | !nzchar(symbols)] <- rownames(counts)[is.na(symbols) | !nzchar(symbols)]
  rownames(counts) <- symbols
  counts <- rowsum(as.matrix(counts), group = rownames(counts), reorder = FALSE)
  counts[nzchar(rownames(counts)), , drop = FALSE]
}

write_manifest <- function(path, fields) {
  stopifnot(!is.null(names(fields)))
  data.table::fwrite(data.table::data.table(field = names(fields), value = unlist(fields)),
                     path, sep = "\t", quote = FALSE)
  invisible(path)
}

matrix_checksum <- function(x) {
  digest::digest(list(dim = dim(x), rownames = rownames(x), colnames = colnames(x),
                      total = sum(x), sample_sums = colSums(x)), algo = "sha256")
}

keep_by_cpm <- function(counts, condition, min_cpm = 1, fraction = 0.25) {
  condition <- droplevels(factor(condition))
  if (length(condition) != ncol(counts) || nlevels(condition) < 2L)
    stop("CPM filter requires one condition per sample and at least two groups")
  min_group_n <- min(table(condition))
  required_n <- max(1L, ceiling(fraction * min_group_n))
  dge <- edgeR::DGEList(counts = counts)
  rowSums(edgeR::cpm(dge) >= min_cpm) >= required_n
}

safe_api_cache <- function(cache_file, fetch, validate = function(x) TRUE,
                           attempts = 3L, retry_seconds = c(1, 3)) {
  invalid_cache <- FALSE
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached) && isTRUE(validate(cached))) return(cached)
    invalid_cache <- TRUE
  }
  attempts <- max(1L, as.integer(attempts))
  value <- NULL
  last_error <- NULL
  for (i in seq_len(attempts)) {
    value <- tryCatch(fetch(), error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })
    if (!is.null(value) && isTRUE(validate(value))) break
    if (i < attempts) Sys.sleep(retry_seconds[min(i, length(retry_seconds))])
  }
  if (is.null(value) || !isTRUE(validate(value))) {
    stop("API unavailable and no valid cache at ", cache_file,
         if (!is.null(last_error)) paste0(": ", last_error) else "")
  }
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(cache_file, ".tmp-", Sys.getpid())
  saveRDS(value, tmp)
  if (invalid_cache && file.exists(cache_file)) {
    backup <- paste0(cache_file, ".invalid-", format(Sys.time(), "%Y%m%d%H%M%S"))
    if (!file.rename(cache_file, backup)) {
      unlink(tmp)
      stop("Could not preserve invalid API cache: ", cache_file)
    }
  }
  if (!file.rename(tmp, cache_file)) {
    unlink(tmp)
    stop("Could not atomically write API cache: ", cache_file)
  }
  value
}

residualize_preserving_condition <- function(expr, metadata, composition) {
  stopifnot(identical(colnames(expr), metadata$sample))
  composition <- as.data.frame(composition[match(metadata$sample, rownames(composition)), , drop = FALSE])
  keep <- vapply(composition, function(x) stats::sd(x, na.rm = TRUE) > 1e-8, logical(1))
  composition <- composition[, keep, drop = FALSE]
  composition[] <- lapply(composition, function(x) {
    x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
    as.numeric(scale(x))
  })
  base <- model.matrix(~ condition, data = metadata)
  full <- cbind(base, as.matrix(composition))
  assert_full_rank(full, "composition-adjusted design")
  fit_full <- limma::lmFit(expr, full)
  fit_base <- limma::lmFit(expr, base)
  fitted_base <- fit_base$coefficients %*% t(base)
  fitted_full <- fit_full$coefficients %*% t(full)
  fitted_base + (expr - fitted_full)
}
