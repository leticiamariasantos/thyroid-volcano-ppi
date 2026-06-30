# ═══════════════════════════════════════════════════════════════════════════════
# R/03d_qc_outliers.R — QC: outlier detection + comprehensive sample composition
# thyroid-volcano-ppi v3.1.0
#
# Outputs: results/tables/T_QC_outlier_candidates.tsv
#          results/tables/T_QC_sample_composition.tsv
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── QC: Outlier detection ──\n")

# ── Sample-level QC metrics ───────────────────────────────────────────────────
# Compute per-sample metrics from expression matrix
sample_metrics <- data.frame(
  sample         = colnames(E),
  mean_expr      = colMeans(E, na.rm = TRUE),
  sd_expr        = apply(E, 2, sd, na.rm = TRUE),
  median_expr    = apply(E, 2, median, na.rm = TRUE),
  iqr_expr       = apply(E, 2, IQR, na.rm = TRUE),
  pct_zero       = colMeans(E == 0, na.rm = TRUE) * 100,
  stringsAsFactors = FALSE
)

# ── Flag outliers: ±3 SD from mean of means ───────────────────────────────────
mu_mean  <- mean(sample_metrics$mean_expr, na.rm = TRUE)
sd_mean  <- sd(sample_metrics$mean_expr, na.rm = TRUE)
mu_sd    <- mean(sample_metrics$sd_expr, na.rm = TRUE)
sd_sd    <- sd(sample_metrics$sd_expr, na.rm = TRUE)

sample_metrics$outlier_mean <- sample_metrics$mean_expr < (mu_mean - 3*sd_mean) |
                                sample_metrics$mean_expr > (mu_mean + 3*sd_mean)
sample_metrics$outlier_sd   <- sample_metrics$sd_expr < (mu_sd - 3*sd_sd) |
                                sample_metrics$sd_expr > (mu_sd + 3*sd_sd)
sample_metrics$is_outlier   <- sample_metrics$outlier_mean | sample_metrics$outlier_sd

# ── Merge with metadata ───────────────────────────────────────────────────────
outliers <- sample_metrics |>
  dplyr::left_join(
    meta |> dplyr::select(sample, condition, study, sample_type),
    by = "sample"
  ) |>
  dplyr::arrange(dplyr::desc(is_outlier), dplyr::desc(abs(mean_expr - mu_mean)))

# ── Export outlier candidates ─────────────────────────────────────────────────
outlier_candidates <- outliers |>
  dplyr::filter(is_outlier) |>
  dplyr::select(sample, condition, study, sample_type, mean_expr, sd_expr,
                median_expr, pct_zero, outlier_mean, outlier_sd)

export_tsv(outlier_candidates,
           here::here(DIRS$tables, "T_QC_outlier_candidates.tsv"),
           "Samples flagged as outliers (±3 SD)")

n_out <- sum(sample_metrics$is_outlier)
cat(sprintf("  Outlier candidates: %d / %d (%.1f%%)\n",
            n_out, nrow(sample_metrics), n_out/nrow(sample_metrics)*100))
if (n_out > 0) {
  cat("  Outlier samples saved to T_QC_outlier_candidates.tsv\n")
}

# ── Comprehensive sample composition ──────────────────────────────────────────
qc_comp <- meta |>
  dplyr::group_by(condition, study, sample_type) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_expr = mean(sample_metrics$mean_expr[sample_metrics$sample %in% sample]),
    sd_expr   = mean(sample_metrics$sd_expr[sample_metrics$sample %in% sample]),
    pct_zero  = mean(sample_metrics$pct_zero[sample_metrics$sample %in% sample]),
    .groups = "drop"
  ) |>
  dplyr::arrange(condition, study)

export_tsv(qc_comp,
           here::here(DIRS$tables, "T_QC_sample_composition.tsv"),
           "Comprehensive sample composition with QC metrics")

cat(sprintf("  QC composition: %d groups\n", nrow(qc_comp)))
