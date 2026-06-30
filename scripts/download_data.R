#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# scripts/download_data.R — Download XENA_THCA.tsv from UCSC Xena Browser
#
# Downloads the TCGA THCA + GTEx Thyroid expression matrix
# used in the thyroid-volcano-ppi analysis.
#
# Source: UCSC Xena Browser
# Bookmark: c486b845ee2e750c3a9d2fc5145c8426
# Data: log2(norm_count + 1) RSEM expected_count from TOIL recompute
#
# Usage:
#   Rscript scripts/download_data.R
# ═══════════════════════════════════════════════════════════════════════════════

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
library(here)
PROJECT_ROOT <- here::here()

XENA_URL <- paste0(
  "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/",
  "TcgaTargetGtex_rsem_gene_tpm.gz"
)

output_path <- file.path(PROJECT_ROOT, "data", "raw", "XENA_THCA.tsv")

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║  Data Download — thyroid-volcano-ppi             ║\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

if (file.exists(output_path)) {
  cat("✓ Data file already exists:", output_path, "\n")
  cat("  Size:", file.info(output_path)$size, "bytes\n")
  cat("  No action needed.\n")
  quit(save = "no", status = 0)
}

cat("Attempting download from Xena Browser...\n")
cat("URL:", XENA_URL, "\n\n")

cat("NOTE: If download fails, please download manually from:\n")
cat("  https://xenabrowser.net/?bookmark=c486b845ee2e750c3a9d2fc5145c8426\n")
cat("  → Click 'Download' → 'Download current visualization data'\n")
cat("  → Save as: data/raw/XENA_THCA.tsv\n\n")

tryCatch({
  options(timeout = 600)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  utils::download.file(url = XENA_URL, destfile = output_path, mode = "wb")

  cat("✓ Download complete:", output_path, "\n")
  cat("  Size:", file.info(output_path)$size, "bytes\n")

  if (file.info(output_path)$size < 1000) {
    stop("Downloaded file is empty or corrupted.")
  }

}, error = function(e) {
  cat("\n✗ Download failed:", conditionMessage(e), "\n\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  MANUAL DOWNLOAD INSTRUCTIONS:                              ║\n")
  cat("║  1. Go to: bit.ly/thyroid-volcano-ppi-data                  ║\n")
  cat("║  2. Click 'Download' button (top-right corner)              ║\n")
  cat("║  3. Select 'Download current visualization data'            ║\n")
  cat("║  4. Save as: data/raw/XENA_THCA.tsv                         ║\n")
  cat("║                                                             ║\n")
  cat("║  Bookmark: c486b845ee2e750c3a9d2fc5145c8426                  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
})
