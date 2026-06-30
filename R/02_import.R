# ═══════════════════════════════════════════════════════════════════════════════
# R/02_import.R — Data import, validation, and preprocessing
# thyroid-volcano-ppi
#
# AUDIT: ✓ NA check  ✓ Duplicate check  ✓ Type coercion guard  ✓ Column validation
# LIMITATION: TCGA vs GTEx comparison without batch effect correction.
#             See AUDIT_REPORT.md §3 and docs/analysis_protocol.md §1.1
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── M1: Data Import & QC ──\n")

INPUT <- here::here("data", "raw", "XENA_THCA.tsv")
if (!file.exists(INPUT)) {
  stop(
    "\n❌ Data file not found: ", INPUT, "\n\n",
    "   ⚠ The file 'data/raw/XENA_THCA.tsv' is required to run the pipeline.\n\n",
    "   How to obtain:\n",
    "   1. Rscript scripts/download_data.R   (automatic download)\n",
    "   2. Or download manually from Xena Browser:\n",
    "      bit.ly/thyroid-volcano-ppi-data\n",
    "      Bookmark: c486b845ee2e750c3a9d2fc5145c8426\n\n",
    "   See README.md for full instructions.\n"
  )
}

raw <- readr::read_tsv(INPUT, show_col_types = FALSE, name_repair = "minimal")
cat(sprintf("  Loaded: %d × %d\n", nrow(raw), ncol(raw)))

# ── Column standardization ────────────────────────────────────────────────────
names(raw) <- sub("^_", "", names(raw))
if ("samples" %in% names(raw)) raw$samples <- NULL

# ── Validate metadata columns ─────────────────────────────────────────────────
META_COLS <- c("sample", "primary_site", "sample_type",
               "study", "TCGA_GTEX_main_category")
missing_meta <- setdiff(META_COLS, names(raw))
if (length(missing_meta) > 0)
  stop("Missing columns: ", paste(missing_meta, collapse = ", "))

# ── Check for duplicate sample IDs ────────────────────────────────────────────
if (anyDuplicated(raw$sample)) {
  warning("Duplicate sample IDs detected — keeping first occurrence")
  raw <- raw[!duplicated(raw$sample), ]
}

# ── Split metadata / expression ───────────────────────────────────────────────
meta <- raw[, META_COLS, drop = FALSE]
expr <- as.matrix(raw[, setdiff(names(raw), META_COLS), drop = FALSE])

if (!all(apply(expr, 2, is.numeric))) {
  stop("Non-numeric columns found in expression matrix.")
}
storage.mode(expr) <- "numeric"

# ── Biological condition ──────────────────────────────────────────────────────
meta$condition <- ifelse(
  meta$TCGA_GTEX_main_category == "GTEX Thyroid", "Normal",
  ifelse(meta$TCGA_GTEX_main_category == "TCGA Thyroid Carcinoma", "THCA", NA)
)
meta$condition <- factor(meta$condition, levels = c("Normal", "THCA"))

keep <- !is.na(meta$condition)
meta <- meta[keep, , drop = FALSE]
expr <- expr[keep, , drop = FALSE]

# ── NA proportion check ───────────────────────────────────────────────────────
na_frac <- mean(is.na(expr))
if (na_frac > 0.05)
  warning(sprintf("%.1f%% missing values in expression data", na_frac * 100))
cat(sprintf("  Missing values: %.2f%%\n", na_frac * 100))

# ── Expression scale QC ───────────────────────────────────────────────────────
validate_expression_scale(expr)

# ── Batch effect caveat ───────────────────────────────────────────────────────
cat("  NOTE: TCGA vs GTEx comparison without batch correction.\n")
cat("  Differential expression may partially reflect technical variation.\n")
cat("  See AUDIT_REPORT.md and docs/analysis_protocol.md\n")

cat("\n  Samples:\n")
print(table(meta$condition))

# Sample composition table
qc_samples <- meta |>
  dplyr::group_by(condition, study, sample_type) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(condition, study)

export_tsv(qc_samples, here::here(DIRS$tables, "T01_sample_composition.tsv"))
cat(sprintf("  Ready: %d samples × %d genes\n", nrow(expr), ncol(expr)))
