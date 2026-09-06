# ═══════════════════════════════════════════════════════════════════════════════
# run_pipeline.R — Master runner | thyroid-volcano-ppi
#
# THCA Transcriptomic Analysis: differential expression (limma),
# Volcano Plot, and STRING PPI network of DEGs.
#
# Outputs (PNG 600 dpi + PDF, Nature Communications / Cell Press standard):
#   1. Fig1_Volcano_THCA_vs_Normal.{png,pdf}
#   2. Fig2_PPI_Network_THCA_DEGs.{png,pdf}
#
# Usage:
#   Rscript run_pipeline.R
#   source("run_pipeline.R")
#
# Repository: https://github.com/santosry/thyroid-volcano-ppi
# ═══════════════════════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════════════════════
# SETUP: project root, internet check, directories, packages, parameters
# ═══════════════════════════════════════════════════════════════════════════════
source(here::here("R", "00_setup.R"), local = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# DATA FILE VERIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

DATA_FILE <- file.path(PROJECT_ROOT, "data", "raw", "XENA_THCA.tsv")
DOWNLOAD_SCRIPT <- file.path(PROJECT_ROOT, "scripts", "download_data.R")

if (!file.exists(DATA_FILE)) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  DATA FILE NOT FOUND                                       ║\n")
  cat("╠══════════════════════════════════════════════════════════════╣\n")
  cat("║                                                            ║\n")
  cat("║  The file 'data/raw/XENA_THCA.tsv' is required.            ║\n")
  cat("║  It contains the gene expression matrix                    ║\n")
  cat("║  (TCGA THCA + GTEx Normal) used in the analysis.           ║\n")
  cat("║                                                            ║\n")
  cat("║  How to obtain the file:                                   ║\n")
  cat("║  1. Go to Xena Browser:                                   ║\n")
  cat("║     bit.ly/thyroid-volcano-ppi-data                        ║\n")
  cat("║  2. Click 'Download' in the top-right corner              ║\n")
  cat("║  3. Select 'Download current visualization data'          ║\n")
  cat("║  4. Save as: data/raw/XENA_THCA.tsv                       ║\n")
  cat("║                                                            ║\n")
  cat("║  Or run the automatic download script:                     ║\n")
  cat("║  Rscript scripts/download_data.R                           ║\n")
  cat("║                                                            ║\n")
  cat("║  Xena Bookmark:                                            ║\n")
  cat("║  c486b845ee2e750c3a9d2fc5145c8426                          ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  if (file.exists(DOWNLOAD_SCRIPT)) {
    cat("Attempting automatic download via scripts/download_data.R...\n\n")
    tryCatch({
      source(DOWNLOAD_SCRIPT, local = FALSE)
    }, error = function(e) {
      cat("\nAutomatic download failed:", conditionMessage(e), "\n")
    })
  }

  if (!file.exists(DATA_FILE)) {
    stop(
      "\nPIPELINE ABORTED: Data file missing.\n",
      "   Expected path: ", DATA_FILE, "\n",
      "   Please download the file manually from Xena Browser and try again.\n",
      "   Full instructions in README.md\n"
    )
  } else {
    cat("Automatic download successful! Proceeding...\n\n")
  }
}

f_info <- file.info(DATA_FILE)
if (f_info$size < 1000) {
  stop("Data file too small (", f_info$size,
       " bytes). File may be corrupted.\n",
       "   Please re-download from Xena Browser.")
}

# ═══════════════════════════════════════════════════════════════════════════════
# INTERNET CONNECTIVITY CHECK (STRING, KEGG require internet)
# ═══════════════════════════════════════════════════════════════════════════════
if (!check_internet()) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║  WARNING: No internet connection detected.                  ║\n")
  cat("║  STRING and KEGG API calls will fail.                       ║\n")
  cat("║  Differential expression and Volcano Plot will still work.  ║\n")
  cat("║  PPI Network requires internet.                             ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")
}

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  thyroid-volcano-ppi : THCA DEG + PPI Analysis      ║\n")
cat("║  Output: Volcano Plot + PPI Network                  ║\n")
cat("║  Standard: Nature Communications / Cell Press        ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PIPELINE ORDER (v3.1.0 — AUDITED):
#   1. Import + validate data
#   2. Filter low expression + limma DEG (gene-by-gene)
#   3. PCA/UMAP QC on filtered expression matrix (diagnostic, NOT DEG replacement)
#   4. Outlier detection
#   5. PPI network (STRING, centrality)
#   6. Volcano plot (with hub labels from PPI)
#   7. Supplementary tables
#
# IMPORTANT: PCA/UMAP are QC only. DEG is done gene-by-gene with limma.
# ═══════════════════════════════════════════════════════════════════════════════

source(here::here("R", "01_functions.R"), local = FALSE)
source(here::here("R", "02_import.R"),    local = FALSE)   # Step 1: Import
source(here::here("R", "03_deg.R"),       local = FALSE)   # Step 2: Filter + DEG (creates E)

# ── QC: PCA + UMAP + Outliers (on filtered E from 03_deg.R) ──────────────────
for (qc_script in c("R/03b_pca.R", "R/03e_umap_qc.R", "R/03d_qc_outliers.R")) {
  if (file.exists(here::here(qc_script))) {
    tryCatch(
      source(here::here(qc_script), local = FALSE),
      error = function(e) cat("  QC skipped (", basename(qc_script), "):",
                               conditionMessage(e), "\n")
    )
  }
}

source(here::here("R", "05_ppi.R"),       local = FALSE)   # Step 5: PPI FIRST
source(here::here("R", "04_volcano.R"),   local = FALSE)   # Step 6: Volcano AFTER PPI

# ── Session info ─────────────────────────────────────────────────────────────
sink(file.path(DIRS$logs, "session_info.txt"))
sessionInfo()
sink()

# ── renv snapshot reminder ────────────────────────────────────────────────────
if (requireNamespace("renv", quietly = TRUE)) {
  cat("\nTip: run renv::snapshot() to update renv.lock with current packages.\n")
}

# ── Final summary ────────────────────────────────────────────────────────────
cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  PIPELINE COMPLETED SUCCESSFULLY                     ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat("Outputs:\n")
cat("  Fig1: Volcano Plot     -> results/figures/Fig1_Volcano_THCA_vs_Normal.{png,pdf}\n")
cat("  Fig2: PPI Network      -> results/figures/Fig2_PPI_Network_THCA_DEGs.{png,pdf}\n")
cat("  Tables (7)             -> results/tables/\n")
cat("  Network metadata (4)   -> results/network/\n")
cat("  Session info           -> logs/session_info.txt\n")
cat("\n")
cat("DISCLAIMER: This pipeline generates hypotheses based on\n")
cat("transcriptomic associations. It does not establish causality,\n")
cat("therapeutic targets, or clinical recommendations. PPI hub\n")
cat("proteins are identified by centrality metrics (exploratory);\n")
cat("they do not imply biological validation or druggability.\n")
cat("\n")
cat("Study type: Exploratory, hypothesis-generating.\n")
cat("TCGA vs GTEx comparison without batch effect correction —\n")
cat("differential expression may partially reflect technical variation.\n")
cat("See README.md#limitações for full scientific limitations.\n")
cat("\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
