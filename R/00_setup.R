# ═══════════════════════════════════════════════════════════════════════════════
# R/00_setup.R — Project setup, parameters, and packages
# thyroid-volcano-ppi — THCA DEG + PPI analysis
#
# Called by run_pipeline.R. Do NOT source standalone without setting PROJECT_ROOT.
# ═══════════════════════════════════════════════════════════════════════════════

# ── Project root detection (unified; no duplication with run_pipeline.R) ──────
if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
library(here)
PROJECT_ROOT <- here::here()

# ── Directories (created if missing) ──────────────────────────────────────────
DIRS <- list(
  tables  = "results/tables",
  figures = "results/figures",
  network = "results/network",
  logs    = "logs"
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Analysis Parameters ───────────────────────────────────────────────────────
THRESHOLD <- list(
  lfc        = 1.0,     # |log2FC| > 1 -> 2-fold biological significance
  fdr        = 0.05,    # Benjamini-Hochberg, 5% false discovery rate
  string     = 700,     # STRING combined_score >= 700 -> high confidence
  expr_min   = 0.5,     # log2(norm_count+1) detection floor
  expr_frac  = 0.1      # expressed in >= 10% of samples
)

KEGG_ID       <- "hsa04919"   # Thyroid hormone signaling pathway
STRING_TAXON  <- 9606         # Homo sapiens
STRING_V      <- "12.0"
PIPELINE_VERSION <- "3.1.0"

# ── Internet connectivity helper ──────────────────────────────────────────────
check_internet <- function(timeout = 5) {
  urls <- c(
    "https://string-db.org",
    "https://rest.kegg.jp"
  )
  for (u in urls) {
    res <- tryCatch(
      httr::HEAD(u, httr::timeout(timeout)),
      error = function(e) NULL
    )
    if (!is.null(res) && res$status_code < 400) return(TRUE)
  }
  # Fallback: generic connectivity
  res <- tryCatch(
    httr::HEAD("https://httpbin.org/get", httr::timeout(timeout)),
    error = function(e) NULL
  )
  return(!is.null(res) && res$status_code < 400)
}

# ── Reproducibility ───────────────────────────────────────────────────────────
set.seed(42)
options(scipen = 999, digits = 4, timeout = 600, stringsAsFactors = FALSE)

# ── Packages ──────────────────────────────────────────────────────────────────
REQUIRED_PKGS <- c(
  "here", "dplyr", "stringr", "tibble", "tidyr", "purrr", "readr",
  "limma", "ggplot2", "ggrepel",
  "igraph", "ggraph",
  "KEGGREST", "org.Hs.eg.db", "AnnotationDbi",
  "httr", "jsonlite"
)

# ── Package configuration ─────────────────────────────────────────────────────
BIOC_PKGS <- c("limma", "KEGGREST", "org.Hs.eg.db", "AnnotationDbi")

for (pkg in REQUIRED_PKGS) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% BIOC_PKGS) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(here); library(dplyr); library(stringr); library(tibble)
  library(tidyr); library(purrr); library(readr)
  library(limma); library(ggplot2); library(ggrepel)
  library(igraph); library(ggraph)
  library(KEGGREST); library(org.Hs.eg.db); library(AnnotationDbi)
  library(httr); library(jsonlite)
})

# ── Publication Color Palette (Nature Communications / Cell Press) ────────────
CELL_COLORS <- list(
  up        = "#4477AA",
  down      = "#AA4488",
  ns        = "#CCCCCC",
  edge      = "#D8D8D8",
  highlight = "#1A1A1A",
  line      = "#4D4D4D",
  bg        = "#FFFFFF",
  grid      = "#E8E8E8",
  outline   = "#333333"
)

# ── Publication Typography ────────────────────────────────────────────────────
FONT_SIZE <- 14
FONT_FAM  <- "sans"

# ── Logging ───────────────────────────────────────────────────────────────────
log_file <- file.path(DIRS$logs, paste0("pipeline_",
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
msg_file <- file.path(DIRS$logs, paste0("msgs_",
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
msg_con  <- file(msg_file, open = "wt")
sink(msg_con, type = "message")

cat("\n═══════════════════════════════════════════\n")
cat("thyroid-volcano-ppi v", PIPELINE_VERSION, "\n")
cat("KEGG:", KEGG_ID, "| STRING v", STRING_V, "\n")
cat("Root:", PROJECT_ROOT, "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("═══════════════════════════════════════════\n")
cat("NOTE: Exploratory study. No batch correction applied.\n")
cat("      TCGA vs GTEx comparison requires cautious interpretation.\n")
cat("      See README.md#limitações for full limitations.\n\n")

# ── renv lockfile reminder ────────────────────────────────────────────────────
if (requireNamespace("renv", quietly = TRUE)) {
  cat("renv detected. Run renv::snapshot() to update renv.lock after package changes.\n")
} else {
  cat("Tip: install.packages('renv') for reproducible package versions.\n")
}
