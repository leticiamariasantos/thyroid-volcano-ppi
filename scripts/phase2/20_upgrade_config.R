# Configuration frozen in docs/METHODOLOGICAL_UPGRADE_2026.md.
suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
source(here("R", "upgrade_utils.R"))

UPGRADE_VERSION <- "5.0.0"
UPGRADE_DATE <- "2026-09-08"
N_BOOTSTRAP <- as.integer(Sys.getenv("THYROID_N_BOOTSTRAP", "100"))
GSEA_BOOTSTRAP_FDR_FRACTION <- 0.70
GSEA_BOOTSTRAP_DIRECTION_FRACTION <- 0.80
COEXPRESSION_RHO <- 0.30
META_I2_MAX <- 75
ITGA2_MIN_LOGFC <- 1.0
ITGA2_MIN_RETENTION <- 0.70
MAX_SURROGATE_VARIABLES <- 10L

DIR_UPGRADE <- file.path(DIR_RES, "upgrade_2026")
DIR_VARIANTS <- file.path(DIR_UPGRADE, "matrices")
DIR_DE_MULTI <- file.path(DIR_UPGRADE, "differential_expression")
DIR_GSEA_MULTI <- file.path(DIR_UPGRADE, "gsea")
DIR_DECONV <- file.path(DIR_UPGRADE, "deconvolution")
DIR_META <- file.path(DIR_UPGRADE, "validation")
DIR_PPI_MULTI <- file.path(DIR_UPGRADE, "ppi")
DIR_CACHE <- here("data", "cache")
DIR_PIPELINE <- file.path(DIR_UPGRADE, "pipeline")

for (d in c(DIR_UPGRADE, DIR_VARIANTS, DIR_DE_MULTI, DIR_GSEA_MULTI, DIR_DECONV,
            DIR_META, DIR_PPI_MULTI, DIR_CACHE, DIR_PIPELINE)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
set.seed(SEED)
