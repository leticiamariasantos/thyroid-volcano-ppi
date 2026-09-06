# ═══════════════════════════════════════════════════════════════════════════════
# 00_config.R — Configuração central da NOVA FASE ANALÍTICA (Phase 2 reboot)
#
# Data da nova execução: 2026-09-06
# Este script define caminhos, thresholds, seeds e helpers usados por todo o
# pipeline. Nenhum resultado da fase anterior é reutilizado.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

# ── Metadados da execução ──────────────────────────────────────────────────────
EXEC_DATE      <- "2026-09-06"
EXEC_PHASE     <- "phase2_reboot"
EXEC_DESC      <- "Reconstrução independente do pipeline com expansão do painel de vias de 10 para 30."
SEED           <- 42
set.seed(SEED)

# ── Thresholds metodológicos (pré-especificados) ───────────────────────────────
LOGFC_THRESH   <- 1.0      # |log2FC| mínimo para classificar DEG
FDR_THRESH     <- 0.05     # FDR máximo (Benjamini-Hochberg)
STRING_SCORE   <- 700      # escore mínimo STRING (alta confiança)
MIN_EXPR_CPM   <- 1        # filtro de baixa expressão (voom/DESeq2): >=1 CPM em >=25% das amostras do menor grupo
MIN_EXPR_TPM   <- 0.25     # filtro de baixa expressão (limma TPM): >=0.25 log2(TPM+0.001) i.e. TPM>=~0.19
EXPR_FRAC      <- 0.25     # fração mínima de amostras expressando acima do limiar

# ── Caminhos (via here::here() → raiz do repositório) ─────────────────────────
DIR_DATA      <- here("data", "global")
DIR_EXTERNAL  <- here("data", "external")
DIR_RES       <- here("results", "phase2")

DIR_AUDIT     <- file.path(DIR_RES, "audit")
DIR_DATAOUT   <- file.path(DIR_RES, "data")
DIR_PREPROC   <- file.path(DIR_RES, "preprocessing")
DIR_DE        <- file.path(DIR_RES, "differential_expression")
DIR_PATH      <- file.path(DIR_RES, "pathways")
DIR_GSEA      <- file.path(DIR_RES, "gsea")
DIR_PPI       <- file.path(DIR_RES, "ppi")
DIR_VAL       <- file.path(DIR_RES, "validation")
DIR_SENS      <- file.path(DIR_RES, "sensitivity")
DIR_FIG       <- file.path(DIR_RES, "figures")
DIR_REP       <- file.path(DIR_RES, "reports")

dirs <- c(DIR_AUDIT, DIR_DATAOUT, DIR_PREPROC, DIR_DE, DIR_PATH, DIR_GSEA,
          DIR_PPI, DIR_VAL, DIR_SENS, DIR_FIG, DIR_REP)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Arquivos de dados
TPM_FILE       <- file.path(DIR_DATA, "TCGA_GTEx_thyroid_tpm.tsv")
COUNTS_FILE    <- file.path(DIR_DATA, "TCGA_GTEx_thyroid_counts.tsv")
RAW_TPM_GZ     <- file.path(DIR_DATA, "TcgaTargetGtex_rsem_gene_tpm.gz")

# ── Helpers de I/O ─────────────────────────────────────────────────────────────
fwrite_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

# Validador cumulativo (para scripts de checagem)
new_checker <- function() {
  n <- 0
  list(
    ck = function(cond, msg) {
      if (isTRUE(cond)) cat("  [OK]  ", msg, "\n")
      else { cat("  [FAIL]", msg, "\n"); n <<- n + 1 }
    },
    total = function() n
  )
}

log_msg("Configuração carregada. Raiz do repositório:", here())
