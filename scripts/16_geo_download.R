# ═══════════════════════════════════════════════════════════════════════════════
# 16_geo_download.R — FASE 18/19: download dos datasets GEO de validação externa
#
# GSE33630 (primária): PTC vs normal thyroid (Affymetrix HG-U133 Plus 2.0)
# GSE60542 (secundária): PTC vs normal + metástase nodal
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(GEOquery); library(data.table)})

cat("══ FASE 18/19 — Download GEO ══\n")
cat("R:", as.character(getRversion()), "| GEOquery:", as.character(packageVersion("GEOquery")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_v <- "10_validation"
dir.create(dir_v, recursive = TRUE, showWarnings = FALSE)

options(download.file.method.GEOquery = "auto")

# ── GSE33630 ──────────────────────────────────────────────────────────────────
cat("── Baixando GSE33630 ──\n")
gse1 <- getGEO("GSE33630", GSEMatrix = TRUE, getGPL = TRUE)
if (is.list(gse1) && length(gse1) >= 1) { gse1 <- gse1[[1]] }
cat("  GSE33630: ", nrow(gse1), "features ×", ncol(gse1), "samples\n")
cat("  GPL:", annotation(gse1), "\n")
cat("  fData colunas:", paste(head(colnames(fData(gse1)), 30), collapse=", "), "\n")
cat("  pData colunas:", paste(head(colnames(pData(gse1)), 30), collapse=", "), "\n")
saveRDS(gse1, file.path(dir_v, "GSE33630_eset.rds"))

# ── GSE60542 ──────────────────────────────────────────────────────────────────
cat("\n── Baixando GSE60542 ──\n")
gse2 <- getGEO("GSE60542", GSEMatrix = TRUE, getGPL = TRUE)
if (is.list(gse2) && length(gse2) >= 1) { gse2 <- gse2[[1]] }
cat("  GSE60542: ", nrow(gse2), "features ×", ncol(gse2), "samples\n")
cat("  GPL:", annotation(gse2), "\n")
cat("  fData colunas:", paste(head(colnames(fData(gse2)), 30), collapse=", "), "\n")
cat("  pData colunas:", paste(head(colnames(pData(gse2)), 30), collapse=", "), "\n")
saveRDS(gse2, file.path(dir_v, "GSE60542_eset.rds"))

cat("\n══ Download GEO CONCLUÍDO ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
