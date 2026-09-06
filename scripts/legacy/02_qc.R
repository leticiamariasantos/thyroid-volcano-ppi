#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 02_qc.R — FASE 2: Controle de Qualidade da matriz GLOBAL (log2(TPM+0.001))
# thyroid-volcano-ppi
#
# Reimplementa o QC de outliers da Fase 1 (que FALHOU por bug de dimensões),
# agora em Python-independente, com log explícito e sem exclusão automática.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(stats)
  library(grDevices)
  library(here)
})

PROJECT_ROOT <- here::here()
dir_qc <- file.path(PROJECT_ROOT, "02_qc")
dir.create(dir_qc, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — QC matriz global ══\n")
cat("R:", as.character(getRversion()), "|", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Carregar matriz (genes × amostras) ─────────────────────────────────────
mat_file <- file.path(PROJECT_ROOT, "data", "global", "TCGA_GTEx_thyroid_tpm.tsv")
X <- as.matrix(data.table::fread(mat_file, showProgress = FALSE), rownames = TRUE)
storage.mode(X) <- "numeric"
cat(sprintf("  Matriz: %d genes × %d amostras\n", nrow(X), ncol(X)))

# ── 2. Metadados (condição derivada do prefixo; cross-check 504/279) ──────────
samples <- colnames(X)
condition <- ifelse(grepl("^GTEX", samples), "Normal",
                    ifelse(grepl("^TCGA", samples), "THCA", NA_character_))
cohort <- ifelse(grepl("^GTEX", samples), "GTEx", "TCGA")
cat(sprintf("  Composição: THCA=%d | Normal=%d | NA=%d\n",
            sum(condition == "THCA", na.rm = TRUE),
            sum(condition == "Normal", na.rm = TRUE), sum(is.na(condition))))
stopifnot(sum(condition == "THCA", na.rm = TRUE) == 504,
          sum(condition == "Normal", na.rm = TRUE) == 279)

meta <- data.frame(sample = samples, condition = condition, cohort = cohort,
                   stringsAsFactors = FALSE)

# ── 3. Métricas por amostra (valores em log2(TPM+0.001)) ─────────────────────
# "expresso" = TPM > 1  <=>  log2(TPM+0.001) > log2(1.001) ≈ 0.0014 → usar > 0
expressed_frac <- colMeans(X > 0)                      # fração de genes com TPM>1
total_tpm      <- colSums(2^X - 0.001)                 # soma de TPM por amostra
median_expr    <- apply(X, 2, median)
mean_expr      <- colMeans(X)

qc <- data.frame(sample = samples, condition = condition, cohort = cohort,
                 expressed_frac = expressed_frac, total_tpm = total_tpm,
                 median_expr = median_expr, mean_expr = mean_expr,
                 stringsAsFactors = FALSE)

cat(sprintf("  Fração expressa: med=%.2f | min=%.3f | max=%.3f\n",
            median(expressed_frac), min(expressed_frac), max(expressed_frac)))
cat(sprintf("  Soma TPM: med=%.1f | min=%.1f | max=%.1f (×1e3)\n",
            median(total_tpm)/1e3, min(total_tpm)/1e3, max(total_tpm)/1e3))

# ── 4. Outliers técnicos (critério PRÉ-definido, dentro da coorte) ────────────
# flag: fração expressa < mediana - 3*IQR  OU  soma TPM < mediana - 3*IQR
flag <- function(v) { m <- median(v); i <- IQR(v); v < (m - 3 * i) }
qc$flag_expr  <- ave(qc$expressed_frac, qc$cohort, FUN = flag) == 1
qc$flag_tpm   <- ave(qc$total_tpm,      qc$cohort, FUN = flag) == 1
qc$is_outlier <- qc$flag_expr | qc$flag_tpm
cat(sprintf("  Amostras flaggeadas (técnico): %d\n", sum(qc$is_outlier)))
if (sum(qc$is_outlier) > 0) print(qc[qc$is_outlier, c("sample","cohort","expressed_frac","total_tpm")])

# ── 5. PCA (top 2000 genes mais variáveis, log2) ──────────────────────────────
set.seed(42)
vars <- matrixStats::rowVars(X)
top <- order(vars, decreasing = TRUE)[1:min(2000, nrow(X))]
pca <- prcomp(t(X[top, ]), center = TRUE, scale. = TRUE)
scores <- as.data.frame(pca$x[, 1:5])
scores$condition <- meta$condition
scores$cohort <- meta$cohort
varexp <- round(100 * summary(pca)$importance[2, 1:5], 1)
cat(sprintf("  PCA: PC1=%.1f%% | PC2=%.1f%% | PC3=%.1f%%\n", varexp[1], varexp[2], varexp[3]))

# ── 6. Figuras QC ─────────────────────────────────────────────────────────────
png(file.path(dir_qc, "QC_PCA_condition.png"), width = 180, height = 150,
    units = "mm", res = 300)
plot(scores$PC1, scores$PC2, col = ifelse(scores$condition == "THCA", "#AA4488", "#4477AA"),
     pch = 19, cex = 0.5, xlab = paste0("PC1 (", varexp[1], "%)"),
     ylab = paste0("PC2 (", varexp[2], "%)"), main = "PCA — condição")
legend("topright", legend = c("THCA","Normal"), col = c("#AA4488","#4477AA"), pch = 19, cex = 0.8)
dev.off()

png(file.path(dir_qc, "QC_PCA_cohort.png"), width = 180, height = 150, units = "mm", res = 300)
plot(scores$PC1, scores$PC2, col = ifelse(scores$cohort == "TCGA", "#CC6677", "#44AA77"),
     pch = 19, cex = 0.5, xlab = paste0("PC1 (", varexp[1], "%)"),
     ylab = paste0("PC2 (", varexp[2], "%)"), main = "PCA — coorte (batch)")
legend("topright", legend = c("TCGA","GTEx"), col = c("#CC6677","#44AA77"), pch = 19, cex = 0.8)
dev.off()

png(file.path(dir_qc, "QC_distribution.png"), width = 180, height = 150, units = "mm", res = 300)
d <- density(X[sample(seq_len(nrow(X)), 50000), sample(seq_len(ncol(X)), 100)])
plot(d, main = "Distribuição de expressão (log2 TPM+0.001)", xlab = "log2(TPM+0.001)", col = "#4477AA")
dev.off()

# ── 7. Salvar QC ──────────────────────────────────────────────────────────────
data.table::fwrite(qc, file.path(dir_qc, "QC_sample_metrics.tsv"), sep = "\t")
data.table::fwrite(as.data.frame(varexp), file.path(dir_qc, "QC_pca_variance.tsv"),
                   sep = "\t", col.names = FALSE, row.names = TRUE)
writeLines(c(
  "=== QC summary ===",
  sprintf("genes=%d samples=%d", nrow(X), ncol(X)),
  sprintf("expressed_frac median=%.3f", median(expressed_frac)),
  sprintf("total_tpm median=%.1f", median(total_tpm)),
  sprintf("outliers_flagged=%d", sum(qc$is_outlier)),
  sprintf("pca_var=PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%", varexp[1], varexp[2], varexp[3]),
  sprintf("cohort=TCGA:tumor, GTEx:normal (CONFUNDIDO por desenho)")
), file.path(dir_qc, "QC_summary.txt"))

cat("\n  QC salvo em 02_qc/\n")
cat("══ FASE 2 — QC CONCLUÍDO ══\n")
