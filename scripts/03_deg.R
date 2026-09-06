#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 03_deg.R — FASE 2: Expressão Diferencial GENOMA-WIDE (limma)
# thyroid-volcano-ppi
#
# Entrada: data/global/TCGA_GTEx_thyroid_tpm.tsv  (log2(TPM+0.001), 58581×783)
# Método: limma (lmFit + contrasts.fit + eBayes, trend=TRUE)
# Design: ~0 + condition ; contraste THCA − Normal ; BH
# Ranking (GSEA): estatística t moderada (preserva magnitude + direção)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(here)
})

PROJECT_ROOT <- here::here()
dir_de <- file.path(PROJECT_ROOT, "04_differential_expression")
dir.create(dir_de, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — DEG genoma-wide (limma) ══\n")
cat("R:", as.character(getRversion()), "| limma:", as.character(packageVersion("limma")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Carregar ────────────────────────────────────────────────────────────────
X <- as.matrix(data.table::fread(file.path(PROJECT_ROOT, "data","global",
                                           "TCGA_GTEx_thyroid_tpm.tsv"),
                                 showProgress = FALSE), rownames = TRUE)
storage.mode(X) <- "numeric"
samples <- colnames(X)
condition <- factor(ifelse(grepl("^GTEX", samples), "Normal", "THCA"),
                    levels = c("Normal", "THCA"))
stopifnot(sum(condition == "THCA") == 504, sum(condition == "Normal") == 279)

# ── 2. Filtro de baixa expressão (critério pré-definido) ──────────────────────
# gene retido se TPM>1 (valor>0) em >=10% das amostras
keep <- rowMeans(X > 0) >= 0.10
Xf <- X[keep, , drop = FALSE]
cat(sprintf("  Filtro baixa expressão: %d mantidos / %d removidos\n",
            sum(keep), sum(!keep)))

# ── 3. limma ───────────────────────────────────────────────────────────────────
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
fit <- lmFit(Xf, design)
contr <- makeContrasts(THCA_vs_Normal = THCA - Normal, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = TRUE)

deg <- topTable(fit2, coef = "THCA_vs_Normal", number = Inf,
                adjust.method = "BH", sort.by = "none")
deg <- as.data.frame(deg)
deg$gene_symbol <- rownames(deg)
deg <- deg[, c("gene_symbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]

deg$regulation <- with(deg, ifelse(adj.P.Val < 0.05 & logFC >  1, "Up",
                            ifelse(adj.P.Val < 0.05 & logFC < -1, "Down", "NS")))
n_up <- sum(deg$regulation == "Up"); n_down <- sum(deg$regulation == "Down")
cat(sprintf("  Genes testados: %d | DEGs: %d (↑%d ↓%d)\n",
            nrow(deg), n_up + n_down, n_up, n_down))

# ── 4. Ranking para GSEA (estatística t moderada; preserva direção+magnitude) ──
# (genes com p=0 recebem valor mínimo de p para evitar -Inf)
pmin_floor <- .Machine$double.xmin
deg$rank_stat <- sign(deg$logFC) * (-log10(pmax(deg$P.Value, pmin_floor)))
# usar t moderada como ranking primário (padrão limma->GSEA)
ranking_t <- setNames(deg$t, deg$gene_symbol)
ranking_t <- ranking_t[is.finite(ranking_t)]
ranking_t <- ranking_t[order(ranking_t, decreasing = TRUE)]

# ── 5. Saídas ─────────────────────────────────────────────────────────────────
deg_out <- deg[order(deg$P.Value), ]
deg_out$logFC <- round(deg_out$logFC, 5); deg_out$AveExpr <- round(deg_out$AveExpr, 4)
deg_out$t <- round(deg_out$t, 3); deg_out$B <- round(deg_out$B, 2)
data.table::fwrite(deg_out, file.path(dir_de, "DEG_full_results.tsv"), sep = "\t")

summ <- data.frame(
  parameter = c("genes_tested","genes_up","genes_down","total_deg",
                "fc_threshold","fdr_threshold","method","contrast","ranking_statistic"),
  value = c(nrow(deg), n_up, n_down, n_up + n_down, 1, 0.05,
            "limma eBayes trend=TRUE", "THCA - Normal", "moderated t-statistic")
)
data.table::fwrite(summ, file.path(dir_de, "DEG_summary.tsv"), sep = "\t")

# ranking para GSEA (t moderada + alternativa signed -log10 p)
data.table::fwrite(data.frame(gene = names(ranking_t), rank = as.numeric(ranking_t)),
                   file.path(dir_de, "ranking_t.rnk"), sep = "\t", col.names = FALSE)
data.table::fwrite(deg[order(-deg$rank_stat), c("gene_symbol","rank_stat")],
                   file.path(dir_de, "ranking_signed_log10p.rnk"), sep = "\t",
                   col.names = FALSE)

# top DEGs para inspeção
top <- deg[deg$regulation != "NS", ][order(-abs(deg[deg$regulation != "NS", ]$logFC)), ]
data.table::fwrite(head(top, 50), file.path(dir_de, "DEG_top50.tsv"), sep = "\t")

cat(sprintf("  Top 5 Up:   %s\n", paste(head(deg$gene_symbol[deg$regulation == "Up"], 5), collapse=", ")))
cat(sprintf("  Top 5 Down: %s\n", paste(head(deg$gene_symbol[deg$regulation == "Down"], 5), collapse=", ")))
cat(sprintf("  Ranking t: %d genes (min=%.2f, max=%.2f)\n",
            length(ranking_t), min(ranking_t), max(ranking_t)))

cat("══ FASE 2 — DEG CONCLUÍDO ══\n")
