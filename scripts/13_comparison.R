# ═══════════════════════════════════════════════════════════════════════════════
# 13_comparison.R — Comparação TPM (limma) vs voom vs DESeq2
#
# Compara logFC, direção, overlap (Jaccard), ranking (Spearman) e classifica
# estabilidade de candidatos.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(data.table))

cat("══════════════════════════════════════════════════════════\n")
cat("FASE 8 — Comparação TPM vs voom vs DESeq2\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("══════════════════════════════════════════════════════════\n\n")

dir_out <- "results/counts_deg"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ── 1. Carrega os três resultados ─────────────────────────────────────────────
tpm  <- fread("04_differential_expression/DEG_full_results.tsv")
voo  <- fread(file.path(dir_out, "voom_full_results.tsv"))
dds  <- fread(file.path(dir_out, "deseq2_full_results.tsv"))

# Normaliza DESeq2 para colunas comuns
dds[, logFC := log2FoldChange]
dds[, adj.P.Val := padj]

setkey(tpm, gene_symbol); setkey(voo, gene_symbol); setkey(dds, gene_symbol)

# Genes testados em comum (símbolo)
common <- Reduce(intersect, list(tpm$gene_symbol, voo$gene_symbol, dds$gene_symbol))
cat(sprintf("Genes testados: TPM=%d voom=%d DESeq2=%d | em comum=%d\n",
            nrow(tpm), nrow(voo), nrow(dds), length(common)))

t <- tpm[common]
v <- voo[common]
d <- dds[common]
stopifnot(identical(t$gene_symbol, v$gene_symbol), identical(t$gene_symbol, d$gene_symbol))

# ── 2. Correlação de logFC ────────────────────────────────────────────────────
cat("\n── Correlação de logFC (Pearson/Spearman) ──\n")
pairs <- list(TPM_voom = c("t","v"), TPM_DESeq2 = c("t","d"), voom_DESeq2 = c("v","d"))
for (nm in names(pairs)) {
  x <- get(pairs[[nm]][1])$logFC; y <- get(pairs[[nm]][2])$logFC
  cat(sprintf("  %-14s Pearson=%.4f Spearman=%.4f\n",
              nm, cor(x, y, use="complete.obs"), cor(x, y, method="spearman", use="complete.obs")))
}

# ── 3. Concordância de direção ────────────────────────────────────────────────
cat("\n── Concordância de direção (signo logFC) ──\n")
dir_agree <- function(x, y) mean(sign(x) == sign(y) & x != 0 & y != 0)
cat(sprintf("  TPM vs voom:   %.3f\n", dir_agree(t$logFC, v$logFC)))
cat(sprintf("  TPM vs DESeq2: %.3f\n", dir_agree(t$logFC, d$logFC)))
cat(sprintf("  voom vs DESeq2: %.3f\n", dir_agree(v$logFC, d$logFC)))

# ── 4. Overlap de DEGs + Jaccard ──────────────────────────────────────────────
cat("\n── Overlap de DEGs (|logFC|>1 & FDR<0.05) ──\n")
deg_t <- t$gene_symbol[t$regulation %in% c("Up","Down")]
deg_v <- v$gene_symbol[v$regulation %in% c("Up","Down")]
deg_d <- d$gene_symbol[d$regulation %in% c("Up","Down")]
jaccard <- function(a, b) length(intersect(a,b)) / length(union(a,b))
cat(sprintf("  DEGs: TPM=%d voom=%d DESeq2=%d\n", length(deg_t), length(deg_v), length(deg_d)))
cat(sprintf("  overlap(TPM∩voom)=%d | (TPM∩DESeq2)=%d | (voom∩DESeq2)=%d | (3-way)=%d\n",
            length(intersect(deg_t, deg_v)), length(intersect(deg_t, deg_d)),
            length(intersect(deg_v, deg_d)), length(Reduce(intersect, list(deg_t, deg_v, deg_d)))))
cat(sprintf("  Jaccard: TPM-voom=%.3f TPM-DESeq2=%.3f voom-DESeq2=%.3f\n",
            jaccard(deg_t, deg_v), jaccard(deg_t, deg_d), jaccard(deg_v, deg_d)))

# ── 5. Correlação de ranking (signed -log10 p) ────────────────────────────────
cat("\n── Correlação de ranking (signed -log10 p) ──\n")
rankstat <- function(df) sign(df$logFC) * (-log10(pmax(df$adj.P.Val, 1e-300)))
rt <- rankstat(t); rv <- rankstat(v); rd <- rankstat(d)
cat(sprintf("  TPM vs voom:   %.4f\n", cor(rt, rv, method="spearman")))
cat(sprintf("  TPM vs DESeq2: %.4f\n", cor(rt, rd, method="spearman")))
cat(sprintf("  voom vs DESeq2: %.4f\n", cor(rv, rd, method="spearman")))

# ── 6. Tabela de concordância por gene + classificação de estabilidade ────────
cat("\n── Tabela de concordância + estabilidade ──\n")
comp <- data.table(
  gene = t$gene_symbol,
  logFC_TPM = t$logFC, logFC_voom = v$logFC, logFC_DESeq2 = d$logFC,
  FDR_TPM = t$adj.P.Val, FDR_voom = v$adj.P.Val, FDR_DESeq2 = d$adj.P.Val
)
comp[, sig_TPM := FDR_TPM < 0.05 & abs(logFC_TPM) > 1]
comp[, sig_voom := FDR_voom < 0.05 & abs(logFC_voom) > 1]
comp[, sig_DESeq2 := FDR_DESeq2 < 0.05 & abs(logFC_DESeq2) > 1]
comp[, n_sig := sig_TPM + sig_voom + sig_DESeq2]
comp[, dir_agree := (sign(logFC_TPM) == sign(logFC_voom)) &
                    (sign(logFC_TPM) == sign(logFC_DESeq2)) &
                    (sign(logFC_TPM) != 0)]
comp[, spread := pmax(abs(logFC_TPM - logFC_voom), abs(logFC_TPM - logFC_DESeq2), abs(logFC_voom - logFC_DESeq2))]

# Classificação de robustez
comp[, robustness := fifelse(
  n_sig == 3 & dir_agree, "ROBUSTO",
  fifelse(n_sig >= 2 & dir_agree, "PARCIALMENTE ROBUSTO",
          fifelse(n_sig == 1, "MÉTODO-DEPENDENTE", "NÃO ROBUSTO")))]
data.table::fwrite(comp, file.path(dir_out, "method_concordance.tsv"), sep = "\t")
cat(sprintf("  ROBUSTO=%d PARCIAL=%d MÉTODO-DEP=%d NÃO-ROBUSTO=%d\n",
            sum(comp$robustness == "ROBUSTO"), sum(comp$robustness == "PARCIALMENTE ROBUSTO"),
            sum(comp$robustness == "MÉTODO-DEPENDENTE"), sum(comp$robustness == "NÃO ROBUSTO")))

# ── 7. Genes-chave (musculares + FN1/ITGA2 + top) ─────────────────────────────
cat("\n── Genes-chave (musculares, FN1/ITGA2) ──\n")
key <- c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM","FN1","ITGA2","CTSS","HLA-DPA1","CCND1","FLNC")
key_out <- comp[gene %in% key]
key_out <- key_out[, .(gene, logFC_TPM, logFC_voom, logFC_DESeq2, FDR_TPM, FDR_voom, FDR_DESeq2, n_sig, robustness)]
print(key_out)
data.table::fwrite(key_out, file.path(dir_out, "key_genes_concordance.tsv"), sep = "\t")

# ── 8. Ranking de estabilidade (top DEGs por método) ──────────────────────────
cat("\n── Ranking de estabilidade (top 30 por TPM) ──\n")
comp[order(-abs(logFC_TPM))][1:30,
  .(gene, logFC_TPM, logFC_voom, logFC_DESeq2, robustness)][, print(.SD)]
data.table::fwrite(comp[order(-abs(logFC_TPM))][1:100],
                   file.path(dir_out, "stability_ranking.tsv"), sep = "\t")

cat("\n══════════════════════════════════════════════════════════\n")
cat("FASE 8 — Comparação CONCLUÍDA.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
