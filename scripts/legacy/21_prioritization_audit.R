# ═══════════════════════════════════════════════════════════════════════════════
# 21_prioritization_audit.R — PARTE 6: auditoria da priorização
#
# Avalia: (a) estabilidade de pesos; (b) remoção de genes musculares;
# (c) troca de método (voom/DESeq2 no lugar de TPM); (d) efeito da composição.
# NÃO altera pesos para favorecer FN1/ITGA2.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table)})

cat("══ PARTE 6 — Auditoria da priorização ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_p <- "results/prioritization"
dir.create(dir_p, recursive = TRUE, showWarnings = FALSE)

pri <- fread("09_target_prioritization/prioritization_scores.tsv")
voo <- fread("results/counts_deg/voom_full_results.tsv")
dds <- fread("results/counts_deg/deseq2_full_results.tsv")

muscle <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
            "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
            "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
            "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")

z <- function(x) { s <- as.numeric(scale(x)); s[is.na(s)] <- 0; s }

# ── 1. Recomputar composite sob pesos variáveis (sensibilidade de pesos) ──────
# composite(w) = w*(0.5*z_lfc + 0.5*z_sig) + (1-w)*(0.6*z_bet + 0.4*z_deg)
rank_for_weight <- function(w, d = pri) {
  de <- 0.5*d$z_lfc + 0.5*d$z_sig
  net <- 0.6*d$z_bet + 0.4*d$z_deg
  w*de + (1-w)*net
}
weight_grid <- c(0, 0.3, 0.5, 0.6, 0.8, 1.0)
top20 <- function(x) head(x[order(-x)], 20)
w_tab <- rbindlist(lapply(weight_grid, function(w) {
  sc <- rank_for_weight(w)
  data.table(DE_weight = w, top20 = paste(pri$gene_symbol[order(-sc)][1:20], collapse = ","))
}))
# Jaccard top20 vs peso original 0.6
base <- pri$gene_symbol[order(-pri$composite)][1:20]
jacc <- function(a,b) length(intersect(a,b))/length(union(a,b))
w_tab[, jaccard_vs_06 := vapply(weight_grid, function(w) {
  jacc(pri$gene_symbol[order(-rank_for_weight(w))][1:20], base)
}, numeric(1))]
cat("── Sensibilidade de pesos (top-20 Jaccard vs peso=0.6) ──\n")
print(w_tab[, .(DE_weight, jaccard_vs_06)])

# ── 2. Remoção de genes musculares ────────────────────────────────────────────
pri_nm <- pri[!(gene_symbol %in% muscle)]
base_nm <- pri_nm$gene_symbol[order(-pri_nm$composite)][1:20]
cat("\n── Top 20 SEM músculo (composite) ──\n")
print(pri_nm$gene_symbol[order(-pri_nm$composite)][1:20])

# ── 3. Troca de método (voom e DESeq2 no lugar de TPM) ────────────────────────
# mantém topologia (degree/betweenness) da PPI, troca DE por voom/DESeq2
build_alt <- function(alt, label) {
  m <- merge(pri[, .(gene_symbol, degree, betweenness, z_bet, z_deg)],
             alt[, .(gene_symbol, logFC, adj.P.Val)], by = "gene_symbol")
  m[, z_lfc := z(abs(logFC))]
  m[, z_sig := z(-log10(pmax(adj.P.Val, .Machine$double.xmin)))]
  m[, composite := 0.30*z_lfc + 0.30*z_sig + 0.25*z_bet + 0.15*z_deg]
  m[order(-composite)]
}
alt_v <- build_alt(voo[, .(gene_symbol, logFC, adj.P.Val)], "voom")
alt_d <- build_alt(dds[, .(gene_symbol, log2FoldChange, padj)][, .(gene_symbol, logFC = log2FoldChange, adj.P.Val = padj)], "DESeq2")

cat("\n── Top 20 por método (composite) ──\n")
cat("  TPM:  ", paste(base, collapse=","), "\n")
cat("  voom: ", paste(alt_v$gene_symbol[1:20], collapse=","), "\n")
cat("  DESeq2:", paste(alt_d$gene_symbol[1:20], collapse=","), "\n")
cat(sprintf("  Jaccard top20 TPM×voom=%.3f TPM×DESeq2=%.3f voom×DESeq2=%.3f\n",
    jacc(base, alt_v$gene_symbol[1:20]), jacc(base, alt_d$gene_symbol[1:20]),
    jacc(alt_v$gene_symbol[1:20], alt_d$gene_symbol[1:20])))

# ── 4. Estabilidade por gene-chave (rank em cada cenário) ─────────────────────
key <- c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")
rank_of <- function(x, gene) { i <- which(x == gene); if (length(i) == 0) NA_integer_ else i[1] }
fin <- rbindlist(lapply(key, function(g) {
  r_tpm <- rank_of(pri$gene_symbol[order(-pri$composite)], g)
  r_nm  <- rank_of(pri_nm$gene_symbol[order(-pri_nm$composite)], g)
  r_v   <- rank_of(alt_v$gene_symbol[order(-alt_v$composite)], g)
  r_d   <- rank_of(alt_d$gene_symbol[order(-alt_d$composite)], g)
  data.table(gene = g, rank_TPM = r_tpm, rank_no_muscle = r_nm,
             rank_voom = r_v, rank_DESeq2 = r_d,
             is_muscle = g %in% muscle)
}))
fin[, rank_shift := rank_TPM - ifelse(is.na(rank_no_muscle), rank_TPM, rank_no_muscle)]
cat("\n── Rank dos genes-chave em cada cenário ──\n")
print(fin)
data.table::fwrite(fin, file.path(dir_p, "final_stability.tsv"), sep = "\t")
data.table::fwrite(w_tab, file.path(dir_p, "weight_sensitivity.tsv"), sep = "\t")

cat("\n══ Auditoria da priorização CONCLUÍDA ══\n")
