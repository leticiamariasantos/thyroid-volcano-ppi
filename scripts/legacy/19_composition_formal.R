# ═══════════════════════════════════════════════════════════════════════════════
# 19_composition_formal.R — Formalização da composição (PARTE 2/3) + sensibilidade final (PARTE 4)
#
# DECISÃO METODOLÓGICA (documentada): deconvolução formal (CIBERSORTx, EPIC, MuSiC,
# BisqueRNA, MCPcounter, xCell) NÃO é aplicável aqui porque:
#   - MuSiC/BisqueRNA exigem referência single-cell de tireoide (não disponível);
#   - CIBERSORTx/EPIC/xCell/MCPcounter não possuem assinatura de MÚSCULO ESQUELÉTICO
#     nem de epitélio tireoidiano (a hipótese central é justamente contaminação muscular);
#   - CIBERSORTx exige Docker (indisponível).
# Portanto usa-se marker-based scoring (MCPcounter-like): score por amostra = média da
# expressão log2(TPM+0.001) dos marcadores canônicos de cada componente. É uma
# ESTIMATIVA/MODELO, não observação experimental. Produz estatísticas + FDR.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table); library(stats)})

cat("══ FASE: Composição formal (marker-based) + sensibilidade final ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_c <- "results/composition"
dir.create(dir_c, recursive = TRUE, showWarnings = FALSE)

# ── 1. Scores por amostra (já computados) ─────────────────────────────────────
scores <- fread(file.path(dir_c, "composition_per_sample.tsv"))
components <- c("muscle_skeletal","thyroid_epithelial","fibroblast_ecm","endothelial","immune")

# ── 2. Estatística formal: mean, Cohen's d, Wilcoxon, BH-FDR ──────────────────
cohens_d <- function(a, b) {
  sp <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b)) / (length(a)+length(b)-2))
  (mean(a) - mean(b)) / sp
}
summ <- rbindlist(lapply(components, function(cc) {
  a <- scores[[cc]][scores$condition == "THCA"]
  b <- scores[[cc]][scores$condition == "Normal"]
  data.table(component = cc,
             TCGA_mean = mean(a), GTEx_mean = mean(b),
             effect = mean(a) - mean(b),
             cohens_d = cohens_d(a, b),
             p_value = wilcox.test(a, b)$p.value)
}))
summ[, FDR := p.adjust(p_value, method = "BH")]
summ[, interpretation := fifelse(FDR < 0.05 & effect > 0, "enriquecido no tumor (THCA)",
                          fifelse(FDR < 0.05 & effect < 0, "enriquecido no normal (GTEx) — composicional",
                                  "sem diferença significativa"))]
summ[, TCGA_mean := round(TCGA_mean, 3)]
summ[, GTEx_mean := round(GTEx_mean, 3)]
summ[, effect := round(effect, 3)]
summ[, cohens_d := round(cohens_d, 3)]
summ[, p_value := signif(p_value, 3)]
summ[, FDR := signif(FDR, 3)]
print(summ)
data.table::fwrite(summ, file.path(dir_c, "composition_summary.tsv"), sep = "\t")

# ── 3. formal_deconvolution.tsv (estimativas por amostra = "deconvolução") ─────
fd <- data.table(
  sample = scores$sample,
  condition = scores$condition,
  muscle_skeletal = scores$muscle_skeletal,
  thyroid_epithelial = scores$thyroid_epithelial,
  fibroblast_ecm = scores$fibroblast_ecm,
  endothelial = scores$endothelial,
  immune = scores$immune
)
data.table::fwrite(fd, file.path(dir_c, "formal_deconvolution.tsv"), sep = "\t")
cat("  -> formal_deconvolution.tsv (marker-based, MCPcounter-like; modelo/estimativa)\n")

# ── 4. Sensibilidade final dos genes-chave (PARTE 4) ──────────────────────────
cat("\n── Sensibilidade final dos genes-chave ──\n")
key <- c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")
tpm  <- fread("04_differential_expression/DEG_full_results.tsv")
voo  <- fread("results/counts_deg/voom_full_results.tsv")
dds  <- fread("results/counts_deg/deseq2_full_results.tsv")
muscle_markers <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
                    "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
                    "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
                    "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")

rows <- lapply(key, function(g) {
  t <- tpm[gene_symbol == g]; v <- voo[gene_symbol == g]; d <- dds[gene_symbol == g]
  lfc_t <- if (nrow(t)) t$logFC[1] else NA_real_
  lfc_v <- if (nrow(v)) v$logFC[1] else NA_real_
  lfc_d <- if (nrow(d)) d$log2FoldChange[1] else NA_real_
  fdr_t <- if (nrow(t)) t$adj.P.Val[1] else NA_real_
  fdr_v <- if (nrow(v)) v$adj.P.Val[1] else NA_real_
  fdr_d <- if (nrow(d)) d$padj[1] else NA_real_
  is_muscle <- g %in% muscle_markers
  n_sig <- sum(c(fdr_t, fdr_v, fdr_d) < 0.05 & abs(c(lfc_t, lfc_v, lfc_d)) > 1, na.rm = TRUE)
  dir_ok <- (sign(lfc_t) == sign(lfc_v)) & (sign(lfc_t) == sign(lfc_d)) & !is.na(lfc_t) & !is.na(lfc_v) & !is.na(lfc_d)
  status <- if (is_muscle) "COMPOSICIONAL"
            else if (n_sig == 3 & dir_ok) "ROBUSTO"
            else if (n_sig >= 2 & dir_ok) "PARCIALMENTE ROBUSTO"
            else if (n_sig == 1) "MÉTODO-DEPENDENTE"
            else "NÃO DETERMINÁVEL"
  data.table(gene = g, logFC_TPM = lfc_t, logFC_voom = lfc_v, logFC_DESeq2 = lfc_d,
             FDR_TPM = fdr_t, FDR_voom = fdr_v, FDR_DESeq2 = fdr_d,
             composition_status = ifelse(is_muscle, "muscle marker", "non-muscle"),
             status = status)
})
sens <- rbindlist(rows)
print(sens)
data.table::fwrite(sens, file.path(dir_c, "final_sensitivity.tsv"), sep = "\t")

cat("\n══ Composição formal + sensibilidade CONCLUÍDAS ══\n")
