#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 09_composition.R — FASE 3: análise de composição celular/tecidual (marcadores)
# thyroid-volcano-ppi
#
# PERGUNTA: quanto do sinal bulk (THCA vs GTEx normal) é compatível com diferenças
# de composição celular/tecidual (ex.: músculo esquelético no tecido normal)?
#
# Marcadores: assinaturas canônicas de literatura (NÃO inventadas; fonte documentada
# como "marcadores canônicos de tipo celular"). Score por amostra = média de expressão
# (log2 TPM+0.001) dos marcadores da categoria.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(stats)
})

PROJECT_ROOT <- here::here()
dir_comp <- file.path(PROJECT_ROOT, "results", "composition")
dir.create(dir_comp, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 3 — Análise de composição (marcadores) ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Carregar matriz TPM ─────────────────────────────────────────────────────
X <- as.matrix(data.table::fread(file.path(PROJECT_ROOT, "data", "global",
                                           "TCGA_GTEx_thyroid_tpm.tsv"),
                                 showProgress = FALSE), rownames = TRUE)
storage.mode(X) <- "numeric"
samples <- colnames(X)
condition <- ifelse(grepl("^GTEX", samples), "Normal", "THCA")
genes <- rownames(X)
cat(sprintf("  Matriz: %d genes × %d amostras\n", nrow(X), ncol(X)))

# ── 2. Assinaturas de marcadores (canônicos de literatura) ─────────────────────
markers <- list(
  muscle_skeletal = c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7",
                      "ACTA1","ACTC1","TNNT1","TNNT3","TNNC1","TNNI1","TNNI2",
                      "TPM1","TPM2","TPM3","CKM","CKMT2","DES","MB","ENO3",
                      "MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
                      "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2"),
  thyroid_epithelial = c("TG","TPO","TSHR","PAX8","NKX2-1","FOXE1","SLC5A5",
                         "DIO1","DIO2","CALCA","TFF3","KRT7","KRT19"),
  fibroblast_ecm = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","DCN",
                     "LUM","VIM","ACTA2","FAP","PDGFRB","POSTN","MMP2","LOX","BGN"),
  endothelial = c("PECAM1","VWF","CDH5","ENG","KDR","FLT1","EMCN","CLDN5","TEK"),
  immune = c("PTPRC","CD3D","CD3E","CD8A","CD4","CD19","MS4A1","CD68","CD14",
             "NCAM1","FCGR3A","ITGAX","NKG7","GNLY","GZMB","PRF1","CD163","CD79A")
)

score_category <- function(mat, genes_avail, mk) {
  g <- intersect(mk, genes_avail)
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  sub <- mat[g, , drop = FALSE]
  colMeans(sub)
}

scores <- data.frame(sample = samples, condition = condition, stringsAsFactors = FALSE)
for (cat in names(markers)) {
  scores[[cat]] <- score_category(X, genes, markers[[cat]])
  cat(sprintf("  %-18s: %d/%d marcadores presentes\n", cat,
              length(intersect(markers[[cat]], genes)), length(markers[[cat]])))
}

# ── 3. Comparação THCA vs Normal (teste de Mann-Whitney por categoria) ─────────
cmp <- lapply(names(markers), function(cat) {
  th <- scores[[cat]][scores$condition == "THCA"]
  no <- scores[[cat]][scores$condition == "Normal"]
  wt <- wilcox.test(th, no)
  data.frame(category = cat,
             mean_THCA = mean(th), mean_Normal = mean(no),
             diff = mean(th) - mean(no),
             p_wilcox = wt$p.value,
             stringsAsFactors = FALSE)
})
cmp <- do.call(rbind, cmp)
cmp$p_wilcox <- signif(cmp$p_wilcox, 4)
print(cmp, row.names = FALSE)
data.table::fwrite(cmp, file.path(dir_comp, "composition_scores.tsv"), sep = "\t")
data.table::fwrite(scores, file.path(dir_comp, "composition_per_sample.tsv"), sep = "\t")

# ── 4. Correlação com PCA (usar scores de composição) ──────────────────────────
set.seed(42)
v <- apply(X, 1, var)
top <- order(v, decreasing = TRUE)[1:2000]
pca <- prcomp(t(X[top, ]), center = TRUE, scale. = TRUE)
pc <- pca$x[, 1:2]
cor_pc <- sapply(names(markers), function(cat) {
  c(cor(scores[[cat]], pc[,1], use="complete.obs"),
    cor(scores[[cat]], pc[,2], use="complete.obs"))
})
colnames(cor_pc) <- names(markers); rownames(cor_pc) <- c("PC1","PC2")
cat("\n  Correlação (score de composição × PC):\n")
print(round(cor_pc, 3))
data.table::fwrite(as.data.frame(cor_pc), file.path(dir_comp, "composition_pca_corr.tsv"),
                   sep = "\t", row.names = TRUE)

# ── 5. Top DEGs que são marcadores musculares ──────────────────────────────────
deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
muscle_genes <- intersect(markers$muscle_skeletal, genes)
deg_muscle <- deg[gene_symbol %in% muscle_genes & regulation != "NS"]
cat(sprintf("\n  DEGs que são marcadores musculares esqueléticos: %d\n", nrow(deg_muscle)))
deg_muscle <- deg_muscle[order(-abs(logFC))]
print(head(deg_muscle[, .(gene_symbol, logFC, adj.P.Val, regulation)], 20), row.names = FALSE)

cat("\n══ FASE 3 — Composição CONCLUÍDA ══\n")
