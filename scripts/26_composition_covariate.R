# ═══════════════════════════════════════════════════════════════════════════════
# 26_composition_covariate.R — ETAPA 3: composição como covariável
#
# Modelo 1: expression ~ condition
# Modelo 2: expression ~ condition + muscle_skeletal (+ outros componentes)
# Compara logFC/FDR/ranking/DEGs para quantificar quanto do sinal é composicional.
# NÃO interpreta o modelo ajustado como causalidade.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table); library(limma)})

cat("══ ETAPA 3 — Composição como covariável ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_c <- "results/composition"
dir.create(dir_c, recursive = TRUE, showWarnings = FALSE)

# ── 1. matriz TPM + filtro (igual ao DEG) ─────────────────────────────────────
X <- as.matrix(fread("data/global/TCGA_GTEx_thyroid_tpm.tsv"), rownames = TRUE)
storage.mode(X) <- "numeric"
samples <- colnames(X)
condition <- factor(ifelse(grepl("^GTEX", samples), "Normal", "THCA"), levels = c("Normal","THCA"))
# filtro TPM>1 em >=10% amostras
keep <- rowSums(X > 1) >= (0.10 * ncol(X))
X <- X[keep, , drop = FALSE]
cat("  Matriz filtrada:", nrow(X), "genes ×", ncol(X), "amostras\n")

# ── 2. scores de composição ───────────────────────────────────────────────────
comp <- fread(file.path(dir_c, "composition_per_sample.tsv"))
comp <- comp[match(samples, comp$sample)]
muscle <- comp$muscle_skeletal
epi <- comp$thyroid_epithelial
endo <- comp$endothelial
fibro <- comp$fibroblast_ecm

# ── 3. Modelo 1: ~ condition ──────────────────────────────────────────────────
design1 <- model.matrix(~ condition)
fit1 <- eBayes(lmFit(X, design1), trend = TRUE)
tt1 <- topTable(fit1, coef = 2, number = Inf, sort.by = "none")

# ── 4. Modelo 2: ~ condition + muscle ─────────────────────────────────────────
design2 <- model.matrix(~ condition + muscle)
fit2 <- eBayes(lmFit(X, design2), trend = TRUE)
tt2 <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")

# ── 5. Modelo 3: ~ condition + 4 componentes (muscle, epi, endo, fibro) ───────
design3 <- model.matrix(~ condition + muscle + epi + endo + fibro)
fit3 <- eBayes(lmFit(X, design3), trend = TRUE)
tt3 <- topTable(fit3, coef = 2, number = Inf, sort.by = "none")

# ── 6. comparação ──────────────────────────────────────────────────────────────
cat("\n── DEGs por modelo (|logFC|>1 & FDR<0.05) ──\n")
cnt <- function(tt) c(up = sum(tt$adj.P.Val<0.05 & tt$logFC>1),
                      down = sum(tt$adj.P.Val<0.05 & tt$logFC< -1),
                      total = sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>1))
cmp <- data.table(model = c("condition","condition+muscle","condition+4components"),
                  rbind(cnt(tt1), cnt(tt2), cnt(tt3)))
print(cmp)
data.table::fwrite(cmp, file.path(dir_c, "composition_covariate_deg_counts.tsv"), sep = "\t")

# ── 7. genes-chave ─────────────────────────────────────────────────────────────
key <- c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")
cat("\n── Genes-chave: logFC/FDR sob cada modelo ──\n")
res <- rbindlist(lapply(key, function(g) {
  a <- tt1[rownames(tt1)==g,]; b <- tt2[rownames(tt2)==g,]; c <- tt3[rownames(tt3)==g,]
  data.table(gene=g,
    logFC_cond = if(nrow(a)) a$logFC else NA_real_,
    logFC_muscle = if(nrow(b)) b$logFC else NA_real_,
    logFC_4comp = if(nrow(c)) c$logFC else NA_real_,
    FDR_cond = if(nrow(a)) a$adj.P.Val else NA_real_,
    FDR_muscle = if(nrow(b)) b$adj.P.Val else NA_real_,
    FDR_4comp = if(nrow(c)) c$adj.P.Val else NA_real_)
}))
res[, logFC_cond := round(logFC_cond, 3)]
res[, logFC_muscle := round(logFC_muscle, 3)]
res[, logFC_4comp := round(logFC_4comp, 3)]
res[, muscle_marker := gene %in% c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")]
print(res)
data.table::fwrite(res, file.path(dir_c, "composition_covariate_genes.tsv"), sep = "\t")

# ── 8. correlação de ranking ──────────────────────────────────────────────────
cat("\n── Correlação de logFC entre modelos (genes comuns) ──\n")
cat("  cond vs +muscle:", cor(tt1$logFC, tt2$logFC, method="spearman"), "\n")
cat("  cond vs +4comp:", cor(tt1$logFC, tt3$logFC, method="spearman"), "\n")

cat("\n══ ETAPA 3 CONCLUÍDA ══\n")
