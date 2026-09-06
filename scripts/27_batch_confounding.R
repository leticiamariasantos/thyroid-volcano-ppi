# ═══════════════════════════════════════════════════════════════════════════════
# 27_batch_confounding.R — Pendência 5: quantificação formal do confundimento TCGA×GTEx
#
# O desenho discovery tem: TCGA = tumor, GTEx = normal. Logo "source" e "condition"
# são perfeitamente colineares (identificabilidade estrutural nula).
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(data.table); library(limma)})

cat("══ Pendência 5 — Confundimento TCGA×GTEx (formal) ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# 1. identidades
samples <- colnames(fread("data/global/TCGA_GTEx_thyroid_tpm.tsv", nrows = 0))
samples <- samples[-1]  # remove coluna de gene ("V1")
condition <- ifelse(grepl("^GTEX", samples), "Normal", "THCA")
source <- ifelse(grepl("^GTEX", samples), "GTEx", "TCGA")
tab <- table(condition, source)
cat("── Tabela condition × source ──\n")
print(tab)
cat("\n  → source e condition são PERFEITAMENTE colineares (cada source = 1 condition).\n")
cat("  → Correlation(source_indicator, condition_indicator) = 1.000\n")
cat("  → VIF = Infinito (não identificável separadamente).\n")
cat("  → ComBat/removeBatchEffect NÃO pode separar source de condition sem destruir o sinal biológico.\n\n")

# 2. demonstração de rank-deficiência: source não é estimável além de condition
design_full <- model.matrix(~ condition + source)
cat("  rank(design ~ condition + source) =", qr(design_full)$rank,
    "| colunas =", ncol(design_full), "→ rank-deficiente (source redundante).\n\n")

# 3. quantificação da variância explicada por condition (proxy do "batch confundido")
X <- as.matrix(fread("data/global/TCGA_GTEx_thyroid_tpm.tsv"), rownames = TRUE)
storage.mode(X) <- "numeric"
keep <- rowSums(X > 1) >= (0.10 * ncol(X))
X <- X[keep, , drop = FALSE]
cond <- factor(condition, levels = c("Normal","THCA"))
d <- model.matrix(~ cond)
# variância total (genes mais variáveis) explicada por condition (R² ajustado)
v <- apply(X, 1, var)
top <- order(v, decreasing = TRUE)[1:2000]
fit <- lmFit(X[top,], d)
fitted_vals <- fit$coefficients %*% t(d)
r2 <- sapply(seq_along(top), function(i) {
  ss_res <- sum((X[top[i],] - fitted_vals[i,])^2)
  ss_tot <- sum((X[top[i],] - mean(X[top[i],]))^2)
  1 - ss_res/ss_tot
})
cat(sprintf("  R² médio de condition nos 2000 genes mais variáveis = %.3f (mediana %.3f)\n",
            mean(r2), median(r2)))

# 4. conclusão formal
cat("\n── CONCLUSÃO ──\n")
cat("  O confundimento TCGA×GTEx é ESTRUTURAL e TOTAL (source == condition).\n")
cat("  Não é modelável separadamente. A única solução metodológica válida é a\n")
cat("  TRIANGULAÇÃO por coortes independentes com desenho pareado/controlado\n")
cat("  (GSE33630 e GSE60542 — já usadas), não correção de batch.\n")
cat("  Registrado como limitação estrutural irremediável do desenho discovery.\n\n")

data.table::fwrite(data.table(
  metric = c("correlation_source_condition","rank_deficiency","mean_R2_condition_top2000","identifiability"),
  value = c("1.000 (perfeita)","sim (rank < colunas)", as.character(round(mean(r2),3)), "não identificável separadamente")
), "results/composition/batch_confounding_summary.tsv", sep = "\t")

cat("══ Pendência 5 CONCLUÍDA ══\n")
