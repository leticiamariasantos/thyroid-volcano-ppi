# ═══════════════════════════════════════════════════════════════════════════════
# 02_qc.R — Controle de qualidade e pré-processamento (reconstruído do zero)
#
# Distribuição de expressão, genes pouco expressos, outliers, PCA, clustering,
# distância entre amostras, composição global, correspondência TCGA/GTEx e
# avaliação de batch/source. Documenta explicitamente: source ≡ condition.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

log_msg("══ QC e pré-processamento ══")

tpm <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
meta <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))

stopifnot(identical(colnames(tpm), meta$sample))
condition <- meta$condition
source_lbl <- meta$source

# ── 1. Distribuição de expressão por amostra ──────────────────────────────────
samp_mean <- colMeans(tpm)
samp_med  <- apply(tpm, 2, median)
samp_sd   <- apply(tpm, 2, sd)
dist_df <- data.frame(sample = colnames(tpm), mean = samp_mean, median = samp_med,
                      sd = samp_sd, condition = condition, source = source_lbl)
fwrite_tsv(dist_df, file.path(DIR_PREPROC, "sample_expression_distribution.tsv"))

p1 <- ggplot(dist_df, aes(x = condition, y = mean, fill = source)) +
  geom_boxplot(outlier.size = 0.6) + theme_minimal(base_size = 11) +
  labs(title = "Distribuição de expressão por amostra (log2 TPM)",
       y = "Média de expressão", x = "Condição") +
  theme(legend.position = "bottom")
ggsave(file.path(DIR_FIG, "QC_expression_distribution.png"), p1, width = 5.5, height = 4.5, dpi = 300)

# densidade global
set.seed(SEED)
idx <- sample(seq_len(nrow(tpm)), 5000)
dens_df <- data.frame(value = as.vector(tpm[idx, ]))
p2 <- ggplot(dens_df, aes(x = value)) + geom_density(fill = "steelblue", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "Densidade de expressão (5.000 genes amostrados)",
       x = "log2(TPM+0.001)", y = "Densidade")
ggsave(file.path(DIR_FIG, "QC_expression_density.png"), p2, width = 5.5, height = 4, dpi = 300)

# ── 2. Genes pouco expressos ───────────────────────────────────────────────────
gene_mean <- rowMeans(tpm)
gene_detected <- rowSums(tpm > log2(0.001 + 1e-9))  # > "zero" artificial
n_low <- sum(gene_detected == 0)
log_msg("Genes nunca detectados em nenhuma amostra:", n_low)
lowexpr_summary <- data.frame(
  parameter = c("genes_total", "genes_never_detected", "genes_detected_in_some"),
  value = c(nrow(tpm), n_low, nrow(tpm) - n_low)
)
fwrite_tsv(lowexpr_summary, file.path(DIR_PREPROC, "low_expression_summary.tsv"))

# ── 3. Seleção de genes variáveis para PCA (top 2000 por MAD) ─────────────────
gene_mad <- apply(tpm, 1, mad)
gene_mad <- gene_mad[gene_mad > 0]
top_var <- names(sort(gene_mad, decreasing = TRUE))[seq_len(min(2000, length(gene_mad)))]
tpm_var <- tpm[top_var, , drop = FALSE]

# ── 4. PCA ─────────────────────────────────────────────────────────────────────
pc <- prcomp(t(tpm_var), center = TRUE, scale. = TRUE)
var_expl <- (pc$sdev^2) / sum(pc$sdev^2) * 100
pca_scores <- as.data.frame(pc$x[, 1:10])
pca_scores$sample <- rownames(pca_scores)
pca_scores$condition <- condition
pca_scores$source <- source_lbl
fwrite_tsv(pca_scores, file.path(DIR_PREPROC, "pca_scores.tsv"))
fwrite_tsv(data.frame(PC = seq_along(var_expl), variance_explained = var_expl),
           file.path(DIR_PREPROC, "pca_variance.tsv"))

p3 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 0.7, alpha = 0.7) + theme_minimal(base_size = 11) +
  labs(title = "PCA por condição (2.000 genes mais variáveis)",
       x = sprintf("PC1 (%.1f%%)", var_expl[1]), y = sprintf("PC2 (%.1f%%)", var_expl[2])) +
  scale_color_manual(values = c("THCA" = "#d95f02", "Normal" = "#1b9e77")) +
  theme(legend.position = "bottom")
ggsave(file.path(DIR_FIG, "QC_PCA_condition.png"), p3, width = 5.5, height = 4.5, dpi = 300)

p4 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = source)) +
  geom_point(size = 0.7, alpha = 0.7) + theme_minimal(base_size = 11) +
  labs(title = "PCA por fonte (source ≡ condition)",
       x = sprintf("PC1 (%.1f%%)", var_expl[1]), y = sprintf("PC2 (%.1f%%)", var_expl[2])) +
  theme(legend.position = "bottom")
ggsave(file.path(DIR_FIG, "QC_PCA_source.png"), p4, width = 5.5, height = 4.5, dpi = 300)

# associação PC vs condição (ANOVA R²) — documenta força do confundimento
pc_r2_cond <- sapply(1:10, function(i) {
  m <- summary(lm(pc$x[, i] ~ condition)); m$r.squared
})
names(pc_r2_cond) <- paste0("PC", 1:10)
log_msg("R² (PC vs condição):", paste(round(pc_r2_cond, 3), collapse = ", "))

# ── 5. Distância entre amostras + clustering hierárquico ──────────────────────
cor_mat <- cor(tpm_var, method = "pearson")
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")

# heatmap de distância (amostragem de 120 amostras para legibilidade)
set.seed(SEED)
samp_idx <- sample(seq_len(ncol(cor_mat)), min(120, ncol(cor_mat)))
cor_sub <- cor_mat[samp_idx, samp_idx]
dist_sub <- as.dist(1 - cor_sub)
ann_col <- data.frame(condition = condition[samp_idx], source = source_lbl[samp_idx],
                      row.names = colnames(cor_mat)[samp_idx])
ann_colors <- list(condition = c(THCA = "#d95f02", Normal = "#1b9e77"),
                   source = c(TCGA = "#7570b3", GTEx = "#e7298a"))
png(file.path(DIR_FIG, "QC_sample_distance_heatmap.png"), width = 2000, height = 1800, res = 200)
pheatmap(cor_sub, clustering_distance_rows = dist_sub,
         clustering_distance_cols = dist_sub, annotation_col = ann_col,
         annotation_colors = ann_colors, show_rownames = FALSE, show_colnames = FALSE,
         main = "Correlação de Pearson entre amostras (120 amostradas)")
dev.off()

# ── 6. Detecção de outliers (distância de Mahalanobis sobre PCs) ──────────────
pc5 <- pc$x[, 1:5]
center <- colMeans(pc5)
cov_inv <- tryCatch(MASS::ginv(cov(pc5)), error = function(e) solve(cov(pc5) + diag(1e-6, 5)))
mah <- apply(pc5, 1, function(x) as.numeric((x - center) %*% cov_inv %*% (x - center)))
q <- qchisq(0.999, df = 5)
outliers <- names(which(mah > q))
log_msg("Outliers (Mahalanobis p<0.001):", length(outliers))
outlier_df <- data.frame(sample = colnames(tpm), mahalanobis = mah,
                         is_outlier = mah > q, condition = condition, source = source_lbl)
fwrite_tsv(outlier_df, file.path(DIR_PREPROC, "sample_outliers.tsv"))

# ── 7. Composição global (proporção de reads em genes marcadores selecionados) ─
markers <- c("EPCAM", "KRT19", "TG", "TPO", "DIO1", "DIO2",          # epitelial/tiroide
             "PTPRC", "CD68", "CD3D", "CD79A",                        # imune
             "COL1A1", "COL1A2", "DCN", "VIM", "ACTA2", "PECAM1",     # estroma/vascular
             "MYH7", "MYL1", "MYL2", "ACTA1", "TNNT3", "CKM")         # músculo
markers <- markers[markers %in% rownames(tpm)]
comp <- t(tpm[markers, , drop = FALSE])
comp_df <- as.data.frame(comp); comp_df$sample <- rownames(comp_df)
comp_df$condition <- condition
comp_long <- reshape2::melt(comp_df, id.vars = c("sample", "condition"),
                            variable.name = "marker", value.name = "expr")
comp_summary <- aggregate(expr ~ condition + marker, data = comp_long, FUN = mean)
fwrite_tsv(comp_summary, file.path(DIR_PREPROC, "marker_composition_summary.tsv"))

# ── 8. Sumário do QC ───────────────────────────────────────────────────────────
qc_lines <- c(
  "═══ QC — FASE 2 REBOOT ═══",
  paste("Data:", EXEC_DATE),
  paste("Amostras:", ncol(tpm), "| Genes:", nrow(tpm)),
  paste("Genes nunca detectados:", n_low),
  paste("Variância explicada PC1..PC5 (%):",
        paste(round(var_expl[1:5], 2), collapse = ", ")),
  paste("R² (PC vs condição) PC1..PC5:",
        paste(round(pc_r2_cond[1:5], 3), collapse = ", ")),
  paste("Outliers detectados (Mahalanobis p<0.001):", length(outliers)),
  "",
  "CONFUNDIMENTO ESTRUTURAL: source ≡ condition.",
  "TCGA ≡ tumor e GTEx ≡ normal estão perfeitamente confundidos. Nenhuma",
  "correção de batch elimina esse problema; a limitação será tratada por",
  "triangulação e análises de sensibilidade."
)
writeLines(qc_lines, file.path(DIR_PREPROC, "QC_summary.txt"))
cat("\n", paste(qc_lines, collapse = "\n"), "\n\n")

# salva a matriz variável para reuso
saveRDS(tpm_var, file.path(DIR_DATAOUT, "tpm_topvar.rds"))
saveRDS(pca_scores, file.path(DIR_DATAOUT, "pca_scores.rds"))

log_msg("══ QC concluído ══")
