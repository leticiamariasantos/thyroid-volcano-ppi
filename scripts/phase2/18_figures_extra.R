# ═══════════════════════════════════════════════════════════════════════════════
# 18_figures_extra.R — Figuras complementares (perguntas científicas adicionais)
#
# 11. Concordância de logFC entre métodos (limma vs voom / limma vs DESeq2)
# 12. Expressão dos candidatos (ITGA2/FN1/CCND1) por condição
# 13. Validação externa (logFC dos candidatos em GSE33630/GSE60542/GSE224356)
# 14. Frequência de mutação (BRAF/RAS/TP53/candidatos)
# 15. Score de composição muscular por condição
# 16. Convergência de genes (leading edge recorrente nas vias robustas)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

log_msg("══ Figuras complementares ══")

limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))
voom  <- fread(file.path(DIR_DE, "voom_full_results.tsv"))
dds   <- fread(file.path(DIR_DE, "deseq2_full_results.tsv"))
tpm   <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
meta  <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))

# ── 11. Concordância de logFC entre métodos ────────────────────────────────────
log_msg("  11. Concordância de métodos...")
cmp <- merge(limma[, .(gene = gene_symbol, limma = logFC)],
             voom[, .(gene = gene_symbol, voom = logFC)], by = "gene", all = FALSE)
cmp <- merge(cmp, dds[, .(gene = gene_symbol, deseq2 = logFC)], by = "gene", all = FALSE)
p11a <- ggplot(cmp, aes(x = limma, y = voom)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 0.3, alpha = 0.4, color = "#377eb8") + theme_minimal(base_size = 10) +
  labs(x = "log2FC limma", y = "log2FC voom", title = "Concordância limma vs voom") +
  coord_fixed()
p11b <- ggplot(cmp, aes(x = limma, y = deseq2)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 0.3, alpha = 0.4, color = "#e41a1c") + theme_minimal(base_size = 10) +
  labs(x = "log2FC limma", y = "log2FC DESeq2", title = "Concordância limma vs DESeq2") +
  coord_fixed()
p11 <- cowplot::plot_grid(p11a, p11b, ncol = 2)
ggsave(file.path(DIR_FIG, "Fig11_Method_Concordance.png"), p11, width = 9, height = 4.5, dpi = 300)

# ── 12. Expressão dos candidatos por condição ──────────────────────────────────
log_msg("  12. Expressão dos candidatos...")
cand <- c("ITGA2", "FN1", "CCND1")
cand_expr <- data.table()
for (g in cand) {
  if (g %in% rownames(tpm)) {
    cand_expr <- rbind(cand_expr, data.table(
      gene = g, sample = colnames(tpm), expr = as.numeric(tpm[g, ]),
      condition = meta$condition[match(colnames(tpm), meta$sample)]
    ))
  }
}
p12 <- ggplot(cand_expr, aes(x = condition, y = expr, fill = condition)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.8) + facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c(THCA = "#d95f02", Normal = "#1b9e77")) +
  theme_minimal(base_size = 10) + labs(x = NULL, y = "log2(TPM+0.001)",
    title = "Expressão dos candidatos por condição") + theme(legend.position = "none")
ggsave(file.path(DIR_FIG, "Fig12_Candidate_Expression.png"), p12, width = 8, height = 4, dpi = 300)

# ── 13. Validação externa ──────────────────────────────────────────────────────
log_msg("  13. Validação externa...")
val <- fread(file.path(DIR_VAL, "validation_candidates.tsv"))
val_g <- val[gene %in% cand]
p13 <- ggplot(val_g, aes(x = gene, y = logFC_external, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 10) + labs(x = NULL, y = "logFC (tumor vs normal)",
    title = "Validação externa (GSE33630 e GSE60542)", fill = NULL)
ggsave(file.path(DIR_FIG, "Fig13_External_Validation.png"), p13, width = 7, height = 4.5, dpi = 300)

# ── 14. Frequência de mutação ──────────────────────────────────────────────────
log_msg("  14. Frequência de mutação...")
mut <- fread(file.path(DIR_VAL, "mutation_frequency.tsv"))
mut <- mut[order(-frequency)]
mut[, gene := factor(gene, levels = gene)]
p14 <- ggplot(mut, aes(x = gene, y = frequency * 100)) +
  geom_col(fill = "#4daf4a") + theme_minimal(base_size = 10) +
  labs(x = NULL, y = "Frequência de mutação (%)",
       title = "Frequência de mutação (THCA PanCan Atlas)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(DIR_FIG, "Fig14_Mutation_Frequency.png"), p14, width = 7, height = 4.5, dpi = 300)

# ── 15. Score de composição muscular ───────────────────────────────────────────
log_msg("  15. Score de composição...")
ms <- fread(file.path(DIR_SENS, "muscle_composition_scores.tsv"))
p15 <- ggplot(ms, aes(x = condition, y = muscle_score, fill = condition)) +
  geom_violin(alpha = 0.7, trim = TRUE) + geom_boxplot(width = 0.15, outlier.size = 0.3) +
  scale_fill_manual(values = c(THCA = "#d95f02", Normal = "#1b9e77")) +
  theme_minimal(base_size = 10) + labs(x = NULL, y = "Score muscular (média de 48 marcadores)",
    title = "Assinatura muscular por condição (artefato de composição)") + theme(legend.position = "none")
ggsave(file.path(DIR_FIG, "Fig15_Composition_Score.png"), p15, width = 6, height = 4.5, dpi = 300)

# ── 16. Convergência de genes (leading edge recorrente) ────────────────────────
log_msg("  16. Convergência de genes...")
conv <- fread(file.path(DIR_VAL, "panel_convergence_genes.tsv"))
conv <- conv[order(-n_robust_pathways_LE, -n_pathways_LE)]
conv <- conv[n_robust_pathways_LE >= 1][1:min(25, .N)]
conv[, gene := factor(gene, levels = rev(gene))]
p16 <- ggplot(conv, aes(x = gene, y = n_robust_pathways_LE)) +
  geom_col(aes(fill = factor(n_robust_pathways_LE))) +
  coord_flip() + theme_minimal(base_size = 9) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1) +
  labs(x = NULL, y = "Nº de vias robustas com o gene no leading edge",
       title = "Convergência molecular (genes recorrentes nas 6 vias robustas)") +
  theme(legend.position = "none")
ggsave(file.path(DIR_FIG, "Fig16_Convergence_Genes.png"), p16, width = 7, height = 6, dpi = 300)

log_msg("══ Figuras complementares concluídas ══")
print(list.files(DIR_FIG, pattern = "Fig1[1-6]"))
