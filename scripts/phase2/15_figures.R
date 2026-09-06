# ═══════════════════════════════════════════════════════════════════════════════
# 15_figures.R — Figuras de alta qualidade (cada figura responde uma pergunta)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(igraph)
  library(RColorBrewer)
})

log_msg("══ Figuras ══")

limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))
tpm <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
panel <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))
gsea_panel <- fread(file.path(DIR_PATH, "PANEL_30_GSEA_RESULTS.tsv"))
gsea_global <- fread(file.path(DIR_GSEA, "GSEA_GLOBAL_ALL.tsv"))
comp <- fread(file.path(DIR_SENS, "gsea_composition_sensitivity.tsv"))
sc_agg <- fread(file.path(DIR_VAL, "singlecell_expression_by_compartment.tsv"))

# ── 1. Volcano plot (limma) ───────────────────────────────────────────────────
log_msg("  Volcano plot...")
vol <- copy(limma)
vol[, sig := fifelse(regulation == "Up", "Up", fifelse(regulation == "Down", "Down", "NS"))]
label_genes <- c("ITGA2","FN1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM","TPO","TG","EPCAM","KRT19")
vol[, label := ifelse(gene_symbol %in% label_genes, gene_symbol, NA_character_)]
p_vol <- ggplot(vol, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.35, alpha = 0.55) +
  scale_color_manual(values = c(Up = "#d95f02", Down = "#1b9e77", NS = "grey80")) +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 20, na.rm = TRUE) +
  theme_minimal(base_size = 11) + labs(x = "log2 fold change (THCA vs Normal)",
    y = "-log10(FDR)", title = "Expressão diferencial (limma, log2 TPM)") +
  guides(color = guide_legend(title = NULL))
ggsave(file.path(DIR_FIG, "Fig2_Volcano.png"), p_vol, width = 7, height = 6, dpi = 300)

# ── 2. Heatmap de genes relevantes ─────────────────────────────────────────────
log_msg("  Heatmap de genes relevantes...")
hm_genes <- c("ITGA2","FN1","CCND1","EPCAM","KRT19","COL1A1","COL3A1","VIM",
              "TPO","TG","DIO1","DIO2","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM",
              "CD3D","CD68","PTPRC","B2M","HLA-DRA")
hm_genes <- hm_genes[hm_genes %in% rownames(tpm)]
meta <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))
set.seed(SEED)
thca_idx <- which(meta$condition == "THCA"); norm_idx <- which(meta$condition == "Normal")
samp_sel <- c(sample(thca_idx, 60), sample(norm_idx, 60))
hm_mat <- tpm[hm_genes, samp_sel]
hm_mat <- t(scale(t(hm_mat)))
hm_ann <- data.frame(condition = meta$condition[samp_sel], row.names = colnames(tpm)[samp_sel])
png(file.path(DIR_FIG, "Fig3_Heatmap.png"), width = 2000, height = 1600, res = 200)
pheatmap(hm_mat, annotation_col = hm_ann,
         annotation_colors = list(condition = c(THCA = "#d95f02", Normal = "#1b9e77")),
         show_colnames = FALSE, cluster_cols = TRUE, fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Genes relevantes (escala z por linha)")
dev.off()

# ── 3. GSEA global (top NES por coleção) ───────────────────────────────────────
log_msg("  GSEA global figure...")
g_top <- gsea_global[!is.na(NES)][order(-abs(NES))]
g_top <- g_top[, .SD[head(order(-abs(NES)), 12)], by = collection]
g_top[, pathway_short := sub("^KEGG_", "", sub("^HALLMARK_", "", sub("^REACTOME_", "", pathway)))]
g_top[, pathway_short := abbreviate(pathway_short, 38)]
g_top[, dir := ifelse(NES > 0, "Up_tumor", "Down_tumor")]
p_gsea <- ggplot(g_top, aes(x = reorder(pathway_short, NES), y = NES, fill = dir)) +
  geom_col() + coord_flip() + facet_wrap(~ collection, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(Up_tumor = "#d95f02", Down_tumor = "#1b9e77")) +
  theme_minimal(base_size = 9) + labs(x = NULL, y = "NES", title = "GSEA global (top 12 por coleção)") +
  theme(legend.position = "none")
ggsave(file.path(DIR_FIG, "Fig4_GSEA_global.png"), p_gsea, width = 8, height = 9, dpi = 300)

# ── 4. Dotplot das 30 vias (limma) ─────────────────────────────────────────────
log_msg("  Dotplot 30 vias...")
gsea_panel[, direction := fifelse(NES_limma > 0, "Up_tumor", "Down_tumor")]
gsea_panel[, sig := ifelse(padj_limma < FDR_THRESH, "yes", "no")]
gsea_panel[, pathway_name := panel$pathway_name[match(database_id, panel$database_id)]]
gsea_panel[, pathway_name := abbreviate(pathway_name, 40)]
p_dot <- ggplot(gsea_panel, aes(x = NES_limma, y = reorder(pathway_name, NES_limma),
                                size = geneset_size, color = direction, alpha = sig)) +
  geom_point() + scale_color_manual(values = c(Up_tumor = "#d95f02", Down_tumor = "#1b9e77")) +
  scale_alpha_manual(values = c(yes = 1, no = 0.35)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_minimal(base_size = 9) + labs(x = "NES (limma)", y = NULL, size = "Gene set size",
    title = "Dotplot das 30 vias (painel)") +
  theme(legend.position = "bottom")
ggsave(file.path(DIR_FIG, "Fig5_Dotplot30.png"), p_dot, width = 8.5, height = 8, dpi = 300)

# ── 5. NES das 30 vias (3 métodos, heatmap) ────────────────────────────────────
log_msg("  NES 30 vias (3 métodos)...")
nes3 <- as.matrix(gsea_panel[, .(limma = NES_limma, voom = NES_voom, DESeq2 = NES_deseq2)])
rownames(nes3) <- gsea_panel$pathway_name
png(file.path(DIR_FIG, "Fig6_NES_heatmap.png"), width = 1200, height = 2200, res = 180)
pheatmap(nes3, cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 7,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "NES das 30 vias (limma | voom | DESeq2)", na_col = "grey90")
dev.off()

# ── 6. FULL vs composition-controlled ──────────────────────────────────────────
log_msg("  Composição figure...")
comp_fig <- comp[!is.na(NES_full) & !is.na(NES_removed)]
comp_fig[, pathway_name := panel$pathway_name[match(pathway, panel$database_id)]]
p_comp <- ggplot(comp_fig, aes(x = NES_full, y = NES_removed, color = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = ifelse(status == "ROBUSTA", pathway_name, "")),
                  size = 2.6, max.overlaps = 20) +
  scale_color_manual(values = c(ROBUSTA = "#d95f02", DESAPARECE = "#1b9e77",
                                SURGE = "#7570b3", NS = "grey75", INDETERMINADO = "grey50")) +
  theme_minimal(base_size = 10) + labs(x = "NES (transcriptoma completo)",
    y = "NES (genes de composição removidos)", color = NULL,
    title = "Sensibilidade de composição muscular (30 vias)")
ggsave(file.path(DIR_FIG, "Fig7_Composition.png"), p_comp, width = 7, height = 6, dpi = 300)

# ── 7. PPI (rede principal) ────────────────────────────────────────────────────
log_msg("  PPI figure...")
g <- readRDS(file.path(DIR_DATAOUT, "ppi_graph.rds"))
if (vcount(g) > 1) {
  set.seed(SEED)
  lay <- layout_with_fr(g)
  deg <- degree(g)
  png(file.path(DIR_FIG, "Fig8_PPI.png"), width = 1800, height = 1800, res = 180)
  plot(g, layout = lay, vertex.size = 3 + 2 * sqrt(deg), vertex.label.cex = 0.7,
       vertex.label.color = "black", vertex.color = "#4daf4a", vertex.frame.color = "grey30",
       edge.width = 0.6, edge.color = "grey60", main = "Rede PPI (STRING >= 700)")
  dev.off()
}

# ── 8. Single-cell localization ITGA2/FN1/CCND1 ────────────────────────────────
log_msg("  Single-cell figure...")
sc_fig <- sc_agg[gene %in% c("ITGA2","FN1","CCND1")]
sc_fig[, cell_type := factor(cell_type, levels = c("Epithelial_tumor","Fibroblast","Endothelial",
        "SMC_pericyte","Myeloid","Mast","T_NK","B_cell"))]
p_sc <- ggplot(sc_fig, aes(x = cell_type, y = mean_expr_across_samples, fill = cell_type)) +
  geom_col() + facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 9) + labs(x = NULL, y = "Expressão média (contagens)",
    title = "Localização celular (GSE232237, marker-based)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none") +
  scale_fill_brewer(palette = "Set2")
ggsave(file.path(DIR_FIG, "Fig9_SingleCell.png"), p_sc, width = 9, height = 4, dpi = 300)

# ── 9. Modelo conceitual final ─────────────────────────────────────────────────
log_msg("  Modelo conceitual...")
model_df <- data.frame(
  stage = factor(c("Biologia molecular","Candidato","Localização celular","Acessibilidade",
                   "Hipótese de targeting","Estratégia nanomedicinal"),
                 levels = c("Biologia molecular","Candidato","Localização celular","Acessibilidade",
                            "Hipótese de targeting","Estratégia nanomedicinal")),
  x = 1:6, y = 1,
  label = c("DE + painel 30 vias\n(ITGA2/FN1/CCND1 Up)",
            "ITGA2 priorizado como\ncandidato translacional",
            "Predominantemente\nepitelial/tumoral (scRNA)",
            "Acessibilidade potencial\n(extracelular/integrina)",
            "Hipótese: investigação de\ndirecionamento molecular",
            "Nanomedicina como\nconsequência (NÃO premissa)")
)
p_model <- ggplot(model_df, aes(x = x, y = y)) +
  geom_segment(aes(xend = c(2:6, 6), yend = y), color = "grey60",
               arrow = arrow(length = unit(0.15, "cm"))) +
  geom_point(size = 6, color = "#377eb8") +
  geom_text(aes(label = label), vjust = -1.3, size = 2.6, lineheight = 0.85) +
  xlim(0.5, 6.5) + ylim(0.8, 2.6) + theme_void() +
  labs(title = "Modelo conceitual: da biologia molecular à hipótese translacional")
ggsave(file.path(DIR_FIG, "Fig10_Conceptual_Model.png"), p_model, width = 10, height = 3.5, dpi = 300)

log_msg("══ Figuras concluídas ══")
cat("\nFiguras geradas em", DIR_FIG, "\n")
print(list.files(DIR_FIG, pattern = "\\.png$"))
