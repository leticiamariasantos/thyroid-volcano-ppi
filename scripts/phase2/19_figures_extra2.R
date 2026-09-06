# ═══════════════════════════════════════════════════════════════════════════════
# 19_figures_extra2.R — Figuras complementares (lote 2)
#
# 17. UMAP das amostras por condição
# 18. Barcode/enrichment plot (fgsea) para via top (Ribosome)
# 19. Volcano plots adicionais (voom e DESeq2)
# 20. Correlação RNA-proteína (CCND1, RPPA vs TPM)
# 21. Heatmap da matriz de robustez (30 vias)
# 22. GSEA global por coleção (NES top 20 KEGG/Reactome/Hallmark)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(uwot)
  library(fgsea)
  library(pheatmap)
  library(RColorBrewer)
  library(httr)
  library(jsonlite)
})

log_msg("══ Figuras complementares (lote 2) ══")

meta <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))
tpm  <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))

# ── 17. UMAP das amostras ──────────────────────────────────────────────────────
log_msg("  17. UMAP...")
tpm_var <- readRDS(file.path(DIR_DATAOUT, "tpm_topvar.rds"))
cond <- meta$condition[match(colnames(tpm_var), meta$sample)]
set.seed(SEED)
um <- uwot::umap(t(tpm_var), n_neighbors = 30, min_dist = 0.3, metric = "correlation",
                 ret_model = FALSE, verbose = FALSE)
um_df <- data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], condition = cond)
p17 <- ggplot(um_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 0.7, alpha = 0.7) + theme_minimal(base_size = 11) +
  scale_color_manual(values = c(THCA = "#d95f02", Normal = "#1b9e77")) +
  labs(title = "UMAP das amostras (2.000 genes mais variáveis)", color = NULL)
ggsave(file.path(DIR_FIG, "Fig17_UMAP.png"), p17, width = 6, height = 5, dpi = 300)

# ── 18. Barcode/enrichment plot (Ribosome) ─────────────────────────────────────
log_msg("  18. Enrichment plot (Ribosome)...")
rank_t <- readRDS(file.path(DIR_DATAOUT, "limma_rank_t.rds"))
kegg_sets <- readRDS(file.path(DIR_DATAOUT, "gsea_kegg_genesets.rds"))
ribo <- kegg_sets[["KEGG_RIBOSOME"]]
ribo <- ribo[ribo %in% names(rank_t)]
p18 <- fgsea::plotEnrichment(ribo, rank_t) + ggplot2::ggtitle("KEGG RIBOSOME (leading edge no topo = Up no tumor)")
ggsave(file.path(DIR_FIG, "Fig18_Enrichment_Ribosome.png"), p18, width = 7, height = 5, dpi = 300)

# ── 19. Volcano plots (voom e DESeq2) ──────────────────────────────────────────
log_msg("  19. Volcano voom/DESeq2...")
voom <- fread(file.path(DIR_DE, "voom_full_results.tsv"))
dds  <- fread(file.path(DIR_DE, "deseq2_full_results.tsv"))
mk_volcano <- function(d, label, lfc_col = "logFC", fdr_col = "adj.P.Val") {
  dd <- d[, .(logFC = get(lfc_col), fdr = get(fdr_col), gene = gene_symbol)]
  dd[, sig := fifelse(fdr < FDR_THRESH & logFC > LOGFC_THRESH, "Up",
               fifelse(fdr < FDR_THRESH & logFC < -LOGFC_THRESH, "Down", "NS"))]
  dd[, lab := ifelse(gene %in% c("ITGA2", "FN1", "CCND1", "MYH7", "TPO", "TG"), gene, NA_character_)]
  ggplot(dd, aes(x = logFC, y = -log10(fdr), color = sig)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_manual(values = c(Up = "#d95f02", Down = "#1b9e77", NS = "grey80")) +
    geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = lab), size = 2.5, max.overlaps = 10, na.rm = TRUE) +
    theme_minimal(base_size = 10) + labs(x = "log2FC", y = "-log10(FDR)", title = label) +
    theme(legend.position = "none")
}
p19v <- mk_volcano(voom, "Volcano (limma-voom)")
p19d <- mk_volcano(dds, "Volcano (DESeq2)")
p19 <- cowplot::plot_grid(p19v, p19d, ncol = 2)
ggsave(file.path(DIR_FIG, "Fig19_Volcano_Voom_DESeq2.png"), p19, width = 11, height = 5.5, dpi = 300)

# ── 20. Correlação RNA-proteína (CCND1, RPPA vs TPM) ───────────────────────────
log_msg("  20. Correlação RNA-proteína (CCND1)...")
rppa <- tryCatch({
  u <- "https://www.cbioportal.org/api/molecular-profiles/thca_tcga_pan_can_atlas_2018_rppa/molecular-data?sampleListId=thca_tcga_pan_can_atlas_2018_all&entrezGeneId=595"
  fromJSON(content(GET(u, timeout(120)), as = "text", encoding = "UTF-8"), flatten = TRUE)
}, error = function(e) NULL)
if (!is.null(rppa) && is.data.frame(rppa) && nrow(rppa) > 0) {
  rppa_dt <- as.data.table(rppa)[, .(sample = sub("-01.*", "", sampleId), rppa = as.numeric(value))]
  rppa_dt <- rppa_dt[!is.na(rppa)]
  ccnd1 <- data.table(sample = sub("-01.*", "", colnames(tpm)), rna = as.numeric(tpm["CCND1", ]))
  m <- merge(rppa_dt, ccnd1, by = "sample", all = FALSE)
  rho <- cor(m$rppa, m$rna, method = "spearman")
  p20 <- ggplot(m, aes(x = rna, y = rppa)) + geom_point(size = 1, alpha = 0.5, color = "#377eb8") +
    geom_smooth(method = "lm", se = FALSE, color = "#e41a1c", linewidth = 0.7) +
    theme_minimal(base_size = 10) +
    labs(x = "CCND1 RNA (log2 TPM)", y = "CCND1 RPPA (Cyclin D1)",
         title = sprintf("Correlação RNA-proteína CCND1 (Spearman %.3f, n=%d)", rho, nrow(m)))
  ggsave(file.path(DIR_FIG, "Fig20_RPPA_Correlation.png"), p20, width = 6, height = 5, dpi = 300)
} else {
  log_msg("    AVISO: RPPA indisponível; figura 20 pulada.")
}

# ── 21. Heatmap da matriz de robustez ──────────────────────────────────────────
log_msg("  21. Heatmap da matriz de robustez...")
rob <- fread(file.path(DIR_PATH, "PATHWAY_ROBUSTNESS_MATRIX.tsv"))
mat <- as.matrix(rob[, .(sig_panel, sig_global, consistent_methods, persists_composition)])
rownames(mat) <- rob$pathway_name
mat <- mat * 1
png(file.path(DIR_FIG, "Fig21_Robustness_Heatmap.png"), width = 1600, height = 2000, res = 180)
pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 8,
         color = c("grey90", "#1b9e77"), breaks = c(-0.5, 0.5, 1.5),
         legend_breaks = c(0, 1), legend_labels = c("não", "sim"),
         main = "Matriz de robustez das 30 vias (4 critérios binários)")
dev.off()

# ── 22. GSEA global por coleção (NES top 20) ───────────────────────────────────
log_msg("  22. GSEA global por coleção (top 20)...")
gsea_all <- fread(file.path(DIR_GSEA, "GSEA_GLOBAL_ALL.tsv"))
gsea_all <- gsea_all[!is.na(NES)]
top <- gsea_all[, .SD[order(-abs(NES))][1:min(20, .N)], by = collection]
top[, short := abbreviate(sub("^KEGG_", "", sub("^HALLMARK_", "", pathway)), 45)]
top[, dir := ifelse(NES > 0, "Up_tumor", "Down_tumor")]
p22 <- ggplot(top, aes(x = reorder(short, NES), y = NES, fill = dir)) +
  geom_col() + coord_flip() + facet_wrap(~ collection, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(Up_tumor = "#d95f02", Down_tumor = "#1b9e77")) +
  theme_minimal(base_size = 8.5) + labs(x = NULL, y = "NES", title = "GSEA global (top 20 por coleção)") +
  theme(legend.position = "none")
ggsave(file.path(DIR_FIG, "Fig22_GSEA_Global_Top20.png"), p22, width = 9, height = 12, dpi = 300)

log_msg("══ Figuras complementares (lote 2) concluídas ══")
print(list.files(DIR_FIG, pattern = "Fig(17|18|19|20|21|22)"))
