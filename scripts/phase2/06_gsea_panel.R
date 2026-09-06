# ═══════════════════════════════════════════════════════════════════════════════
# 06_gsea_panel.R — GSEA do painel de 30 vias (pré-especificado), por método
#
# Para cada via calcula: NES, pval, padj, nº de genes sobrepostos, leading edge
# e direção. Executa de forma independente para limma (t), voom (t) e DESeq2
# (stat). Pergunta: "Como se comportam programas biologicamente pré-especificados?"
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
})

log_msg("══ GSEA do painel de 30 vias ══")

panel <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))
genesets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))

# ── Rankings nomeados por símbolo gênico ──────────────────────────────────────
limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))
voom  <- fread(file.path(DIR_DE, "voom_full_results.tsv"))
dds   <- fread(file.path(DIR_DE, "deseq2_full_results.tsv"))

rank_limma <- setNames(limma$t, limma$gene_symbol)
rank_limma <- sort(rank_limma[!duplicated(names(rank_limma))], decreasing = TRUE)

rank_voom <- setNames(voom$t, voom$gene_symbol)
rank_voom <- sort(rank_voom[!duplicated(names(rank_voom))], decreasing = TRUE)

rank_dds <- setNames(dds$stat, dds$gene_symbol)
rank_dds <- sort(rank_dds[!duplicated(names(rank_dds))], decreasing = TRUE)

# ── GSEA para um método ────────────────────────────────────────────────────────
gsea_one <- function(rankvec, label) {
  gs <- lapply(genesets, function(g) unique(g[g %in% names(rankvec)]))
  gs <- gs[lengths(gs) > 0]
  set.seed(SEED)
  res <- fgseaMultilevel(pathways = gs, stats = rankvec,
                         minSize = 5, maxSize = 2000, eps = 0, nPermSimple = 10000)
  res <- as.data.table(res)
  res$method <- label
  res
}

log_msg("  fgsea: limma...")
g_limma <- gsea_one(rank_limma, "limma")
log_msg("  fgsea: voom...")
g_voom  <- gsea_one(rank_voom, "limma-voom")
log_msg("  fgsea: DESeq2...")
g_dds   <- gsea_one(rank_dds, "DESeq2")

# ── Formato longo (para figuras e auditoria) ──────────────────────────────────
g_long <- rbindlist(list(
  copy(g_limma)[, method := "limma"],
  copy(g_voom)[, method := "limma-voom"],
  copy(g_dds)[, method := "DESeq2"]
), fill = TRUE)
g_long <- merge(g_long, panel[, .(database_id, pathway_name, panel_group)],
                by.x = "pathway", by.y = "database_id", all.x = TRUE)
fwrite_tsv(g_long, file.path(DIR_PATH, "PANEL_30_GSEA_LONG.tsv"))

# ── Consolidação em formato largo (uma linha por via) ─────────────────────────
build_wide <- function(g, lab) {
  g <- copy(g)
  g[, overlap := sapply(g$leadingEdge, length)]
  g[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ",")]
  out <- data.table(
    pathway = g$pathway,
    NES = g$NES, pval = g$pval, padj = g$padj,
    overlap = g$overlap, leadingEdge = g$leadingEdge
  )
  setnames(out, c("NES","pval","padj","overlap","leadingEdge"),
           paste0(c("NES","pval","padj","overlap","leadingEdge"), "_", lab))
  out
}
w_limma <- build_wide(g_limma, "limma")
w_voom  <- build_wide(g_voom, "voom")
w_dds   <- build_wide(g_dds, "deseq2")

panel_wide <- merge(panel[, .(panel_group, database_id, pathway_name, geneset_size)],
                    w_limma, by.x = "database_id", by.y = "pathway", all.x = TRUE)
panel_wide <- merge(panel_wide, w_voom, by.x = "database_id", by.y = "pathway", all.x = TRUE)
panel_wide <- merge(panel_wide, w_dds, by.x = "database_id", by.y = "pathway", all.x = TRUE)

# direção do enriquecimento (limma, método principal)
panel_wide$direction_limma <- fifelse(panel_wide$NES_limma > 0, "Up_tumor",
                               fifelse(panel_wide$NES_limma < 0, "Down_tumor", "NS"))
panel_wide$direction_limma[is.na(panel_wide$NES_limma)] <- "NS"

fwrite_tsv(panel_wide, file.path(DIR_PATH, "PANEL_30_GSEA_RESULTS.tsv"))

# ── Significância por painel e método (multiple testing já embutido no fgsea) ──
sig_summary <- g_long[, .(
  n_pathways = uniqueN(pathway),
  n_sig_limma = sum(method == "limma" & padj < FDR_THRESH, na.rm = TRUE),
  n_sig_voom  = sum(method == "limma-voom" & padj < FDR_THRESH, na.rm = TRUE),
  n_sig_deseq2 = sum(method == "DESeq2" & padj < FDR_THRESH, na.rm = TRUE)
), by = panel_group]
fwrite_tsv(sig_summary, file.path(DIR_PATH, "PANEL_30_SIGNIFICANCE_SUMMARY.tsv"))

log_msg("══ GSEA do painel concluída ══")
cat("\n=== Resultado do painel (limma) — ordenado por |NES| ===\n")
show <- panel_wide[order(-abs(NES_limma))]
print(show[, .(database_id, pathway_name, geneset_size,
               NES_limma, padj_limma, direction_limma)], topn = 30)
