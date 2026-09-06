# ═══════════════════════════════════════════════════════════════════════════════
# 07_redundancy.R — Redundância entre as 30 vias do painel
#
# Calcula overlap de genes (Jaccard), similaridade de gene sets, correlação de
# NES (entre métodos) e clustering das vias. Gera PATHWAY_REDUNDANCY.tsv e figura.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

log_msg("══ Redundância entre vias ══")

panel <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))
genesets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))
gsea <- fread(file.path(DIR_PATH, "PANEL_30_GSEA_RESULTS.tsv"))

ids <- panel$database_id
names(ids) <- panel$database_id
n <- length(ids)

# ── 1. Matriz de Jaccard entre gene sets ──────────────────────────────────────
jaccard_mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
overlap_mat <- matrix(NA_integer_, n, n, dimnames = list(ids, ids))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    gi <- genesets[[ids[i]]]; gj <- genesets[[ids[j]]]
    inter <- length(intersect(gi, gj)); uni <- length(union(gi, gj))
    overlap_mat[i, j] <- inter
    jaccard_mat[i, j] <- if (uni > 0) inter / uni else NA_real_
  }
}

# ── 2. Correlação de NES entre métodos (para pares de vias) ───────────────────
nes_limma <- setNames(gsea$NES_limma, gsea$database_id)
nes_voom  <- setNames(gsea$NES_voom, gsea$database_id)
nes_dds   <- setNames(gsea$NES_deseq2, gsea$database_id)
nes_limma <- nes_limma[ids]; nes_voom <- nes_voom[ids]; nes_dds <- nes_dds[ids]

# ── 3. Tabela de redundância (pares) ──────────────────────────────────────────
pair_rows <- list()
k <- 1
for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    pair_rows[[k]] <- data.table(
      pathway_A = ids[i], pathway_B = ids[j],
      name_A = panel$pathway_name[match(ids[i], panel$database_id)],
      name_B = panel$pathway_name[match(ids[j], panel$database_id)],
      overlap_genes = overlap_mat[i, j],
      jaccard = jaccard_mat[i, j]
    )
    k <- k + 1
  }
}
pairs <- rbindlist(pair_rows)
pairs <- pairs[order(-jaccard)]
fwrite_tsv(pairs, file.path(DIR_PATH, "PATHWAY_REDUNDANCY.tsv"))

# ── 4. Clustering hierárquico das vias por Jaccard ─────────────────────────────
jd <- as.dist(1 - jaccard_mat)
diag(jaccard_mat) <- 1
hc <- hclust(as.dist(1 - jaccard_mat), method = "average")

# módulos por corte de árvore (altura 0.85 → grupos com Jaccard >= 0.15)
modules <- cutree(hc, h = 0.85)
module_df <- data.table(database_id = ids, module = modules,
                        pathway_name = panel$pathway_name[match(ids, panel$database_id)],
                        panel_group = panel$panel_group[match(ids, panel$database_id)])
module_df <- module_df[order(module, database_id)]
fwrite_tsv(module_df, file.path(DIR_PATH, "PATHWAY_MODULES.tsv"))
log_msg("Nº de módulos de redundância (corte h=0.85):", max(modules))

# ── 5. Figura: heatmap de Jaccard com clustering ──────────────────────────────
ann <- data.frame(panel = panel$panel_group[match(ids, panel$database_id)],
                  row.names = ids)
ann_colors <- list(panel = c(ORIGINAL_10 = "#7570b3", ADDITIONAL_20 = "#1b9e77"))
png(file.path(DIR_FIG, "PATHWAY_REDUNDANCY_heatmap.png"), width = 2200, height = 2000, res = 200)
pheatmap(jaccard_mat, clustering_distance_rows = as.dist(1 - jaccard_mat),
         clustering_distance_cols = as.dist(1 - jaccard_mat),
         annotation_row = ann, annotation_col = ann, annotation_colors = ann_colors,
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         main = "Similaridade de Jaccard entre as 30 vias do painel",
         fontsize_row = 7, fontsize_col = 7)
dev.off()

log_msg("══ Redundância concluída ══")
cat("\n=== Top 15 pares de vias mais redundantes ===\n")
print(head(pairs[, .(name_A, name_B, overlap_genes, jaccard)], 15))
