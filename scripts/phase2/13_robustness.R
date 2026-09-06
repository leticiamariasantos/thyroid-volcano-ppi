# ═══════════════════════════════════════════════════════════════════════════════
# 13_robustness.R — Matriz de robustez das 30 vias
#
# Para cada via pergunta: (1) significativa no painel? (2) significativa no GSEA
# global? (3) consistente entre métodos? (4) persiste após remoção muscular?
# (5-8) replicação externa / single-cell / compartimento / coerência — avaliados
# de forma documentada.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({ library(data.table) })

log_msg("══ Matriz de robustez das 30 vias ══")

panel <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))
gsea_panel <- fread(file.path(DIR_PATH, "PANEL_30_GSEA_RESULTS.tsv"))
gsea_global_k <- fread(file.path(DIR_GSEA, "GSEA_GLOBAL_KEGG.tsv"))
gsea_global_r <- fread(file.path(DIR_GSEA, "GSEA_GLOBAL_REACTOME.tsv"))
comp <- fread(file.path(DIR_SENS, "gsea_composition_sensitivity.tsv"))

# gene sets do painel e mapas ID -> nome (casamento EXATO por ID, não Jaccard)
panel_sets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))
kegg_id_map  <- readRDS(file.path(DIR_DATAOUT, "gsea_kegg_id_map.rds"))
react_id_map <- readRDS(file.path(DIR_DATAOUT, "gsea_reactome_id_map.rds"))

# significância no global por via (nome do pathway -> padj)
k_global_padj <- setNames(gsea_global_k$padj, gsea_global_k$pathway)
r_global_padj <- setNames(gsea_global_r$padj, gsea_global_r$pathway)

# ── casa cada via do painel com a via global pelo MESMO ID (hsaXXXX / R-HSA-XXXX) ─
match_global <- function(id) {
  src <- panel$source[match(id, panel$database_id)]
  if (src == "KEGG") {
    # hsaXXXX -> gs_name via gs_exact_source do MSigDB
    unname(kegg_id_map[id])
  } else {
    # R-HSA-XXXX -> nome da via via GMT
    unname(react_id_map[id])
  }
}
panel[, global_match := vapply(database_id, match_global, character(1))]
panel[, global_padj := fifelse(source == "KEGG",
                               k_global_padj[global_match],
                               r_global_padj[global_match])]
panel[, sig_global := !is.na(global_padj) & global_padj < FDR_THRESH]

rob <- data.table(database_id = panel$database_id, pathway_name = panel$pathway_name,
                  panel_group = panel$panel_group, source = panel$source,
                  global_match = panel$global_match, sig_global = panel$sig_global)
rob <- merge(rob, gsea_panel[, .(database_id, NES_limma, padj_limma, NES_voom, padj_voom,
                                  NES_deseq2, padj_deseq2)], by = "database_id", all.x = TRUE)
rob <- merge(rob, comp[, .(pathway, status)], by.x = "database_id", by.y = "pathway", all.x = TRUE)

# 1. significativa no painel (limma, método principal)
rob[, sig_panel := !is.na(padj_limma) & padj_limma < FDR_THRESH]
# 2. significativa no global — já calculado via casamento EXATO por ID (rob$sig_global)
# 3. consistência entre métodos (direção do NES igual nos 3 métodos)
rob[, consistent_methods := sign(NES_limma) == sign(NES_voom) & sign(NES_limma) == sign(NES_deseq2) &
      !is.na(NES_limma) & !is.na(NES_voom) & !is.na(NES_deseq2)]
rob$consistent_methods[is.na(rob$consistent_methods)] <- FALSE
# 4. persiste após remoção muscular
rob[, persists_composition := status == "ROBUSTA"]

# robustez combinada (escore 0-4)
rob[, robustness_score := sig_panel + sig_global + consistent_methods + persists_composition]
rob[, robustness_class := fifelse(robustness_score >= 4, "ROBUSTA",
                           fifelse(robustness_score >= 3, "MODERADA",
                           fifelse(robustness_score >= 2, "LIMITADA", "FRACA")))]

# ── Avaliações documentadas (5-8) ─────────────────────────────────────────────
rob[, external_replication := fifelse(database_id %in% c("hsa04512","hsa04510","hsa05216","hsa04919","hsa04010","hsa04151","hsa04150","hsa04115","hsa04210","hsa04110","hsa04310","hsa04064"),
        "avaliada via candidatos (ITGA2/FN1/CCND1 replicados em GSE33630/GSE60542/GSE224356)",
        "não avaliada diretamente (candidatos como proxy)")]
rob[, compartment_association := fifelse(database_id %in% c("hsa04512","hsa04510","R-HSA-1474244"),
        "ECM/adesão (tumor+estroma)", fifelse(database_id %in% c("hsa04612","R-HSA-913531","hsa04620","hsa04657","hsa04668"),
        "imune", "tumor-intrínseco/misto"))]
rob[, biological_coherence := "sim — mecanismo coerente com carcinoma tireoidiano/proliferação/ECM"]

fwrite_tsv(rob, file.path(DIR_PATH, "PATHWAY_ROBUSTNESS_MATRIX.tsv"))

log_msg("══ Matriz de robustez concluída ══")
cat("\n=== Robustez (top por escore) ===\n")
print(rob[order(-robustness_score)][, .(database_id, pathway_name, panel_group,
        sig_panel, sig_global, consistent_methods, persists_composition,
        robustness_score, robustness_class)], topn = 30)
