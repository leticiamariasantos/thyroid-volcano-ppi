# ═══════════════════════════════════════════════════════════════════════════════
# 05_gsea_global.R — GSEA global (sem restrição ao painel de 30 vias)
#
# Pergunta respondida: "Quais programas biológicos emergem do transcriptoma sem
# restringir previamente o espaço de hipóteses?" Coleções: KEGG (CP:KEGG_LEGACY),
# Reactome (GMT oficial) e Hallmark (MSigDB H). Ranking: estatística t do limma
# (análise principal).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(msigdbr)
})

log_msg("══ GSEA global ══")

rank_t <- readRDS(file.path(DIR_DATAOUT, "limma_rank_t.rds"))
# rank_t já é ordenado decrescente; fgsea aceita vetor nomeado ordenado

run_fgsea <- function(genesets, rankvec, label) {
  gs <- lapply(genesets, function(g) unique(g[g %in% names(rankvec)]))
  gs <- gs[lengths(gs) >= 5]
  if (length(gs) == 0) return(data.table())
  set.seed(SEED)
  res <- fgseaMultilevel(pathways = gs, stats = rankvec,
                         minSize = 5, maxSize = 500, eps = 0, nPermSimple = 10000)
  res <- as.data.table(res)
  res$collection <- label
  res$pathway_name <- res$pathway
  res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ",")]
  res
}

# ── 1. KEGG (CP:KEGG_LEGACY via MSigDB) ───────────────────────────────────────
log_msg("Carregando gene sets KEGG (MSigDB CP:KEGG_LEGACY)...")
kegg_ms <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_sets <- split(kegg_ms$gene_symbol, kegg_ms$gs_name)
log_msg("  KEGG:", length(kegg_sets), "gene sets")
saveRDS(kegg_sets, file.path(DIR_DATAOUT, "gsea_kegg_genesets.rds"))
# mapa ID -> nome (para casamento exato na matriz de robustez)
kegg_id_map <- setNames(kegg_ms$gs_name, kegg_ms$gs_exact_source)
kegg_id_map <- kegg_id_map[!duplicated(names(kegg_id_map)) & !is.na(names(kegg_id_map))]
saveRDS(kegg_id_map, file.path(DIR_DATAOUT, "gsea_kegg_id_map.rds"))
gsea_kegg <- run_fgsea(kegg_sets, rank_t, "KEGG")
fwrite_tsv(gsea_kegg, file.path(DIR_GSEA, "GSEA_GLOBAL_KEGG.tsv"))

# ── 2. Reactome (GMT oficial) ─────────────────────────────────────────────────
log_msg("Carregando gene sets Reactome (GMT)...")
gmt <- readLines(file.path(DIR_EXTERNAL, "reactome", "ReactomePathways.gmt"))
react_sets <- lapply(gmt, function(line) {
  p <- strsplit(line, "\t")[[1]]
  if (length(p) < 3) return(NULL)
  setNames(list(unique(p[-c(1,2)])), p[1])
})
react_sets <- unlist(react_sets, recursive = FALSE)
log_msg("  Reactome:", length(react_sets), "gene sets")
saveRDS(react_sets, file.path(DIR_DATAOUT, "gsea_reactome_genesets.rds"))
# mapa ID -> nome (para casamento exato na matriz de robustez)
gmt_parts <- strsplit(gmt, "\t")
react_id_map <- setNames(vapply(gmt_parts, function(p) p[1], character(1)),
                         vapply(gmt_parts, function(p) p[2], character(1)))
react_id_map <- react_id_map[!duplicated(names(react_id_map)) & !is.na(names(react_id_map))]
saveRDS(react_id_map, file.path(DIR_DATAOUT, "gsea_reactome_id_map.rds"))
gsea_react <- run_fgsea(react_sets, rank_t, "Reactome")
fwrite_tsv(gsea_react, file.path(DIR_GSEA, "GSEA_GLOBAL_REACTOME.tsv"))

# ── 3. Hallmark (MSigDB H) ────────────────────────────────────────────────────
log_msg("Carregando gene sets Hallmark (MSigDB H)...")
hall_ms <- msigdbr(species = "Homo sapiens", collection = "H")
hall_sets <- split(hall_ms$gene_symbol, hall_ms$gs_name)
log_msg("  Hallmark:", length(hall_sets), "gene sets")
saveRDS(hall_sets, file.path(DIR_DATAOUT, "gsea_hallmark_genesets.rds"))
gsea_hall <- run_fgsea(hall_sets, rank_t, "Hallmark")
fwrite_tsv(gsea_hall, file.path(DIR_GSEA, "GSEA_GLOBAL_HALLMARK.tsv"))

# ── 4. Combinação e sumário ───────────────────────────────────────────────────
gsea_all <- rbindlist(list(gsea_kegg, gsea_react, gsea_hall), fill = TRUE)
fwrite_tsv(gsea_all, file.path(DIR_GSEA, "GSEA_GLOBAL_ALL.tsv"))

# top por coleção (significativas FDR<0.05)
sig <- gsea_all[padj < FDR_THRESH & !is.na(padj)]
log_msg(sprintf("GSEA global: %d gene sets testados, %d significativos (FDR<0.05)",
                nrow(gsea_all), nrow(sig)))

top_summary <- sig[, .(collection, pathway, size, NES, pval, padj)]
top_summary[, abs_NES := abs(NES)]
setorder(top_summary, collection, -abs_NES)
fwrite_tsv(top_summary, file.path(DIR_GSEA, "GSEA_GLOBAL_significant_summary.tsv"))

log_msg("══ GSEA global concluída ══")
cat("\n=== TOP por coleção (FDR<0.05) ===\n")
for (cc in c("KEGG", "Reactome", "Hallmark")) {
  sub <- top_summary[collection == cc]
  cat("\n--", cc, sprintf("(%d sig)", nrow(sub)), "--\n")
  print(head(sub[order(-abs(NES))], 15))
}
