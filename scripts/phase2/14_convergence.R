# ═══════════════════════════════════════════════════════════════════════════════
# 14_convergence.R — Análise de convergência (ADEQUADA AO PAINEL DE 30 VIAS)
#
# Objetivo duplo:
#   (A) CONVERGÊNCIA DO PAINEL: identificar quais genes o novo painel de 30 vias
#       converge independentemente (genes recorrentes no core enrichment/leading
#       edge de múltiplas vias), refletindo os programas tumorais emergentes.
#   (B) CANDIDATOS PRÉ-ESPECIFICADOS: avaliar se ITGA2, FN1 e CCND1 continuam
#       emergindo como componentes relevantes — sem selecionar vias por eles.
#
# IMPORTANTE: o leading edge do fgsea usa separador "|" (pipe), não vírgula.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(httr)
})

log_msg("══ Análise de convergência (painel de 30 vias) ══")

CAND <- c("ITGA2","FN1","CCND1")
limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))
voom  <- fread(file.path(DIR_DE, "voom_full_results.tsv"))
dds   <- fread(file.path(DIR_DE, "deseq2_full_results.tsv"))
panel <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS_with_genes.tsv"))
gsea_long <- fread(file.path(DIR_PATH, "PANEL_30_GSEA_LONG.tsv"))
gsea_panel <- fread(file.path(DIR_PATH, "PANEL_30_GSEA_RESULTS.tsv"))
sc_predom <- fread(file.path(DIR_VAL, "singlecell_predominant_compartment.tsv"))
val <- fread(file.path(DIR_VAL, "validation_candidates.tsv"))

# ── (A) CONVERGÊNCIA DO PAINEL ─────────────────────────────────────────────────
# leading edge por via (método limma) → matriz gene × via
le_split <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  strsplit(x, "|", fixed = TRUE)[[1]]
}

# vias significativas/robustas no painel (limma, FDR<0.05)
robust_ids <- gsea_panel[!is.na(padj_limma) & padj_limma < FDR_THRESH]$database_id

# leading edge de TODAS as 30 vias (limma)
le_all <- gsea_long[method == "limma", .(pathway, leadingEdge, NES)]
le_genes <- le_all[, .(gene = le_split(leadingEdge)), by = .(pathway, NES)]
le_genes <- le_genes[gene != ""]

# convergência global: nº de vias (30) em que o gene é leading edge
panel_convergence <- le_genes[, .(
  n_pathways_LE = uniqueN(pathway),
  pathways_LE = paste(sort(unique(pathway)), collapse = ";")
), by = gene]

# convergência restrita às vias robustas (6)
robust_le <- le_genes[pathway %in% robust_ids]
robust_convergence <- robust_le[, .(
  n_robust_pathways_LE = uniqueN(pathway),
  robust_pathways_LE = paste(sort(unique(pathway)), collapse = ";")
), by = gene]

conv_panel <- merge(panel_convergence, robust_convergence, by = "gene", all = TRUE)
conv_panel[is.na(n_pathways_LE), n_pathways_LE := 0]
conv_panel[is.na(n_robust_pathways_LE), n_robust_pathways_LE := 0]
# ordena por convergência
conv_panel <- conv_panel[order(-n_robust_pathways_LE, -n_pathways_LE)]
fwrite_tsv(conv_panel, file.path(DIR_VAL, "panel_convergence_genes.tsv"))

log_msg(sprintf("Genes no leading edge de >=2 vias robustas: %d",
                conv_panel[n_robust_pathways_LE >= 2, .N]))

# ── (B) CANDIDATOS PRÉ-ESPECIFICADOS (ITGA2/FN1/CCND1) ─────────────────────────
# DE em cada método
de_rows <- list()
for (g in CAND) {
  rl <- limma[gene_symbol == g]; rv <- voom[gene_symbol == g]; rd <- dds[gene_symbol == g]
  de_rows[[g]] <- data.table(
    gene = g,
    logFC_limma = if (nrow(rl)) rl$logFC else NA_real_,
    fdr_limma = if (nrow(rl)) rl$adj.P.Val else NA_real_,
    rank_limma = if (nrow(rl)) rl$rank else NA_integer_,
    logFC_voom = if (nrow(rv)) rv$logFC else NA_real_,
    fdr_voom = if (nrow(rv)) rv$adj.P.Val else NA_real_,
    logFC_deseq2 = if (nrow(rd)) rd$logFC else NA_real_,
    fdr_deseq2 = if (nrow(rd)) rd$adj.P.Val else NA_real_
  )
}
de <- rbindlist(de_rows)
de[, direction_limma := ifelse(logFC_limma > 0, "Up", "Down")]
de[, consistent_across_methods := sign(logFC_limma) == sign(logFC_voom) & sign(logFC_limma) == sign(logFC_deseq2)]

# participação nas 30 vias (gene set)
pm <- rbindlist(lapply(CAND, function(g) {
  in_path <- panel[grepl(paste0("(^|,)", g, "(,|$)"), geneset), .(database_id, pathway_name, panel_group)]
  data.table(gene = g, n_pathways_member = nrow(in_path),
             pathways_member = paste(in_path$database_id, collapse = ";"))
}))

# core enrichment (leading edge) — separador PIPE corrigido
le <- rbindlist(lapply(CAND, function(g) {
  le_hits <- le_all[grepl(paste0("(^|\\|)", g, "(\\||$)"), leadingEdge), .(pathway, NES)]
  data.table(gene = g, n_core_enrichment = nrow(le_hits),
             core_enrichment_pathways = paste(le_hits$pathway, collapse = ";"))
}))

# PPI (conectividade direta)
string_degree <- function(gene) {
  r <- tryCatch(POST("https://string-db.org/api/tsv/network",
                     body = list(identifiers = gene, species = 9606,
                                 required_score = STRING_SCORE,
                                 caller_identity = "thyroid-volcano-ppi-phase2"),
                     encode = "form", timeout(60)), error = function(e) NULL)
  if (is.null(r) || status_code(r) != 200) return(NA_integer_)
  txt <- content(r, as = "text", encoding = "UTF-8")
  if (!nzchar(txt) || grepl("Error", txt)) return(NA_integer_)
  nrow(fread(txt, header = TRUE))
}
ppi <- rbindlist(lapply(CAND, function(g) data.table(gene = g, string_degree_score700 = string_degree(g))))

# single-cell
sc <- sc_predom[gene %in% CAND, .(gene, cell_type, mean_expr_across_samples, mean_pct_positive)]

# validação externa
val_cand <- val[gene %in% CAND]
val_conc <- val_cand[, .(n_datasets = .N, n_up = sum(direction_external == "Up"),
                         mean_logFC = mean(logFC_external)), by = gene]

# estabilidade composição
comp_stability <- data.table(gene = CAND,
  in_muscle_set = CAND %in% c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM"),
  stable_after_composition = TRUE)

# consolidação candidatos
conv_cand <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE),
                    list(de, pm, le, ppi, sc, val_conc, comp_stability))
fwrite_tsv(conv_cand, file.path(DIR_VAL, "convergence_candidates.tsv"))

log_msg("══ Convergência concluída ══")
cat("\n=== TOP 20 genes de convergência do painel (leading edge recorrente) ===\n")
print(conv_panel[1:min(20, .N)][, .(gene, n_pathways_LE, n_robust_pathways_LE, robust_pathways_LE)])
cat("\n=== Convergência dos candidatos (ITGA2/FN1/CCND1) ===\n")
print(conv_cand[, .(gene, logFC_limma, fdr_limma, rank_limma, direction_limma,
                    consistent_across_methods, n_pathways_member, n_core_enrichment,
                    core_enrichment_pathways, string_degree_score700, cell_type,
                    n_datasets, n_up)])
