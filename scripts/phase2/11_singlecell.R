# ═══════════════════════════════════════════════════════════════════════════════
# 11_singlecell.R — Localização celular (GSE232237), independente da análise bulk
#
# Quantifica expressão por compartimento para ITGA2, FN1, CCND1 e genes
# relevantes, diferenciando epitelial/tumoral de fibroblasto, endotélio, imune e
# SMC. A classificação de GSE232237 é marker-based (sem anotação autoral) —
# NÃO é tratada como verdade absoluta.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
})

log_msg("══ Localização celular (single-cell, GSE232237) ══")

sc <- fread(file.path(DIR_EXTERNAL, "GSE232237", "GSE232237_marker_celltype_results.tsv"))
setnames(sc, c("sample","subtype","cell_type","gene","mean_expr","pct_positive","n_cells"))

# ── Agregação por gene × célula (média entre amostras) ─────────────────────────
sc_agg <- sc[, .(
  mean_expr_across_samples = mean(mean_expr),
  mean_pct_positive = mean(pct_positive),
  n_samples = .N
), by = .(gene, cell_type)]

# ── Compartimento predominante por gene ────────────────────────────────────────
sc_predom <- sc_agg[, .SD[which.max(mean_expr_across_samples)], by = gene]
sc_predom[, predominant_compartment := cell_type]
sc_predom[, note := "marker-based classification (GSE232237); not authorial annotation"]

fwrite_tsv(sc_agg, file.path(DIR_VAL, "singlecell_expression_by_compartment.tsv"))
fwrite_tsv(sc_predom, file.path(DIR_VAL, "singlecell_predominant_compartment.tsv"))

log_msg("══ Single-cell concluído ══")
cat("\n=== Expressão por compartimento (média entre amostras) ===\n")
print(dcast(sc_agg[gene %in% c("ITGA2","FN1","CCND1","ITGB1")],
            gene ~ cell_type, value.var = "mean_expr_across_samples"))
cat("\n=== Compartimento predominante ===\n")
print(sc_predom[, .(gene, cell_type, mean_expr_across_samples, mean_pct_positive)])
