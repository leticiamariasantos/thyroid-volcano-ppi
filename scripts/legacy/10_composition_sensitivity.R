#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 10_composition_sensitivity.R — FASE 3: sensibilidade à composição
# thyroid-volcano-ppi
#
# A. re-ranking dos DEGs excluindo marcadores musculares esqueléticos;
# B. re-GSEA (fgsea) excluindo marcadores musculares do ranking;
# C. comparação antes/depois (vias preservadas/perdidas; status FN1/ITGA2).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table); library(here); library(fgsea)
})

PROJECT_ROOT <- here::here()
dir_sens <- file.path(PROJECT_ROOT, "results", "sensitivity")
dir.create(dir_sens, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 3 — Sensibilidade à composição ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
kegg <- readRDS(file.path(PROJECT_ROOT, "07_pathway_redundancy", "kegg_pathways.rds"))

muscle <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
            "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
            "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
            "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")

# ── A. Re-ranking sem marcadores musculares ─────────────────────────────────────
deg_nomuscle <- deg[!(gene_symbol %in% muscle), ]
deg_sig_nm <- deg_nomuscle[regulation != "NS"]
top_nm <- deg_sig_nm[order(-abs(logFC))]
cat("A) Top 25 DEGs APÓS remover marcadores musculares:\n")
print(head(top_nm[, .(gene_symbol, logFC, adj.P.Val, regulation)], 25), row.names = FALSE)

# status FN1 / ITGA2
for (g in c("FN1","ITGA2","CTSS","CCND1","HLA-DPA1","LAMA2","HSPG2","ITGA2B")) {
  r <- deg[gene_symbol == g]
  if (nrow(r)) cat(sprintf("  %-8s logFC=%+.3f FDR=%.1e (%s)\n",
                           g, r$logFC, r$adj.P.Val, r$regulation))
}

# ── B. Re-GSEA sem marcadores musculares ───────────────────────────────────────
rnk <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "ranking_t.rnk"), header = FALSE)
ranking <- setNames(rnk$V2, rnk$V1)
ranking_nm <- ranking[!(names(ranking) %in% muscle)]
ranking_nm <- ranking_nm[order(ranking_nm, decreasing = TRUE)]
cat(sprintf("\nB) Ranking original: %d genes | sem músculo: %d genes\n",
            length(ranking), length(ranking_nm)))

set.seed(42)
gsea_orig <- fgsea(pathways = kegg, stats = ranking, minSize = 15, maxSize = 500, nPermSimple = 10000)
gsea_nm   <- fgsea(pathways = kegg, stats = ranking_nm, minSize = 15, maxSize = 500, nPermSimple = 10000)

merge_g <- merge(gsea_orig[, .(pathway, NES_orig = NES, padj_orig = padj)],
                 gsea_nm[, .(pathway, NES_nm = NES, padj_nm = padj)], by = "pathway")
merge_g[, pathway_short := sub("__.*$", "", pathway)]

# vias significativas antes vs depois
sig_orig <- merge_g[padj_orig < 0.05, pathway_short]
sig_nm   <- merge_g[padj_nm < 0.05, pathway_short]
lost <- setdiff(sig_orig, sig_nm)
kept <- intersect(sig_orig, sig_nm)
gained <- setdiff(sig_nm, sig_orig)
cat(sprintf("  Vias significativas: antes=%d | depois=%d\n", length(sig_orig), length(sig_nm)))
cat(sprintf("  PRESERVADAS: %d | PERDIDAS: %d | NOVAS: %d\n", length(kept), length(lost), length(gained)))
cat("  Perdidas após remover músculo:\n"); print(lost)
cat("  Preservadas:\n"); print(kept)

out <- merge_g[, .(pathway = pathway_short, NES_orig, padj_orig, NES_nm, padj_nm)]
out <- out[order(padj_orig)]
data.table::fwrite(out, file.path(dir_sens, "gsea_composition_sensitivity.tsv"), sep = "\t")
data.table::fwrite(top_nm, file.path(dir_sens, "DEG_no_muscle_top.tsv"), sep = "\t")

# painel a priori: NES antes/depois
A_PRIORI <- c("05216","04919","04010","04151","04150","04115","04210","04110","04310","04064")
cat("\n  Painel a priori (NES antes → depois):\n")
for (h in A_PRIORI) {
  r <- out[pathway == h]
  if (nrow(r)) cat(sprintf("  hsa%s: NES %.2f -> %.2f (padj %.3g -> %.3g)\n",
                           h, r$NES_orig, r$NES_nm, r$padj_orig, r$padj_nm))
  else cat(sprintf("  hsa%s: NA\n", h))
}

cat("\n══ FASE 3 — Sensibilidade à composição CONCLUÍDA ══\n")
