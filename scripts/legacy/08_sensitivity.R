#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 08_sensitivity.R — FASE 2: análise de sensibilidade (thresholds e parâmetros)
# thyroid-volcano-ppi
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({ library(data.table); library(here); library(igraph) })
PROJECT_ROOT <- here::here()
dir_sens <- file.path(PROJECT_ROOT, "11_sensitivity")
dir.create(dir_sens, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — Sensibilidade ══\n")
res <- list()

# ── 1. Sensibilidade de DEG (limiares FDR e |logFC|) ──────────────────────────
deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
deg[, logFC_mag := abs(logFC)]
for (fdr in c(0.01, 0.05, 0.10)) {
  for (lfc in c(0.5, 1, 1.5)) {
    n <- deg[adj.P.Val < fdr & logFC_mag > lfc, .N]
    res[[paste0("DEG_FDR", fdr, "_LFC", lfc)]] <- n
  }
}
cat("  DEG counts (FDR × |logFC|):\n")
print(as.data.frame(res), row.names = FALSE)

# ── 2. Sensibilidade do score STRING (rede) ───────────────────────────────────
edges <- data.table::fread(file.path(PROJECT_ROOT, "08_ppi", "PPI_edges.tsv"))
mapped <- data.table::fread(file.path(PROJECT_ROOT, "08_ppi", "PPI_string_mapping.tsv"))
for (sc in c(400, 700, 900)) {
  e <- edges[score >= sc]
  if (nrow(e) > 0) {
    g <- graph_from_data_frame(e[, .(from, to)], directed = FALSE)
    cc <- components(g)
    gc <- induced_subgraph(g, V(g)[cc$membership == which.max(cc$csize)])
    gc <- delete_vertices(gc, V(gc)[degree(gc) == 0])
    cat(sprintf("  STRING>=%d: arestas=%d, nós(componente gigante)=%d\n", sc, nrow(e), vcount(gc)))
  } else cat(sprintf("  STRING>=%d: 0 arestas\n", sc))
}

# ── 3. Correlação entre rankings alternativos (t vs signed -log10 p) ──────────
r1 <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "ranking_t.rnk"), header = FALSE)
r2 <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "ranking_signed_log10p.rnk"), header = FALSE)
m <- merge(r1, r2, by = "V1")
cat(sprintf("  Spearman(rating t, signed -log10p) = %.4f\n", cor(m$V2.x, m$V2.y, method = "spearman")))

cat("══ Sensibilidade CONCLUÍDA ══\n")
