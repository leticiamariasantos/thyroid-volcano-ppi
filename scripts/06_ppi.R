#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 06_ppi.R — FASE 2: PPI (STRING) para genes DEG das vias selecionadas
# thyroid-volcano-ppi
#
# Conjunto de entrada (pré-definido): genes que são (i) DEG (FDR<0.05, |logFC|>1)
# E (ii) pertencem a pelo menos uma das 7 vias não-redundantes selecionadas.
# STRING v12.0 REST, score combinado >= 700. Centralidade + comunidades (igraph).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(httr)
  library(jsonlite)
  library(igraph)
})

PROJECT_ROOT <- here::here()
dir_ppi <- file.path(PROJECT_ROOT, "08_ppi")
dir.create(dir_ppi, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — PPI (STRING) ══\n")
cat("R:", as.character(getRversion()), "| igraph:", as.character(packageVersion("igraph")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

TAXON <- 9606
SCORE <- 700

# ── 1. Conjunto de entrada ─────────────────────────────────────────────────────
sel <- data.table::fread(file.path(PROJECT_ROOT, "07_pathway_redundancy", "pathways_selected.tsv"),
                         colClasses = c(id = "character"))
deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
deg_sig <- deg[regulation %in% c("Up", "Down"), ]

kegg <- readRDS(file.path(PROJECT_ROOT, "07_pathway_redundancy", "kegg_pathways.rds"))
sel_genes <- unique(unlist(kegg[sel$id]))
input_genes <- intersect(sel_genes, deg_sig$gene_symbol)
cat(sprintf("  Genes de entrada (DEG ∩ vias selecionadas): %d\n", length(input_genes)))

# ── 2. STRING mapping ─────────────────────────────────────────────────────────
cat("  Mapeando IDs STRING...\n")
map_resp <- httr::POST("https://string-db.org/api/json/get_string_ids",
  body = list(identifiers = paste(input_genes, collapse = "\r\n"),
              species = as.character(TAXON), limit = "1",
              caller_identity = "thyroid_global_ppi"),
  encode = "form", httr::timeout(120))
httr::stop_for_status(map_resp)
mapped <- jsonlite::fromJSON(httr::content(map_resp, as = "text", encoding = "UTF-8"), flatten = TRUE)
mapped <- data.table(mapped)[, .(gene_symbol = preferredName, STRING_id = stringId)]
mapped <- mapped[!is.na(STRING_id) & STRING_id != ""]
mapped <- merge(mapped, deg_sig[, .(gene_symbol, logFC, adj.P.Val, regulation)], by = "gene_symbol", all.x = TRUE)
cat(sprintf("  Mapeados: %d / %d\n", nrow(mapped), length(input_genes)))

# ── 3. STRING interações ──────────────────────────────────────────────────────
cat("  Buscando interações (score>=700)...\n")
net_resp <- httr::POST("https://string-db.org/api/tsv/network",
  body = list(identifiers = paste(mapped$STRING_id, collapse = "\r\n"),
              species = as.character(TAXON), required_score = "0",
              caller_identity = "thyroid_global_ppi"),
  encode = "form", httr::timeout(180))
httr::stop_for_status(net_resp)
raw <- data.table::fread(text = httr::content(net_resp, as = "text", encoding = "UTF-8"),
                         header = TRUE, colClasses = "character")
if (nrow(raw) == 0) stop("Sem interações retornadas.")
edges <- raw[, .(from = stringId_A, to = stringId_B, score = as.numeric(score) * 1000)]
edges <- edges[score >= SCORE & from %in% mapped$STRING_id & to %in% mapped$STRING_id]
cat(sprintf("  Arestas (>=%d): %d\n", SCORE, nrow(edges)))

# ── 4. Grafo + componente gigante ─────────────────────────────────────────────
g <- graph_from_data_frame(edges[, .(from, to, weight = score)], directed = FALSE,
                           vertices = mapped[, .(name = STRING_id, gene_symbol, logFC, adj.P.Val, regulation)])
comp <- components(g)
gcc <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
gcc <- delete_vertices(gcc, V(gcc)[degree(gcc) == 0])
cat(sprintf("  Rede: %d nós | %d arestas\n", vcount(gcc), ecount(gcc)))

# ── 5. Centralidade + comunidades ─────────────────────────────────────────────
cen <- data.frame(
  gene_symbol = V(gcc)$gene_symbol,
  logFC = V(gcc)$logFC,
  regulation = V(gcc)$regulation,
  degree = degree(gcc),
  betweenness = round(betweenness(gcc, normalized = TRUE), 4),
  closeness = round(closeness(gcc, normalized = TRUE), 6),
  stringsAsFactors = FALSE
)
cen <- cen[order(-cen$betweenness, -cen$degree), ]
set.seed(42)
wc <- if (ecount(gcc) >= 2) cluster_walktrap(gcc) else NULL
if (!is.null(wc)) V(gcc)$community <- wc$membership else V(gcc)$community <- 1L
cen$community <- V(gcc)$community[match(cen$gene_symbol, V(gcc)$gene_symbol)]

cat(sprintf("  Comunidades: %d | Modularidade: %.4f\n",
            length(unique(cen$community)), if (!is.null(wc)) modularity(wc) else NA))
cat(sprintf("  Densidade: %.4f | Clustering: %.4f\n",
            edge_density(gcc), transitivity(gcc, type = "global")))

# ── 6. Salvar ─────────────────────────────────────────────────────────────────
data.table::fwrite(mapped, file.path(dir_ppi, "PPI_string_mapping.tsv"), sep = "\t")
data.table::fwrite(edges, file.path(dir_ppi, "PPI_edges.tsv"), sep = "\t")
data.table::fwrite(cen, file.path(dir_ppi, "PPI_centrality.tsv"), sep = "\t")
data.table::fwrite(data.table::data.table(
  metric = c("input_genes","mapped","edges","nodes","communities","density",
             "clustering","STRING_score","STRING_version"),
  value = c(length(input_genes), nrow(mapped), nrow(edges), vcount(gcc),
            length(unique(cen$community)), round(edge_density(gcc), 5),
            round(transitivity(gcc, type = "global"), 5), SCORE, "12.0")),
  file.path(dir_ppi, "PPI_summary.tsv"), sep = "\t")

cat("\n  Top 15 por betweenness:\n")
print(head(cen[, c("gene_symbol","regulation","logFC","degree","betweenness","community")], 15), row.names = FALSE)
cat("══ FASE 2 — PPI CONCLUÍDO ══\n")
