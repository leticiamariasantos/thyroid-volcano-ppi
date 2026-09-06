# ═══════════════════════════════════════════════════════════════════════════════
# 20_ppi_audit.R — PARTE 5: auditoria quantitativa da estabilidade do PPI
#
# Compara 4 redes (subgrafos induzidos da rede STRING discovery):
#   1. original (177 genes -> 148 nós)
#   2. count-based  (genes DEG em voom E DESeq2)
#   3. no-composition (exclui marcadores musculares)
#   4. robust candidates (não-musculares ROBUSTOS nos 3 métodos)
# Métricas: nós, arestas, Jaccard nós/arestas, correlação de grau, hubs, comunidades.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table); library(igraph)})

cat("══ PARTE 5 — Auditoria do PPI ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_p <- "results/ppi_audit"
dir.create(dir_p, recursive = TRUE, showWarnings = FALSE)

# ── carrega ───────────────────────────────────────────────────────────────────
edges <- fread("08_ppi/PPI_edges.tsv")
map   <- fread("08_ppi/PPI_string_mapping.tsv")
cen   <- fread("08_ppi/PPI_centrality.tsv")
# mapping ENSP -> symbol
id2sym <- setNames(map$gene_symbol, map$STRING_id)

# gene-level graph
e <- data.table(from = id2sym[edges$from], to = id2sym[edges$to], score = edges$score)
e <- e[!is.na(from) & !is.na(to)]
g0 <- graph_from_data_frame(e[, .(from, to)], directed = FALSE)

# muscle markers
muscle <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
            "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
            "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
            "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")

# genes DEG em voom E DESeq2
voo <- fread("results/counts_deg/voom_full_results.tsv")
dds <- fread("results/counts_deg/deseq2_full_results.tsv")
deg_v <- voo[regulation %in% c("Up","Down"), gene_symbol]
deg_d <- dds[regulation %in% c("Up","Down"), gene_symbol]
count_robust <- intersect(deg_v, deg_d)

# genes ROBUSTOS não-musculares (dos 3 métodos)
mc <- fread("results/counts_deg/method_concordance.tsv")
robust_genes <- mc[robustness == "ROBUSTO" & !(gene %in% muscle), gene]

subgraph_metrics <- function(genes_keep, label) {
  genes_keep <- intersect(genes_keep, V(g0)$name)
  g <- induced_subgraph(g0, genes_keep)
  comps <- components(g)
  gcc <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  deg <- degree(g)
  hub_rank <- names(sort(deg, decreasing = TRUE))[1:min(5, length(deg))]
  comm <- tryCatch(cluster_louvain(gcc), error = function(e) NULL)
  list(label = label,
       nodes = vcount(g),
       edges = ecount(g),
       gcc_nodes = vcount(gcc),
       communities = if (!is.null(comm)) length(unique(membership(comm))) else NA_integer_,
       top_hubs = paste(hub_rank, collapse = ","))
}

variants <- list(
  original = subgraph_metrics(V(g0)$name, "original"),
  count_based = subgraph_metrics(count_robust, "count_based"),
  no_composition = subgraph_metrics(setdiff(V(g0)$name, muscle), "no_composition"),
  robust_candidates = subgraph_metrics(robust_genes, "robust_candidates")
)
out <- rbindlist(lapply(variants, as.data.table))
print(out)
data.table::fwrite(out, file.path(dir_p, "ppi_variants.tsv"), sep = "\t")

# ── Jaccard de nós e arestas vs original ──────────────────────────────────────
jaccard <- function(a, b) length(intersect(a,b)) / length(union(a,b))
cat("\n── Jaccard vs original ──\n")
V0 <- V(g0)$name
E0 <- paste(pmin(ends(g0, E(g0))[,1], ends(g0, E(g0))[,2]), pmax(ends(g0, E(g0))[,1], ends(g0, E(g0))[,2]))
for (nm in c("count_based","no_composition","robust_candidates")) {
  keep <- switch(nm,
    count_based = intersect(count_robust, V0),
    no_composition = setdiff(V0, muscle),
    robust_candidates = intersect(robust_genes, V0))
  g <- induced_subgraph(g0, keep)
  Vv <- V(g)$name
  Ev <- paste(pmin(ends(g, E(g))[,1], ends(g, E(g))[,2]), pmax(ends(g, E(g))[,1], ends(g, E(g))[,2]))
  cat(sprintf("  %-18s Jaccard nós=%.3f Jaccard arestas=%.3f\n", nm, jaccard(Vv, V0), jaccard(Ev, E0)))
}

# ── Correlação de grau (original vs variantes, nós compartilhados) ────────────
cat("\n── Correlação de grau (nós compartilhados) ──\n")
d0 <- degree(g0)
for (nm in c("count_based","no_composition","robust_candidates")) {
  keep <- switch(nm,
    count_based = intersect(count_robust, V0),
    no_composition = setdiff(V0, muscle),
    robust_candidates = intersect(robust_genes, V0))
  g <- induced_subgraph(g0, keep)
  shared <- intersect(V(g)$name, names(d0))
  if (length(shared) > 3) {
    cat(sprintf("  %-18s Spearman grau=%.3f\n", nm,
        cor(d0[shared], degree(g)[shared], method = "spearman")))
  }
}

# ── Preservação de hubs ───────────────────────────────────────────────────────
cat("\n── Preservação de hubs (top 5 por grau, por rede) ──\n")
for (v in variants) cat(sprintf("  %-18s: %s\n", v$label, v$top_hubs))

cat("\n══ Auditoria do PPI CONCLUÍDA ══\n")
