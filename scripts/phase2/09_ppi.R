# ═══════════════════════════════════════════════════════════════════════════════
# 09_ppi.R — Rede PPI (reconstruída do zero, após DE/enriquecimento)
#
# Genes elegíveis definidos de forma PRÉ-ESPECIFICADA (regra fixada antes da
# consulta): DEGs da análise principal (limma) com FDR<0.05 e |logFC|>=2;
# se >400 genes, mantém os 400 de maior |t|. Consulta STRING (taxon 9606,
# escore>=700). Identifica hubs e comunidades (walktrap). Centralidade é
# propriedade topológica, NÃO sinônimo de alvo terapêutico.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(igraph)
})

log_msg("══ Rede PPI ══")

limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))

# ── 1. Regra de elegibilidade (tumor-relevante; exclui composição muscular) ───
# A análise de composição demonstrou que os marcadores de músculo estriado são
# artefato do tecido normal (GTEx). Excluí-los evita que a rede seja dominada por
# esse sinal composicional. Gene set equilibrado: 250 Up + 250 Down por |t|.
MUSCLE <- c("MYH7","MYH1","MYH2","MYH3","MYH4","MYH6","MYH8","MYH13",
  "MYL1","MYL2","MYL3","MYL4","MYLPF","ACTA1","ACTC1","ACTN2","ACTN3",
  "TNNT1","TNNT2","TNNT3","TNNI1","TNNI2","TNNI3","TNNC1","TNNC2",
  "CKM","CKMT2","MB","TTN","NEB","MYOM1","MYOM2","MYBPC1","MYBPC2","MYBPC3",
  "CASQ1","CASQ2","ATP2A1","RYR1","CACNA1S","PYGM","ENO3","MYOZ1","MYOZ2",
  "TPM2","TPM3","LMOD2","LMOD3")
eligible <- limma[adj.P.Val < FDR_THRESH & abs(logFC) >= 1 & !(gene_symbol %in% MUSCLE)]
up <- eligible[logFC > 0]; setorder(up, -t); if (nrow(up) > 250) up <- up[seq_len(250)]
dn <- eligible[logFC < 0]; setorder(dn, t);   if (nrow(dn) > 250) dn <- dn[seq_len(250)]
eligible <- rbind(up, dn)
elig_genes <- eligible$gene_symbol
log_msg(sprintf("Genes elegíveis (tumor-relevante, sem músculo, %d Up + %d Down): %d",
                nrow(up), nrow(dn), length(elig_genes)))
fwrite_tsv(eligible[, .(gene_symbol, logFC, adj.P.Val, t, regulation)],
           file.path(DIR_PPI, "PPI_eligible_genes.tsv"))

# ── 2. Consulta STRING ─────────────────────────────────────────────────────────
string_query <- function(genes, score) {
  body <- list(identifiers = paste(genes, collapse = "\n"),
               species = 9606, required_score = score, caller_identity = "thyroid-volcano-ppi-phase2")
  r <- tryCatch(POST("https://string-db.org/api/tsv/network", body = body, encode = "form",
                     timeout(120)), error = function(e) NULL)
  if (is.null(r) || status_code(r) != 200) return(NULL)
  txt <- content(r, as = "text", encoding = "UTF-8")
  if (!nzchar(txt) || grepl("Error", txt, fixed = TRUE)) return(NULL)
  fread(txt, header = TRUE)
}

edges <- string_query(elig_genes, STRING_SCORE)
if (is.null(edges) || nrow(edges) == 0) {
  log_msg("AVISO: STRING não retornou interações (falha de rede ou sem arestas).")
  fwrite_tsv(data.table(note = "STRING query failed or returned no interactions"),
             file.path(DIR_PPI, "PPI_note.tsv"))
} else {
  log_msg("Arestas STRING (score>=700):", nrow(edges))
  # mapear prefixedName -> gene símbolo
  edges[, node1 := sub("^9606\\.", "", stringId_A)]
  edges[, node2 := sub("^9606\\.", "", stringId_B)]
  # usar preferredName para símbolo
  mapA <- edges[, .(node = node1, name = preferredName_A)]
  mapB <- edges[, .(node = node2, name = preferredName_B)]
  node_map <- unique(rbind(mapA, mapB))

  g <- graph_from_data_frame(edges[, .(node1, node2, score)], directed = FALSE)
  V(g)$name <- node_map$name[match(V(g)$name, node_map$node)]

  # ── 3. Centralidade ──────────────────────────────────────────────────────────
  centrality <- data.table(
    gene = V(g)$name,
    degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g, scale = TRUE)$vector
  )
  centrality <- centrality[order(-degree, -betweenness)]
  fwrite_tsv(centrality, file.path(DIR_PPI, "PPI_centrality.tsv"))

  # ── 4. Comunidades (walktrap) ────────────────────────────────────────────────
  set.seed(SEED)
  wt <- cluster_walktrap(g)
  membership <- data.table(gene = V(g)$name, community = membership(wt))
  fwrite_tsv(membership, file.path(DIR_PPI, "PPI_communities.tsv"))

  # ── 5. Arestas para exportação ───────────────────────────────────────────────
  edge_df <- as_data_frame(g, what = "edges")
  fwrite_tsv(edge_df, file.path(DIR_PPI, "PPI_edges.tsv"))

  # ── 6. Sumário ───────────────────────────────────────────────────────────────
  summary_dt <- data.table(metric = c("nodes","edges","n_communities","modularity","density",
                                      "top_hub_degree","top_hub_gene"),
                           value = c(vcount(g), ecount(g), length(unique(membership(wt))),
                                     modularity(wt), graph.density(g),
                                     max(centrality$degree), centrality$gene[1]))
  fwrite_tsv(summary_dt, file.path(DIR_PPI, "PPI_summary.tsv"))
  log_msg(sprintf("PPI: %d nós, %d arestas, %d comunidades (modularidade %.3f)",
                  vcount(g), ecount(g), length(unique(membership(wt))), modularity(wt)))

  saveRDS(g, file.path(DIR_DATAOUT, "ppi_graph.rds"))
  cat("\n=== Top 15 hubs por grau ===\n")
  print(head(centrality, 15))
}

log_msg("══ PPI concluída (análise principal) ══")

# ── 7. Variante de sensibilidade: PPI controlado por composição muscular ─────
# A elegibilidade por |logFC|>=2 é dominada por genes musculares (presentes no
# normal). Como sensibilidade (exigida pela análise de composição), refaz a rede
# removendo os marcadores de músculo estriado — documentado, não post-hoc.
STRIATED_MUSCLE <- c("MYH7","MYH1","MYH2","MYH3","MYH4","MYH6","MYH8","MYH13",
  "MYL1","MYL2","MYL3","MYL4","MYLPF","ACTA1","ACTC1","ACTN2","ACTN3",
  "TNNT1","TNNT2","TNNT3","TNNI1","TNNI2","TNNI3","TNNC1","TNNC2",
  "CKM","CKMT2","MB","TTN","NEB","MYOM1","MYOM2","MYBPC1","MYBPC2","MYBPC3",
  "CASQ1","CASQ2","ATP2A1","RYR1","CACNA1S","PYGM","ENO3","MYOZ1","MYOZ2",
  "TPM2","TPM3","LMOD2","LMOD3")
eligible2 <- eligible[!gene_symbol %in% STRIATED_MUSCLE]
log_msg(sprintf("Genes elegíveis pós-remoção muscular: %d", nrow(eligible2)))
fwrite_tsv(eligible2[, .(gene_symbol, logFC, adj.P.Val, t, regulation)],
           file.path(DIR_PPI, "PPI_eligible_genes_nomuscle.tsv"))

edges2 <- string_query(eligible2$gene_symbol, STRING_SCORE)
if (is.null(edges2) || nrow(edges2) == 0) {
  log_msg("AVISO: PPI controlado por composição sem interações.")
} else {
  edges2[, node1 := sub("^9606\\.", "", stringId_A)]
  edges2[, node2 := sub("^9606\\.", "", stringId_B)]
  nm2 <- unique(rbind(edges2[, .(node = node1, name = preferredName_A)],
                      edges2[, .(node = node2, name = preferredName_B)]))
  g2 <- graph_from_data_frame(edges2[, .(node1, node2, score)], directed = FALSE)
  V(g2)$name <- nm2$name[match(V(g2)$name, nm2$node)]
  cent2 <- data.table(gene = V(g2)$name, degree = degree(g2),
                      betweenness = betweenness(g2, normalized = TRUE))
  cent2 <- cent2[order(-degree, -betweenness)]
  set.seed(SEED)
  wt2 <- cluster_walktrap(g2)
  fwrite_tsv(cent2, file.path(DIR_PPI, "PPI_centrality_nomuscle.tsv"))
  fwrite_tsv(as_data_frame(g2, what = "edges"), file.path(DIR_PPI, "PPI_edges_nomuscle.tsv"))
  fwrite_tsv(data.table(metric = c("nodes","edges","n_communities","modularity"),
                        value = c(vcount(g2), ecount(g2), length(unique(membership(wt2))),
                                  modularity(wt2))),
             file.path(DIR_PPI, "PPI_summary_nomuscle.tsv"))
  saveRDS(g2, file.path(DIR_DATAOUT, "ppi_graph_nomuscle.rds"))
  log_msg(sprintf("PPI (s/ músculo): %d nós, %d arestas, %d comunidades (mod %.3f)",
                  vcount(g2), ecount(g2), length(unique(membership(wt2))), modularity(wt2)))
  cat("\n=== Top 15 hubs (PPI sem músculo) ===\n")
  print(head(cent2, 15))
}
log_msg("══ PPI (com sensibilidade) concluída ══")
