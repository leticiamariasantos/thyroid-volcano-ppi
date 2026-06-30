# ═══════════════════════════════════════════════════════════════════════════════
# R/05_ppi.R — PPI Network (Nature Communications / Frontiers standard)
# thyroid-volcano-ppi — MAIN OUTPUT 2/2
#
# Style: Professional Fruchterman-Reingold layout,
#        functional community colors, edge thickness ∝ STRING score,
#        hub genes highlighted, labels only for key nodes,
#        excellent spatial separation, minimal edge crossings.
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── M3: PPI Network ──\n")

if (length(deg_genes) < 2) stop("Need ≥2 DEGs for PPI network.")
cat(sprintf("  DEGs: %d (↑%d ↓%d)\n", length(deg_genes), length(up_genes), length(down_genes)))

nodes <- tibble(
  gene_symbol = deg_genes,
  regulation  = ifelse(gene_symbol %in% up_genes, "Up",
                       ifelse(gene_symbol %in% down_genes, "Down", "NS"))
)

# ── STRING mapping ────────────────────────────────────────────────────────────
cat("  STRING mapping...\n")
mapped <- map_string_ids(nodes, STRING_TAXON)
cat(sprintf("  Mapped: %d/%d\n", nrow(mapped), nrow(nodes)))
if (nrow(mapped) < 2) stop("Need ≥2 mapped proteins.")
mapped_ids <- mapped$STRING_id

# ── Fetch interactions ────────────────────────────────────────────────────────
cat("  Fetching edges...\n")
edges <- fetch_string_edges(mapped_ids, STRING_TAXON) |>
  dplyr::filter(from %in% mapped_ids, to %in% mapped_ids,
                combined_score >= THRESHOLD$string)
cat(sprintf("  Edges: %d (≥%d)\n", nrow(edges), THRESHOLD$string))
if (nrow(edges) == 0) stop("No edges. Lower THRESHOLD$string to 400.")

# ── Build igraph ──────────────────────────────────────────────────────────────
g <- build_igraph(edges, mapped)
comp  <- igraph::components(g)
giant <- which.max(comp$csize)
g_cc  <- igraph::induced_subgraph(g, igraph::V(g)[comp$membership == giant])
g_cc  <- igraph::delete_vertices(g_cc, igraph::V(g_cc)[igraph::degree(g_cc) == 0])

n_nodes <- igraph::vcount(g_cc)
n_edges <- igraph::ecount(g_cc)
cat(sprintf("  Network: %d nodes, %d edges\n", n_nodes, n_edges))

# ── Centrality metrics ────────────────────────────────────────────────────────
cent <- compute_centrality(g_cc)
# Hub criterion: betweenness >= 75th percentile of non-zero values AND degree >= 2,
# OR degree >= 4 (highly connected nodes in small networks may have low betweenness)
hub_th_bet <- quantile(cent$betweenness[cent$betweenness > 0], 0.75, na.rm = TRUE)
if (is.na(hub_th_bet) || hub_th_bet <= 0) hub_th_bet <- 0.01
cent$is_hub <- (cent$betweenness >= hub_th_bet & cent$degree >= 2) | cent$degree >= 4
cent$is_hub[is.na(cent$is_hub)] <- FALSE
cat(sprintf("  Hubs: %d\n", sum(cent$is_hub)))

# ── Community detection (walktrap) ────────────────────────────────────────────
set.seed(42)
if (n_edges >= 2 && n_nodes >= 3) {
  wc <- igraph::cluster_walktrap(g_cc, steps = 4)
  V(g_cc)$community <- wc$membership
} else {
  V(g_cc)$community <- 1L
}
n_comm <- length(unique(V(g_cc)$community))
cat(sprintf("  Communities: %d (modularity %.4f)\n",
            n_comm, if (exists("wc")) igraph::modularity(wc) else NA_real_))

# ── Community color palette ───────────────────────────────────────────────────
COMM_PALETTE <- c("#4477AA", "#AA4488", "#44AA77", "#DDCC77", "#88CCEE", "#CC6677")
comm_colors <- setNames(COMM_PALETTE[seq_len(n_comm)], seq_len(n_comm))
V(g_cc)$comm_fill <- comm_colors[as.character(V(g_cc)$community)]

# ── Fruchterman-Reingold layout ───────────────────────────────────────────────
cat("  Layout: Fruchterman-Reingold...\n")
set.seed(42)
lay <- tryCatch(
  igraph::layout_with_fr(g_cc, niter = 3000, area = n_nodes^2.5,
                         repulserad = n_nodes^2.5 * 0.15,
                         weights = igraph::E(g_cc)$weight / 1000),
  error = function(e) igraph::layout_with_fr(g_cc, niter = 1000)
)

if (any(is.na(lay)) || nrow(lay) != n_nodes) {
  set.seed(42)
  lay <- tryCatch(
    igraph::layout_with_graphopt(g_cc, charge = 0.005),
    error = function(e) igraph::layout_nicely(g_cc)
  )
}
if (any(is.na(lay)) || nrow(lay) != n_nodes) {
  lay <- igraph::layout_nicely(g_cc)
}

lay <- scale(lay, scale = TRUE)

V(g_cc)$x       <- lay[, 1]
V(g_cc)$y       <- lay[, 2]
V(g_cc)$nodedeg <- igraph::degree(g_cc)

# ── Hub and label settings ────────────────────────────────────────────────────
hub_ppi    <- cent |> dplyr::filter(is_hub) |> dplyr::pull(gene_symbol)
kegg_net   <- intersect(V(g_cc)$gene_symbol, kegg_list)
V(g_cc)$show_label <- V(g_cc)$gene_symbol %in% union(hub_ppi, kegg_net)

V(g_cc)$is_hub      <- V(g_cc)$gene_symbol %in% hub_ppi
V(g_cc)$node_border <- ifelse(V(g_cc)$is_hub, CELL_COLORS$highlight, "grey55")
V(g_cc)$bstroke     <- ifelse(V(g_cc)$is_hub, 0.55, 0.18)
V(g_cc)$node_size   <- ifelse(V(g_cc)$is_hub,
  scales::rescale(V(g_cc)$nodedeg, to = c(4.5, 10)),
  scales::rescale(V(g_cc)$nodedeg, to = c(3.0, 6.5)))

# ── Edge attributes ───────────────────────────────────────────────────────────
E(g_cc)$score_norm <- igraph::E(g_cc)$weight / 1000
E(g_cc)$is_intra   <- igraph::ends(g_cc, igraph::E(g_cc)) |>
  apply(1, function(e) {
    V(g_cc)[e[1]]$community == V(g_cc)[e[2]]$community
  })
E(g_cc)$edge_color <- ifelse(E(g_cc)$is_intra, "#777777", "#CCBBBB")
E(g_cc)$edge_alpha <- ifelse(E(g_cc)$is_intra, 0.30, 0.14)

cat(sprintf("  Intra-community edges: %d | Inter-community: %d\n",
            sum(E(g_cc)$is_intra), sum(!E(g_cc)$is_intra)))

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN PPI NETWORK
# ═══════════════════════════════════════════════════════════════════════════════

p_ppi <- ggraph(g_cc, layout = "manual", x = V(g_cc)$x, y = V(g_cc)$y) +

  # Edges: thickness ∝ STRING score
  geom_edge_link(
    aes(width = score_norm, color = edge_color, alpha = edge_alpha),
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.20, 1.4), guide = "none") +
  scale_edge_color_identity() +
  scale_edge_alpha_identity() +

  # Nodes: community-colored, degree-sized, hub-highlighted
  geom_node_point(
    aes(size = node_size, fill = comm_fill,
        color = node_border, stroke = bstroke),
    shape = 21
  ) +
  scale_size_continuous(range = c(2.5, 10), guide = "none") +
  scale_fill_identity() +
  scale_color_identity() +

  # Labels: large, italic, white halo
  ggrepel::geom_text_repel(
    data = . %>% dplyr::filter(show_label),
    aes(label = gene_symbol, x = x, y = y),
    size = 4.8, fontface = "italic", color = "#1A1A1A",
    bg.color = "white", bg.r = 0.16,
    segment.color = "grey50", segment.size = 0.25,
    segment.alpha = 0.6,
    box.padding = 0.50, point.padding = 0.35,
    max.overlaps = 60, force = 6, force_pull = 0.6,
    min.segment.length = 0.005,
    seed = 42
  ) +

  theme_void(base_size = FONT_SIZE) +
  theme(
    text             = element_text(family = FONT_FAM, color = "black"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(6, 6, 6, 6)
  )

# ── Inset legend (top-left) ───────────────────────────────────────────────────
xr <- range(lay[, 1]); yr <- range(lay[, 2])
lx <- xr[1] + diff(xr) * 0.03
ly <- yr[2] - diff(yr) * 0.025
gap_y <- diff(yr) * 0.05

n_up_net   <- sum(V(g_cc)$regulation == "Up")
n_down_net <- sum(V(g_cc)$regulation == "Down")
n_hub_net  <- sum(V(g_cc)$is_hub)

comm_counts <- table(V(g_cc)$community)
comm_labels <- sapply(names(comm_counts), function(ci) {
  paste0("Módulo ", ci, " (", comm_counts[ci], ")")
})

p_ppi <- p_ppi +
  # Module 1
  annotate("point", x = lx, y = ly, size = 4.5,
           color = "grey55", fill = comm_colors["1"], stroke = 0.30, shape = 21) +
  annotate("text",  x = lx + diff(xr) * 0.075, y = ly,
           label = comm_labels[1],
           size = 3.5, hjust = 0, color = "black", family = FONT_FAM) +
  # Module 2+
  annotate("point", x = lx, y = ly - gap_y, size = 4.5,
           color = "grey55", fill = comm_colors["2"], stroke = 0.30, shape = 21) +
  annotate("text",  x = lx + diff(xr) * 0.075, y = ly - gap_y,
           label = if (n_comm >= 2) comm_labels[2] else "",
           size = 3.5, hjust = 0, color = "black", family = FONT_FAM)

# Hub marker
hub_y_off <- if (n_comm >= 2) gap_y * 2.8 else gap_y * 1.5
p_ppi <- p_ppi +
  annotate("point", x = lx, y = ly - hub_y_off, size = 4.5,
           color = CELL_COLORS$highlight, fill = "white",
           stroke = 0.65, shape = 21) +
  annotate("text",  x = lx + diff(xr) * 0.075, y = ly - hub_y_off,
           label = paste0("Gene Hub (", n_hub_net, ")"),
           size = 3.5, hjust = 0, color = "black", family = FONT_FAM)

# ── Export (PNG 600 dpi + PDF) ──────────────────────────────────────────────
ggsave(here::here(DIRS$figures, "Fig2_PPI_Network_THCA_DEGs.png"),
       p_ppi, width = 180, height = 180, units = "mm", dpi = 600, bg = "white")
ggsave(here::here(DIRS$figures, "Fig2_PPI_Network_THCA_DEGs.pdf"),
       p_ppi, width = 180, height = 180, units = "mm", device = "pdf", bg = "white")
cat(sprintf("  PPI: %d nodes | %d edges | %d communities | %d hubs\n",
            n_nodes, n_edges, n_comm, sum(cent$is_hub)))
cat("  NOTE: Hubs are exploratory centrality metrics; they do not imply\n")
cat("  therapeutic targets, druggability, or causal mechanisms.\n")

# ── Tables ────────────────────────────────────────────────────────────────────
export_tsv(mapped, here::here(DIRS$network, "N01_string_mapping.tsv"))
export_tsv(edges,  here::here(DIRS$network, "N02_string_interactions.tsv"))
export_tsv(cent,   here::here(DIRS$network, "N03_centrality_metrics.tsv"))

ppi_summary <- tibble(
  metric = c("Input DEGs","Mapped STRING","Mapping rate",
             "Edges (≥threshold)","Nodes","Edges",
             "Hub proteins","Top hub","Top hub degree",
             "Top hub betweenness","Communities","Modularity",
             "Clustering coeff","Density",
             "STRING score","STRING version","Layout"),
  value = c(
    length(deg_genes), nrow(mapped),
    sprintf("%.1f%%", nrow(mapped)/length(deg_genes)*100),
    nrow(edges), n_nodes, n_edges,
    sum(cent$is_hub), cent$gene_symbol[1], cent$degree[1],
    round(cent$betweenness[1], 4),
    n_comm, if (exists("wc")) round(igraph::modularity(wc), 4) else NA,
    round(igraph::transitivity(g_cc, type = "global"), 4),
    round(igraph::edge_density(g_cc), 6),
    THRESHOLD$string, STRING_V, "Fruchterman-Reingold"
  )
)
export_tsv(ppi_summary, here::here(DIRS$network, "N04_network_summary.tsv"))

hub_table <- cent |> dplyr::filter(is_hub) |>
  dplyr::select(gene_symbol, regulation, degree, betweenness, closeness, hub_score) |>
  dplyr::arrange(dplyr::desc(betweenness)) |>
  dplyr::mutate(betweenness = round(betweenness, 4),
                closeness   = round(closeness, 6),
                hub_score   = round(hub_score, 4))
export_tsv(hub_table, here::here(DIRS$tables, "T06_hub_proteins.tsv"))

kegg_cent <- kegg_degs |>
  dplyr::left_join(cent |> dplyr::select(gene_symbol, degree, betweenness, is_hub),
                   by = "gene_symbol") |>
  dplyr::arrange(dplyr::desc(abs(logFC))) |>
  dplyr::mutate(logFC = round(logFC, 3),
                adj.P.Val = format(adj.P.Val, digits = 2, scientific = TRUE),
                betweenness = round(betweenness, 4))
export_tsv(kegg_cent, here::here(DIRS$tables, "T07_kegg_degs_ppi.tsv"))
cat("  PPI done.\n")
