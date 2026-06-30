# ═══════════════════════════════════════════════════════════════════════════════
# R/04_volcano.R — Volcano Plot (Nature Communications / Cell Press standard)
# thyroid-volcano-ppi — MAIN OUTPUT 1/2
#
# Style: White background, clean axes, dashed significance thresholds,
#        points with thin black outline + colored fill,
#        blue = upregulated, magenta = downregulated,
#        large ggrepel labels with leader lines, compact right-side legend.
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── M4: Volcano Plot (hubs from PPI) ──\n")

deg$pval_plot <- pmax(deg$adj.P.Val, .Machine$double.xmin)
deg$log_pval  <- -log10(deg$pval_plot)

# ── Axis limits ───────────────────────────────────────────────────────────────
x_max   <- max(abs(deg$logFC),  na.rm = TRUE) * 1.12
y_max   <- max(deg$log_pval,    na.rm = TRUE) * 1.10
x_breaks <- pretty(c(-x_max, x_max), 7)
x_breaks <- x_breaks[x_breaks >= -x_max & x_breaks <= x_max]
y_breaks <- Filter(function(y) y <= y_max, pretty(c(0, y_max), 6))

# ── Point sizing (proportional to significance, no external scales dep) ──────
deg$pt_size <- ifelse(deg$regulation == "NS", 0.65,
                      scales::rescale(deg$log_pval, to = c(1.0, 3.2)))

# ── Label selection (most biologically relevant genes) ────────────────────────
top_fdr     <- deg |> dplyr::filter(regulation != "NS") |>
  dplyr::arrange(pval_plot) |> dplyr::slice_head(n = 8) |> dplyr::pull(gene_symbol)
top_lfc_up   <- deg |> dplyr::filter(regulation == "Up") |>
  dplyr::arrange(dplyr::desc(logFC)) |> dplyr::slice_head(n = 5) |> dplyr::pull(gene_symbol)
top_lfc_down <- deg |> dplyr::filter(regulation == "Down") |>
  dplyr::arrange(logFC) |> dplyr::slice_head(n = 5) |> dplyr::pull(gene_symbol)
top_kegg     <- deg |> dplyr::filter(regulation != "NS", in_kegg) |>
  dplyr::arrange(pval_plot) |> dplyr::slice_head(n = 6) |> dplyr::pull(gene_symbol)

label_set <- Reduce(union, list(top_fdr, top_lfc_up, top_lfc_down, top_kegg))

# ── Hub genes from PPI network (cent must exist; pipeline now runs PPI first) ─
if (exists("cent") && is.data.frame(cent)) {
  hubs <- cent |> dplyr::filter(is_hub) |> dplyr::pull(gene_symbol)
  label_set <- union(label_set, intersect(hubs, deg$gene_symbol))
}

# ── FALLBACK: force critical hub genes even if cent unavailable ───────────
FORCE_LABEL <- c("PRKCA", "PLCG1", "PLCD4", "PRKCG", "PLCD3")
label_set <- union(label_set, intersect(FORCE_LABEL, deg$gene_symbol))
deg$label <- ifelse(deg$gene_symbol %in% label_set, deg$gene_symbol, "")

cat(sprintf("  Labels: %d genes\n", length(label_set)))

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN VOLCANO PLOT
# ═══════════════════════════════════════════════════════════════════════════════

p <- ggplot(deg, aes(x = logFC, y = log_pval)) +

  # ── Minimal grid ─────────────────────────────────────────────────────────
  geom_hline(yintercept = y_breaks,
             linewidth = 0.10, color = CELL_COLORS$grid, linetype = "dotted") +
  geom_vline(xintercept = x_breaks,
             linewidth = 0.10, color = CELL_COLORS$grid, linetype = "dotted") +

  # ── NS points: small, grey fill, thin dark outline ───────────────────────
  geom_point(
    data = subset(deg, regulation == "NS"),
    aes(size = pt_size),
    shape = 21, fill = CELL_COLORS$ns,
    color = CELL_COLORS$outline, stroke = 0.08, alpha = 0.40
  ) +

  # ── DEGs: colored fill with thin dark outline ────────────────────────────
  geom_point(
    data = subset(deg, regulation != "NS"),
    aes(size = pt_size, fill = regulation),
    shape = 21, color = CELL_COLORS$outline, stroke = 0.12, alpha = 0.92
  ) +

  # ── KEGG ring: subtle open circle around pathway genes ───────────────────
  geom_point(
    data = subset(deg, regulation != "NS" & in_kegg),
    shape = 1, size = 3.2, stroke = 0.25,
    color = "#555555", show.legend = FALSE
  ) +

  # ── Dashed significance thresholds ───────────────────────────────────────
  geom_vline(xintercept = c(-THRESHOLD$lfc, THRESHOLD$lfc),
             linetype = "dashed", linewidth = 0.40, color = "#555555") +
  geom_hline(yintercept = -log10(THRESHOLD$fdr),
             linetype = "dashed", linewidth = 0.40, color = "#555555") +

  # ── Gene labels: large, italic, with leader lines ────────────────────────
  ggrepel::geom_text_repel(
    data = subset(deg, label != ""),
    aes(label = label),
    size = 5.2, fontface = "italic", color = "#1A1A1A",
    max.overlaps = 50,
    box.padding = 0.60, point.padding = 0.28,
    segment.size = 0.25, segment.color = "grey45",
    segment.alpha = 0.7, segment.linetype = "solid",
    min.segment.length = 0.02,
    force = 3, force_pull = 0.7,
    nudge_x = 0.25, nudge_y = 0.25,
    direction = "both",
    seed = 42, show.legend = FALSE
  ) +

  # ── Color scale (legend on right) ────────────────────────────────────────
  scale_fill_manual(
    values = c("Up"   = CELL_COLORS$up,
               "Down" = CELL_COLORS$down),
    labels = c(
      "Up"   = paste0("Superexpresso  (", n_up,  ")"),
      "Down" = paste0("Subexpresso    (", n_down, ")")
    ),
    name = "THCA vs Normal"
  ) +

  scale_size(range = c(0.5, 3.5), guide = "none") +

  # ── Axis scales ──────────────────────────────────────────────────────────
  scale_x_continuous(
    limits = c(-x_max, x_max),
    breaks = x_breaks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = expansion(mult = c(0.02, 0.05))
  ) +

  # ── Axis labels ──────────────────────────────────────────────────────────
  labs(
    x = expression(log[2] ~ "(fold change)  —  THCA / Normal"),
    y = expression(-log[10] ~ "(valor-p ajustado)")
  ) +

  # ═══════════════════════════════════════════════════════════════════════════
  # THEME: Nature Communications / Cell Press style
  # ═══════════════════════════════════════════════════════════════════════════
  theme_minimal(base_size = FONT_SIZE, base_family = FONT_FAM) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.line.x       = element_line(color = "black", linewidth = 0.45),
    axis.line.y       = element_line(color = "black", linewidth = 0.45),
    axis.ticks        = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(2.0, "mm"),

    axis.text  = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = FONT_SIZE, face = "plain"),

    legend.position      = "right",
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", color = "grey75",
                                        linewidth = 0.25),
    legend.key           = element_rect(fill = "white", color = NA),
    legend.key.size      = unit(6, "mm"),
    legend.text          = element_text(size = 11, color = "black"),
    legend.title         = element_text(size = 12, face = "bold", color = "black"),
    legend.margin        = margin(6, 8, 6, 6),

    plot.margin = margin(8, 10, 8, 8)
  )

# ── Export (PNG 600 dpi + PDF) ──────────────────────────────────────────────
ggsave(here::here(DIRS$figures, "Fig1_Volcano_THCA_vs_Normal.png"),
       p, width = 180, height = 150, units = "mm", dpi = 600, bg = "white")
ggsave(here::here(DIRS$figures, "Fig1_Volcano_THCA_vs_Normal.pdf"),
       p, width = 180, height = 150, units = "mm", device = "pdf", bg = "white")
cat(sprintf("  Volcano: %d genes | %d DEGs (up %d, down %d)\n", nrow(deg), n_deg, n_up, n_down))
