# ═══════════════════════════════════════════════════════════════════════════════
# R/03c_heatmap.R — Heatmap of top DEGs
# thyroid-volcano-ppi v3.1.0
#
# Generates heatmap of top DEGs (by |logFC|) across all samples.
# Uses pheatmap with row-scaling (z-score).
#
# Output: results/figures/FigS2_Heatmap_TopDEGs.{png,pdf}
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── S2: Heatmap (Top DEGs) ──\n")

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}
library(pheatmap)

tryCatch({

  # ── Select top 20 DEGs by |logFC| ───────────────────────────────────────────
  top_deg <- deg |>
    dplyr::filter(regulation != "NS") |>
    dplyr::arrange(dplyr::desc(abs(logFC))) |>
    dplyr::slice_head(n = 20)

  # ── Get genes present in expression matrix ──────────────────────────────────
  genes_present <- intersect(top_deg$gene_symbol, rownames(E))
  
  if (length(genes_present) < 5) {
    stop("Fewer than 5 top genes found in expression matrix")
  }

  # ── Subset expression ───────────────────────────────────────────────────────
  expr_sub <- E[genes_present, , drop = FALSE]

  # ── Sample annotation ───────────────────────────────────────────────────────
  ann_col <- data.frame(
    Condition = meta$condition,
    row.names = colnames(expr_sub)
  )

  ann_colors <- list(
    Condition = c(Normal = "#4477AA", THCA = "#AA4488")
  )

  # ── Row-scaling (z-score) ───────────────────────────────────────────────────
  expr_scaled <- t(scale(t(expr_sub)))

  n_genes <- nrow(expr_scaled)
  n_samples <- ncol(expr_scaled)

  # ── Generate heatmap ────────────────────────────────────────────────────────
  p_hm <- pheatmap(
    expr_scaled,
    color = colorRampPalette(c("#4477AA", "white", "#AA4488"))(100),
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = if (n_genes > 20) 7 else 9,
    fontsize_col = 5,
    border_color = NA,
    main = paste0("Top ", n_genes, " DEGs — THCA vs Normal"),
    treeheight_row = 15,
    treeheight_col = 15,
    silent = TRUE
  )

  # ── Export ──────────────────────────────────────────────────────────────────
  fig_h <- max(100, n_genes * 8)
  png(here::here(DIRS$figures, "FigS2_Heatmap_TopDEGs.png"),
      width = 180, height = fig_h, units = "mm", res = 600, bg = "white")
  grid::grid.draw(p_hm$gtable)
  dev.off()

  pdf(here::here(DIRS$figures, "FigS2_Heatmap_TopDEGs.pdf"),
      width = 180/25.4, height = fig_h/25.4)
  grid::grid.draw(p_hm$gtable)
  dev.off()

  cat(sprintf("  Heatmap: %d genes x %d samples\n", n_genes, n_samples))

}, error = function(e) {
  cat("  Heatmap skipped:", conditionMessage(e), "\n")
})
