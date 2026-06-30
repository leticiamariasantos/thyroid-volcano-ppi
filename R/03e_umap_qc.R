# ═══════════════════════════════════════════════════════════════════════════════
# R/03e_umap_qc.R — UMAP QC plots (Nature Comm / Cell Press typography)
# thyroid-volcano-ppi v3.1.0
#
# UMAP is QC/diagnostic — does NOT replace gene-level DEG with limma.
# Font sizes match volcano plot (FONT_SIZE=14, axis.text=12, axis.title=14).
#
# Output: results/figures/Fig_QC_UMAP_condition.{png,pdf}
#         results/figures/Fig_QC_UMAP_source.{png,pdf}
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── QC: UMAP ──\n")

if (!requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot", repos = "https://cloud.r-project.org")
}
library(uwot)

# ═══════════════════════════════════════════════════════════════════════════════
# UMAP theme — matches volcano plot (FONT_SIZE=14, axis.text=12, axis.title=14)
# ═══════════════════════════════════════════════════════════════════════════════
umap_theme <- theme_minimal(base_size = FONT_SIZE, base_family = FONT_FAM) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = CELL_COLORS$grid, linewidth = 0.15),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.45),
    axis.ticks        = element_line(color = "black", linewidth = 0.45),
    axis.ticks.length = unit(2.0, "mm"),
    axis.text         = element_text(color = "black", size = 12),
    axis.title        = element_text(color = "black", size = FONT_SIZE, face = "plain"),
    legend.position      = "right",
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", color = "grey75", linewidth = 0.25),
    legend.key           = element_rect(fill = "white", color = NA),
    legend.key.size      = unit(6, "mm"),
    legend.text          = element_text(size = 11, color = "black"),
    legend.title         = element_text(size = 12, face = "bold", color = "black"),
    legend.margin        = margin(6, 8, 6, 6),
    plot.margin          = margin(8, 10, 8, 8)
  )

tryCatch({

  set.seed(42)
  X <- unname(as.matrix(t(E)))
  mode(X) <- "numeric"

  umap_res <- uwot::umap(
    X,
    n_neighbors  = 30,
    min_dist     = 0.3,
    n_components = 2,
    metric       = "euclidean",
    verbose      = FALSE
  )

  colnames(umap_res) <- c("UMAP1", "UMAP2")
  umap_df <- data.frame(
    UMAP1     = umap_res[, 1],
    UMAP2     = umap_res[, 2],
    condition = meta$condition,
    study     = meta$study,
    stringsAsFactors = FALSE
  )

  # ── UMAP by CONDITION ─────────────────────────────────────────────────────
  p_umap_cond <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
    geom_point(size = 1.8, alpha = 0.75, stroke = 0.3) +
    scale_color_manual(
      values = c("Normal" = CELL_COLORS$up, "THCA" = CELL_COLORS$down),
      labels = c("Normal" = "Normal (GTEx)", "THCA"  = "THCA (TCGA)")
    ) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Condicao") +
    umap_theme

  ggsave(here::here(DIRS$figures, "Fig_QC_UMAP_condition.png"),
         p_umap_cond, width = 180, height = 150, units = "mm",
         dpi = 600, bg = "white")
  ggsave(here::here(DIRS$figures, "Fig_QC_UMAP_condition.pdf"),
         p_umap_cond, width = 180, height = 150, units = "mm",
         device = "pdf", bg = "white")

  # ── UMAP by SOURCE (batch effect assessment) ──────────────────────────────
  p_umap_src <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = study)) +
    geom_point(size = 1.8, alpha = 0.75, stroke = 0.3) +
    scale_color_manual(
      values = c("TCGA" = "#CC6677", "GTEX" = "#44AA77"),
      labels = c("TCGA" = "TCGA",     "GTEX" = "GTEx")
    ) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Estudo") +
    umap_theme

  ggsave(here::here(DIRS$figures, "Fig_QC_UMAP_source.png"),
         p_umap_src, width = 180, height = 150, units = "mm",
         dpi = 600, bg = "white")
  ggsave(here::here(DIRS$figures, "Fig_QC_UMAP_source.pdf"),
         p_umap_src, width = 180, height = 150, units = "mm",
         device = "pdf", bg = "white")

  cat("  UMAP QC plots exported.\n")

}, error = function(e) {
  cat("  UMAP skipped:", conditionMessage(e), "\n")
})
