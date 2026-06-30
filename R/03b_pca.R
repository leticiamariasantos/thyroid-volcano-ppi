# ═══════════════════════════════════════════════════════════════════════════════
# R/03b_pca.R — PCA/MDS QC plots (Nature Comm / Cell Press typography)
# thyroid-volcano-ppi v3.1.0
#
# PCA is QC/diagnostic — does NOT replace gene-level DEG with limma.
# Font sizes match volcano plot (FONT_SIZE=14, axis.text=12, axis.title=14).
#
# Output: results/figures/Fig_QC_PCA_condition.{png,pdf}
#         results/figures/Fig_QC_PCA_source.{png,pdf}
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── QC: PCA / MDS ──\n")

mds <- limma::plotMDS(E, top = 500, plot = FALSE)

pca_df <- data.frame(
  PC1 = mds$x,
  PC2 = mds$y,
  condition = meta$condition,
  study     = meta$study,
  stringsAsFactors = FALSE
)

var_explained <- round(mds$var.explained * 100, 1)

# ═══════════════════════════════════════════════════════════════════════════════
# PCA theme — matches volcano plot (FONT_SIZE=14, axis.text=12, axis.title=14)
# ═══════════════════════════════════════════════════════════════════════════════
pca_theme <- theme_minimal(base_size = FONT_SIZE, base_family = FONT_FAM) +
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

# ── PCA by CONDITION ──────────────────────────────────────────────────────────
p_pca_cond <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 1.8, alpha = 0.75, stroke = 0.3) +
  scale_color_manual(
    values = c("Normal" = CELL_COLORS$up, "THCA" = CELL_COLORS$down),
    labels = c("Normal" = "Normal (GTEx)",  "THCA"  = "THCA (TCGA)")
  ) +
  labs(
    x     = paste0("PC1 (", var_explained[1], "%)"),
    y     = paste0("PC2 (", var_explained[2], "%)"),
    color = "Condicao"
  ) +
  pca_theme

ggsave(here::here(DIRS$figures, "Fig_QC_PCA_condition.png"),
       p_pca_cond, width = 180, height = 150, units = "mm", dpi = 600, bg = "white")
ggsave(here::here(DIRS$figures, "Fig_QC_PCA_condition.pdf"),
       p_pca_cond, width = 180, height = 150, units = "mm", device = "pdf", bg = "white")

# ── PCA by SOURCE (batch effect assessment) ───────────────────────────────────
p_pca_src <- ggplot(pca_df, aes(x = PC1, y = PC2, color = study)) +
  geom_point(size = 1.8, alpha = 0.75, stroke = 0.3) +
  scale_color_manual(
    values = c("TCGA" = "#CC6677", "GTEX" = "#44AA77"),
    labels = c("TCGA" = "TCGA",     "GTEX" = "GTEx")
  ) +
  labs(
    x     = paste0("PC1 (", var_explained[1], "%)"),
    y     = paste0("PC2 (", var_explained[2], "%)"),
    color = "Estudo"
  ) +
  pca_theme

ggsave(here::here(DIRS$figures, "Fig_QC_PCA_source.png"),
       p_pca_src, width = 180, height = 150, units = "mm", dpi = 600, bg = "white")
ggsave(here::here(DIRS$figures, "Fig_QC_PCA_source.pdf"),
       p_pca_src, width = 180, height = 150, units = "mm", device = "pdf", bg = "white")

cat(sprintf("  PCA: %d samples | PC1: %.1f%% | PC2: %.1f%%\n",
            nrow(pca_df), var_explained[1], var_explained[2]))
cat("  NOTE: PCA is QC only — does not replace gene-level DEG with limma.\n")
cat("  WARNING: TCGA=cancer, GTEx=normal -> batch confounded with condition.\n")
cat("  ComBat/removeBatchEffect NOT applied — would erase biological signal.\n")
