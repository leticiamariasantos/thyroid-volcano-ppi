# ═══════════════════════════════════════════════════════════════════════════════
# 18_figures.R — Figuras da Fase 2 estendida (contagens + validação)
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

dir_f <- "14_figures"
dir.create(dir_f, recursive = TRUE, showWarnings = FALSE)

cat("══ Figuras (Fase 2 estendida) ══\n")

# ── Volcano voom + DESeq2 ─────────────────────────────────────────────────────
voom <- fread("results/counts_deg/voom_full_results.tsv")
dds  <- fread("results/counts_deg/deseq2_full_results.tsv")

volcano <- function(d, logfc_col, padj_col, title, out) {
  dd <- data.frame(logFC = d[[logfc_col]], padj = d[[padj_col]], regulation = d$regulation)
  dd$sig <- ifelse(is.na(dd$padj), "NS",
            ifelse(dd$padj < 0.05 & dd$logFC > 1, "Up",
            ifelse(dd$padj < 0.05 & dd$logFC < -1, "Down", "NS")))
  dd$nlogp <- -log10(pmax(dd$padj, 1e-300))
  dd$nlogp <- pmin(dd$nlogp, 300)
  p <- ggplot(dd, aes(logFC, nlogp, color = sig)) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "grey80")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.4) +
    labs(title = title, x = "log2 fold change", y = "-log10(adj. p)") +
    theme_bw() + theme(legend.position = "none")
  ggsave(out, p, width = 6, height = 5, dpi = 150)
  cat("  ->", out, "\n")
}
volcano(voom, "logFC", "adj.P.Val", "limma-voom (counts)", file.path(dir_f, "Fig_voom_volcano.png"))
volcano(dds, "log2FoldChange", "padj", "DESeq2 (counts)", file.path(dir_f, "Fig_deseq2_volcano.png"))

# ── Concordância de logFC (TPM × voom / TPM × DESeq2) ─────────────────────────
mc <- fread("results/counts_deg/method_concordance.tsv")
p1 <- ggplot(mc, aes(logFC_TPM, logFC_voom)) + geom_point(size = 0.3, alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "logFC TPM vs voom", x = "logFC (TPM)", y = "logFC (voom)") + theme_bw()
p2 <- ggplot(mc, aes(logFC_TPM, logFC_DESeq2)) + geom_point(size = 0.3, alpha = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "logFC TPM vs DESeq2", x = "logFC (TPM)", y = "logFC (DESeq2)") + theme_bw()
ggsave(file.path(dir_f, "Fig_logFC_concordance.png"),
       gridExtra::grid.arrange(p1, p2, nrow = 1), width = 11, height = 5, dpi = 150)
cat("  -> Fig_logFC_concordance.png\n")

# ── Validação GEO: genes-chave (discovery vs GSE33630 vs GSE60542) ────────────
vk <- fread("10_validation/validation_key_genes.tsv")
vl <- melt(vk[, .(gene, discovery_logFC, GSE33630_logFC, GSE60542_logFC)],
           id.vars = "gene", variable.name = "dataset", value.name = "logFC")
vl[dataset == "discovery_logFC", dataset := "discovery"]
vl[dataset == "GSE33630_logFC", dataset := "GSE33630"]
vl[dataset == "GSE60542_logFC", dataset := "GSE60542"]
vl[, dataset := factor(dataset, levels = c("discovery", "GSE33630", "GSE60542"))]
vl[, gene := factor(gene, levels = rev(unique(vk$gene)))]
p3 <- ggplot(vl, aes(dataset, gene, fill = logFC)) + geom_tile() +
  scale_fill_gradient2(low = "#1f77b4", mid = "white", high = "#d62728") +
  labs(title = "Validação GEO — genes congelados", x = NULL, y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(dir_f, "Fig_validation_key_genes.png"), p3, width = 7, height = 5, dpi = 150)
cat("  -> Fig_validation_key_genes.png\n")

cat("══ Figuras CONCLUÍDAS ══\n")
