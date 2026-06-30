# ═══════════════════════════════════════════════════════════════════════════════
# R/03_deg.R — Differential expression analysis (limma)
# thyroid-volcano-ppi
#
# AUDIT: ✓ Variable naming  ✓ P-value formatting  ✓ dplyr:: qualified
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── M2: Differential Expression (limma) ──\n")

# ── Transpose to gene × sample ────────────────────────────────────────────────
E <- t(expr)
storage.mode(E) <- "numeric"

# ── Expression filter ─────────────────────────────────────────────────────────
keep <- rowMeans(E > THRESHOLD$expr_min) >= THRESHOLD$expr_frac
E <- E[keep, , drop = FALSE]
cat(sprintf("  Filter: %d kept / %d removed (%.1f%%)\n",
            sum(keep), sum(!keep), sum(keep)/length(keep)*100))

# ── Design matrix ─────────────────────────────────────────────────────────────
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
stopifnot(ncol(E) == nrow(design))

# ── Fit + contrast + empirical Bayes ──────────────────────────────────────────
fit      <- lmFit(E, design)
contrast <- makeContrasts(THCA_vs_Normal = THCA - Normal, levels = design)
fit2     <- contrasts.fit(fit, contrast)
fit2     <- eBayes(fit2)

# ── Results table ─────────────────────────────────────────────────────────────
deg <- topTable(fit2, coef = "THCA_vs_Normal", number = Inf,
                adjust.method = "BH", sort.by = "P") |>
  tibble::rownames_to_column(var = "gene_symbol") |>
  dplyr::mutate(gene_symbol = as.character(gene_symbol))

# ── DEG classification ────────────────────────────────────────────────────────
deg$regulation <- with(deg, ifelse(
  adj.P.Val < THRESHOLD$fdr & logFC >  THRESHOLD$lfc, "Up",
  ifelse(adj.P.Val < THRESHOLD$fdr & logFC < -THRESHOLD$lfc, "Down", "NS")
))

cat("\n  DEGs (|log2FC| >", THRESHOLD$lfc, ", FDR <", THRESHOLD$fdr, "):\n", sep = "")
print(table(deg$regulation))

n_up   <- sum(deg$regulation == "Up")
n_down <- sum(deg$regulation == "Down")
n_deg  <- n_up + n_down

# ── KEGG pathway annotation ───────────────────────────────────────────────────
kegg_genes <- fetch_kegg_genes(KEGG_ID)
kegg_list  <- kegg_genes$gene_symbol
deg$in_kegg <- deg$gene_symbol %in% kegg_list

kegg_degs <- deg |>
  dplyr::filter(in_kegg, regulation != "NS") |>
  dplyr::select(gene_symbol, logFC, AveExpr, adj.P.Val, regulation, B) |>
  dplyr::arrange(dplyr::desc(abs(logFC)))

cat(sprintf("\n  KEGG: %d genes | %d in data | %d DEGs (%d↑ %d↓)\n",
            length(kegg_list), sum(deg$in_kegg), nrow(kegg_degs),
            sum(kegg_degs$regulation == "Up"),
            sum(kegg_degs$regulation == "Down")))

# Missing KEGG genes
kegg_missing <- setdiff(kegg_list, deg$gene_symbol)
if (length(kegg_missing) > 0) {
  cat("  Missing from expression:", length(kegg_missing), "\n")
  export_tsv(tibble(gene_symbol = sort(kegg_missing),
                    reason = "Not detected"),
             here::here(DIRS$tables, "T05_kegg_missing_genes.tsv"))
}

# ── DEG gene lists ────────────────────────────────────────────────────────────
up_genes   <- deg |> dplyr::filter(regulation == "Up")   |> dplyr::pull(gene_symbol)
down_genes <- deg |> dplyr::filter(regulation == "Down") |> dplyr::pull(gene_symbol)
deg_genes  <- sort(unique(c(up_genes, down_genes)))

# ── QC summary table ──────────────────────────────────────────────────────────
deg_summary <- tibble(
  parameter = c(
    "Total genes tested", "Expression filter retained",
    "Upregulated (THCA)", "Downregulated (THCA)", "Total DEGs",
    "|log2FC| threshold", "FDR threshold (BH)", "Expression filter min",
    "Expression filter fraction", "KEGG pathway",
    "KEGG genes in data", "KEGG DEGs", "KEGG Up", "KEGG Down"
  ),
  value = c(
    nrow(deg), sum(keep), n_up, n_down, n_deg,
    THRESHOLD$lfc, THRESHOLD$fdr, THRESHOLD$expr_min,
    paste0(THRESHOLD$expr_frac * 100, "%"),
    KEGG_ID, sum(deg$in_kegg), nrow(kegg_degs),
    sum(kegg_degs$regulation == "Up"), sum(kegg_degs$regulation == "Down")
  )
)
export_tsv(deg_summary, here::here(DIRS$tables, "T02_deg_summary.tsv"))

# ── Full results (clean format) ───────────────────────────────────────────────
deg_out <- deg |>
  dplyr::mutate(
    logFC     = round(logFC, 4),
    AveExpr   = round(AveExpr, 2),
    t         = round(t, 2),
    P.Value   = format(P.Value, digits = 3, scientific = TRUE),
    adj.P.Val = format(adj.P.Val, digits = 3, scientific = TRUE),
    B         = round(B, 2)
  ) |>
  dplyr::select(gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B,
                regulation, in_kegg) |>
  dplyr::arrange(adj.P.Val)

export_tsv(deg_out, here::here(DIRS$tables, "T03_deg_full_results.tsv"))

# ── Top 20 DEGs ───────────────────────────────────────────────────────────────
top20 <- deg |>
  dplyr::filter(regulation != "NS") |>
  dplyr::arrange(dplyr::desc(abs(logFC))) |>
  dplyr::slice_head(n = 20) |>
  dplyr::select(gene_symbol, logFC, AveExpr, adj.P.Val, regulation, B, in_kegg) |>
  dplyr::mutate(
    logFC = round(logFC, 3),
    AveExpr = round(AveExpr, 2),
    adj.P.Val = format(adj.P.Val, digits = 2, scientific = TRUE),
    B = round(B, 1)
  )

export_tsv(top20, here::here(DIRS$tables, "T04_top20_degs.tsv"))

cat("\n  DEG done: ↑", n_up, "| ↓", n_down, "| Σ", n_deg, "\n")
