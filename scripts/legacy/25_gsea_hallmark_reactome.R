# ═══════════════════════════════════════════════════════════════════════════════
# 25_gsea_hallmark_reactome.R — ETAPA 8/9: GSEA MSigDB Hallmark + Reactome
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(data.table); library(fgsea); library(msigdbr)})

cat("══ ETAPA 8/9 — GSEA Hallmark + Reactome ══\n")
cat("R:", as.character(getRversion()), "| msigdbr:", as.character(packageVersion("msigdbr")),
    "| fgsea:", as.character(packageVersion("fgsea")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_g <- "results/gsea"
dir.create(dir_g, recursive = TRUE, showWarnings = FALSE)

# ── ranking (t moderada do discovery TPM) ─────────────────────────────────────
rnk <- fread("04_differential_expression/ranking_t.rnk")
r <- rnk[[2]]; names(r) <- rnk[[1]]
r <- r[!is.na(r) & !duplicated(names(r))]
cat("  Ranking:", length(r), "genes\n")

# ── gene sets ─────────────────────────────────────────────────────────────────
msig <- msigdbr(species = "Homo sapiens")
hall <- msig[msig$gs_collection == "H", ]
react <- msig[grepl("REACTOME", msig$gs_subcollection), ]
hall_sets <- split(hall$gene_symbol, hall$gs_name)
react_sets <- split(react$gene_symbol, react$gs_name)
cat("  Hallmark sets:", length(hall_sets), "| Reactome sets:", length(react_sets), "\n")
cat("  MSigDB version:", unique(msig$db_version), "| access:", Sys.Date(), "\n")

run_fgsea <- function(sets, stats, label) {
  set.seed(42)
  res <- fgsea(pathways = sets, stats = stats, minSize = 15, maxSize = 500, nPermSimple = 10000)
  out <- as.data.table(res)[, .(pathway, pval, padj, NES, ES, size)]
  data.table::fwrite(out, file.path(dir_g, paste0("GSEA_", label, ".tsv")), sep = "\t")
  cat(sprintf("  %s: %d sets | padj<0.05: %d | padj<0.25: %d\n",
              label, nrow(out), sum(out$padj < 0.05, na.rm=TRUE), sum(out$padj < 0.25, na.rm=TRUE)))
  out
}

cat("\n── Hallmark ──\n")
hall_out <- run_fgsea(hall_sets, r, "Hallmark")
cat("\n── Reactome ──\n")
react_out <- run_fgsea(react_sets, r, "Reactome")

# ── Hallmark top (reportar TODOS, mas destacar top) ────────────────────────────
cat("\n── Hallmark (top 15 por |NES|) ──\n")
print(hall_out[order(-abs(NES))][1:15, .(pathway, NES, padj)])
cat("\n── Reactome (top 15 por |NES|) ──\n")
print(react_out[order(-abs(NES))][1:15, .(pathway, NES, padj)])

cat("\n══ ETAPA 8/9 CONCLUÍDA ══\n")
