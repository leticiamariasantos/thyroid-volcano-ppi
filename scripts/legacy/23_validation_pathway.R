# ═══════════════════════════════════════════════════════════════════════════════
# 23_validation_pathway.R — PARTE 15: validação de vias (GSEA) nos GEO
#
# GSEA pre-ranked (t moderada) nos GEO, usando os MESMOS gene sets KEGG do discovery.
# Compara NES/direção entre discovery e cada GEO.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table); library(fgsea)})

cat("══ PARTE 15 — Validação de vias (GSEA GEO) ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_v <- "results/validation"
dir.create(dir_v, recursive = TRUE, showWarnings = FALSE)

kegg <- readRDS("07_pathway_redundancy/kegg_pathways.rds")  # named list pathway -> genes
disco_gsea <- fread("05_gsea/GSEA_KEGG_all.tsv")

run_gsea <- function(deg_file) {
  deg <- fread(deg_file)
  # ranking por t moderada (descarta NA)
  rnk <- deg[!is.na(t) & !is.na(gene_symbol)]
  rnk <- rnk[order(-t)]
  r <- rnk$t; names(r) <- rnk$gene_symbol
  r <- r[!duplicated(names(r))]
  set.seed(42)
  fgsea(pathways = kegg, stats = r, minSize = 15, maxSize = 500, nPermSimple = 10000)
}

g1 <- run_gsea("10_validation/GSE33630_deg.tsv")
g2 <- run_gsea("10_validation/GSE60542_deg.tsv")

g1_out <- as.data.table(g1)[, .(pathway, pval, padj, NES, ES, size)]
g2_out <- as.data.table(g2)[, .(pathway, pval, padj, NES, ES, size)]
data.table::fwrite(g1_out, file.path(dir_v, "GSE33630_GSEA.tsv"), sep = "\t")
data.table::fwrite(g2_out, file.path(dir_v, "GSE60542_GSEA.tsv"), sep = "\t")

# ── compara vias não-redundantes do discovery ─────────────────────────────────
sel <- fread("07_pathway_redundancy/pathways_selected.tsv", colClasses = c(id = "character"))
disco_gsea <- fread("05_gsea/GSEA_KEGG_all.tsv", colClasses = c(pathway = "character"))
cat("\n── Vias discovery não-redundantes: NES discovery vs GEO ──\n")
cmp <- rbindlist(lapply(sel$id, function(p) {
  d <- disco_gsea[pathway == p]
  a <- g1_out[pathway == p]; b <- g2_out[pathway == p]
  data.table(pathway = p,
             discovery_NES = if(nrow(d)) d$NES[1] else NA_real_,
             discovery_padj = if(nrow(d)) d$padj[1] else NA_real_,
             GSE33630_NES = if(nrow(a)) a$NES[1] else NA_real_,
             GSE33630_padj = if(nrow(a)) a$padj[1] else NA_real_,
             GSE60542_NES = if(nrow(b)) b$NES[1] else NA_real_,
             GSE60542_padj = if(nrow(b)) b$padj[1] else NA_real_)
}))
cmp[, replicated := (sign(GSE33630_NES) == sign(discovery_NES) & GSE33630_padj < 0.25) |
                     (sign(GSE60542_NES) == sign(discovery_NES) & GSE60542_padj < 0.25)]
cmp[, discovery_NES := round(discovery_NES, 2)]
cmp[, GSE33630_NES := round(GSE33630_NES, 2)]
cmp[, GSE60542_NES := round(GSE60542_NES, 2)]
print(cmp)
data.table::fwrite(cmp, file.path(dir_v, "pathway_validation.tsv"), sep = "\t")

cat("\n══ PARTE 15 CONCLUÍDA ══\n")
