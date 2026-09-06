#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 04_gsea.R — FASE 2: GSEA global (exploratório, todas as vias KEGG)
#            + painel a priori de 10 vias KEGG (pré-especificado, imutável)
# thyroid-volcano-ppi
#
# Ranking: estatística t moderada (limma eBayes) — preserva magnitude e direção.
# Método: fgsea::fgsea (pre-ranked), nPermSimple=10000, seed=42.
# Gene sets: KEGG (keggLink pathway/hsa + org.Hs.eg.db Entrez->símbolo).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(fgsea)
  library(KEGGREST)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

PROJECT_ROOT <- here::here()
dir_gsea <- file.path(PROJECT_ROOT, "05_gsea")
dir_panel <- file.path(PROJECT_ROOT, "06_kegg_panel")
dir.create(dir_gsea, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_panel, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — GSEA global + painel KEGG ══\n")
cat("R:", as.character(getRversion()), "| fgsea:", as.character(packageVersion("fgsea")),
    "| KEGGREST:", as.character(packageVersion("KEGGREST")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Ranking (t moderada) ───────────────────────────────────────────────────
rnk <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "ranking_t.rnk"),
                         header = FALSE)
ranking <- setNames(rnk$V2, rnk$V1)
ranking <- ranking[order(ranking, decreasing = TRUE)]
ranking <- ranking[!duplicated(names(ranking))]
cat(sprintf("  Ranking: %d genes (min=%.2f max=%.2f)\n", length(ranking), min(ranking), max(ranking)))

# ── 2. Gene sets KEGG (todas as vias hsa) ─────────────────────────────────────
cat("  Baixando gene sets KEGG (link pathway/hsa)...\n")
link <- tryCatch(keggLink("pathway", "hsa"), error = function(e) stop("KEGG link falhou: ", e$message))
# link: nomes = "hsa:ENTREZ", valores = "path:hsaXXXXX"
path_map <- sub("path:", "", unname(link))
entrez   <- sub("hsa:", "", names(link))
df <- data.frame(entrez = entrez, pathway = sub("hsa", "", path_map), stringsAsFactors = FALSE)

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(df$entrez),
                             keytype = "ENTREZID", column = "SYMBOL",
                             multiVals = "first")
df$symbol <- sym[df$entrez]
df <- df[!is.na(df$symbol), ]

pathway_names <- tryCatch(keggList("pathway", "hsa"), error = function(e) character(0))
pathway_name_map <- sub(" - Homo sapiens.*$", "", pathway_names)
names(pathway_name_map) <- sub("hsa", "", names(pathway_names))

pathways <- split(df$symbol, df$pathway)
pathways <- lapply(pathways, unique)
# nomear via: hsaID -> descrição curta
for (p in names(pathways)) {
  desc <- if (p %in% names(pathway_name_map)) pathway_name_map[[p]] else p
  names(pathways)[names(pathways) == p] <- paste0(p, "__", desc)
}
cat(sprintf("  Vias KEGG com genes: %d (total genes únicos %d)\n",
            length(pathways), length(unique(df$symbol))))

# ── 3. fgsea (global exploratório) ────────────────────────────────────────────
cat("  Executando fgsea (nPermSimple=10000, seed=42)...\n")
set.seed(42)
gsea <- fgsea(pathways = pathways, stats = ranking,
              minSize = 15, maxSize = 500, nPermSimple = 10000)
gsea <- gsea[order(gsea$pval), ]
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))
gsea_out <- as.data.frame(gsea)
gsea_out$description <- sub("^.*__", "", gsea$pathway)
gsea_out$pathway <- sub("__.*$", "", gsea$pathway)
data.table::fwrite(gsea_out[, c("pathway","description","pval","padj","NES","ES","size","leadingEdge")],
                   file.path(dir_gsea, "GSEA_KEGG_all.tsv"), sep = "\t")
cat(sprintf("  Vias significativas (padj<0.25): %d | (padj<0.05): %d\n",
            sum(gsea$padj < 0.25, na.rm = TRUE), sum(gsea$padj < 0.05, na.rm = TRUE)))

# ── 4. Painel a priori (10 vias, IMUTÁVEL) ────────────────────────────────────
A_PRIORI <- c(
  "hsa05216","hsa04919","hsa04010","hsa04151","hsa04150",
  "hsa04115","hsa04210","hsa04110","hsa04310","hsa04064"
)
cat("\n── Painel a priori (10 vias KEGG) ──\n")
# carrega DEG para contagens up/down no painel
deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
up_genes_sym   <- deg$gene_symbol[deg$regulation == "Up"]
down_genes_sym <- deg$gene_symbol[deg$regulation == "Down"]

panel <- lapply(A_PRIORI, function(h) {
  h_short <- sub("^hsa", "", h)
  hit <- which(sub("__.*$", "", names(pathways)) == h_short)
  if (length(hit) == 0) return(data.frame(pathway=h, present=FALSE))
  p <- pathways[[hit[1]]]
  gin <- intersect(p, names(ranking))
  g <- gsea[sub("__.*$", "", gsea$pathway) == h_short, ]
  data.frame(
    pathway = h,
    description = if (h_short %in% names(pathway_name_map)) pathway_name_map[[h_short]] else h_short,
    genes_in_pathway = length(p),
    genes_in_matrix = length(gin),
    n_up = sum(gin %in% up_genes_sym),
    n_down = sum(gin %in% down_genes_sym),
    NES = if (nrow(g) > 0) round(g$NES[1], 4) else NA_real_,
    pval = if (nrow(g) > 0) signif(g$pval[1], 4) else NA_real_,
    padj = if (nrow(g) > 0) signif(g$padj[1], 4) else NA_real_,
    present = TRUE
  )
})
panel <- do.call(rbind, panel)
print(panel, row.names = FALSE)
data.table::fwrite(panel, file.path(dir_panel, "KEGG_a_priori_panel.tsv"), sep = "\t")

# ── 5. Top vias (para inspeção, NÃO seleção) ──────────────────────────────────
cat("\n  Top 20 vias por pval (inspeção, NÃO = seleção):\n")
print(head(gsea[, c("pathway","NES","pval","padj","size")], 20), row.names = FALSE)

cat("\n══ FASE 2 — GSEA CONCLUÍDO ══\n")
