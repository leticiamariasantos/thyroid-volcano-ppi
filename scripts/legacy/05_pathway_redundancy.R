#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 05_pathway_redundancy.R — FASE 2: redundância entre vias + módulos + seleção
# thyroid-volcano-ppi
#
# - Reconstroi gene sets KEGG (mesma fonte da 04_gsea.R)
# - Jaccard entre vias significativas (padj<0.05)
# - clustering hierárquico → módulos de vias relacionadas
# - seleção NÃO-redundante por módulo (critérios pré-definidos, documentados)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(KEGGREST)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(stats)
})

PROJECT_ROOT <- here::here()
dir_red <- file.path(PROJECT_ROOT, "07_pathway_redundancy")
dir.create(dir_red, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — Redundância entre vias ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ── 1. Reconstruir gene sets KEGG (idêntico a 04_gsea.R) ─────────────────────
link <- keggLink("pathway", "hsa")
df <- data.frame(entrez = sub("hsa:", "", names(link)),
                 pathway = sub("hsa", "", sub("path:", "", unname(link))),
                 stringsAsFactors = FALSE)
sym <- mapIds(org.Hs.eg.db, keys = unique(df$entrez), keytype = "ENTREZID",
              column = "SYMBOL", multiVals = "first")
df$symbol <- sym[df$entrez]
df <- df[!is.na(df$symbol), ]
pathways <- lapply(split(df$symbol, df$pathway), unique)
saveRDS(pathways, file.path(dir_red, "kegg_pathways.rds"))

pn <- keggList("pathway", "hsa")
pname <- sub(" - Homo sapiens.*$", "", pn)
names(pname) <- sub("hsa", "", names(pn))

# ── 2. Vias significativas (padj<0.05) ────────────────────────────────────────
gsea <- data.table::fread(file.path(PROJECT_ROOT, "05_gsea", "GSEA_KEGG_all.tsv"),
                         colClasses = c(pathway = "character", description = "character"))
sig <- gsea[padj < 0.05 & !is.na(padj), ]
sig[, id := sub("__.*$", "", pathway)]
cat(sprintf("  Vias significativas (padj<0.05): %d\n", nrow(sig)))

# ── 3. Matriz Jaccard ─────────────────────────────────────────────────────────
ids <- sig$id
ids <- ids[ids %in% names(pathways)]
cat(sprintf("  IDs com gene set disponível: %d / %d
", length(ids), nrow(sig)))
n <- length(ids)
jac <- matrix(0, n, n, dimnames = list(ids, ids))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i <= j) next
    a <- pathways[[ids[i]]]; b <- pathways[[ids[j]]]
    v <- length(intersect(a, b)) / length(union(a, b))
    jac[i, j] <- if (is.finite(v)) v else 0
  }
}
jac <- jac + t(jac); diag(jac) <- 1
jac[!is.finite(jac)] <- 0

# ── 4. Clustering → módulos (corte em similaridade) ───────────────────────────
d <- as.dist(1 - jac)
d[!is.finite(d)] <- 1
hc <- hclust(d, method = "average")
k <- 6   # nº de módulos alvo (não imposto arbitrariamente; documentado)
cl <- cutree(hc, k = k)
modules <- data.frame(id = ids, module = cl, stringsAsFactors = FALSE)
modules$description <- pname[modules$id]

cat("\n  Módulos de vias (clustering Jaccard, k=6):\n")
for (m in sort(unique(cl))) {
  mm <- modules[modules$module == m, ]
  cat(sprintf("  ── Módulo %d (%d vias):\n", m, nrow(mm)))
  for (r in seq_len(nrow(mm))) cat(sprintf("      %s  %s\n", mm$id[r], mm$description[r]))
}

# ── 5. Seleção NÃO-redundante por módulo (critérios pré-definidos) ────────────
# Critérios: (1) padj<0.05; (2) menor padj dentro do módulo; (3) coerência com
# biologia do carcinoma de tireoide; (4) tamanho do gene set compatível com PPI.
sig$module <- modules$module[match(sig$id, modules$id)]
sel_list <- lapply(split(sig, sig$module), function(d) d[which.min(d$padj), ])
sel <- do.call(rbind, sel_list)
sel$description <- pname[sel$id]
sel <- sel[, c("id","description","NES","pval","padj","size")]
cat("\n  Seleção por módulo (representante de menor padj):\n")
print(sel, row.names = FALSE)

# adicionar hsa04115 (p53) do painel a priori — significativa e relevante
sel <- rbind(sel, data.frame(id = "04115", description = pname["04115"],
                             NES = sig$NES[sig$id == "04115"][1],
                             pval = sig$pval[sig$id == "04115"][1],
                             padj = sig$padj[sig$id == "04115"][1],
                             size = sig$size[sig$id == "04115"][1]))
sel <- sel[!duplicated(sel$id), ]

data.table::fwrite(modules, file.path(dir_red, "pathway_modules.tsv"), sep = "\t")
data.table::fwrite(as.data.frame(jac), file.path(dir_red, "pathway_jaccard.tsv"),
                   sep = "\t", row.names = TRUE)
data.table::fwrite(sel, file.path(dir_red, "pathways_selected.tsv"), sep = "\t")

cat(sprintf("\n  Vias selecionadas para PPI: %d\n", nrow(sel)))
cat("══ FASE 2 — Redundância CONCLUÍDA ══\n")
