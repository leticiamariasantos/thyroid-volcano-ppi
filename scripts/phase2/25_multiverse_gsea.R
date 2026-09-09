# fgsea + clusterProfiler + GSVA across matrix variants and five collections.
suppressPackageStartupMessages({
  library(here); library(data.table); library(fgsea); library(msigdbr)
  library(clusterProfiler); library(GSVA)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

msig_cache <- file.path(DIR_CACHE, "msigdb", as.character(packageVersion("msigdbr")), "gene_sets.rds")
sets <- safe_api_cache(msig_cache, function() {
  mk <- function(collection, subcollection = NULL) {
    z <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    split(z$gene_symbol, z$gs_name)
  }
  list(Hallmark = mk("H"), KEGG = mk("C2", "CP:KEGG_LEGACY"),
       GO_BP = mk("C5", "GO:BP"))
}, validate = function(x) is.list(x) && all(c("Hallmark", "KEGG", "GO_BP") %in% names(x)))

react_lines <- readLines(file.path(DIR_EXTERNAL, "reactome", "ReactomePathways.gmt"))
parse_gmt <- function(z) { p <- strsplit(z, "\t", fixed = TRUE); setNames(lapply(p, function(x) x[-c(1, 2)]), vapply(p, `[`, "", 1)) }
sets$Reactome <- parse_gmt(react_lines)
sets$Panel30 <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))

variants <- c("tcga_matched", "tcga_paired", "raw", "batch_corrected", "composition_adjusted")
methods <- c(voom_qw = "voom_qw", DESeq2 = "DESeq2")
all_fg <- list()
for (variant in variants) for (method in names(methods)) {
  f <- file.path(DIR_DE_MULTI, variant, paste0(methods[[method]], "_full_results.tsv"))
  tab <- fread(f)
  ranks <- setNames(tab$statistic, tab$gene_symbol)
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  for (collection in names(sets)) {
    gs <- lapply(sets[[collection]], intersect, y = names(ranks))
    gs <- gs[lengths(gs) >= 10L & lengths(gs) <= 1000L]
    set.seed(SEED)
    fg <- as.data.table(fgseaMultilevel(gs, ranks, minSize = 10, maxSize = 1000,
                                        eps = 0, nPermSimple = 10000))
    if (nrow(fg)) {
      fg[, leadingEdge := vapply(leadingEdge, paste, collapse = "|", FUN.VALUE = character(1))]
      fg[, `:=`(variant = variant, method = method, engine = "fgsea", collection = collection)]
      all_fg[[paste(variant, method, collection)]] <- fg
    }
  }
}
fg_all <- rbindlist(all_fg, fill = TRUE)
fwrite_tsv(fg_all, file.path(DIR_GSEA_MULTI, "GSEA_fgsea_all_variants.tsv"))

# Independent engine replication on the pre-specified primary ranking.
primary <- fread(file.path(DIR_DE_MULTI, "tcga_matched", "voom_qw_full_results.tsv"))
primary_rank <- sort(setNames(primary$statistic, primary$gene_symbol), decreasing = TRUE)
cp_all <- rbindlist(lapply(names(sets), function(collection) {
  term2gene <- rbindlist(lapply(names(sets[[collection]]), function(nm)
    data.table(term = nm, gene = sets[[collection]][[nm]])))
  z <- suppressMessages(clusterProfiler::GSEA(primary_rank, TERM2GENE = term2gene,
    minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE, verbose = FALSE))
  out <- as.data.table(as.data.frame(z))
  if (nrow(out)) out[, `:=`(collection = collection, engine = "clusterProfiler",
                             variant = "tcga_matched", method = "voom_qw")]
  out
}), fill = TRUE)
fwrite_tsv(cp_all, file.path(DIR_GSEA_MULTI, "GSEA_clusterProfiler_primary.tsv"))

# Sample-level pathway scores (Hallmark + frozen panel) for each continuous matrix.
matrix_files <- c(raw = "raw_logcpm.rds", batch_corrected = "batch_corrected_logcpm.rds",
                  composition_adjusted = "composition_adjusted_logcpm.rds")
gsva_sets <- c(sets$Hallmark, sets$Panel30)
for (variant in names(matrix_files)) {
  x <- readRDS(file.path(DIR_VARIANTS, matrix_files[[variant]]))
  gs <- lapply(gsva_sets, intersect, y = rownames(x)); gs <- gs[lengths(gs) >= 10L]
  if (packageVersion("GSVA") >= "2.0.0") score <- gsva(gsvaParam(x, gs), verbose = FALSE)
  else score <- gsva(x, gs, method = "gsva", verbose = FALSE)
  saveRDS(score, file.path(DIR_GSEA_MULTI, paste0("GSVA_", variant, ".rds")))
}

comp_genes <- unique(unlist(list(
  muscle = c("ACTA1", "MYH7", "MYL1", "MYL2", "TNNT3", "CKM", "DES"),
  fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2"),
  immune = c("PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "LYZ", "FCER1G"))))
le_overlap <- fg_all[, .(leading_edge_n = lengths(strsplit(leadingEdge, "|", fixed = TRUE)),
  composition_overlap_n = lengths(lapply(strsplit(leadingEdge, "|", fixed = TRUE), intersect, y = comp_genes)),
  composition_overlap = vapply(lapply(strsplit(leadingEdge, "|", fixed = TRUE), intersect, y = comp_genes), paste, collapse = "|", FUN.VALUE = character(1))),
  by = .(variant, method, collection, pathway, NES, padj)]
fwrite_tsv(le_overlap, file.path(DIR_GSEA_MULTI, "leading_edge_composition_overlap.tsv"))
log_msg("Multiverse GSEA complete:", nrow(fg_all), "fgsea results")
