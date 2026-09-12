# Three evidence-filtered PPI networks; centrality is descriptive only.
suppressPackageStartupMessages({
  library(here); library(data.table); library(httr); library(igraph); library(clusterProfiler); library(msigdbr)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

de <- fread(file.path(DIR_DE_MULTI, "batch_corrected", "voom_qw_full_results.tsv"))
up <- head(de[logFC > 0][order(-abs(statistic)), gene_symbol], 250L)
down <- head(de[logFC < 0][order(-abs(statistic)), gene_symbol], 250L)
panel_sets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))
gsea <- fread(file.path(DIR_GSEA_MULTI, "GSEA_fgsea_all_variants.tsv"))
robust_pathways <- gsea[method == "voom_qw" & variant %in% c("tcga_matched", "batch_corrected", "composition_adjusted"),
  .(n_sig = sum(padj < FDR_THRESH, na.rm = TRUE), direction_consistent = uniqueN(sign(NES)) == 1L), by = pathway][n_sig == 3L & direction_consistent, pathway]
leading <- gsea[variant == "tcga_matched" & method == "voom_qw" & pathway %in% robust_pathways,
                unique(unlist(strsplit(leadingEdge, "|", fixed = TRUE)))]
universes <- list(genome_wide_top500 = unique(c(up, down)), panel30 = unique(unlist(panel_sets)),
                  robust_leading_edge = leading)

expr <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
tumor_expr <- expr[, meta$condition == "Tumor", drop = FALSE]

query_string <- function(genes, network_name) {
  key <- digest::digest(sort(unique(genes)), algo = "sha256")
  cache <- file.path(DIR_CACHE, "STRING", "v12.0", paste0(network_name, "_", key, ".rds"))
  safe_api_cache(cache, function() {
    ids <- unique(genes)
    # STRING network endpoint caps requests at 2000 proteins; chunk and merge to
    # preserve the a-priori universe definitions (panel30/leading edge can exceed 2000).
    chunks <- split(ids, ceiling(seq_along(ids) / 1900))
    rbindlist(lapply(chunks, function(chunk) {
      r <- POST("https://string-db.org/api/tsv/network", body = list(
        identifiers = paste(chunk, collapse = "\n"), species = 9606,
        required_score = STRING_SCORE, network_type = "physical",
        caller_identity = "thyroid-volcano-ppi-v5"), encode = "form", timeout(180))
      stop_for_status(r)
      fread(content(r, as = "text", encoding = "UTF-8"))
    }), fill = TRUE)
  }, validate = function(x) is.data.frame(x) && all(c("preferredName_A", "preferredName_B", "score") %in% names(x)))
}

go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
term2gene <- unique(go[, c("gs_name", "gene_symbol")])
setnames(term2gene, c("term", "gene"))

network_summary <- list()
for (nm in names(universes)) {
  genes <- intersect(universes[[nm]], rownames(tumor_expr))
  if (length(genes) < 10L) {
    network_summary[[nm]] <- data.table(network = nm, input_genes = length(genes),
      experimental_edges = 0L, coexpressed_edges = 0L, nodes = 0L,
      status = "too_few_genes_for_network")
    next
  }
  edges <- query_string(genes, nm)
  exp_col <- intersect(c("escore", "experimental", "experiments"), names(edges))
  if (length(exp_col)) exp_col <- exp_col[1]
  if (!length(exp_col)) stop("STRING response lacks experimental evidence channel")
  edges <- edges[get(exp_col) > 0]
  experimental_edges <- nrow(edges)
  if (!experimental_edges) {
    network_summary[[nm]] <- data.table(network = nm, input_genes = length(genes),
      experimental_edges = 0L, coexpressed_edges = 0L, nodes = 0L,
      status = "no_experimental_edges")
    next
  }
  edges[, `:=`(gene_A = preferredName_A, gene_B = preferredName_B)]
  edges <- edges[gene_A %in% rownames(tumor_expr) & gene_B %in% rownames(tumor_expr)]
  cors <- rbindlist(lapply(seq_len(nrow(edges)), function(i) {
    z <- tryCatch(suppressWarnings(cor.test(tumor_expr[edges$gene_A[i], ],
      tumor_expr[edges$gene_B[i], ], method = "spearman", exact = FALSE)),
      error = function(e) NULL)
    if (is.null(z)) data.table(rho = NA_real_, p = NA_real_)
    else data.table(rho = unname(z$estimate), p = z$p.value)
  }))
  edges <- cbind(edges, cors); edges[, coexpression_fdr := p.adjust(p, "BH")]
  edges <- edges[is.finite(rho) & abs(rho) >= COEXPRESSION_RHO & coexpression_fdr < FDR_THRESH]
  if (!nrow(edges)) {
    fwrite_tsv(data.table(network = nm, status = "no edges after experimental+coexpression filters"),
               file.path(DIR_PPI_MULTI, paste0(nm, "_status.tsv")))
    network_summary[[nm]] <- data.table(network = nm, input_genes = length(genes),
      experimental_edges = experimental_edges, coexpressed_edges = 0L, nodes = 0L,
      status = "no_coexpressed_edges")
    next
  }
  g <- igraph::simplify(graph_from_data_frame(edges[, .(gene_A, gene_B, score, rho)], directed = FALSE),
    remove.multiple = TRUE, remove.loops = TRUE,
    edge.attr.comb = list(score = "max", rho = "mean", .default = "ignore"))
  cent <- data.table(network = nm, gene = V(g)$name, degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE), closeness = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g)$vector,
    interpretation = "topological descriptor; not therapeutic evidence")
  set.seed(SEED)
  algorithms <- list(walktrap = cluster_walktrap(g), louvain = cluster_louvain(g),
                     leiden = cluster_leiden(g, objective_function = "modularity"))
  modules <- rbindlist(lapply(names(algorithms), function(a)
    data.table(network = nm, algorithm = a, gene = V(g)$name,
               module = as.integer(membership(algorithms[[a]])))))
  enrich <- modules[, {
    z <- tryCatch(enricher(gene, TERM2GENE = term2gene, universe = genes,
                           pvalueCutoff = 1, qvalueCutoff = 1), error = function(e) NULL)
    if (is.null(z)) data.table() else head(as.data.table(as.data.frame(z)), 10L)
  }, by = .(network, algorithm, module)]
  fwrite_tsv(edges, file.path(DIR_PPI_MULTI, paste0(nm, "_edges.tsv")))
  fwrite_tsv(cent, file.path(DIR_PPI_MULTI, paste0(nm, "_centrality.tsv")))
  fwrite_tsv(modules, file.path(DIR_PPI_MULTI, paste0(nm, "_modules.tsv")))
  fwrite_tsv(enrich, file.path(DIR_PPI_MULTI, paste0(nm, "_module_enrichment.tsv")))
  saveRDS(g, file.path(DIR_PPI_MULTI, paste0(nm, "_graph.rds")))
  network_summary[[nm]] <- data.table(network = nm, input_genes = length(genes),
    experimental_edges = experimental_edges, coexpressed_edges = ecount(g),
    nodes = vcount(g), status = "completed")
}
fwrite_tsv(rbindlist(network_summary, fill = TRUE),
  file.path(DIR_PPI_MULTI, "PPI_upgrade_summary.tsv"))
log_msg("Three PPI universes processed with experimental and tumor-only coexpression filters")
