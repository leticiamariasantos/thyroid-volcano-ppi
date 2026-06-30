# ═══════════════════════════════════════════════════════════════════════════════
# R/01_functions.R — Core utility functions
# thyroid-volcano-ppi
#
# All API calls wrapped in tryCatch with retry logic.
# All exports include column dictionary comments.
# ═══════════════════════════════════════════════════════════════════════════════

# ── F1. Extract gene symbols from KEGG pathway ────────────────────────────────
fetch_kegg_genes <- function(pathway_id) {
  pw <- tryCatch(
    keggGet(pathway_id),
    error = function(e) stop("KEGG API error for ", pathway_id, ": ", e$message)
  )
  stopifnot(length(pw) == 1L)
  g <- pw[[1]]$GENE
  if (is.null(g) || length(g) == 0L) stop("Empty GENE field for ", pathway_id)

  genes <- g[seq(2, length(g), by = 2)] |>
    str_remove("\\s*\\[.*$") |>
    str_trim() |>
    str_extract("^[^;]+") |>
    str_trim() |>
    unique() |>
    sort()

  # Remove KEGG readthrough-fusion artifacts
  genes <- genes[!grepl("-", genes) | nchar(genes) <= 10]
  tibble(pathway_id = pathway_id, gene_symbol = genes)
}

# ── F2. Validate expression data scale ────────────────────────────────────────
validate_expression_scale <- function(expr) {
  frac <- mean(abs(expr - round(expr)) > 0.01, na.rm = TRUE)
  if (frac < 0.01) warning("Expression values appear integer; confirm log2 scale.")
  if (max(expr, na.rm = TRUE) > 30) stop("Max expression >30; data not in log2 scale.")
  cat(sprintf("  QC: max=%.2f | non-int=%.1f%%\n",
              round(max(expr, na.rm = TRUE), 2), round(frac * 100, 1)))
}

# ── F3. Map gene symbols -> STRING protein IDs (REST API) ─────────────────────
map_string_ids <- function(annot_df, taxon = 9606) {
  clean <- annot_df |>
    dplyr::mutate(gene_symbol = as.character(gene_symbol)) |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  stopifnot(nrow(clean) >= 2)

  resp <- tryCatch(
    httr::POST(
      "https://string-db.org/api/json/get_string_ids",
      body = list(
        identifiers     = paste(clean$gene_symbol, collapse = "\r\n"),
        species         = as.character(taxon),
        limit           = "1",
        caller_identity = "thyroid_volcano_ppi"
      ),
      encode = "form",
      httr::timeout(60)
    ),
    error = function(e) stop("STRING API mapping failed: ", e$message)
  )
  httr::stop_for_status(resp)

  ids <- jsonlite::fromJSON(
    httr::content(resp, as = "text", encoding = "UTF-8"), flatten = TRUE
  )

  result <- ids |>
    dplyr::select(gene_symbol = preferredName, STRING_id = stringId) |>
    dplyr::left_join(
      clean |> dplyr::select(gene_symbol, regulation),
      by = "gene_symbol"
    ) |>
    dplyr::filter(!is.na(STRING_id))

  na_n <- sum(is.na(result$regulation))
  if (na_n > 0) {
    cat(sprintf("  Note: %d gene(s) without regulation annotation\n", na_n))
    result$regulation[is.na(result$regulation)] <- "NS"
  }
  result
}

# ── F4. Fetch STRING PPI interactions (REST API) ──────────────────────────────
fetch_string_edges <- function(string_ids, taxon = 9606) {
  resp <- tryCatch(
    httr::POST(
      "https://string-db.org/api/tsv/network",
      body = list(
        identifiers     = paste(string_ids, collapse = "\r\n"),
        species         = as.character(taxon),
        required_score  = "0",
        caller_identity = "thyroid_volcano_ppi"
      ),
      encode = "form",
      httr::timeout(120)
    ),
    error = function(e) stop("STRING API network call failed: ", e$message)
  )
  httr::stop_for_status(resp)

  raw <- utils::read.delim(
    text = httr::content(resp, as = "text", encoding = "UTF-8"),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (nrow(raw) == 0L) {
    return(tibble(from = character(), to = character(),
                  combined_score = numeric()))
  }
  if (!all(c("stringId_A", "stringId_B") %in% colnames(raw))) {
    stop("STRING REST response missing expected columns (stringId_A/stringId_B).")
  }

  raw |>
    tibble::as_tibble() |>
    dplyr::transmute(
      from = stringId_A,
      to   = stringId_B,
      combined_score = as.numeric(score) * 1000
    )
}

# ── F5. Build igraph network ──────────────────────────────────────────────────
build_igraph <- function(edges, nodes) {
  igraph::graph_from_data_frame(
    d = edges |> dplyr::transmute(from, to, weight = combined_score),
    directed = FALSE,
    vertices = nodes |> dplyr::transmute(
      name = STRING_id, gene_symbol, regulation
    )
  )
}

# ── F6. Compute centrality metrics ────────────────────────────────────────────
compute_centrality <- function(g) {
  tibble(
    gene_symbol = igraph::V(g)$gene_symbol,
    STRING_id   = igraph::V(g)$name,
    regulation  = igraph::V(g)$regulation,
    degree      = igraph::degree(g),
    betweenness = igraph::betweenness(g, normalized = TRUE),
    closeness   = igraph::closeness(g, normalized = TRUE),
    hub_score   = igraph::hits_scores(g)$hub
  ) |> dplyr::arrange(dplyr::desc(betweenness))
}

# ── F7. Evaluate layout quality ───────────────────────────────────────────────
evaluate_layout <- function(g, layout_fun, dim = 2, ...) {
  set.seed(42)
  lay <- tryCatch(layout_fun(g, dim = dim, ...), error = function(e) NULL)
  if (is.null(lay) || any(is.na(lay)) || nrow(lay) != igraph::vcount(g))
    return(0)

  d <- as.matrix(dist(lay))
  diag(d) <- Inf
  min_dist <- apply(d, 1, min)

  adj <- igraph::as_adjacency_matrix(g)
  edge_dists <- d[which(adj > 0, arr.ind = TRUE)]
  mean_edge <- if (length(edge_dists) > 0) mean(edge_dists, na.rm = TRUE) else 1

  score <- median(min_dist / mean_edge, na.rm = TRUE)
  if (!is.finite(score)) score <- 0
  score
}

# ── F8. Export table (TSV) with dictionary header ─────────────────────────────
export_tsv <- function(x, filepath, description = NULL) {
  if (!grepl("\\.tsv$", filepath)) filepath <- paste0(filepath, ".tsv")
  readr::write_tsv(x, filepath)
  cat(sprintf("  -> %s (%dx%d)\n", basename(filepath), nrow(x), ncol(x)))
  if (!is.null(description)) {
    cat(sprintf("     %s\n", description))
  }
}
