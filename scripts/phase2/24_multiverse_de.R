# Differential expression across tcga_matched/raw/batch/composition variants.
suppressPackageStartupMessages({
  library(here); library(data.table); library(edgeR); library(limma); library(DESeq2)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))
source(here("R", "quality_weights_equivalent.R"))

standardize <- function(tt, method, variant) {
  tt <- as.data.table(tt, keep.rownames = "gene_symbol")
  aliases <- c(log2FoldChange = "logFC", pvalue = "P.Value", padj = "adj.P.Val")
  for (old in names(aliases)) if (old %in% names(tt)) setnames(tt, old, aliases[[old]])
  if (!"statistic" %in% names(tt)) {
    stat_col <- intersect(c("t", "F", "stat"), names(tt))[1]
    tt[, statistic := if (length(stat_col)) get(stat_col) else logFC / pmax(1e-12, abs(logFC))]
  }
  tt[, `:=`(method = method, variant = variant,
            regulation = fifelse(!is.na(adj.P.Val) & adj.P.Val < FDR_THRESH & logFC >= LOGFC_THRESH, "Up",
                          fifelse(!is.na(adj.P.Val) & adj.P.Val < FDR_THRESH & logFC <= -LOGFC_THRESH, "Down", "NS")))]
  setorder(tt, P.Value)
  tt[, rank := .I]
  tt[]
}

run_variant <- function(counts, meta, design, variant, continuous = NULL, cache_key) {
  stopifnot(identical(colnames(counts), meta$sample))
  assert_full_rank(design, paste(variant, "design"))
  coef_name <- grep("^conditionTumor$", colnames(design), value = TRUE)
  if (length(coef_name) != 1L) stop("conditionTumor coefficient unavailable for ", variant)
  dge <- DGEList(counts = round(counts))
  keep <- keep_by_cpm(counts, meta$condition, min_cpm = MIN_EXPR_CPM, fraction = EXPR_FRAC)
  dge <- normLibSizes(dge[keep, , keep.lib.sizes = FALSE])
  valid_table <- function(x) is.data.frame(x) && nrow(x) == nrow(dge) &&
    all(c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "statistic", "variant") %in% names(x)) &&
    !anyDuplicated(x$gene_symbol) && setequal(x$gene_symbol, rownames(dge)) &&
    all(x$variant == variant) && all(is.finite(x$logFC)) &&
    all(is.na(x$P.Value) | (x$P.Value >= 0 & x$P.Value <= 1)) &&
    all(is.na(x$adj.P.Val) | (x$adj.P.Val >= 0 & x$adj.P.Val <= 1))
  checkpoint <- function(method, calculate, validate = valid_table) {
    path <- file.path(DIR_CACHE, "multiverse_de", paste0("v", UPGRADE_VERSION),
                      "methods", paste0(variant, "-", cache_key), paste0(method, ".rds"))
    began <- Sys.time()
    hit <- file.exists(path)
    value <- safe_api_cache(path, calculate, validate = validate, attempts = 1L)
    elapsed <- as.numeric(difftime(Sys.time(), began, units = "secs"))
    log_msg("DE method ready:", variant, method, sprintf("%.1f seconds", elapsed))
    fwrite_tsv(data.table(variant = variant, method = method, cache_key = cache_key,
      cache_present_at_start = hit, elapsed_seconds = elapsed,
      completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      sha256 = digest::digest(path, algo = "sha256", file = TRUE)),
      file.path(DIR_PIPELINE, paste0("DE_", variant, "_", method, "_benchmark.tsv")))
    value
  }

  log_msg("DE", variant, "voom")
  v <- voom(dge, design, plot = FALSE)
  voom_res <- checkpoint("voom", function() {
    fit <- eBayes(lmFit(v, design), robust = TRUE)
    standardize(topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), "voom", variant)
  })

  log_msg("DE", variant, "voom quality weights")
  qw <- checkpoint("voom_qw", function() {
    set.seed(SEED)
    vqw <- voom_quality_weights_equivalent(dge, design,
      progress = function(i, n) log_msg("DE", variant, "quality weights genes", i, "/", n))
    fit_qw <- eBayes(lmFit(vqw, design), robust = TRUE)
    list(table = standardize(topTable(fit_qw, coef = coef_name, number = Inf, sort.by = "none"),
                              "voom_quality_weights", variant),
         weights = setNames(vqw$targets$sample.weights, colnames(counts)))
  }, validate = function(x) is.list(x) && isTRUE(valid_table(x$table)) &&
       identical(names(x$weights), colnames(counts)) &&
       all(is.finite(x$weights) & x$weights > 0))
  vqw_res <- qw$table

  log_msg("DE", variant, "edgeR quasi-likelihood")
  ql_res <- checkpoint("edgeR_QL", function() {
  y <- estimateDisp(dge, design, robust = TRUE)
  qlf <- glmQLFit(y, design, robust = TRUE)
  qlt <- glmQLFTest(qlf, coef = match(coef_name, colnames(design)))
  ql_tab <- topTags(qlt, n = Inf, sort.by = "none")$table
  ql_tab$P.Value <- ql_tab$PValue
  ql_tab$adj.P.Val <- p.adjust(ql_tab$PValue, "BH")
  ql_tab$statistic <- sign(ql_tab$logFC) * sqrt(pmax(ql_tab$F, 0))
  standardize(ql_tab, "edgeR_QL", variant)
  })

  log_msg("DE", variant, "DESeq2")
  dds_res <- checkpoint("DESeq2", function() {
  coldata <- as.data.frame(meta)
  rownames(coldata) <- coldata$sample
  form_terms <- setdiff(colnames(design), "(Intercept)")
  # DESeq2 gets the same covariates through a formula reconstructed from metadata.
  rhs <- c(if (any(grepl("^participant", form_terms))) "participant" else NULL,
           if ("sourceTCGA" %in% form_terms) "source" else NULL,
           grep("^(comp_|sv_)", names(coldata), value = TRUE), "condition")
  dds <- DESeqDataSetFromMatrix(round(counts[keep, , drop = FALSE]), coldata,
                                as.formula(paste("~", paste(rhs, collapse = "+"))))
  dds <- DESeq(dds, quiet = TRUE)
  dr <- as.data.frame(results(dds, contrast = c("condition", "Tumor", "Normal")))
  dr$statistic <- dr$stat
  standardize(dr, "DESeq2", variant)
  })

  log_msg("DE", variant, "limma continuous")
  if (is.null(continuous)) continuous <- v$E
  limma_res <- checkpoint("limma", function() {
    fit_l <- eBayes(lmFit(continuous[rownames(v$E), , drop = FALSE], design), robust = TRUE)
    standardize(topTable(fit_l, coef = coef_name, number = Inf, sort.by = "none"), "limma_continuous", variant)
  })
  sample_weights <- qw$weights
  list(voom = voom_res, voom_qw = vqw_res, edgeR_QL = ql_res, DESeq2 = dds_res,
       limma = limma_res, sample_weights = setNames(sample_weights, colnames(counts)))
}

hash_memo <- new.env(parent = emptyenv())
file_sha256 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!exists(path, envir = hash_memo, inherits = FALSE)) {
    assign(path, digest::digest(path, algo = "sha256", file = TRUE), envir = hash_memo)
  }
  get(path, envir = hash_memo, inherits = FALSE)
}

cached_variant <- function(variant, counts, meta, design, input_paths, continuous = NULL) {
  package_versions <- vapply(c("edgeR", "limma", "DESeq2"), function(pkg) {
    as.character(packageVersion(pkg))
  }, character(1))
  contract <- list(
    variant = variant,
    inputs = setNames(vapply(input_paths, file_sha256, character(1)), basename(input_paths)),
    design = design,
    metadata = meta,
    script = file_sha256(here("scripts", "phase2", "24_multiverse_de.R")),
    package_versions = package_versions,
    seed = SEED,
    min_cpm = MIN_EXPR_CPM,
    expression_fraction = EXPR_FRAC
  )
  legacy_contract <- contract
  legacy_contract$script <- "f1ed30b423655d0541f6017cb36842a4ad8b7728b3e1eabd34533ac94e51d6af"
  legacy_key <- digest::digest(legacy_contract, algo = "sha256")
  contract$helpers <- vapply(c(here("R", "quality_weights_equivalent.R"),
    here("R", "upgrade_utils.R"), here("scripts", "phase2", "00_config.R"),
    here("scripts", "phase2", "20_upgrade_config.R"), here("renv.lock")), file_sha256, character(1))
  cache_key <- digest::digest(contract, algo = "sha256")
  cache_file <- file.path(DIR_CACHE, "multiverse_de", paste0("v", UPGRADE_VERSION),
                          paste0(variant, "-", cache_key, ".rds"))
  validate <- function(x) {
    required <- c("voom", "voom_qw", "edgeR_QL", "DESeq2", "limma", "sample_weights")
    is.list(x) && all(required %in% names(x)) &&
      all(vapply(x[required[1:5]], function(tab) is.data.frame(tab) && nrow(tab) > 0L, logical(1))) &&
      all(vapply(x[required[1:5]], function(tab) !anyDuplicated(tab$gene_symbol) &&
        all(tab$gene_symbol %in% rownames(counts)) && all(tab$variant == variant), logical(1))) &&
      identical(names(x$sample_weights), colnames(counts)) &&
      all(is.finite(x$sample_weights) & x$sample_weights > 0)
  }
  legacy_file <- file.path(dirname(cache_file), paste0(variant, "-", legacy_key, ".rds"))
  result <- safe_api_cache(cache_file, function() {
    if (file.exists(legacy_file)) {
      original <- readRDS(legacy_file)
      if (!isTRUE(validate(original))) stop("Invalid legacy DE cache: ", legacy_file)
      log_msg("DE validated reference cache reused:", variant, legacy_key)
      return(original)
    }
    log_msg("DE cache miss:", variant, cache_key)
    run_variant(counts, meta, design, variant, continuous, cache_key)
  }, validate = validate, attempts = 1L)
  fwrite_tsv(data.table(variant = variant, cache_key = cache_key,
    cache_sha256 = digest::digest(cache_file, algo = "sha256", file = TRUE),
    legacy_key = if (file.exists(legacy_file)) legacy_key else NA_character_,
    engine = if (file.exists(legacy_file)) "official_reference_cache" else "equivalent_linear_algebra"),
    file.path(DIR_PIPELINE, paste0("DE_", variant, "_cache_provenance.tsv")))
  log_msg("DE variant ready:", variant, cache_key)
  result
}

tcga_counts <- readRDS(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"))
tcga_meta <- fread(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv"))
tcga_meta[, condition := factor(condition, levels = c("Normal", "Tumor"))]

analyses <- list()
analyses$tcga_matched <- cached_variant("tcga_matched", tcga_counts, tcga_meta,
  model.matrix(~ condition, tcga_meta),
  c(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"),
    file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv")))
rm(tcga_counts)
gc(verbose = FALSE)

# Load only the matrices needed by the active variant on memory-limited hosts.
raw <- readRDS(file.path(DIR_VARIANTS, "raw_counts.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
meta[, `:=`(condition = factor(condition, levels = c("Normal", "Tumor")),
            source = factor(source, levels = c("GTEx", "TCGA")))]
analyses$raw <- cached_variant("raw", raw, meta, model.matrix(~ condition + source, meta),
  c(file.path(DIR_VARIANTS, "raw_counts.rds"), file.path(DIR_VARIANTS, "analysis_metadata.tsv")))
rm(raw)
gc(verbose = FALSE)
combat <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_counts.rds"))
batch_log <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
analyses$batch_corrected <- cached_variant("batch_corrected", combat, meta,
  model.matrix(~ condition, meta),
  c(file.path(DIR_VARIANTS, "batch_corrected_counts.rds"),
    file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"),
    file.path(DIR_VARIANTS, "analysis_metadata.tsv")), batch_log)
rm(combat, batch_log)
gc(verbose = FALSE)

comp <- readRDS(file.path(DIR_DECONV, "composition_consensus.rds"))
sv <- readRDS(file.path(DIR_VARIANTS, "surrogate_variables.rds"))
comp2 <- as.data.frame(comp[match(meta$sample, rownames(comp)), , drop = FALSE])
names(comp2) <- paste0("comp_", make.names(names(comp2)))
meta_comp <- cbind(as.data.frame(meta), comp2)
if (ncol(sv)) {
  sv2 <- as.data.frame(sv[match(meta$sample, rownames(sv)), , drop = FALSE])
  names(sv2) <- paste0("sv_", seq_len(ncol(sv2)))
  meta_comp <- cbind(meta_comp, sv2)
}
design_comp <- model.matrix(as.formula(paste("~ condition + source +",
  paste(c(names(comp2), grep("^sv_", names(meta_comp), value = TRUE)), collapse = " + "))), meta_comp)
raw <- readRDS(file.path(DIR_VARIANTS, "raw_counts.rds"))
composition_log <- readRDS(file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"))
analyses$composition_adjusted <- cached_variant("composition_adjusted", raw, meta_comp, design_comp,
  c(file.path(DIR_VARIANTS, "raw_counts.rds"),
    file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"),
    file.path(DIR_DECONV, "composition_consensus.rds"),
    file.path(DIR_VARIANTS, "surrogate_variables.rds"),
    file.path(DIR_VARIANTS, "analysis_metadata.tsv")), composition_log)
rm(raw, composition_log)
gc(verbose = FALSE)

# Strict paired sensitivity: only participants with both 01 and 11.
paired <- tcga_meta[paired == TRUE]
tcga_counts <- readRDS(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"))
paired_counts <- tcga_counts[, paired$sample, drop = FALSE]
rm(tcga_counts)
gc(verbose = FALSE)
paired[, participant := factor(participant)]
analyses$tcga_paired <- cached_variant("tcga_paired", paired_counts, paired,
  model.matrix(~ participant + condition, paired),
  c(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"),
    file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv")))

summary <- rbindlist(lapply(names(analyses), function(vn) {
  rbindlist(lapply(analyses[[vn]][1:5], function(x)
    data.table(variant = vn, method = unique(x$method), genes_tested = nrow(x),
               total_deg = sum(x$regulation != "NS"), up = sum(x$regulation == "Up"),
               down = sum(x$regulation == "Down"))))
}))
fwrite_tsv(summary, file.path(DIR_DE_MULTI, "DEG_summary_all_variants.tsv"))

for (vn in names(analyses)) {
  dir.create(file.path(DIR_DE_MULTI, vn), recursive = TRUE, showWarnings = FALSE)
  for (mn in names(analyses[[vn]][1:5])) {
    fwrite_tsv(analyses[[vn]][[mn]], file.path(DIR_DE_MULTI, vn, paste0(mn, "_full_results.tsv")))
  }
  saveRDS(analyses[[vn]]$sample_weights, file.path(DIR_DE_MULTI, vn, "sample_weights.rds"))
}
log_msg("Multiverse DE complete:", nrow(summary), "variant-method combinations")
