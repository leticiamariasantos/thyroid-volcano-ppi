# Differential expression across tcga_matched/raw/batch/composition variants.
suppressPackageStartupMessages({
  library(here); library(data.table); library(edgeR); library(limma); library(DESeq2)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

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

run_variant <- function(counts, meta, design, variant, continuous = NULL) {
  stopifnot(identical(colnames(counts), meta$sample))
  assert_full_rank(design, paste(variant, "design"))
  coef_name <- grep("^conditionTumor$", colnames(design), value = TRUE)
  if (length(coef_name) != 1L) stop("conditionTumor coefficient unavailable for ", variant)
  dge <- DGEList(counts = round(counts))
  keep <- keep_by_cpm(counts, meta$condition, min_cpm = MIN_EXPR_CPM, fraction = EXPR_FRAC)
  dge <- normLibSizes(dge[keep, , keep.lib.sizes = FALSE])

  log_msg("DE", variant, "voom")
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  voom_res <- standardize(topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), "voom", variant)

  log_msg("DE", variant, "voom quality weights")
  set.seed(SEED)
  vqw <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit_qw <- eBayes(lmFit(vqw, design), robust = TRUE)
  vqw_res <- standardize(topTable(fit_qw, coef = coef_name, number = Inf, sort.by = "none"), "voom_quality_weights", variant)

  log_msg("DE", variant, "edgeR quasi-likelihood")
  y <- estimateDisp(dge, design, robust = TRUE)
  qlf <- glmQLFit(y, design, robust = TRUE)
  qlt <- glmQLFTest(qlf, coef = match(coef_name, colnames(design)))
  ql_tab <- topTags(qlt, n = Inf, sort.by = "none")$table
  ql_tab$P.Value <- ql_tab$PValue
  ql_tab$adj.P.Val <- p.adjust(ql_tab$PValue, "BH")
  ql_tab$statistic <- sign(ql_tab$logFC) * sqrt(pmax(ql_tab$F, 0))
  ql_res <- standardize(ql_tab, "edgeR_QL", variant)

  log_msg("DE", variant, "DESeq2")
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
  dds_res <- standardize(dr, "DESeq2", variant)

  log_msg("DE", variant, "limma continuous")
  if (is.null(continuous)) continuous <- v$E
  fit_l <- eBayes(lmFit(continuous[rownames(v$E), , drop = FALSE], design), robust = TRUE)
  limma_res <- standardize(topTable(fit_l, coef = coef_name, number = Inf, sort.by = "none"), "limma_continuous", variant)
  sample_weights <- vqw$targets$sample.weights
  if (is.null(sample_weights)) sample_weights <- rep(NA_real_, ncol(counts))
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
  cache_key <- digest::digest(list(
    variant = variant,
    inputs = setNames(vapply(input_paths, file_sha256, character(1)), basename(input_paths)),
    design = design,
    metadata = meta,
    script = file_sha256(here("scripts", "phase2", "24_multiverse_de.R")),
    package_versions = package_versions,
    seed = SEED,
    min_cpm = MIN_EXPR_CPM,
    expression_fraction = EXPR_FRAC
  ), algo = "sha256")
  cache_file <- file.path(DIR_CACHE, "multiverse_de", paste0("v", UPGRADE_VERSION),
                          paste0(variant, "-", cache_key, ".rds"))
  result <- safe_api_cache(cache_file, function() {
    log_msg("DE cache miss:", variant, cache_key)
    run_variant(counts, meta, design, variant, continuous)
  }, validate = function(x) {
    required <- c("voom", "voom_qw", "edgeR_QL", "DESeq2", "limma", "sample_weights")
    is.list(x) && all(required %in% names(x)) &&
      all(vapply(x[required[1:5]], function(tab) is.data.frame(tab) && nrow(tab) > 0L, logical(1))) &&
      length(x$sample_weights) == ncol(counts)
  }, attempts = 1L)
  log_msg("DE variant ready:", variant, cache_key)
  result
}

tcga_counts <- readRDS(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"))
tcga_meta <- fread(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv"))
tcga_meta[, condition := factor(condition, levels = c("Normal", "Tumor"))]
raw <- readRDS(file.path(DIR_VARIANTS, "raw_counts.rds"))
combat <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_counts.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
meta[, `:=`(condition = factor(condition, levels = c("Normal", "Tumor")),
            source = factor(source, levels = c("GTEx", "TCGA")))]
batch_log <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
composition_log <- readRDS(file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"))
comp <- readRDS(file.path(DIR_DECONV, "composition_consensus.rds"))
sv <- readRDS(file.path(DIR_VARIANTS, "surrogate_variables.rds"))

analyses <- list()
analyses$tcga_matched <- cached_variant("tcga_matched", tcga_counts, tcga_meta,
  model.matrix(~ condition, tcga_meta),
  c(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"),
    file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv")))
analyses$raw <- cached_variant("raw", raw, meta, model.matrix(~ condition + source, meta),
  c(file.path(DIR_VARIANTS, "raw_counts.rds"), file.path(DIR_VARIANTS, "analysis_metadata.tsv")))
analyses$batch_corrected <- cached_variant("batch_corrected", combat, meta,
  model.matrix(~ condition, meta),
  c(file.path(DIR_VARIANTS, "batch_corrected_counts.rds"),
    file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"),
    file.path(DIR_VARIANTS, "analysis_metadata.tsv")), batch_log)

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
analyses$composition_adjusted <- cached_variant("composition_adjusted", raw, meta_comp, design_comp,
  c(file.path(DIR_VARIANTS, "raw_counts.rds"),
    file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"),
    file.path(DIR_DECONV, "composition_consensus.rds"),
    file.path(DIR_VARIANTS, "surrogate_variables.rds"),
    file.path(DIR_VARIANTS, "analysis_metadata.tsv")), composition_log)

# Strict paired sensitivity: only participants with both 01 and 11.
paired <- tcga_meta[paired == TRUE]
paired_counts <- tcga_counts[, paired$sample, drop = FALSE]
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
