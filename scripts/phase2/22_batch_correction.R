# Build raw and batch-corrected matrices with an identifiable condition/source design.
suppressPackageStartupMessages({
  library(here); library(data.table); library(edgeR); library(limma); library(sva); library(BiocParallel)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

tcga <- readRDS(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"))
tcga_meta <- fread(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv"))
old <- readRDS(file.path(DIR_DATAOUT, "counts_matrix.rds"))
gtex <- old[, grepl("^GTEX", colnames(old)), drop = FALSE]
common <- intersect(rownames(tcga), rownames(gtex))
if (length(common) < 15000L) stop("Insufficient common recount3 genes: ", length(common))

raw <- cbind(tcga[common, , drop = FALSE], gtex[common, , drop = FALSE])
meta <- rbind(
  tcga_meta[, .(sample, condition = as.character(condition), source, cohort)],
  data.table(sample = colnames(gtex), condition = "Normal", source = "GTEx", cohort = "GTEx_recount3")
)
meta[, condition := factor(condition, levels = c("Normal", "Tumor"))]
meta[, source := factor(source, levels = c("GTEx", "TCGA"))]
stopifnot(identical(colnames(raw), meta$sample))
design <- model.matrix(~ condition + source, meta)
assert_full_rank(design, "condition + source design")

dge <- DGEList(raw)
keep <- keep_by_cpm(raw, meta$condition, min_cpm = MIN_EXPR_CPM, fraction = EXPR_FRAC)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
raw_f <- dge$counts
logcpm <- cpm(dge, log = TRUE, prior.count = 2)

# removeBatchEffect is a continuous-scale sensitivity; condition is preserved.
batch_logcpm <- removeBatchEffect(logcpm, batch = meta$source,
                                  design = model.matrix(~ condition, meta))

# ComBat-seq supplies a count-scale sensitivity. It is valid here only because
# TCGA normals make source and condition partially separable.
set.seed(SEED)
combat_workers <- suppressWarnings(as.integer(Sys.getenv("THYROID_COMBAT_WORKERS", "2")))
if (!is.finite(combat_workers) || combat_workers < 1L) combat_workers <- 1L
combat_workers <- min(combat_workers, nrow(raw_f))
gene_cost <- pmax(1, rowSums(raw_f > 0))
assignment <- integer(nrow(raw_f))
estimated_load <- numeric(combat_workers)
for (i in order(gene_cost, decreasing = TRUE, method = "radix")) {
  worker <- which.min(estimated_load)
  assignment[i] <- worker
  estimated_load[worker] <- estimated_load[worker] + gene_cost[i]
}
chunks <- split(seq_len(nrow(raw_f)), assignment)
combat_key <- digest::digest(list(raw = matrix_checksum(raw_f), batch = meta$source,
  group = meta$condition, package = as.character(packageVersion("sva")),
  partition = "greedy_nonzero_v1", shrink = FALSE), algo = "sha256")
combat_cache <- file.path(DIR_CACHE, "ComBat_seq", paste0("sva-", packageVersion("sva")),
                          paste0(combat_key, ".rds"))
combat_counts <- safe_api_cache(combat_cache, function() {
  count_parts <- lapply(chunks, function(i) as.matrix(round(raw_f[i, , drop = FALSE])))
  bp <- if (combat_workers > 1L && .Platform$OS.type == "windows") {
    SnowParam(combat_workers, type = "SOCK", progressbar = TRUE, RNGseed = SEED)
  } else if (combat_workers > 1L) {
    MulticoreParam(combat_workers, progressbar = TRUE, RNGseed = SEED)
  } else {
    SerialParam(progressbar = TRUE, RNGseed = SEED)
  }
  parts <- bplapply(count_parts, function(part, batch, group) {
    sva::ComBat_seq(part, batch = batch, group = group, full_mod = TRUE)
  }, batch = meta$source, group = meta$condition, BPPARAM = bp)
  out <- do.call(rbind, parts)
  out[rownames(raw_f), , drop = FALSE]
}, validate = function(x) is.matrix(x) && identical(dimnames(x), dimnames(raw_f)), attempts = 1L)
storage.mode(combat_counts) <- "integer"

# Surrogate variables are saved for model-based sensitivity, not subtracted blindly.
mod <- model.matrix(~ condition + source, meta)
mod0 <- model.matrix(~ source, meta)
n_sv_estimated <- tryCatch(sva::num.sv(logcpm, mod, method = "be", B = 20, seed = SEED),
  error = function(e) 0L)
if (!is.finite(n_sv_estimated) || n_sv_estimated < 0L) n_sv_estimated <- 0L
n_sv <- min(n_sv_estimated, MAX_SURROGATE_VARIABLES)
sv <- if (n_sv > 0L) sva::sva(logcpm, mod, mod0, n.sv = n_sv)$sv else matrix(nrow = ncol(logcpm), ncol = 0L)
rownames(sv) <- meta$sample
if (ncol(sv)) colnames(sv) <- sprintf("SV%02d", seq_len(ncol(sv)))

saveRDS(raw_f, file.path(DIR_VARIANTS, "raw_counts.rds"))
saveRDS(logcpm, file.path(DIR_VARIANTS, "raw_logcpm.rds"))
saveRDS(combat_counts, file.path(DIR_VARIANTS, "batch_corrected_counts.rds"))
saveRDS(batch_logcpm, file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
saveRDS(sv, file.path(DIR_VARIANTS, "surrogate_variables.rds"))
fwrite_tsv(meta, file.path(DIR_VARIANTS, "analysis_metadata.tsv"))

pca_r2 <- function(x, md) {
  pc <- prcomp(t(x[order(matrixStats::rowVars(x), decreasing = TRUE)[seq_len(min(2000L, nrow(x)))], ]), scale. = TRUE)$x[, 1:10, drop = FALSE]
  rbindlist(lapply(seq_len(ncol(pc)), function(i) {
    data.table(PC = i,
      condition_R2 = summary(lm(pc[, i] ~ md$condition))$r.squared,
      source_R2 = summary(lm(pc[, i] ~ md$source))$r.squared)
  }))
}
diag <- rbind(data.table(matrix = "raw", pca_r2(logcpm, meta)),
              data.table(matrix = "batch_corrected", pca_r2(batch_logcpm, meta)))
fwrite_tsv(diag, file.path(DIR_VARIANTS, "batch_diagnostics.tsv"))
write_manifest(file.path(DIR_VARIANTS, "batch_correction_manifest.tsv"), list(
  filter = sprintf("CPM >= %s in >= %s samples (25%% of smaller condition group)",
                   MIN_EXPR_CPM, ceiling(EXPR_FRAC * min(table(meta$condition)))),
  genes_retained = nrow(raw_f), combat_workers = combat_workers,
  partition = "deterministic greedy balance by nonzero sample count",
  estimated_worker_load = paste(estimated_load, collapse = "|"),
  combat_cache_key = combat_key,
  combat_shrinkage = FALSE, condition_preserved = TRUE,
  n_sv_method = "Buja-Eyuboglu B=20 seed=42", n_sv_estimated = n_sv_estimated,
  n_sv_retained = n_sv, n_sv_safety_cap = MAX_SURROGATE_VARIABLES,
  raw_matrix_fingerprint_sha256 = matrix_checksum(raw_f)
))
log_msg("Batch correction complete;", nrow(raw_f), "genes; design rank",
        qr(design)$rank, "/", ncol(design), "; n.sv=", n_sv,
        "of", n_sv_estimated, "estimated")
