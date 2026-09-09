# Stratified bootstrap and pre-defined outlier-removal GSEA robustness.
suppressPackageStartupMessages({ library(here); library(data.table); library(limma); library(fgsea) })
source(here("scripts", "phase2", "20_upgrade_config.R"))

expr <- readRDS(file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
sets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))
sets <- lapply(sets, intersect, y = rownames(expr)); sets <- sets[lengths(sets) >= 10L]

rank_samples <- function(idx) {
  md <- meta[idx]; md[, condition := factor(condition, levels = c("Normal", "Tumor"))]
  design <- model.matrix(~ condition + source, md)
  assert_full_rank(design, "bootstrap design")
  fit <- eBayes(lmFit(expr[, idx, drop = FALSE], design), robust = TRUE)
  tt <- topTable(fit, coef = "conditionTumor", number = Inf, sort.by = "none")
  sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
}
run_gsea <- function(rank, replicate, analysis) {
  set.seed(SEED + replicate)
  z <- as.data.table(fgseaMultilevel(sets, rank, minSize = 10, maxSize = 1000,
                                     eps = 1e-10, nPermSimple = 1000))
  z[, `:=`(replicate = replicate, analysis = analysis)]
  z
}

set.seed(SEED)
boot <- rbindlist(lapply(seq_len(N_BOOTSTRAP), function(b) {
  idx <- unlist(lapply(split(seq_len(nrow(meta)), interaction(meta$condition, meta$source, drop = TRUE)),
                       function(i) sample(i, length(i), replace = TRUE)), use.names = FALSE)
  run_gsea(rank_samples(idx), b, "stratified_bootstrap")
}))
fwrite_tsv(boot, file.path(DIR_GSEA_MULTI, "GSEA_bootstrap_panel30.tsv"))

out <- fread(file.path(DIR_PREPROC, "sample_outliers.tsv"))
remove <- out[is_outlier == TRUE, sample]
idx_clean <- which(!meta$sample %in% remove)
clean <- run_gsea(rank_samples(idx_clean), 0L, "outliers_removed")
fwrite_tsv(clean, file.path(DIR_GSEA_MULTI, "GSEA_outliers_removed_panel30.tsv"))

rob <- boot[, .(direction_fraction = max(mean(NES > 0), mean(NES < 0)),
                significant_fraction = mean(padj < FDR_THRESH, na.rm = TRUE),
                median_NES = median(NES, na.rm = TRUE)), by = pathway]
rob[, robust := direction_fraction >= GSEA_BOOTSTRAP_DIRECTION_FRACTION &
                 significant_fraction >= GSEA_BOOTSTRAP_FDR_FRACTION]
rob <- merge(rob, clean[, .(pathway, outlier_removed_NES = NES, outlier_removed_FDR = padj)], by = "pathway", all.x = TRUE)
fwrite_tsv(rob, file.path(DIR_GSEA_MULTI, "GSEA_robustness_summary.tsv"))
log_msg("GSEA robustness complete:", N_BOOTSTRAP, "bootstraps;", sum(rob$robust), "robust pathways")
