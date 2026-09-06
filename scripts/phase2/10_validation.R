# ═══════════════════════════════════════════════════════════════════════════════
# 10_validation.R — Validação externa (GSE33630, GSE60542, GSE224356)
#
# GSE33630 (PTC vs normal pareado) e GSE60542 (PTC vs normal thyroid) usam GPL570
# (hgu133plus2.db). GSE224356 usa listas de DEGs pré-computadas (xlsx).
# GSE224357 não é tratada como coorte independente (mesmo experimento/SuperSeries).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(limma)
  library(readxl)
})

log_msg("══ Validação externa ══")

CANDIDATES <- c("ITGA2","FN1","CCND1","MYH7","TPO","TG","EPCAM","KRT19","CDH1","VIM","COL1A1","COL3A1")

# ── helpers ────────────────────────────────────────────────────────────────────
parse_geo_matrix <- function(path) {
  lines <- readLines(gzfile(path))
  hdr_idx <- which(grepl("ID_REF", lines, fixed = TRUE) & !grepl("^!", lines))[1]
  hdr <- gsub('"', '', lines[hdr_idx])
  sample_ids <- strsplit(hdr, "\t")[[1]][-1]
  dat_idx <- which(!grepl("^!", lines) & nzchar(lines))
  dat_idx <- dat_idx[dat_idx > hdr_idx]
  dat_lines <- lines[dat_idx]
  # source_name
  src_lines <- lines[grep("^!Sample_source_name_ch1", lines)]
  source <- strsplit(gsub('"', '', src_lines), "\t")[[1]][-1]
  # characteristics (todas as linhas)
  ch_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
  char <- lapply(ch_lines, function(l) strsplit(gsub('"', '', l), "\t")[[1]][-1])
  names(source) <- sample_ids
  char <- lapply(char, function(x) { names(x) <- sample_ids; x })
  list(sample_ids = sample_ids, source = source, char = char, dat_lines = dat_lines)
}

compute_candidate_de <- function(mat, sym, cond, dataset) {
  gm <- rowMeans(mat, na.rm = TRUE)
  dt <- data.table(sym = sym, mean = gm, idx = seq_along(sym))
  setorder(dt, -mean)
  dt <- dt[!duplicated(sym) & !is.na(sym)]
  mat_g <- mat[dt$idx, , drop = FALSE]
  rownames(mat_g) <- dt$sym
  design <- model.matrix(~ cond)
  fit <- eBayes(lmFit(mat_g, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[tt$gene %in% CANDIDATES, ]
  data.table(gene = tt$gene, dataset = dataset,
             logFC_external = tt$logFC, pval = tt$P.Value,
             fdr = p.adjust(tt$P.Value, method = "BH"),
             direction_external = ifelse(tt$logFC > 0, "Up", "Down"),
             n_normal = sum(cond == levels(cond)[1]),
             n_tumor = sum(cond == levels(cond)[2]))
}

build_matrix <- function(g) {
  mat <- do.call(rbind, lapply(g$dat_lines, function(l) {
    x <- strsplit(gsub('"', '', l), "\t")[[1]]
    suppressWarnings(as.numeric(x[-1]))
  }))
  probe <- vapply(g$dat_lines, function(l) strsplit(gsub('"', '', l), "\t")[[1]][1], character(1))
  rownames(mat) <- probe
  colnames(mat) <- g$sample_ids
  list(mat = mat, probe = probe)
}

# ── GSE33630 (PTC vs normal) ──────────────────────────────────────────────────
log_msg("  Processando GSE33630...")
g1 <- parse_geo_matrix(file.path(DIR_EXTERNAL, "GSE33630", "GSE33630_matrix.txt.gz"))
b1 <- build_matrix(g1)
src1 <- g1$source
# condição a partir de characteristics (pathological diagnostic) + source
diag1 <- g1$char[[which(vapply(g1$char, function(x) any(grepl("diagnostic|pathological", x, ignore.case = TRUE)), logical(1)))[1]]]
diag1 <- if (is.null(diag1)) src1 else diag1
is_normal1 <- grepl("non-tumor|normal", src1, ignore.case = TRUE)
is_ptc1 <- grepl("papillary", diag1, ignore.case = TRUE) & !is_normal1
cond1 <- ifelse(is_ptc1, "PTC", ifelse(is_normal1, "Normal", "Other"))
log_msg(sprintf("    GSE33630: PTC=%d Normal=%d Other=%d", sum(is_ptc1), sum(is_normal1), sum(!is_ptc1 & !is_normal1)))
sym1 <- mapIds(hgu133plus2.db, keys = rownames(b1$mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
keep1 <- !is.na(sym1)
samp1 <- cond1 %in% c("PTC", "Normal")
mat1k <- b1$mat[keep1, samp1]
sym1k <- sym1[keep1]
cond1k <- factor(cond1[samp1], levels = c("Normal", "PTC"))
gse33630 <- compute_candidate_de(mat1k, sym1k, cond1k, "GSE33630")

# ── GSE60542 (PTC vs normal thyroid) ──────────────────────────────────────────
log_msg("  Processando GSE60542...")
g2 <- parse_geo_matrix(file.path(DIR_EXTERNAL, "GSE60542", "GSE60542_matrix.txt.gz"))
b2 <- build_matrix(g2)
src2 <- g2$source
is_normal2 <- grepl("Normal thyroid", src2, fixed = TRUE)
is_tumor2 <- grepl("Papillary thyroid carcinoma", src2, fixed = TRUE)
cond2 <- ifelse(is_tumor2, "Tumor", ifelse(is_normal2, "Normal", "Other"))
log_msg(sprintf("    GSE60542: Tumor=%d Normal=%d Other=%d", sum(is_tumor2), sum(is_normal2), sum(!is_tumor2 & !is_normal2)))
sym2 <- mapIds(hgu133plus2.db, keys = rownames(b2$mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
keep2 <- !is.na(sym2)
samp2 <- cond2 %in% c("Tumor", "Normal")
mat2k <- b2$mat[keep2, samp2]
sym2k <- sym2[keep2]
cond2k <- factor(cond2[samp2], levels = c("Normal", "Tumor"))
gse60542 <- compute_candidate_de(mat2k, sym2k, cond2k, "GSE60542")

# ── GSE224356 (listas pré-computadas) ─────────────────────────────────────────
log_msg("  Processando GSE224356 (listas de DEGs)...")
gse224356_files <- list.files(file.path(DIR_EXTERNAL, "GSE224356"), pattern = "\\.xlsx$", full.names = TRUE)
gse224356_rows <- list()
for (f in gse224356_files) {
  nm <- sub("GSE224356_", "", basename(f)); nm <- sub("\\.xlsx$", "", nm)
  parts <- strsplit(nm, "_")[[1]]
  comparison <- paste0(parts[1], "_", parts[2])
  direction <- parts[3]
  genes <- suppressMessages(read_excel(f))
  gene_col <- intersect(c("Other Gene ID","Gene","Gene Symbol","Symbol","gene","GENE","SYMBOL","GeneSymbol"), colnames(genes))[1]
  if (is.na(gene_col)) gene_col <- colnames(genes)[1]
  glist <- toupper(unique(as.character(genes[[gene_col]])))
  for (cg in CANDIDATES) {
    if (toupper(cg) %in% glist) {
      gse224356_rows[[length(gse224356_rows) + 1]] <-
        data.table(gene = cg, comparison = comparison, direction = direction, present = TRUE)
    }
  }
}
gse224356 <- if (length(gse224356_rows) > 0) rbindlist(gse224356_rows) else
  data.table(gene = character(), comparison = character(), direction = character(), present = logical())

# ── Consolidação ──────────────────────────────────────────────────────────────
val_all <- rbindlist(list(gse33630, gse60542), fill = TRUE)
fwrite_tsv(val_all, file.path(DIR_VAL, "validation_candidates.tsv"))
fwrite_tsv(gse224356, file.path(DIR_VAL, "GSE224356_candidates.tsv"))

limma <- fread(file.path(DIR_DE, "limma_full_results.tsv"))
dir_main <- setNames(sign(limma$logFC), limma$gene_symbol)
concord <- val_all[, .(gene, dataset, logFC_external, direction_external, pval, fdr)]
concord[, direction_main := dir_main[gene]]
concord[, concordant := sign(logFC_external) == direction_main]
fwrite_tsv(concord, file.path(DIR_VAL, "validation_concordance.tsv"))

log_msg("══ Validação externa concluída ══")
cat("\n=== Validação externa (candidatos) ===\n")
print(val_all[, .(gene, dataset, n_normal, n_tumor, logFC_external, pval, fdr, direction_external)])
cat("\n=== GSE224356 presença dos candidatos ===\n")
print(gse224356)
