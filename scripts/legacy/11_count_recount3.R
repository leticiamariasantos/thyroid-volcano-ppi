# ═══════════════════════════════════════════════════════════════════════════════
# 11_count_recount3.R — Count matrix (raw gene-level counts) via recount3
#
# MOTIVAÇÃO: o arquivo Xena `TcgaTargetGtex_rsem_gene_count.gz` retornou
# HTTP 403 (AccessDenied). Fonte alternativa legítima para AS MESMAS amostras:
# recount3 (https://github.com/LieberInstitute/recount3) fornece gene-level
# raw counts (STAR, annotation G026) para TCGA THCA e GTEx THYROID.
#
# PROVENIÊNCIA (registrada):
#   - TCGA THCA : project "THCA", file_source "tcga"  -> 572 amostras
#   - GTEx THYROID : project "THYROID", file_source "gtex" -> 706 amostras
#   - Assay: "raw_counts" (somas de contagens por gene, inteiros)
#   - Anotação: GENCODE G026 (63856 features)
#
# SAÍDA:
#   data/global/TCGA_GTEx_thyroid_counts.tsv  (genes × amostras matched)
#   data/global/MANIFEST_counts.tsv
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(recount3)
  library(SummarizedExperiment)
  library(data.table)
})

cat("══════════════════════════════════════════════════════════\n")
cat("FASE 4/7 — Matriz de contagens (recount3)\n")
cat("R:", as.character(getRversion()), "| recount3:", as.character(packageVersion("recount3")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("══════════════════════════════════════════════════════════\n\n")

dir_global <- "data/global"
if (!dir.exists(dir_global)) dir.create(dir_global, recursive = TRUE)

# ── 1. IDs de amostras da matriz TPM (referência das 783) ────────────────────
hdr <- fread(file.path(dir_global, "TCGA_GTEx_thyroid_tpm.tsv"), nrows = 0)
sids <- colnames(hdr)
sids <- sids[nzchar(sids) & !grepl("^V[0-9]+$", sids) & !is.na(sids)]
tcga_tpm <- sids[grepl("^TCGA", sids)]
gtex_tpm <- sids[grepl("^GTEX", sids)]
cat(sprintf("  Amostras TPM de referência: %d (TCGA=%d, GTEx=%d)\n",
            length(sids), length(tcga_tpm), length(gtex_tpm)))

# ── 2. Carrega recount3 ───────────────────────────────────────────────────────
ap <- available_projects()
rse_t <- create_rse(ap[ap$project == "THCA" & ap$file_source == "tcga", ])
rse_g <- create_rse(ap[ap$project == "THYROID" & ap$file_source == "gtex", ])

# TCGA: barcode completo -> chave "participante-sample(2 dígitos)"
tb <- colData(rse_t)$tcga.tcga_barcode
tb_key <- paste0(substr(tb, 1, 12), "-", substr(tb, 14, 15))
# GTEx: sampid direto
gs <- colData(rse_g)$gtex.sampid

# ── 3. Matching ───────────────────────────────────────────────────────────────
mt_t <- match(tcga_tpm, tb_key)
mt_g <- match(gtex_tpm, gs)
cat(sprintf("  TCGA matched: %d / %d\n", sum(!is.na(mt_t)), length(tcga_tpm)))
cat(sprintf("  GTEx matched: %d / %d\n", sum(!is.na(mt_g)), length(gtex_tpm)))
cat("  GTEx ausentes no recount3:", paste(gtex_tpm[is.na(mt_g)], collapse = ", "), "\n")

# ── 4. Extrai contagens ───────────────────────────────────────────────────────
counts_t <- assay(rse_t, "raw_counts")[, mt_t[!is.na(mt_t)], drop = FALSE]
counts_g <- assay(rse_g, "raw_counts")[, mt_g[!is.na(mt_g)], drop = FALSE]

# Nomes de genes (símbolo). Mantém Ensembl quando símbolo ausente.
gname_t <- rowData(rse_t)$gene_name
gname_g <- rowData(rse_g)$gene_name
stopifnot(identical(gname_t, gname_g))   # mesma anotação G026
gn <- ifelse(is.na(gname_t) | gname_t == "", rowData(rse_t)$gene_id, gname_t)

# Colunas com nomes de amostra originais (TPM)
colnames(counts_t) <- tcga_tpm[!is.na(mt_t)]
colnames(counts_g) <- gtex_tpm[!is.na(mt_g)]

# ── 5. Junta TCGA + GTEx e colapsa símbolos duplicados (soma) ────────────────
C <- cbind(counts_t, counts_g)
stopifnot(identical(rownames(C), rownames(counts_t)))
rownames(C) <- gn

# Colapsa duplicatas de símbolo somando contagens (apropriado para count data)
if (anyDuplicated(rownames(C))) {
  dup <- rownames(C)[duplicated(rownames(C))]
  cat(sprintf("  Símbolos duplicados: %d únicos -> colapsando por SOMA\n", length(unique(dup))))
  C <- rowsum(C, group = rownames(C), reorder = FALSE)
}
C <- C[rownames(C) != "", , drop = FALSE]

# Validação de integridade
cat(sprintf("  Matriz final: %d genes × %d amostras\n", nrow(C), ncol(C)))
cat(sprintf("  Valores: min=%g max=%g | inteiros=%s | NAs=%d | negativos=%d\n",
            min(C), max(C),
            all(C == round(C)), sum(is.na(C)), sum(C < 0)))

# ── 6. Grava ──────────────────────────────────────────────────────────────────
out <- file.path(dir_global, "TCGA_GTEx_thyroid_counts.tsv")
dt <- data.table(gene = rownames(C), C)
data.table::fwrite(dt, out, sep = "\t", quote = FALSE)
cat("  -> ", out, "\n")

manifest <- data.table(
  field = c("dataset", "source", "project_tcga", "project_gtex", "assay",
            "annotation", "genes", "samples", "samples_matched_tcga",
            "samples_matched_gtex", "samples_missing_gtex", "unit",
            "R_version", "recount3_version", "access_date"),
  value = c("recount3 gene-level raw counts (STAR)", "recount3 (Lieber Institute)",
            "THCA/tcga (572)", "THYROID/gtex (706)", "raw_counts",
            "GENCODE G026 (63856 features)", as.character(nrow(C)),
            as.character(ncol(C)), as.character(sum(!is.na(mt_t))),
            as.character(sum(!is.na(mt_g))), paste(gtex_tpm[is.na(mt_g)], collapse=","),
            "raw counts (inteiros, não normalizados)",
            as.character(getRversion()), as.character(packageVersion("recount3")),
            format(Sys.Date(), "%Y-%m-%d"))
)
data.table::fwrite(manifest, file.path(dir_global, "MANIFEST_counts.tsv"), sep = "\t", quote = FALSE)

cat("\n══════════════════════════════════════════════════════════\n")
cat("FASE 4/7 — Matriz de contagens CONCLUÍDA.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
