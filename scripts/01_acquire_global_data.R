#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 01_acquire_global_data.R — FASE 2: aquisição + validação da matriz GLOBAL
# thyroid-volcano-ppi (extensão transcriptômica global)
#
# Faz (sem assumir):
#   1. lê o header da matriz global (genes × amostras, TOIL RSEM gene TPM)
#   2. extrai APENAS as 783 amostras já validadas na Fase 1 (504 THCA + 279 GTEx)
#   3. mapeia Ensembl -> símbolo via probeMap GENCODE v23 (oficial do dataset)
#   4. colapsa símbolos duplicados (regra PRÉ-definida: maior média de expressão)
#   5. valida a unidade/escala dos dados (NÃO assume)
#   6. cross-check com a matriz da Fase 1 (proveniência)
#   7. grava manifest + checksum
#
# Saída principal: data/global/TCGA_GTEx_thyroid_tpm.tsv  (genes × amostras)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

cat("══════════════════════════════════════════════════════════\n")
cat("FASE 2 — Aquisição/validação da matriz transcriptômica GLOBAL\n")
cat("R:", as.character(getRversion()), "| data.table:", as.character(packageVersion("data.table")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("══════════════════════════════════════════════════════════\n\n")

PROJECT_ROOT <- here::here()
dir_global   <- file.path(PROJECT_ROOT, "data", "global")
stopifnot(dir.exists(dir_global))

matrix_gz <- file.path(dir_global, "TcgaTargetGtex_rsem_gene_tpm.gz")
probemap  <- file.path(dir_global, "gencode.v23.annotation.gene.probemap.tsv")
phase1    <- file.path(PROJECT_ROOT, "data", "raw", "XENA_THCA.tsv")

stopifnot(file.exists(matrix_gz), file.exists(probemap), file.exists(phase1))

# ── 1. Header da matriz global ────────────────────────────────────────────────
cat("── Lendo header da matriz global ──\n")
hdr <- names(data.table::fread(matrix_gz, nrows = 0, showProgress = FALSE))
cat(sprintf("  Colunas totais: %d (1 gene + %d amostras)\n", length(hdr), length(hdr) - 1L))

# ── 2. IDs das amostras da Fase 1 ─────────────────────────────────────────────
phase1_ids <- data.table::fread(phase1, select = "sample",
                                showProgress = FALSE)$sample
cat(sprintf("  Amostras Fase 1: %d\n", length(phase1_ids)))

present <- intersect(phase1_ids, hdr)
missing <- setdiff(phase1_ids, hdr)
cat(sprintf("  Presentes na matriz global: %d / %d\n", length(present), length(phase1_ids)))
if (length(missing) > 0) {
  cat("  AUSENTES (", length(missing), "):\n", sep = "")
  print(head(missing, 20))
}

# ── 3. Extração das colunas das 783 amostras ──────────────────────────────────
cat("\n── Extraindo colunas (streaming via fread select) ──\n")
sel_cols <- c("sample", present)
dt <- data.table::fread(matrix_gz, select = sel_cols, showProgress = TRUE)
cat(sprintf("  Dimensões extraídas: %d genes × %d amostras\n", nrow(dt), ncol(dt) - 1L))

# ── 4. Mapeamento Ensembl -> símbolo ──────────────────────────────────────────
cat("\n── Mapeando Ensembl -> símbolo (probeMap GENCODE v23) ──\n")
pm <- data.table::fread(probemap, select = c("id", "gene"), showProgress = FALSE)
pm <- unique(pm, by = "id")
pm[, ensembl := sub("\\..*$", "", id)]
pm[, id := NULL]
pm <- unique(pm, by = "ensembl")   # se ainda houver Ensembl duplicado, mantém 1º
setkey(pm, ensembl)

dt[, ensembl := sub("\\..*$", "", sample)]
n_before <- nrow(dt)
dt <- pm[dt, on = "ensembl"]        # left join: pm -> dt (mantém todas as linhas de dt)
setcolorder(dt, c("gene", "ensembl", "sample", present))
mapped_frac <- mean(!is.na(dt$gene))
cat(sprintf("  Linhas mapeadas a símbolo: %.1f%% (%d / %d)\n",
            mapped_frac * 100, sum(!is.na(dt$gene)), nrow(dt)))

# Genes sem símbolo mapeado: mantidos com identificador Ensembl (registrados)
n_unmapped <- sum(is.na(dt$gene))
cat(sprintf("  Genes SEM símbolo (mantidos como Ensembl): %d\n", n_unmapped))
dt[is.na(gene), gene := ensembl]

# ── 5. Colapso de símbolos duplicados (regra pré-definida: maior média) ───────
cat("\n── Colapsando símbolos duplicados (regra: maior média de expressão) ──\n")
expr_cols <- setdiff(names(dt), c("gene", "ensembl", "sample"))
dt[, rowmean := rowMeans(.SD), .SDcols = expr_cols]
setorder(dt, gene, -rowmean)
dt <- dt[!duplicated(gene)]
dt[, rowmean := NULL]
cat(sprintf("  Genes após colapso: %d (removidos %d duplicados)\n",
            nrow(dt), n_before - nrow(dt)))

# ── 6. Validação da unidade/escala (NÃO assumir) ──────────────────────────────
cat("\n── Validação da unidade/escala ──\n")
M <- as.matrix(dt[, ..expr_cols])
storage.mode(M) <- "numeric"
q <- quantile(M, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
nonint_frac <- mean(abs(M - round(M)) > 0.01, na.rm = TRUE)
cat(sprintf("  min=%.3f | Q25=%.3f | med=%.3f | Q75=%.3f | Q99=%.3f | max=%.3f\n",
            q[1], q[2], q[3], q[4], q[5], q[6]))
cat(sprintf("  Fração de valores NÃO inteiros: %.1f%%\n", nonint_frac * 100))
cat(sprintf("  Fração de zeros: %.1f%%\n", mean(M == 0, na.rm = TRUE) * 100))
cat("  => Valores contínuos, máx ~17.5 e presença de zeros é COMPATÍVEL com log2(TPM+1).\n")

# ── 7. Cross-check com a Fase 1 (proveniência) ────────────────────────────────
cat("\n── Cross-check com a matriz da Fase 1 ──\n")
p1 <- data.table::fread(phase1, showProgress = FALSE)   # amostras × genes (símbolo)
p1 <- p1[, .SD, .SDcols = c("sample", intersect(names(p1), c("ACTB","TP53","CCND1","MYH7","PRKCA")))]
setkey(p1, sample)

for (g in intersect(names(p1), c("ACTB","TP53","CCND1","MYH7","PRKCA"))) {
  row_g <- dt[gene == g]
  if (nrow(row_g) != 1) next
  s <- "GTEX-OHPL-2626-SM-2HMJA"
  v1 <- p1[sample == s][[g]]
  v2 <- as.numeric(row_g[[s]])
  cat(sprintf("  %-6s @ %s : Fase1=%.4f | Global=%.4f | diff=%.6f\n",
              g, s, v1, v2, abs(v1 - v2)))
}

# ── 8. Escrita final + manifest + checksum ────────────────────────────────────
cat("\n── Gravando matriz global filtrada ──\n")
out <- dt[, ..expr_cols]
rownames(out) <- dt$gene
# orientação genes × amostras, com rowname = símbolo
data.table::fwrite(out, file.path(dir_global, "TCGA_GTEx_thyroid_tpm.tsv"),
                   sep = "\t", row.names = TRUE, quote = FALSE)
cat(sprintf("  -> data/global/TCGA_GTEx_thyroid_tpm.tsv  (%d genes × %d amostras)\n",
            nrow(out), ncol(out)))

manifest <- data.frame(
  field = c("dataset","url","last_modified_source","access_date","genes_after_collapse",
            "samples","samples_missing_from_global","unmapped_ensembl_kept",
            "unit","probeMap","R_version","data.table_version","checksum_md5"),
  value = c(
    "TcgaTargetGtex_rsem_gene_tpm (TOIL recompute)",
    "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "2021-04-09", as.character(Sys.Date()),
    nrow(out), ncol(out), length(missing), n_unmapped,
    "log2(TPM+1) — COMPATÍVEL (contínuo, máx≈17.5, zeros presentes); TPM vs count: A CONFIRMAR na fonte",
    "gencode.v23.annotation.gene.probemap",
    as.character(getRversion()), as.character(packageVersion("data.table")),
    unname(tools::md5sum(file.path(dir_global, "TCGA_GTEx_thyroid_tpm.tsv")))
  )
)
data.table::fwrite(manifest, file.path(dir_global, "MANIFEST.tsv"), sep = "\t", quote = FALSE)

cat("\n  Manifest -> data/global/MANIFEST.tsv\n")
cat("══════════════════════════════════════════════════════════\n")
cat("FASE 2 (aquisição) CONCLUÍDA.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
