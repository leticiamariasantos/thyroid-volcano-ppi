# ═══════════════════════════════════════════════════════════════════════════════
# 01_audit_data.R — Auditoria integral do dataset de entrada (reconstruída do zero)
#
# Verifica integridade, dimensões, IDs, correspondência matriz↔metadados,
# duplicatas, valores ausentes e distribuição. Nenhuma conclusão biológica é
# derivada aqui — apenas caracterização técnica do input.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))

log_msg("══ Auditoria dos dados de entrada ══")
chk <- new_checker()

audit_log <- list()

# ── 1. Presença e tamanho dos arquivos ────────────────────────────────────────
for (f in c(TPM_FILE, COUNTS_FILE, RAW_TPM_GZ)) {
  chk$ck(file.exists(f) && file.info(f)$size > 0,
     sprintf("existe e não-vazio: %s (%.1f MB)", basename(f), file.info(f)$size/1e6))
}

# ── 2. Leitura dos cabeçalhos (não carrega matriz inteira ainda) ──────────────
tpm_hdr <- fread(TPM_FILE, nrows = 0)
tpm_ncol <- ncol(tpm_hdr)
log_msg("TPM header colunas:", tpm_ncol)

# Leitura eficiente da matriz TPM (símbolos nas linhas)
tpm <- fread(TPM_FILE, header = TRUE, sep = "\t", showProgress = FALSE)
tpm_genes <- tpm[[1]]
tpm_samples <- names(tpm)[-1]
tpm_mat <- as.matrix(tpm[, -1, with = FALSE])
rownames(tpm_mat) <- tpm_genes
rm(tpm); invisible(gc())

log_msg("TPM matriz:", nrow(tpm_mat), "genes ×", ncol(tpm_mat), "amostras")

# ── 3. Dimensões esperadas ─────────────────────────────────────────────────────
chk$ck(nrow(tpm_mat) == 58581, sprintf("TPM genes=%d (esperado 58581)", nrow(tpm_mat)))
chk$ck(ncol(tpm_mat) == 783, sprintf("TPM amostras=%d (esperado 783)", ncol(tpm_mat)))

# ── 4. IDs das amostras ────────────────────────────────────────────────────────
pref <- sub("-.*", "", tpm_samples)
ntcga <- sum(pref == "TCGA"); ngt <- sum(pref == "GTEX")
chk$ck(ntcga == 504 && ngt == 279, sprintf("TPM grupos: TCGA=%d, GTEx=%d", ntcga, ngt))
chk$ck(all(pref %in% c("TCGA", "GTEX")), "TPM: todos os prefixos são TCGA/GTEX")
audit_log$tpm_ntcga <- ntcga; audit_log$tpm_ngtex <- ngt

# ── 5. Duplicatas de genes ─────────────────────────────────────────────────────
chk$ck(!anyDuplicated(tpm_genes), "TPM: sem genes duplicados")
chk$ck(!anyDuplicated(tpm_samples), "TPM: sem amostras duplicadas")

# ── 6. Valores ausentes / infinitos ────────────────────────────────────────────
n_na <- sum(is.na(tpm_mat)); n_inf <- sum(is.infinite(tpm_mat))
chk$ck(n_na == 0 && n_inf == 0, sprintf("TPM: sem NA/Inf (NA=%d, Inf=%d)", n_na, n_inf))
audit_log$tpm_na <- n_na; audit_log$tpm_inf <- n_inf

# ── 7. Distribuição da expressão ───────────────────────────────────────────────
mn <- mean(tpm_mat); med <- median(tpm_mat)
q0 <- min(tpm_mat); q100 <- max(tpm_mat)
frac_zero <- mean(tpm_mat <= log2(0.001 + 1e-9))
log_msg(sprintf("TPM distribuição: mean=%.3f med=%.3f min=%.3f max=%.3f frac_zero=%.3f",
                mn, med, q0, q100, frac_zero))
audit_log$tpm_mean <- mn; audit_log$tpm_median <- med
audit_log$tpm_min <- q0; audit_log$tpm_max <- q100
audit_log$tpm_frac_zero <- frac_zero

# faixa esperada de log2(TPM+0.001)
chk$ck(q100 < 18 && q100 > 10, sprintf("TPM max=%.2f (faixa esperada log2(TPM+0.001))", q100))

# ── 8. Verificação de que é log2(TPM+0.001): mínimo ≈ log2(0.001) = -9.966 ─────
lo_expected <- log2(0.001)
chk$ck(abs(q0 - lo_expected) < 0.05, sprintf("TPM min=%.3f ≈ log2(0.001)=%.3f", q0, lo_expected))

# ── 9. Matriz de contagens ─────────────────────────────────────────────────────
cnt <- fread(COUNTS_FILE, header = TRUE, sep = "\t", showProgress = FALSE)
cnt_genes <- cnt[[1]]; cnt_samples <- names(cnt)[-1]
cnt_mat <- as.matrix(cnt[, -1, with = FALSE]); rownames(cnt_mat) <- cnt_genes
rm(cnt); invisible(gc())

log_msg("Counts matriz:", nrow(cnt_mat), "genes ×", ncol(cnt_mat), "amostras")
chk$ck(nrow(cnt_mat) == 56937, sprintf("Counts genes=%d (esperado 56937)", nrow(cnt_mat)))
chk$ck(ncol(cnt_mat) == 782, sprintf("Counts amostras=%d (esperado 782)", ncol(cnt_mat)))

cpref <- sub("-.*", "", cnt_samples)
cnt_tcga <- sum(cpref == "TCGA"); cnt_gt <- sum(cpref == "GTEX")
chk$ck(cnt_tcga == 504 && cnt_gt == 278, sprintf("Counts grupos: TCGA=%d, GTEx=%d", cnt_tcga, cnt_gt))
audit_log$counts_ntcga <- cnt_tcga; audit_log$counts_ngtex <- cnt_gt

chk$ck(!anyDuplicated(cnt_genes), "Counts: sem genes duplicados")
chk$ck(!anyDuplicated(cnt_samples), "Counts: sem amostras duplicadas")
chk$ck(all(is.finite(cnt_mat)), "Counts: sem NA/Inf")
chk$ck(all(cnt_mat >= 0), "Counts: não-negativas")
chk$ck(all(cnt_mat == round(cnt_mat)), "Counts: inteiras")

# ── 10. Sobreposição de amostras TPM × Counts ─────────────────────────────────
common_samples <- intersect(tpm_samples, cnt_samples)
log_msg("Amostras comuns TPM×Counts:", length(common_samples))
audit_log$common_samples <- length(common_samples)
missing_in_counts <- setdiff(tpm_samples, cnt_samples)
audit_log$missing_in_counts <- missing_in_counts

# ── 11. Metadados de condição (derivado do próprio ID, documentado) ────────────
meta <- data.table(
  sample = tpm_samples,
  condition = ifelse(grepl("^TCGA", tpm_samples), "THCA", "Normal"),
  source = ifelse(grepl("^TCGA", tpm_samples), "TCGA", "GTEx")
)
fwrite_tsv(meta, file.path(DIR_AUDIT, "sample_metadata.tsv"))
chk$ck(nrow(meta) == 783, "metadados gerados para 783 amostras")

# ── 12. Duplicata de amostras entre TCGA e GTEx (contaminação cruzada) ────────
chk$ck(length(intersect(tpm_samples[grepl("^TCGA", tpm_samples)],
                    tpm_samples[grepl("^GTEX", tpm_samples)])) == 0,
   "nenhuma amostra presente simultaneamente em TCGA e GTEx")

# ── 13. Consistência com o arquivo .gz bruto ───────────────────────────────
# O .gz é a matriz GLOBAL TOIL (todos os tecidos); o arquivo processado é o
# subconjunto tireoide (TCGA-THCA + GTEx thyroid). Espera-se que o .gz seja um
# SUPERCONJUNTO em número de amostras, não igual.
if (file.exists(RAW_TPM_GZ)) {
  md5 <- tools::md5sum(RAW_TPM_GZ)
  log_msg("MD5 do .gz bruto:", unname(md5))
  audit_log$raw_gz_md5 <- unname(md5)
  con <- gzfile(RAW_TPM_GZ, "rt")
  hdr_line <- readLines(con, n = 1)
  close(con)
  raw_ncol <- length(strsplit(hdr_line, "\t")[[1]])
  log_msg("Raw .gz header colunas (global):", raw_ncol)
  audit_log$raw_gz_ncol <- raw_ncol
  audit_log$raw_gz_nsamples <- raw_ncol - 1
  chk$ck(raw_ncol - 1 > ncol(tpm_mat),
     sprintf(".gz global tem %d amostras > processado tireoide=%d (superconjunto)",
             raw_ncol - 1, ncol(tpm_mat)))
  # nº de genes: o .gz global contém 60498 genes; após colapso/remoção de
  # unmapped Ensembl o processado fica com 58581 (documentado em MANIFEST.tsv).
  raw_genecount <- nrow(fread(RAW_TPM_GZ, select = 1, header = TRUE, showProgress = FALSE))
  log_msg("Raw .gz nº de genes:", raw_genecount)
  audit_log$raw_gz_ngenes <- raw_genecount
  chk$ck(raw_genecount >= nrow(tpm_mat),
     sprintf(".gz global genes=%d >= processado=%d (colapso documentado)", raw_genecount, nrow(tpm_mat)))
} else {
  log_msg("AVISO: arquivo .gz bruto não encontrado (mantido processado).")
}

# ── 14. Gravação do relatório de auditoria ────────────────────────────────────
audit_log$exec_date <- EXEC_DATE
audit_log$r_version <- as.character(getRversion())
audit_log$n_tpm_genes <- nrow(tpm_mat)
audit_log$n_tpm_samples <- ncol(tpm_mat)
audit_log$n_counts_genes <- nrow(cnt_mat)
audit_log$n_counts_samples <- ncol(cnt_mat)

log_lines <- c(
  "═══ AUDITORIA DE DADOS — FASE 2 REBOOT ═══",
  paste("Data:", EXEC_DATE),
  paste("R:", R.version.string),
  "",
  paste("TPM genes:", nrow(tpm_mat), "| amostras:", ncol(tpm_mat)),
  paste("  TCGA:", ntcga, "| GTEx:", ngt),
  paste("TPM distribuição: mean", round(mn,3), "median", round(med,3),
        "min", round(q0,3), "max", round(q100,3)),
  paste("TPM fração de zeros:", round(frac_zero,4)),
  paste("TPM NA:", n_na, "| Inf:", n_inf),
  "",
  paste("Counts genes:", nrow(cnt_mat), "| amostras:", ncol(cnt_mat)),
  paste("  TCGA:", cnt_tcga, "| GTEx:", cnt_gt),
  paste("Amostras comuns TPM×Counts:", length(common_samples)),
  paste("Amostras ausentes em Counts:", paste(missing_in_counts, collapse=", ")),
  "",
  "NOTA source≡condition: TCGA≡tumor e GTEx≡normal estão perfeitamente confundidos.",
  "Nenhuma correção de batch é alegada como eliminadora desse confundimento."
)
writeLines(log_lines, file.path(DIR_AUDIT, "DATA_AUDIT.txt"))
cat("\n", paste(log_lines, collapse = "\n"), "\n\n")

# salva matriz TPM/contagens em formato leve para etapas seguintes (RDS)
saveRDS(tpm_mat, file.path(DIR_DATAOUT, "tpm_matrix.rds"))
saveRDS(cnt_mat, file.path(DIR_DATAOUT, "counts_matrix.rds"))

fail <- chk$total()
log_msg("══ Auditoria de dados:", fail, "falha(s) ══")
if (fail > 0) quit(status = 1)
