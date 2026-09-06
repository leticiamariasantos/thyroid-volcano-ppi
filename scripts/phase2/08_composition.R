# ═══════════════════════════════════════════════════════════════════════════════
# 08_composition.R — Análise de composição celular e sensibilidade muscular
#
# Avalia sinais decorrentes de composição celular, com atenção aos marcadores de
# músculo estriado/esquelético (MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM e afins).
# Compara FULL TRANSCRIPTOME vs COMPOSITION-SENSITIVE GENES REMOVED e verifica
# quais vias desaparecem, permanecem ou mudam de magnitude.
# source≡condition é documentado; não há alegação de correção completa.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(data.table)
})

log_msg("══ Análise de composição celular ══")

tpm  <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
meta <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))
genesets <- readRDS(file.path(DIR_DATAOUT, "panel30_genesets.rds"))
condition <- factor(meta$condition[match(colnames(tpm), meta$sample)], levels = c("Normal", "THCA"))

# ── Conjuntos de genes de composição muscular (pré-especificados) ──────────────
CORE_MUSCLE <- c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")   # nomeados no protocolo
STRIATED_MUSCLE <- c(
  "MYH7","MYH1","MYH2","MYH3","MYH4","MYH6","MYH8","MYH13",
  "MYL1","MYL2","MYL3","MYL4","MYLPF","ACTA1","ACTC1","ACTN2","ACTN3",
  "TNNT1","TNNT2","TNNT3","TNNI1","TNNI2","TNNI3","TNNC1","TNNC2",
  "CKM","CKMT2","MB","TTN","NEB","MYOM1","MYOM2","MYBPC1","MYBPC2","MYBPC3",
  "CASQ1","CASQ2","ATP2A1","RYR1","CACNA1S","PYGM","ENO3","MYOZ1","MYOZ2",
  "TPM2","TPM3","LMOD2","LMOD3"
)
muscle_present <- STRIATED_MUSCLE[STRIATED_MUSCLE %in% rownames(tpm)]
log_msg("Marcadores de músculo estriado presentes na matriz:", length(muscle_present))

# ── 1. Score muscular por amostra (média de expressão dos marcadores) ──────────
muscle_score <- colMeans(tpm[muscle_present, , drop = FALSE])
score_df <- data.table(sample = colnames(tpm), muscle_score = muscle_score,
                       condition = condition)
fwrite_tsv(score_df, file.path(DIR_SENS, "muscle_composition_scores.tsv"))
log_msg(sprintf("Score muscular: média tumor=%.2f normal=%.2f",
        mean(muscle_score[condition == "THCA"]),
        mean(muscle_score[condition == "Normal"])))

# ── 2. Função DE (limma) ───────────────────────────────────────────────────────
run_limma_rank <- function(mat) {
  expr_thr <- log2(0.1 + 0.001)
  keep <- rowSums(mat > expr_thr) >= (EXPR_FRAC * ncol(mat))
  m <- mat[keep, , drop = FALSE]
  design <- model.matrix(~ condition)
  fit <- eBayes(lmFit(m, design))
  tt <- topTable(fit, coef = "conditionTHCA", number = Inf, sort.by = "none")
  r <- setNames(tt$t, rownames(tt))
  sort(r[!duplicated(names(r))], decreasing = TRUE)
}

# ── 3. FULL transcriptome ──────────────────────────────────────────────────────
rank_full <- run_limma_rank(tpm)
set.seed(SEED)
gsea_full <- fgseaMultilevel(
  pathways = lapply(genesets, function(g) unique(g[g %in% names(rank_full)])),
  stats = rank_full, minSize = 1, maxSize = 2000, eps = 0, nPermSimple = 10000)
setnames(gsea_full, c("NES","pval","padj"), c("NES_full","pval_full","padj_full"))

# ── 4. COMPOSITION-SENSITIVE GENES REMOVED ─────────────────────────────────────
tpm_nomuscle <- tpm[setdiff(rownames(tpm), muscle_present), , drop = FALSE]
rank_removed <- run_limma_rank(tpm_nomuscle)
set.seed(SEED)
gsea_removed <- fgseaMultilevel(
  pathways = lapply(genesets, function(g) unique(g[g %in% names(rank_removed)])),
  stats = rank_removed, minSize = 1, maxSize = 2000, eps = 0, nPermSimple = 10000)
setnames(gsea_removed, c("NES","pval","padj"), c("NES_removed","pval_removed","padj_removed"))

# ── 5. Comparação FULL vs REMOVED ──────────────────────────────────────────────
comp <- merge(as.data.table(gsea_full)[, .(pathway, NES_full, padj_full)],
              as.data.table(gsea_removed)[, .(pathway, NES_removed, padj_removed)],
              by = "pathway")
comp[, delta_NES := NES_removed - NES_full]
comp[, status := "NS"]
comp[!is.na(padj_full) & !is.na(padj_removed) & padj_full < FDR_THRESH & padj_removed < FDR_THRESH, status := "ROBUSTA"]
comp[!is.na(padj_full) & !is.na(padj_removed) & padj_full < FDR_THRESH & padj_removed >= FDR_THRESH, status := "DESAPARECE"]
comp[!is.na(padj_full) & !is.na(padj_removed) & padj_full >= FDR_THRESH & padj_removed < FDR_THRESH, status := "SURGE"]
comp[is.na(NES_full) | is.na(NES_removed), status := "INDETERMINADO"]
# adiciona nomes reais das vias (a partir do painel)
panel_nm <- fread(file.path(DIR_PATH, "PANEL_30_PATHWAYS.tsv"))[, .(database_id, pathway_name)]
comp <- merge(comp, panel_nm, by.x = "pathway", by.y = "database_id", all.x = TRUE)
fwrite_tsv(comp, file.path(DIR_SENS, "gsea_composition_sensitivity.tsv"))
log_msg("══ Sensibilidade de composição concluída ══")
cat("\n=== Comparação FULL vs composition-removed (30 vias) ===\n")
print(comp[order(-abs(delta_NES))][, .(pathway, NES_full, NES_removed, delta_NES, padj_full, padj_removed, status)], topn = 30)
