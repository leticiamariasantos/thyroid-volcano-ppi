# ═══════════════════════════════════════════════════════════════════════════════
# 03_de.R — Expressão diferencial (reconstruída do zero)
#
# Métodos: limma (log2 TPM, análise principal), limma-voom (contagens) e
# DESeq2 (contagens). Gera logFC, p-value, FDR, ranking completo e UP/DOWN.
# Compara métodos (correlação logFC, direção, ranking, Jaccard).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
})

log_msg("══ Expressão diferencial ══")

tpm  <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
cnt  <- readRDS(file.path(DIR_DATAOUT, "counts_matrix.rds"))
meta <- fread(file.path(DIR_AUDIT, "sample_metadata.tsv"))

# Alinhar metadados às colunas
tpm_meta <- meta[match(colnames(tpm), meta$sample)]
cnt_meta <- meta[match(colnames(cnt), meta$sample)]
condition_tpm <- factor(tpm_meta$condition, levels = c("Normal", "THCA"))
condition_cnt <- factor(cnt_meta$condition, levels = c("Normal", "THCA"))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. limma (análise PRINCIPAL — log2 TPM)
# ═══════════════════════════════════════════════════════════════════════════════
log_msg("── limma (log2 TPM) ──")
# filtro de baixa expressão: expresso (TPM>0.1) em >=25% das amostras
expr_thr <- log2(0.1 + 0.001)
keep_tpm <- rowSums(tpm > expr_thr) >= (EXPR_FRAC * ncol(tpm))
log_msg("Genes retidos p/ limma:", sum(keep_tpm), "de", nrow(tpm))
tpm_f <- tpm[keep_tpm, , drop = FALSE]

design <- model.matrix(~ condition_tpm)
fit <- lmFit(tpm_f, design)
fit <- eBayes(fit)
limma_res <- topTable(fit, coef = "condition_tpmTHCA", number = Inf, sort.by = "none")
limma_res$gene_symbol <- rownames(limma_res)
limma_res$method <- "limma"
limma_res$regulation <- ifelse(limma_res$adj.P.Val < FDR_THRESH & limma_res$logFC > LOGFC_THRESH, "Up",
                        ifelse(limma_res$adj.P.Val < FDR_THRESH & limma_res$logFC < -LOGFC_THRESH, "Down", "NS"))
limma_res <- limma_res[order(limma_res$P.Value), ]
limma_res$rank <- seq_len(nrow(limma_res))
fwrite_tsv(limma_res, file.path(DIR_DE, "limma_full_results.tsv"))

n_limma <- sum(limma_res$regulation != "NS")
n_limma_up <- sum(limma_res$regulation == "Up")
n_limma_dn <- sum(limma_res$regulation == "Down")
log_msg(sprintf("limma: %d testados | DEG=%d (up=%d, down=%d)", nrow(limma_res), n_limma, n_limma_up, n_limma_dn))

# ranking pré-classificado (signed -log10(p) * sign(logFC) e t-statistic)
limma_rank_p <- setNames(-log10(limma_res$P.Value) * sign(limma_res$logFC), limma_res$gene_symbol)
limma_rank_t <- setNames(limma_res$t, limma_res$gene_symbol)
limma_rank_p <- sort(limma_rank_p, decreasing = TRUE)
limma_rank_t <- sort(limma_rank_t, decreasing = TRUE)
saveRDS(limma_rank_t, file.path(DIR_DATAOUT, "limma_rank_t.rds"))
saveRDS(limma_rank_p, file.path(DIR_DATAOUT, "limma_rank_p.rds"))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. limma-voom (contagens recount3)
# ═══════════════════════════════════════════════════════════════════════════════
log_msg("── limma-voom (contagens) ──")
cnt_f <- cnt[rowSums(cnt) > 0, , drop = FALSE]
dge <- DGEList(counts = cnt_f, group = condition_cnt)
keep_cnt <- filterByExpr(dge, group = condition_cnt, min.count = MIN_EXPR_CPM)
log_msg("Genes retidos p/ voom/DESeq2:", sum(keep_cnt), "de", nrow(cnt_f))
dge <- dge[keep_cnt, ]
dge <- calcNormFactors(dge)
design_cnt <- model.matrix(~ condition_cnt)
v <- voom(dge, design_cnt, plot = FALSE)
fit_v <- lmFit(v, design_cnt)
fit_v <- eBayes(fit_v)
voom_res <- topTable(fit_v, coef = "condition_cntTHCA", number = Inf, sort.by = "none")
# nota: design usa condition_tpm (mesma ordem das colunas do tpm); corrigir p/ counts
# (a ordem de colunas de cnt == ordem de tpm? verificar via common samples)
voom_res$gene_symbol <- rownames(voom_res)
voom_res$method <- "limma-voom"
voom_res$regulation <- ifelse(voom_res$adj.P.Val < FDR_THRESH & voom_res$logFC > LOGFC_THRESH, "Up",
                       ifelse(voom_res$adj.P.Val < FDR_THRESH & voom_res$logFC < -LOGFC_THRESH, "Down", "NS"))
voom_res <- voom_res[order(voom_res$P.Value), ]
voom_res$rank <- seq_len(nrow(voom_res))
fwrite_tsv(voom_res, file.path(DIR_DE, "voom_full_results.tsv"))
voom_rank_t <- setNames(voom_res$t, voom_res$gene_symbol)
voom_rank_t <- sort(voom_rank_t[!duplicated(names(voom_rank_t))], decreasing = TRUE)
saveRDS(voom_rank_t, file.path(DIR_DATAOUT, "voom_rank_t.rds"))

n_voom <- sum(voom_res$regulation != "NS")
log_msg(sprintf("voom: %d testados | DEG=%d (up=%d, down=%d)", nrow(voom_res), n_voom,
                sum(voom_res$regulation == "Up"), sum(voom_res$regulation == "Down")))

# ═══════════════════════════════════════════════════════════════════════════════
# 3. DESeq2 (contagens recount3)
# ═══════════════════════════════════════════════════════════════════════════════
log_msg("── DESeq2 (contagens) ──")
cnt_dds <- cnt_f[keep_cnt, , drop = FALSE]
# metadados na ordem das colunas de counts
cnt_meta2 <- meta[match(colnames(cnt_dds), meta$sample)]
condition_dds <- factor(cnt_meta2$condition, levels = c("Normal", "THCA"))
coldata <- data.frame(condition = condition_dds, row.names = colnames(cnt_dds))
dds <- DESeqDataSetFromMatrix(countData = round(cnt_dds), colData = coldata, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
dds_res <- results(dds, contrast = c("condition", "THCA", "Normal"), independentFiltering = TRUE)
dds_res <- as.data.frame(dds_res)
dds_res$gene_symbol <- rownames(dds_res)
dds_res$method <- "DESeq2"
dds_res$regulation <- ifelse(!is.na(dds_res$padj) & dds_res$padj < FDR_THRESH & dds_res$log2FoldChange > LOGFC_THRESH, "Up",
                      ifelse(!is.na(dds_res$padj) & dds_res$padj < FDR_THRESH & dds_res$log2FoldChange < -LOGFC_THRESH, "Down", "NS"))
dds_res$regulation[is.na(dds_res$padj)] <- "NS"
dds_res <- dds_res[order(dds_res$pvalue), ]
dds_res$rank <- seq_len(nrow(dds_res))
# renomear para consistência
names(dds_res)[names(dds_res) == "log2FoldChange"] <- "logFC"
names(dds_res)[names(dds_res) == "padj"] <- "adj.P.Val"
names(dds_res)[names(dds_res) == "pvalue"] <- "P.Value"
fwrite_tsv(dds_res, file.path(DIR_DE, "deseq2_full_results.tsv"))
deseq2_rank_stat <- setNames(dds_res$stat, dds_res$gene_symbol)
deseq2_rank_stat <- sort(deseq2_rank_stat[!duplicated(names(deseq2_rank_stat))], decreasing = TRUE)
saveRDS(deseq2_rank_stat, file.path(DIR_DATAOUT, "deseq2_rank_stat.rds"))

n_dds <- sum(dds_res$regulation != "NS")
log_msg(sprintf("DESeq2: %d testados | DEG=%d (up=%d, down=%d)", nrow(dds_res), n_dds,
                sum(dds_res$regulation == "Up"), sum(dds_res$regulation == "Down")))

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Sumário de DEGs por método
# ═══════════════════════════════════════════════════════════════════════════════
deg_summary <- data.table(
  method = c("limma", "limma-voom", "DESeq2"),
  genes_tested = c(nrow(limma_res), nrow(voom_res), nrow(dds_res)),
  total_deg = c(n_limma, n_voom, n_dds),
  up = c(n_limma_up, sum(voom_res$regulation == "Up"), sum(dds_res$regulation == "Up")),
  down = c(n_limma_dn, sum(voom_res$regulation == "Down"), sum(dds_res$regulation == "Down")),
  threshold = sprintf("|log2FC|>%.1f & FDR<%.2f", LOGFC_THRESH, FDR_THRESH)
)
fwrite_tsv(deg_summary, file.path(DIR_DE, "DEG_summary.tsv"))
log_msg("DEG summary gravado.")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Comparação entre métodos
# ═══════════════════════════════════════════════════════════════════════════════
log_msg("── Comparação entre métodos ──")
# voom e DESeq2 compartilham IDs de genes (contagens); limma usa símbolos
voom_key <- data.table(gene = voom_res$gene_symbol, logFC_v = voom_res$logFC,
                       fdr_v = voom_res$adj.P.Val, reg_v = voom_res$regulation,
                       p_v = voom_res$P.Value)
dds_key  <- data.table(gene = dds_res$gene_symbol, logFC_d = dds_res$logFC,
                       fdr_d = dds_res$adj.P.Val, reg_d = dds_res$regulation,
                       p_d = dds_res$P.Value)
limma_key <- data.table(gene = limma_res$gene_symbol, logFC_l = limma_res$logFC,
                        fdr_l = limma_res$adj.P.Val, reg_l = limma_res$regulation,
                        p_l = limma_res$P.Value)

# pares de métodos
pairs <- list(
  list(nm = "limma_vs_voom", a = limma_key, b = voom_key, na_ = c("logFC_l", "logFC_v"),
       na_fdr = c("fdr_l", "fdr_v"), na_reg = c("reg_l", "reg_v")),
  list(nm = "limma_vs_deseq2", a = limma_key, b = dds_key, na_ = c("logFC_l", "logFC_d"),
       na_fdr = c("fdr_l", "fdr_d"), na_reg = c("reg_l", "reg_d")),
  list(nm = "voom_vs_deseq2", a = voom_key, b = dds_key, na_ = c("logFC_v", "logFC_d"),
       na_fdr = c("fdr_v", "fdr_d"), na_reg = c("reg_v", "reg_d"))
)

compare_methods <- function(a, b, na_, na_fdr, na_reg) {
  m <- merge(a, b, by = "gene", all = FALSE)
  rho_lfc <- cor(m[[na_[1]]], m[[na_[2]]], method = "spearman", use = "complete.obs")
  # concordância de direção entre DEGs significativos em ambos
  sig_both <- m[[na_fdr[1]]] < FDR_THRESH & m[[na_fdr[2]]] < FDR_THRESH
  dir_conc <- if (sum(sig_both) > 0) {
    mean(sign(m[[na_[1]]][sig_both]) == sign(m[[na_[2]]][sig_both]))
  } else NA_real_
  # ranking concordance (top 500)
  r1 <- m[[na_[1]]]; r2 <- m[[na_[2]]]
  top500 <- union(head(m$gene[order(-abs(r1))], 500), head(m$gene[order(-abs(r2))], 500))
  # Jaccard de DEGs
  setA <- m$gene[m[[na_reg[1]]] != "NS"]; setB <- m$gene[m[[na_reg[2]]] != "NS"]
  jacc <- length(intersect(setA, setB)) / length(union(setA, setB))
  data.table(pair = NA, spearman_logFC = rho_lfc, direction_concordance = dir_conc,
             jaccard_deg = jacc, n_shared_deg = length(intersect(setA, setB)),
             n_deg_a = length(setA), n_deg_b = length(setB))
}
cmp_list <- lapply(pairs, function(p) {
  r <- compare_methods(p$a, p$b, p$na_, p$na_fdr, p$na_reg)
  r$pair <- p$nm; r
})
method_concordance <- rbindlist(cmp_list)
fwrite_tsv(method_concordance, file.path(DIR_DE, "method_concordance.tsv"))
log_msg("Concordância entre métodos gravada.")

# gene-level robustness (para genes-chave e futuramente)
allg <- merge(merge(limma_key, voom_key, by = "gene", all = TRUE),
              dds_key, by = "gene", all = TRUE)
allg[is.na(logFC_l), logFC_l := 0]; allg[is.na(logFC_v), logFC_v := 0]; allg[is.na(logFC_d), logFC_d := 0]
allg$sign_l <- sign(allg$logFC_l); allg$sign_v <- sign(allg$logFC_v); allg$sign_d <- sign(allg$logFC_d)
allg$n_same_dir <- (allg$sign_l == allg$sign_v) + (allg$sign_l == allg$sign_d) + (allg$sign_v == allg$sign_d)
allg$robustness <- fifelse(allg$n_same_dir >= 3 & (allg$sign_l != 0), "ROBUSTO",
                    fifelse(allg$n_same_dir == 2, "PARCIALMENTE ROBUSTO", "INSTAVEL"))
fwrite_tsv(allg, file.path(DIR_DE, "gene_concordance_all.tsv"))

log_msg("══ Expressão diferencial concluída ══")
cat("\n=== RESUMO DEG ===\n")
print(deg_summary)
cat("\n=== CONCORDÂNCIA ===\n")
print(method_concordance)
