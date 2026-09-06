# ═══════════════════════════════════════════════════════════════════════════════
# 12_counts_deg.R — DE baseada em contagens (limma-voom + DESeq2)
#
# Entrada: data/global/TCGA_GTEx_thyroid_counts.tsv (recount3 raw counts)
#          56937 genes × 782 amostras (504 THCA + 278 GTEx)
#
# Método A (limma-voom): edgeR DGEList + TMM + voom + lmFit + eBayes
# Método B (DESeq2):     raw counts, design ~ condition
# Contraste: THCA - Normal (igual ao TPM)
#
# Filtro (pré-definido, comparável ao TPM): manter genes com CPM > 1 em ≥10% amostras.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(DESeq2)
})

cat("══════════════════════════════════════════════════════════\n")
cat("FASE 7/8 — DE por contagem (voom + DESeq2)\n")
cat("R:", as.character(getRversion()),
    "| edgeR:", as.character(packageVersion("edgeR")),
    "| limma:", as.character(packageVersion("limma")),
    "| DESeq2:", as.character(packageVersion("DESeq2")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("══════════════════════════════════════════════════════════\n\n")

dir_out <- "results/counts_deg"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ── 1. Leitura ────────────────────────────────────────────────────────────────
cat("── Lendo matriz de contagens ──\n")
raw <- fread("data/global/TCGA_GTEx_thyroid_counts.tsv")
gn <- raw$gene
nm <- names(raw)[-1]
M <- as.matrix(raw[, -1, with = FALSE])
rownames(M) <- gn
rm(raw); invisible(gc())

condition <- factor(ifelse(grepl("^TCGA", nm), "THCA", "Normal"),
                    levels = c("Normal", "THCA"))
cat(sprintf("  Matriz: %d genes × %d amostras (THCA=%d, Normal=%d)\n",
            nrow(M), ncol(M), sum(condition == "THCA"), sum(condition == "Normal")))

# ── 2. QC de contagens ────────────────────────────────────────────────────────
cat("\n── QC de contagens ──\n")
lib <- colSums(M)
cat(sprintf("  Library sizes: min=%.2e med=%.2e max=%.2e\n",
            min(lib), median(lib), max(lib)))
# remover genes sem nenhuma contagem (soma 0)
nz <- rowSums(M) > 0
cat(sprintf("  Genes com soma de contagens > 0: %d / %d\n", sum(nz), length(nz)))
M <- M[nz, , drop = FALSE]

# ── 3. Filtro de baixa expressão (CPM > 1 em ≥10%) ────────────────────────────
cat("\n── Filtro de baixa expressão (CPM > 1 em ≥10% amostras) ──\n")
dge <- DGEList(counts = M, group = condition)
cpm_m <- cpm(dge)
keep <- rowSums(cpm_m > 1) >= (0.10 * ncol(M))
cat(sprintf("  Genes mantidos: %d / %d (removidos %d)\n",
            sum(keep), nrow(M), nrow(M) - sum(keep)))
dge <- dge[keep, , keep.lib.sizes = FALSE]

# ── 4. limma-voom ─────────────────────────────────────────────────────────────
cat("\n── Método A: limma-voom ──\n")
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ condition)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)
voom_res <- topTable(fit, coef = "conditionTHCA", number = Inf, sort.by = "none")
voom_res$gene_symbol <- rownames(voom_res)
voom_res$regulation <- ifelse(voom_res$adj.P.Val < 0.05 & voom_res$logFC > 1, "Up",
                       ifelse(voom_res$adj.P.Val < 0.05 & voom_res$logFC < -1, "Down", "NS"))
cat(sprintf("  Genes testados: %d\n", nrow(voom_res)))
cat(sprintf("  DEGs (|logFC|>1, FDR<0.05): up=%d down=%d total=%d\n",
            sum(voom_res$regulation == "Up"), sum(voom_res$regulation == "Down"),
            sum(voom_res$regulation != "NS")))
data.table::fwrite(as.data.table(voom_res), file.path(dir_out, "voom_full_results.tsv"), sep = "\t")

# ── 5. DESeq2 ─────────────────────────────────────────────────────────────────
cat("\n── Método B: DESeq2 ──\n")
colData <- data.frame(condition = condition, row.names = nm)
cds <- DESeqDataSetFromMatrix(countData = M, colData = colData, design = ~ condition)
# aplicar o mesmo filtro (manter genes mantidos no voom para comparabilidade)
cds <- cds[rownames(dge$counts), ]
set.seed(42)
cds <- DESeq(cds, quiet = TRUE)
dds_res <- results(cds, contrast = c("condition", "THCA", "Normal"))
dds_res <- as.data.frame(dds_res)
dds_res$gene_symbol <- rownames(dds_res)
dds_res$regulation <- ifelse(!is.na(dds_res$padj) & dds_res$padj < 0.05 & dds_res$log2FoldChange > 1, "Up",
                      ifelse(!is.na(dds_res$padj) & dds_res$padj < 0.05 & dds_res$log2FoldChange < -1, "Down", "NS"))
cat(sprintf("  Genes testados: %d\n", nrow(dds_res)))
cat(sprintf("  DEGs (|log2FC|>1, padj<0.05): up=%d down=%d total=%d\n",
            sum(dds_res$regulation == "Up", na.rm = TRUE),
            sum(dds_res$regulation == "Down", na.rm = TRUE),
            sum(dds_res$regulation != "NS", na.rm = TRUE)))
data.table::fwrite(as.data.table(dds_res), file.path(dir_out, "deseq2_full_results.tsv"), sep = "\t")

# ── 6. Sumário ────────────────────────────────────────────────────────────────
cat("\n── Sumários ──\n")
summ <- data.table(
  method = c("voom", "DESeq2"),
  genes_tested = c(nrow(voom_res), nrow(dds_res)),
  up = c(sum(voom_res$regulation == "Up"), sum(dds_res$regulation == "Up", na.rm = TRUE)),
  down = c(sum(voom_res$regulation == "Down"), sum(dds_res$regulation == "Down", na.rm = TRUE)),
  total_deg = c(sum(voom_res$regulation != "NS"), sum(dds_res$regulation != "NS", na.rm = TRUE))
)
data.table::fwrite(summ, file.path(dir_out, "counts_deg_summary.tsv"), sep = "\t")
print(summ)

cat("\n══════════════════════════════════════════════════════════\n")
cat("FASE 7/8 — DE por contagem CONCLUÍDA.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
