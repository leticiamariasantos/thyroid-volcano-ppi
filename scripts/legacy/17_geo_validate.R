# ═══════════════════════════════════════════════════════════════════════════════
# 17_geo_validate.R — Validação externa (GEO): GSE33630 (primária) e GSE60542
#
# Compara DIREÇÃO e SIGNIFICÂNCIA dos genes discovery (congelados) em coortes
# independentes de microarray (Affymetrix HG-U133 Plus 2.0, GPL570).
# NÃO compara expressão absoluta entre plataformas.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

cat("══ FASE 18/19 — Validação externa GEO ══\n")
cat("R:", as.character(getRversion()), "| limma:", as.character(packageVersion("limma")), "\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_v <- "10_validation"
dir.create(dir_v, recursive = TRUE, showWarnings = FALSE)

# ── helpers ───────────────────────────────────────────────────────────────────
read_series_matrix <- function(path) {
  # lê a tabela de expressão de um series matrix GEO (plain text)
  lines <- readLines(gzfile(path), warn = FALSE)
  tb <- which(lines == "!series_matrix_table_begin")
  te <- which(lines == "!series_matrix_table_end")
  title_line <- lines[grep("^!Sample_title", lines)][1]
  title_fields <- strsplit(title_line, "\t")[[1]][-1]
  title_fields <- gsub('"', '', title_fields)
  # matriz
  con <- textConnection(lines[(tb+1):(te-1)])
  mat <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
  close(con)
  list(mat = mat, titles = title_fields)
}

map_probes <- function(probes) {
  # probe -> symbol (HG-U133 Plus 2.0)
  m <- AnnotationDbi::select(hgu133plus2.db, keys = probes, keytype = "PROBEID",
                             columns = "SYMBOL")
  m <- m[!is.na(m$SYMBOL) & m$SYMBOL != "", ]
  m
}

collapse_by_symbol <- function(mat, idcol = "ID_REF") {
  # colapsa probes para gene (representante = maior média)
  probes <- mat[[idcol]]
  map <- map_probes(probes)
  keep <- mat[[idcol]] %in% map$PROBEID
  mat <- mat[keep, , drop = FALSE]
  sym <- map$SYMBOL[match(mat[[idcol]], map$PROBEID)]
  X <- as.matrix(mat[, -which(colnames(mat) == idcol), drop = FALSE])
  storage.mode(X) <- "numeric"
  # média por símbolo; escolher probe de maior média
  gm <- rowMeans(X, na.rm = TRUE)
  ord <- order(-gm)
  X <- X[ord, , drop = FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym)
  X <- X[keep, , drop = FALSE]; rownames(X) <- sym[keep]
  X
}

run_limma <- function(X, groups, g1, g2) {
  sel <- groups %in% c(g1, g2)
  X <- X[, sel, drop = FALSE]
  gr <- factor(ifelse(groups[sel] == g1, g1, g2), levels = c(g2, g1))
  design <- model.matrix(~ gr)
  fit <- lmFit(X, design)
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt
}

# ── GSE33630 ──────────────────────────────────────────────────────────────────
cat("── GSE33630 (PTC vs normal) ──\n")
r1 <- read_series_matrix(file.path(dir_v, "raw", "GSE33630_matrix.txt.gz"))
char_line <- readLines(gzfile(file.path(dir_v, "raw", "GSE33630_matrix.txt.gz")))[
  grep("^!Sample_characteristics_ch1", readLines(gzfile(file.path(dir_v, "raw", "GSE33630_matrix.txt.gz"))))][1]
char_fields <- gsub('"', '', strsplit(char_line, "\t")[[1]][-1])
grp1 <- ifelse(grepl("papillary", char_fields), "PTC",
        ifelse(grepl("non-tumor", char_fields), "Normal",
        ifelse(grepl("anaplastic", char_fields), "ATC", "Other")))
cat("  grupos:", paste(names(table(grp1)), table(grp1), collapse=", "), "\n")
X1 <- collapse_by_symbol(r1$mat)
cat("  genes após colapso:", nrow(X1), "\n")
# alinhar colunas de X1 com grp1 (X1 colunas = GSM ids; mat colunas = titles? verificar)
# mat colnames são GSM ids; grp1 ordenado como char_fields que corresponde às colunas
stopifnot(ncol(X1) == length(grp1))
deg1 <- run_limma(X1, grp1, "PTC", "Normal")
cat("  DEGs GSE33630 (|logFC|>1, FDR<0.05): up=", sum(deg1$adj.P.Val<0.05 & deg1$logFC>1),
    " down=", sum(deg1$adj.P.Val<0.05 & deg1$logFC< -1), "\n")
data.table::fwrite(as.data.table(deg1), file.path(dir_v, "GSE33630_deg.tsv"), sep = "\t")

# ── GSE60542 ──────────────────────────────────────────────────────────────────
cat("\n── GSE60542 (PTC vs normal) ──\n")
r2 <- read_series_matrix(file.path(dir_v, "raw", "GSE60542_matrix.txt.gz"))
title2 <- r2$titles
# 2º campo separado por vírgula descreve o tipo
grp2 <- vapply(strsplit(title2, ","), function(x) trimws(x[2]), character(1))
grp2 <- ifelse(grp2 == "Papillary thyroid carcinoma", "PTC",
        ifelse(grp2 == "Normal thyroid", "Normal", grp2))
cat("  grupos:", paste(names(table(grp2)), table(grp2), collapse=", "), "\n")
X2 <- collapse_by_symbol(r2$mat)
cat("  genes após colapso:", nrow(X2), "\n")
stopifnot(ncol(X2) == length(grp2))
deg2 <- run_limma(X2, grp2, "PTC", "Normal")
cat("  DEGs GSE60542 (|logFC|>1, FDR<0.05): up=", sum(deg2$adj.P.Val<0.05 & deg2$logFC>1),
    " down=", sum(deg2$adj.P.Val<0.05 & deg2$logFC< -1), "\n")
data.table::fwrite(as.data.table(deg2), file.path(dir_v, "GSE60542_deg.tsv"), sep = "\t")

# ── Comparação com discovery (genes congelados) ───────────────────────────────
cat("\n── Comparação com discovery (genes congelados) ──\n")
disco <- fread("04_differential_expression/DEG_full_results.tsv")
key <- c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM","SOX10","TFF1","DIO3")

comp_rows <- list()
for (g in key) {
  d <- disco[gene_symbol == g]
  a <- deg1[deg1$gene_symbol == g, ]
  b <- deg2[deg2$gene_symbol == g, ]
  comp_rows[[g]] <- data.table(
    gene = g,
    discovery_logFC = if (nrow(d)) d$logFC[1] else NA_real_,
    discovery_FDR = if (nrow(d)) d$adj.P.Val[1] else NA_real_,
    GSE33630_logFC = if (nrow(a)) a$logFC[1] else NA_real_,
    GSE33630_FDR = if (nrow(a)) a$adj.P.Val[1] else NA_real_,
    GSE60542_logFC = if (nrow(b)) b$logFC[1] else NA_real_,
    GSE60542_FDR = if (nrow(b)) b$adj.P.Val[1] else NA_real_
  )
}
comp <- rbindlist(comp_rows)
comp[, direction_33630 := sign(GSE33630_logFC) == sign(discovery_logFC)]
comp[, direction_60542 := sign(GSE60542_logFC) == sign(discovery_logFC)]
comp[, rep_33630 := direction_33630 & GSE33630_FDR < 0.05 & !is.na(direction_33630)]
comp[, rep_60542 := direction_60542 & GSE60542_FDR < 0.05 & !is.na(direction_60542)]
comp[, classification := fifelse(rep_33630 & rep_60542, "INDEPENDENTLY REPLICATED",
    fifelse(rep_33630 | rep_60542, "PARTIALLY REPLICATED",
            fifelse(direction_33630 & direction_60542, "DIRECTION CONSISTENT (NS)", "NOT REPLICATED/CONTRADICTED")))]
print(comp)
data.table::fwrite(comp, file.path(dir_v, "validation_key_genes.tsv"), sep = "\t")

# ── Concordância global de direção (todos os genes discovery DEG) ─────────────
cat("\n── Concordância global de direção ──\n")
disco_sig <- disco[regulation %in% c("Up","Down")]
d1 <- as.data.table(deg1)[, .(gene_symbol, logFC, adj.P.Val)]
d2 <- as.data.table(deg2)[, .(gene_symbol, logFC, adj.P.Val)]
m1 <- merge(disco_sig[, .(gene_symbol, logFC, adj.P.Val, regulation)], d1, by = "gene_symbol", suffixes = c(".disc", ".geo"))
cat("  GSE33630: ", nrow(m1), "genes discovery-DEG testados | direção concordante=",
    mean(sign(m1$logFC.disc) == sign(m1$logFC.geo)), "\n")
m2 <- merge(disco_sig[, .(gene_symbol, logFC, adj.P.Val, regulation)], d2, by = "gene_symbol", suffixes = c(".disc", ".geo"))
cat("  GSE60542: ", nrow(m2), "genes discovery-DEG testados | direção concordante=",
    mean(sign(m2$logFC.disc) == sign(m2$logFC.geo)), "\n")
cat("  Spearman logFC discovery×GSE33630 =", cor(m1$logFC.disc, m1$logFC.geo, method="spearman"), "\n")
cat("  Spearman logFC discovery×GSE60542 =", cor(m2$logFC.disc, m2$logFC.geo, method="spearman"), "\n")

summ <- data.table(
  dataset = c("GSE33630", "GSE60542"),
  n_discovery_deg_tested = c(nrow(m1), nrow(m2)),
  direction_concordance = c(mean(sign(m1$logFC.disc)==sign(m1$logFC.geo)),
                            mean(sign(m2$logFC.disc)==sign(m2$logFC.geo))),
  spearman_logFC = c(cor(m1$logFC.disc, m1$logFC.geo, method="spearman"),
                     cor(m2$logFC.disc, m2$logFC.geo, method="spearman"))
)
data.table::fwrite(summ, file.path(dir_v, "validation_concordance_summary.tsv"), sep = "\t")

cat("\n══ FASE 18/19 — Validação externa CONCLUÍDA ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
