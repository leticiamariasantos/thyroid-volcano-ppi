# ═══════════════════════════════════════════════════════════════════════════════
# 22_validation_enhance.R — PARTE 8–22: metadata, composição (stress test),
#                            candidatos integrados, proveniência.
#
# Reutiliza as matrizes já processadas; produz:
#   results/validation/GSE33630/metadata.tsv
#   results/validation/GSE60542/metadata.tsv
#   results/validation/composition_stress.tsv
#   results/validation/integrated_candidates.tsv
#   results/validation/provenance.tsv
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table); library(limma); library(hgu133plus2.db); library(AnnotationDbi)})

cat("══ PARTE 8–22 — Validação GEO (metadata, composição, integração) ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_v <- "results/validation"
dir.create(dir_v, recursive = TRUE, showWarnings = FALSE)
for (s in c("GSE33630","GSE60542")) dir.create(file.path(dir_v, s, "metadata"), recursive = TRUE, showWarnings = FALSE)

# ── helpers (idênticas ao script 17) ──────────────────────────────────────────
read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  tb <- which(lines == "!series_matrix_table_begin")
  te <- which(lines == "!series_matrix_table_end")
  title_line <- lines[grep("^!Sample_title", lines)][1]
  title_fields <- gsub('"', '', strsplit(title_line, "\t")[[1]][-1])
  mat <- read.delim(textConnection(lines[(tb+1):(te-1)]), check.names = FALSE, stringsAsFactors = FALSE)
  list(mat = mat, titles = title_fields, lines = lines)
}
collapse_by_symbol <- function(mat, idcol = "ID_REF") {
  probes <- mat[[idcol]]
  map <- AnnotationDbi::select(hgu133plus2.db, keys = probes, keytype = "PROBEID", columns = "SYMBOL")
  keep <- mat[[idcol]] %in% map$PROBEID
  mat <- mat[keep, , drop = FALSE]
  sym <- map$SYMBOL[match(mat[[idcol]], map$PROBEID)]
  X <- as.matrix(mat[, -which(colnames(mat) == idcol), drop = FALSE])
  storage.mode(X) <- "numeric"
  gm <- rowMeans(X, na.rm = TRUE)
  ord <- order(-gm); X <- X[ord, , drop = FALSE]; sym <- sym[ord]
  keep2 <- !duplicated(sym)
  X <- X[keep2, , drop = FALSE]; rownames(X) <- sym[keep2]
  X
}

muscle <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
            "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
            "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
            "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")
immune <- c("PTPRC","CD3D","CD3E","CD8A","CD4","CD19","MS4A1","CD68","CD14","NCAM1",
            "FCGR3A","ITGAX","NKG7","GNLY","GZMB","PRF1","CD163","CD79A")

marker_score <- function(X, mk) {
  g <- intersect(mk, rownames(X))
  if (length(g) == 0) return(rep(NA_real_, ncol(X)))
  colMeans(X[g, , drop = FALSE])
}

# ── GSE33630 ──────────────────────────────────────────────────────────────────
cat("── GSE33630 ──\n")
r1 <- read_series_matrix("10_validation/raw/GSE33630_matrix.txt.gz")
char1 <- gsub('"', '', strsplit(r1$lines[grep("^!Sample_characteristics_ch1", r1$lines)][1], "\t")[[1]][-1])
grp1 <- ifelse(grepl("papillary", char1), "PTC",
        ifelse(grepl("non-tumor", char1), "Normal", "ATC"))
# patient id (pareado): prefixo numérico do título
pat1 <- sub("^(ATC|[0-9]+).*", "\\1", r1$titles)
meta1 <- data.table(sample_id = r1$titles, group = grp1,
                    tissue = ifelse(grp1=="Normal","non-tumor thyroid","tumor"),
                    subtype = ifelse(grp1=="PTC","papillary thyroid carcinoma (PTC)",
                              ifelse(grp1=="ATC","anaplastic thyroid carcinoma (ATC)","normal")),
                    paired = ifelse(grp1 %in% c("PTC","Normal"), "patient-matched", "no"),
                    patient_id = ifelse(grepl("^ATC", r1$titles), NA_character_, pat1),
                    platform = "GPL570", source = "GEO GSE33630")
data.table::fwrite(meta1, file.path(dir_v, "GSE33630", "metadata", "metadata.tsv"), sep = "\t")

X1 <- collapse_by_symbol(r1$mat)
# stress test de composição
sel1 <- grp1 %in% c("PTC","Normal")
sc_m1 <- marker_score(X1[, sel1, drop=FALSE], muscle)
sc_i1 <- marker_score(X1[, sel1, drop=FALSE], immune)
comp1 <- data.table(dataset = "GSE33630",
  muscle_PTC = mean(sc_m1[grp1[sel1]=="PTC"]), muscle_Normal = mean(sc_m1[grp1[sel1]=="Normal"]),
  muscle_diff = mean(sc_m1[grp1[sel1]=="PTC"]) - mean(sc_m1[grp1[sel1]=="Normal"]),
  immune_PTC = mean(sc_i1[grp1[sel1]=="PTC"]), immune_Normal = mean(sc_i1[grp1[sel1]=="Normal"]))
rm(X1); invisible(gc())

# ── GSE60542 ──────────────────────────────────────────────────────────────────
cat("── GSE60542 ──\n")
r2 <- read_series_matrix("10_validation/raw/GSE60542_matrix.txt.gz")
t2 <- r2$titles
parts <- strsplit(t2, ",")
grp2 <- vapply(parts, function(x) trimws(x[2]), character(1))
grp2 <- ifelse(grp2 == "Papillary thyroid carcinoma", "PTC",
        ifelse(grp2 == "Normal thyroid", "Normal", grp2))
pat2 <- vapply(parts, function(x) trimws(x[1]), character(1))
meta2 <- data.table(sample_id = t2, patient_id = pat2, group = grp2,
                    tumor_status = ifelse(grp2=="PTC","tumor","normal"),
                    nodal_status = vapply(parts, function(x) trimws(x[3]), character(1)),
                    tissue = ifelse(grp2 %in% c("PTC","Normal"), "thyroid", "other"),
                    metastasis = grp2 %in% c("Lymph node metastasis","Pleural metastasis"),
                    paired = ifelse(grp2 %in% c("PTC","Normal"), "patient-matched (parcial)", "no"),
                    platform = "GPL570", source = "GEO GSE60542")
data.table::fwrite(meta2, file.path(dir_v, "GSE60542", "metadata", "metadata.tsv"), sep = "\t")

X2 <- collapse_by_symbol(r2$mat)
sel2 <- grp2 %in% c("PTC","Normal")
sc_m2 <- marker_score(X2[, sel2, drop=FALSE], muscle)
sc_i2 <- marker_score(X2[, sel2, drop=FALSE], immune)
comp2 <- data.table(dataset = "GSE60542",
  muscle_PTC = mean(sc_m2[grp2[sel2]=="PTC"]), muscle_Normal = mean(sc_m2[grp2[sel2]=="Normal"]),
  muscle_diff = mean(sc_m2[grp2[sel2]=="PTC"]) - mean(sc_m2[grp2[sel2]=="Normal"]),
  immune_PTC = mean(sc_i2[grp2[sel2]=="PTC"]), immune_Normal = mean(sc_i2[grp2[sel2]=="Normal"]))
rm(X2); invisible(gc())

comp <- rbind(comp1, comp2)
comp[, muscle_diff := round(muscle_diff, 3)]
cat("\n── Composição (stress test) nos GEO ──\n")
print(comp)
data.table::fwrite(comp, file.path(dir_v, "composition_stress.tsv"), sep = "\t")

# ── Candidatos integrados ─────────────────────────────────────────────────────
cat("\n── Candidatos integrados ──\n")
disco <- fread("04_differential_expression/DEG_full_results.tsv")
voo <- fread("results/counts_deg/voom_full_results.tsv")
dds <- fread("results/counts_deg/deseq2_full_results.tsv")
g1 <- fread("10_validation/GSE33630_deg.tsv")
g2 <- fread("10_validation/GSE60542_deg.tsv")
fin <- fread("results/composition/final_sensitivity.tsv")
key <- fin$gene

integ <- rbindlist(lapply(key, function(g) {
  d <- disco[gene_symbol==g]; v <- voo[gene_symbol==g]; dd <- dds[gene_symbol==g]
  a <- g1[gene_symbol==g]; b <- g2[gene_symbol==g]
  d_lfc <- if(nrow(d)) d$logFC[1] else NA_real_
  a_lfc <- if(nrow(a)) a$logFC[1] else NA_real_; a_fdr <- if(nrow(a)) a$adj.P.Val[1] else NA_real_
  b_lfc <- if(nrow(b)) b$logFC[1] else NA_real_; b_fdr <- if(nrow(b)) b$adj.P.Val[1] else NA_real_
  dir_a <- sign(a_lfc)==sign(d_lfc); dir_b <- sign(b_lfc)==sign(d_lfc)
  rep_a <- !is.na(dir_a) & dir_a & a_fdr<0.05; rep_b <- !is.na(dir_b) & dir_b & b_fdr<0.05
  cl <- if (rep_a & rep_b) "INDEPENDENTLY REPLICATED"
        else if (rep_a | rep_b) "PARTIALLY REPLICATED"
        else if (!is.na(dir_a) & !is.na(dir_b) & dir_a & dir_b) "DIRECTION CONSISTENT (NS)"
        else if (!is.na(dir_a) & dir_a & !dir_b | !is.na(dir_b) & dir_b & !dir_a) "PARTIALLY REPLICATED"
        else "NOT REPLICATED"
  data.table(gene=g, discovery_logFC=d_lfc,
             voom_logFC=if(nrow(v)) v$logFC[1] else NA_real_,
             DESeq2_logFC=if(nrow(dd)) dd$log2FoldChange[1] else NA_real_,
             GSE33630_effect=a_lfc, GSE33630_direction=dir_a,
             GSE60542_effect=b_lfc, GSE60542_direction=dir_b,
             composition_status=fin[gene==g]$status,
             PPI_status=ifelse(g %in% c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1"), "hub/ECM-immune", "—"),
             final_status=cl)
}))
integ[, discovery_logFC := round(discovery_logFC, 3)]
integ[, voom_logFC := round(voom_logFC, 3)]
integ[, DESeq2_logFC := round(DESeq2_logFC, 3)]
integ[, GSE33630_effect := round(GSE33630_effect, 3)]
integ[, GSE60542_effect := round(GSE60542_effect, 3)]
print(integ)
data.table::fwrite(integ, file.path(dir_v, "integrated_candidates.tsv"), sep = "\t")

# ── Proveniência (SHA256 dos raw) ─────────────────────────────────────────────
cat("\n── Proveniência (SHA256 raw) ──\n")
prov <- rbindlist(lapply(c("GSE33630","GSE60542"), function(g) {
  f <- sprintf("10_validation/raw/%s_matrix.txt.gz", g)
  data.table(accession = g, file = f,
             size_bytes = file.info(f)$size,
             sha256 = system2("sha256sum", f, stdout = TRUE))
}))
data.table::fwrite(prov, file.path(dir_v, "provenance.tsv"), sep = "\t")
print(prov)

cat("\n══ PARTE 8–22 CONCLUÍDA ══\n")
