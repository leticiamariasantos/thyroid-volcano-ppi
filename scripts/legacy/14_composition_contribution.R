# ═══════════════════════════════════════════════════════════════════════════════
# 14_composition_contribution.R — Quantifica contribuição da composição
#
# Responde (FASE 9/10): "Quanto do padrão transcriptômico pode ser compatível
# com diferenças de composição tecidual?" usando marcadores canônicos.
#
# NOTA METODOLÓGICA: marker-based scoring é uma ESTIMATIVA/MODELO de composição,
# não observação experimental direta (sem CIBERSORTx/assinatura de referência).
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(data.table))

cat("══ FASE 9/10 — Contribuição da composição ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_out <- "results/composition"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

markers <- list(
  muscle_skeletal = c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7",
                      "ACTA1","ACTC1","TNNT1","TNNT3","TNNC1","TNNI1","TNNI2",
                      "TPM1","TPM2","TPM3","CKM","CKMT2","DES","MB","ENO3",
                      "MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
                      "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2"),
  thyroid_epithelial = c("TG","TPO","TSHR","PAX8","NKX2-1","FOXE1","SLC5A5",
                         "DIO1","DIO2","CALCA","TFF3","KRT7","KRT19"),
  fibroblast_ecm = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","DCN",
                     "LUM","VIM","ACTA2","FAP","PDGFRB","POSTN","MMP2","LOX","BGN"),
  endothelial = c("PECAM1","VWF","CDH5","ENG","KDR","FLT1","EMCN","CLDN5","TEK"),
  immune = c("PTPRC","CD3D","CD3E","CD8A","CD4","CD19","MS4A1","CD68","CD14",
             "NCAM1","FCGR3A","ITGAX","NKG7","GNLY","GZMB","PRF1","CD163","CD79A")
)

deg <- fread("04_differential_expression/DEG_full_results.tsv")
deg_sig <- deg[regulation %in% c("Up","Down")]

out <- list()
for (cat in names(markers)) {
  mk <- intersect(markers[[cat]], deg$gene_symbol)
  mk_sig <- deg_sig[gene_symbol %in% mk]
  n_up <- sum(mk_sig$regulation == "Up")
  n_down <- sum(mk_sig$regulation == "Down")
  out[[cat]] <- data.table(
    category = cat,
    markers_total = length(markers[[cat]]),
    markers_tested = length(mk),
    markers_deg = nrow(mk_sig),
    up = n_up,
    down = n_down,
    frac_of_DEGs = nrow(mk_sig) / nrow(deg_sig)
  )
}
res <- rbindlist(out)
res[, frac_of_DEGs := round(frac_of_DEGs, 4)]
print(res)
data.table::fwrite(res, file.path(dir_out, "composition_contribution.tsv"), sep = "\t")

# Muscle genes específicos levantados anteriormente
muscle_of_interest <- c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")
cat("\n  Genes musculares de interesse (logFC TPM + voom + DESeq2):\n")
comp <- fread("results/counts_deg/method_concordance.tsv")
print(comp[gene %in% muscle_of_interest,
           .(gene, logFC_TPM, logFC_voom, logFC_DESeq2, robustness)])

# Direção do sinal muscular: concentrado em GTEx normal?
scores <- fread(file.path(dir_out, "composition_per_sample.tsv"))
cat("\n  Score muscular: THCA vs Normal (log2 TPM)\n")
print(scores[, .(mean_thca = mean(muscle_skeletal[condition=="THCA"]),
                 mean_normal = mean(muscle_skeletal[condition=="Normal"]))])

cat("\n══ FASE 9/10 — Contribuição da composição CONCLUÍDA ══\n")
