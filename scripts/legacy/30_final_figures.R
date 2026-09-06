# ═══════════════════════════════════════════════════════════════════════════════
# 30_final_figures.R — Figuras finais: single-cell localization + evidence ladder
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

dir_f <- "14_figures"
dir.create(dir_f, recursive = TRUE, showWarnings = FALSE)
cat("══ Figuras finais (single-cell + evidence ladder) ══\n")

# ── 1. Single-cell: FN1 e ITGA2 por cell type ─────────────────────────────────
sc <- fread("results/single_cell/single_cell_candidate_localization.tsv")
sc <- sc[gene %in% c("FN1","ITGA2") & dataset=="GSE182416 (normal thyroid)"]
sc <- sc[!grepl("T/NK/B/Myeloid", cell_type)]  # remove agrupado
sc[, cell_type := factor(cell_type, levels=c("Epithelial cells","Fibroblasts","Endothelial cells","SMCs/Pericytes","Myeloid cells","T/NK cells","B cells"))]
p1 <- ggplot(sc, aes(cell_type, mean_expr, fill=gene)) +
  geom_col(position="dodge") +
  labs(title="FN1 e ITGA2 em tireoide normal (GSE182416, 54.726 células)",
       x=NULL, y="expressão média (raw counts)") +
  theme_minimal() + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(dir_f,"Fig_single_cell_localization.png"), p1, width=8, height=5, dpi=150)

# ── 2. Evidence ladder (conceitual) ───────────────────────────────────────────
steps <- data.frame(
  step = c("TCGA/GTEx DE","Múltiplos métodos (limma/voom/DESeq2)","Validação microarray (GSE33630/GSE60542)",
           "Validação RNA-seq (GSE224356)","Proteína (RPPA)","Multiômica (CNV/mutação)","Composição controlada",
           "Single-cell localization","Hipótese translacional","Hipótese nanomédica"),
  level = c("EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","EVIDÊNCIA","HIPÓTESE","HIPÓTESE"),
  y = 10:1
)
steps$level <- factor(steps$level, levels=c("EVIDÊNCIA","HIPÓTESE"))
p2 <- ggplot(steps, aes(1, y)) +
  geom_segment(aes(x=1, xend=1, y=1, yend=10), color="grey80") +
  geom_point(aes(color=level), size=5) +
  geom_text(aes(x=1.12, label=step), hjust=0, size=3.2) +
  geom_hline(yintercept=2.5, linetype="dashed", color="#d62728") +
  annotate("text", x=1, y=2.2, label="← evidência termina | hipótese começa →", color="#d62728", size=3) +
  scale_color_manual(values=c("EVIDÊNCIA"="#1f77b4","HIPÓTESE"="#ff7f0e")) +
  xlim(1, 3.2) + labs(title="Escada de evidência: da observação à hipótese", x=NULL, y=NULL) +
  theme_void() + theme(legend.position="none")
ggsave(file.path(dir_f,"Fig_evidence_ladder.png"), p2, width=10, height=6, dpi=150)

cat("══ Figuras finais CONCLUÍDAS ══\n")
