# ═══════════════════════════════════════════════════════════════════════════════
# 29_ranking_figures.R — Figuras do ranking composto
# ═══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

dir_f <- "results/prioritization/figures"
dir.create(dir_f, recursive = TRUE, showWarnings = FALSE)
cat("══ Figuras do ranking ══\n")

U <- fread("results/prioritization/integrated_gene_ranking.tsv")

# 1. Top 20 molecular
top20 <- U[order(-FINAL_MOLECULAR)][1:20]
top20[, gene := factor(gene, levels = rev(top20$gene))]
p1 <- ggplot(top20, aes(FINAL_MOLECULAR, gene, fill = TIER)) +
  geom_col() + labs(title="Ranking molecular integrado (top 20)", x="Molecular evidence score", y=NULL) +
  scale_fill_manual(values=c("TIER 1"="#d62728","TIER 2"="#ff7f0e","TIER 3"="#2ca02c","TIER 4"="#1f77b4","TIER 5 (composicional)"="grey70")) +
  theme_minimal()
ggsave(file.path(dir_f,"Fig_ranking_molecular_top20.png"), p1, width=7, height=6, dpi=150)

# 2. Scatter molecular x translacional
p2 <- ggplot(U, aes(FINAL_MOLECULAR, FINAL_TRANSLATIONAL, color=is_muscle)) +
  geom_point(alpha=0.6) +
  geom_text(data=U[gene %in% c("FN1","ITGA2","CCND1","CTSS","HLA-DPA1","LAMA2","ITGA2B")], aes(label=gene), vjust=-1, size=3) +
  scale_color_manual(values=c("FALSE"="#1f77b4","TRUE"="grey70"), labels=c("não-muscular","muscular")) +
  labs(title="Robustez molecular × prioridade translacional", x="Molecular", y="Translacional") +
  theme_minimal()
ggsave(file.path(dir_f,"Fig_scatter_molecular_translational.png"), p2, width=7, height=5, dpi=150)

# 3. Scatter molecular x nanomédica
p3 <- ggplot(U, aes(FINAL_MOLECULAR, FINAL_NANOMEDICINE, color=is_muscle)) +
  geom_point(alpha=0.6) +
  geom_text(data=U[gene %in% c("FN1","ITGA2","LAMA2","ITGA2B","SDC4")], aes(label=gene), vjust=-1, size=3) +
  scale_color_manual(values=c("FALSE"="#1f77b4","TRUE"="grey70")) +
  labs(title="Robustez molecular × plausibilidade nanomédica", x="Molecular", y="Nanomedicina (hipótese)") +
  theme_minimal()
ggsave(file.path(dir_f,"Fig_scatter_molecular_nano.png"), p3, width=7, height=5, dpi=150)

# 4. Heatmap dimensões (top 15 × dimensões)
top15 <- U[order(-FINAL_MOLECULAR)][1:15]
hm <- melt(top15[, .(gene, Transcriptômica=TRANSCRIPTOMIC_ROBUSTNESS, Replicação=ifelse(is.na(EXTERNAL_REPLICATION),0.5,EXTERNAL_REPLICATION),
                     Proteína=ifelse(is.na(PROTEIN_SUPPORT),0.5,PROTEIN_SUPPORT), Via=PATHWAY_CONTEXT, PPI=PPI_SUPPORT, Acessibilidade=TARGET_ACCESSIBILITY)],
           id.vars="gene", variable.name="dimension", value.name="score")
hm[, gene := factor(gene, levels = rev(top15$gene))]
p4 <- ggplot(hm, aes(dimension, gene, fill=score)) + geom_tile() +
  scale_fill_gradient2(low="#1f77b4", mid="white", high="#d62728", midpoint=0.5) +
  labs(title="Matriz de evidências (top 15)", x=NULL, y=NULL) +
  theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(dir_f,"Fig_evidence_heatmap.png"), p4, width=9, height=6, dpi=150)

cat("══ Figuras do ranking CONCLUÍDAS ══\n")
