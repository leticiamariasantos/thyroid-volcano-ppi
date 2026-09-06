# ═══════════════════════════════════════════════════════════════════════════════
# 17_validate.R — Validação técnica final da FASE 2 (reboot)
#
# Verifica existência, dimensões, consistência e integridade de todos os outputs
# da nova fase. Objetivo: 0 failures. Substitui o antigo 15_validate_outputs.R.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({ library(data.table) })

cat("══ VALIDAÇÃO TÉCNICA — FASE 2 REBOOT ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
chk <- new_checker()
ck <- chk$ck

# ── 1. Arquivos críticos existem e não são vazios ──────────────────────────────
required <- c(
  "results/phase2/audit/DATA_AUDIT.txt",
  "results/phase2/audit/sample_metadata.tsv",
  "results/phase2/preprocessing/QC_summary.txt",
  "results/phase2/preprocessing/pca_scores.tsv",
  "results/phase2/differential_expression/DEG_summary.tsv",
  "results/phase2/differential_expression/limma_full_results.tsv",
  "results/phase2/differential_expression/voom_full_results.tsv",
  "results/phase2/differential_expression/deseq2_full_results.tsv",
  "results/phase2/differential_expression/method_concordance.tsv",
  "results/phase2/pathways/PANEL_30_PATHWAYS.tsv",
  "results/phase2/pathways/PANEL_30_PATHWAYS_with_genes.tsv",
  "results/phase2/pathways/PANEL_SELECTION_RATIONALE.md",
  "results/phase2/pathways/PANEL_30_GSEA_RESULTS.tsv",
  "results/phase2/pathways/PATHWAY_REDUNDANCY.tsv",
  "results/phase2/pathways/PATHWAY_MODULES.tsv",
  "results/phase2/pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv",
  "results/phase2/gsea/GSEA_GLOBAL_KEGG.tsv",
  "results/phase2/gsea/GSEA_GLOBAL_REACTOME.tsv",
  "results/phase2/gsea/GSEA_GLOBAL_HALLMARK.tsv",
  "results/phase2/sensitivity/gsea_composition_sensitivity.tsv",
  "results/phase2/sensitivity/muscle_composition_scores.tsv",
  "results/phase2/ppi/PPI_centrality.tsv",
  "results/phase2/ppi/PPI_edges.tsv",
  "results/phase2/ppi/PPI_summary.tsv",
  "results/phase2/ppi/PPI_eligible_genes.tsv",
  "results/phase2/validation/validation_candidates.tsv",
  "results/phase2/validation/GSE224356_candidates.tsv",
  "results/phase2/validation/singlecell_predominant_compartment.tsv",
  "results/phase2/validation/RPPA_summary.tsv",
  "results/phase2/validation/mutation_frequency.tsv",
  "results/phase2/validation/convergence_candidates.tsv",
  "results/phase2/validation/panel_convergence_genes.tsv"
)
for (f in required) ck(file.exists(f) && file.info(f)$size > 0, paste("existe e não-vazio:", f))

# ── 2. Painel: exatamente 30 vias (10 + 20) ────────────────────────────────────
panel <- fread("results/phase2/pathways/PANEL_30_PATHWAYS.tsv")
ck(nrow(panel) == 30, sprintf("painel = 30 vias (atual: %d)", nrow(panel)))
ck(sum(panel$panel_group == "ORIGINAL_10") == 10, "10 vias ORIGINAIS")
ck(sum(panel$panel_group == "ADDITIONAL_20") == 20, "20 vias ADICIONAIS")
expected_orig <- c("hsa05216","hsa04919","hsa04010","hsa04151","hsa04150",
                   "hsa04115","hsa04210","hsa04110","hsa04310","hsa04064")
ck(setequal(panel[panel_group == "ORIGINAL_10"]$database_id, expected_orig),
   "painel original = 10 vias KEGG corretas (inalterado)")
ck(!anyDuplicated(panel$database_id), "sem IDs duplicados no painel")
ck(all(panel$geneset_size > 0), "todos os gene sets não-vazios")

# ── 3. DEG consistency ─────────────────────────────────────────────────────────
deg_sum <- fread("results/phase2/differential_expression/DEG_summary.tsv")
limma <- fread("results/phase2/differential_expression/limma_full_results.tsv")
voom  <- fread("results/phase2/differential_expression/voom_full_results.tsv")
dds   <- fread("results/phase2/differential_expression/deseq2_full_results.tsv")
ck(deg_sum[method == "limma", genes_tested] == nrow(limma), "limma genes_tested coerente")
ck(deg_sum[method == "limma", total_deg] == sum(limma$regulation != "NS"), "limma DEG count coerente")
ck(sum(limma$regulation == "Up") == deg_sum[method == "limma", up] &&
   sum(limma$regulation == "Down") == deg_sum[method == "limma", down], "limma up/down coerentes")
ck(all(is.finite(limma$logFC)) && all(is.finite(limma$adj.P.Val)), "limma logFC/FDR finitos")
ck(all(is.finite(voom$logFC)) && all(is.finite(voom$adj.P.Val)), "voom logFC/FDR finitos")
ck(all(is.finite(dds$logFC)) && all(is.finite(dds$adj.P.Val)), "DESeq2 logFC/padj finitos")
ck(!anyDuplicated(limma$gene_symbol), "limma sem genes duplicados")

# ── 4. GSEA painel (30 vias × 3 métodos) ───────────────────────────────────────
gsea_panel <- fread("results/phase2/pathways/PANEL_30_GSEA_RESULTS.tsv")
ck(nrow(gsea_panel) == 30, sprintf("GSEA painel = 30 vias (atual: %d)", nrow(gsea_panel)))
ck(all(c("NES_limma","NES_voom","NES_deseq2","padj_limma") %in% names(gsea_panel)),
   "GSEA painel com colunas NES/padj por método")

# ── 5. GSEA global com coleções ────────────────────────────────────────────────
g_kegg <- fread("results/phase2/gsea/GSEA_GLOBAL_KEGG.tsv")
g_react <- fread("results/phase2/gsea/GSEA_GLOBAL_REACTOME.tsv")
g_hall <- fread("results/phase2/gsea/GSEA_GLOBAL_HALLMARK.tsv")
ck(nrow(g_kegg) >= 100, sprintf("GSEA KEGG global >= 100 vias (atual %d)", nrow(g_kegg)))
ck(nrow(g_react) >= 1000, sprintf("GSEA Reactome global >= 1000 (atual %d)", nrow(g_react)))
ck(nrow(g_hall) == 50, sprintf("GSEA Hallmark = 50 (atual %d)", nrow(g_hall)))
ck(sum(g_kegg$padj < 0.05, na.rm = TRUE) >= 10, "GSEA KEGG: >= 10 vias significativas")

# ── 6. Redundância e robustez ──────────────────────────────────────────────────
red <- fread("results/phase2/pathways/PATHWAY_REDUNDANCY.tsv")
ck(nrow(red) == choose(30, 2), sprintf("pares de redundância = 435 (atual %d)", nrow(red)))
rob <- fread("results/phase2/pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv")
ck(nrow(rob) == 30, "matriz de robustez = 30 vias")
ck(sum(rob$robustness_class == "ROBUSTA") >= 1, ">= 1 via ROBUSTA")

# ── 7. Sensibilidade de composição ─────────────────────────────────────────────
comp <- fread("results/phase2/sensitivity/gsea_composition_sensitivity.tsv")
ck(nrow(comp) == 30, "sensibilidade de composição = 30 vias")
ck(sum(comp$status == "ROBUSTA", na.rm = TRUE) >= 1, ">= 1 via ROBUSTA pós-composição")

# ── 8. PPI ─────────────────────────────────────────────────────────────────────
ppi_sum <- fread("results/phase2/ppi/PPI_summary.tsv")
cent <- fread("results/phase2/ppi/PPI_centrality.tsv")
ck(as.integer(ppi_sum[metric == "nodes", value]) == nrow(cent), "PPI nós == centralidade")

# ── 9. Candidatos na convergência ──────────────────────────────────────────────
conv <- fread("results/phase2/validation/convergence_candidates.tsv")
ck(all(c("ITGA2","FN1","CCND1") %in% conv$gene), "ITGA2/FN1/CCND1 na convergência")
val <- fread("results/phase2/validation/validation_candidates.tsv")
ck(all(c("ITGA2","FN1","CCND1") %in% val$gene), "ITGA2/FN1/CCND1 na validação externa")
# convergência do painel: CCND1 como gene de core enrichment (correção do separador |)
pconv <- fread("results/phase2/validation/panel_convergence_genes.tsv")
ck("CCND1" %in% pconv$gene, "CCND1 presente na convergência do painel")
ck(conv[gene == "CCND1", n_core_enrichment] >= 1, "CCND1 é core enrichment (>=1 via)")

# ── 10. Figuras ────────────────────────────────────────────────────────────────
figs <- c("Fig2_Volcano.png","Fig3_Heatmap.png","Fig4_GSEA_global.png",
          "Fig5_Dotplot30.png","Fig6_NES_heatmap.png","Fig7_Composition.png",
          "Fig8_PPI.png","Fig9_SingleCell.png","Fig10_Conceptual_Model.png",
          "QC_PCA_condition.png","PATHWAY_REDUNDANCY_heatmap.png")
for (f in figs) ck(file.exists(file.path("results/phase2/figures", f)), paste("figura:", f))

# ── 11. Relatórios ─────────────────────────────────────────────────────────────
reps <- c("PHASE2_FINAL_REPORT.md","PHASE2_METHODS.md","PHASE2_RESULTS.md",
          "PHASE2_LIMITATIONS.md","PHASE2_FALSIFICATION_AUDIT.md","PHASE2_REPRODUCIBILITY.md")
for (f in reps) ck(file.exists(file.path("results/phase2/reports", f)), paste("relatório:", f))

# ── 12. Dados brutos preservados ───────────────────────────────────────────────
ck(file.exists("data/global/TCGA_GTEx_thyroid_tpm.tsv"), "dados brutos TPM preservados")
ck(file.exists("data/global/TCGA_GTEx_thyroid_counts.tsv"), "dados brutos counts preservados")
ck(file.exists("renv.lock"), "renv.lock preservado")
ck(file.exists("Dockerfile"), "Dockerfile preservado")

fail <- chk$total()
cat(sprintf("\n══ VALIDAÇÃO: %d falha(s) ══\n", fail))
if (fail > 0) { cat("  → FALHAS PRESENTES\n"); quit(status = 1) } else cat("  → VALIDAÇÃO OK (0 failures)\n")
