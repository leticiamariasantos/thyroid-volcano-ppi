# PHASE 2 (reboot) — Métodos

**Data da execução:** 2026-09-06
**Natureza:** reconstrução independente do pipeline com expansão do painel de vias de 10 → 30.
**Ambiente:** R 4.6.1 (ucrt); pacotes-chave: limma, edgeR, DESeq2, fgsea (1.37.4),
msigdbr (2026.1.Hs), KEGGREST, hgu133plus2.db, igraph, data.table, ggplot2.

## 1. Dados

| Camada | Arquivo | Genes | Amostras | Unidade |
|---|---|---|---|---|
| Principal (TPM) | `data/global/TCGA_GTEx_thyroid_tpm.tsv` | 58.581 | 783 (504 TCGA-THCA, 279 GTEx thyroid) | log2(TPM+0.001) |
| Contagens (sensibilidade) | `data/global/TCGA_GTEx_thyroid_counts.tsv` | 56.937 | 782 (504 TCGA, 278 GTEx) | raw counts (recount3, GENCODE) |

Auditoria (`results/phase2/audit/DATA_AUDIT.txt`): sem genes/amostras duplicadas, sem
NA/Inf, contagens inteiras e não-negativas, 782 amostras comuns entre matrizes. O `.gz`
global (60.498 genes × 19.131 amostras) é superconjunto; o processado é o subconjunto
tireoide (colapso documentado em `01_data/MANIFEST_TPM.tsv`).

## 2. QC e pré-processamento

- **Distribuição:** média −4,85; mediana −9,97; mínimo −9,97 (≈log2(0.001)); máximo 17,54;
  53,1% de zeros. Confirma escala log2(TPM+0.001).
- **Genes pouco expressos:** 9.037 genes nunca detectados.
- **PCA:** 2.000 genes mais variáveis (MAD). PC1 = 21,6% da variância; R²(PC1~condição) = 0,862.
- **Outliers:** 1 amostra (Mahalanobis p<0,001 sobre PC1–5).
- **Clustering/distância:** heatmap de correlação de Pearson (amostragem 120).

**Confundimento estrutural documentado:** `source ≡ condition` (TCGA ≡ tumor, GTEx ≡ normal).
Nenhuma correção de batch é alegada como eliminadora; a limitação é tratada por triangulação
e análises de sensibilidade (não por alegação de correção completa).

## 3. Expressão diferencial

- **limma (principal):** `lmFit`/`eBayes` sobre log2 TPM, filtro de baixa expressão
  (TPM>0,1 em ≥25% das amostras) → 26.011 genes testados. Contraste THCA vs Normal.
- **limma-voom:** `filterByExpr` (CPM>1) → 39.702 genes; `voom` + `lmFit` + `eBayes`.
- **DESeq2:** mesmas contagens filtradas; `DESeq` (Wald), `independentFiltering=TRUE`.
- **Thresholds (pré-especificados):** |log2FC| > 1,0 e FDR (Benjamini-Hochberg) < 0,05.

## 4. Painel de 30 vias (definido a priori)

- 10 vias originais (INALTERADAS) + 20 vias adicionais selecionadas **antes** da
  interpretação dos resultados (critérios em `pathways/PANEL_SELECTION_RATIONALE.md`).
- Gene sets: KEGG via `KEGGREST::keggGet` (símbolos parseados do campo GENE); Reactome via
  GMT oficial (`ReactomePathways.gmt`).
- `pathways/PANEL_30_PATHWAYS.tsv` (painel), `PANEL_30_PATHWAYS_with_genes.tsv` (com genes).

## 5. GSEA

- **Painel (30 vias):** `fgseaMultilevel` (minSize=1, maxSize=2000, nPermSimple=10.000,
  seed=42) sobre rankings de t (limma, voom) e stat (DESeq2). padj = BH dentro das 30 vias.
- **Global:** `fgseaMultilevel` sobre coleções KEGG (MSigDB CP:KEGG_LEGACY, 186 sets),
  Reactome (GMT, 2.868 sets) e Hallmark (50 sets). padj = BH dentro de cada coleção.

## 6. Múltiplos testes (separados)

- A. painel original (10 vias) — hipótese a priori;
- B. 20 vias adicionais — expansão pré-especificada (linguagem exploratória);
- C. GSEA global — descoberta exploratória.
Correção: Benjamini-Hochberg (padj do fgsea), reportado o nº de hipóteses por grupo.

## 7. Redundância entre vias

Jaccard de gene sets (matriz 30×30), clustering hierárquico (average, corte h=0,85) e
heatmap. `pathways/PATHWAY_REDUNDANCY.tsv` e `PATHWAY_MODULES.tsv`.

## 8. Composição celular

Score muscular (média de 48 marcadores estriados: MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM, …).
Sensibilidade: remoção dos 48 marcadores → re-limma + re-fgsea. Comparação
FULL vs REMOVED (`sensitivity/gsea_composition_sensitivity.tsv`).

## 9. PPI

Genes elegíveis (tumor-relevante, pré-especificado): DEGs limma com FDR<0,05 e |logFC|≥1,
excluindo os 48 marcadores de músculo estriado (artefato de composição do tecido normal),
equilibrados em 250 Up + 250 Down por |t|. STRING (taxon 9606, escore≥700). Centralidade
(degree, betweenness, closeness, eigenvector) e comunidades (walktrap).

## 10. Validação externa

GSE33630 (PTC vs normal pareado, GPL570) e GSE60542 (PTC vs normal, GPL570) — mapeamento
probe→símbolo via `hgu133plus2.db`; DE dos candidatos por limma. GSE224356 (listas de DEGs
pré-computadas). GSE224357 NÃO tratada como coorte independente (mesma SuperSeries).

## 11. Single-cell

GSE232237 — resultados marker-based pré-computados (`GSE232237_marker_celltype_results.tsv`),
por compartimento (Epithelial_tumor, Fibroblast, Endothelial, T_NK, B_cell, Myeloid,
SMC_pericyte, Mast). Classificação marker-based documentada (sem anotação autoral).

## 12. Proteína e mutação

cBioPortal (THCA PanCan Atlas 2018): mutações (BRAF, NRAS, HRAS, KRAS, TP53, TERT, ITGA2,
FN1, CCND1) e RPPA (CCND1/Cyclin D1). ITGA2 e FN1 documentados como fora do painel RPPA.
