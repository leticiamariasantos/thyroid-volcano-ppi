# DISCOVERY FREEZE — Congelamento do Discovery

**Data do congelamento:** 2026-09-06
**Versão do estudo:** v3.1.0 (Fase 1) → Fase 2 global (discovery) → validação avançada
**Commit Git (discovery):** f4c0e1f (base versionada; Fase 2 em working tree não commitada)
**Regra:** a partir deste ponto, **É PROIBIDO** alterar candidatos discovery com base nos
resultados GEO. Os datasets externos (GSE33630, GSE60542) serão usados **exclusivamente**
para validação independente.

---

## 1. Dataset discovery

| Item | Valor |
|---|---|
| Dataset | `TcgaTargetGtex_rsem_gene_tpm.gz` (UCSC Xena TOIL recompute) |
| Matriz | 58.581 genes × 783 amostras (504 THCA + 279 GTEx) |
| Unidade | log2(TPM + 0.001) |
| Fonte counts | recount3 (56.937 genes × 782 amostras; STAR G026) |

## 2. Contraste

- **Principal:** THCA (tumor) − Normal (GTEx).
- **Métodos:** limma sobre TPM (trend=TRUE); limma-voom (TMM); DESeq2 (raw counts).

## 3. Métodos e parâmetros (congelados)

- Filtro: TPM>1 (ou CPM>1) em ≥10% amostras.
- FDR: Benjamini–Hochberg < 0,05; |log2FC| > 1 para DEG.
- GSEA: fgsea pre-ranked (t moderada), 335 vias KEGG.
- PPI: STRING v12.0, score ≥ 700.
- Priorização: composite = 0,30·z(|logFC|) + 0,30·z(−log10 FDR) + 0,25·z(betweenness) + 0,15·z(degree).

## 4. Genes prioritários (congelados ANTES da validação)

**Candidatos tumor-intrínsecos (não composicionais):**
FN1, ITGA2, CTSS, HLA-DPA1, CCND1.

**Genes musculares (classificados COMPOSICIONAL — excluídos de candidatura):**
MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM.

## 5. Vias prioritárias (congeladas)

- Painel a priori (10): hsa05216, hsa04919, hsa04010, hsa04151, hsa04150, hsa04115, hsa04210, hsa04110, hsa04310, hsa04064 → **apenas p53 (hsa04115) significativa**.
- Vias não-redundantes (7): hsa04612, hsa03050, hsa04820 (composicional), hsa03030, hsa00510, hsa01232, hsa04115.
- Vias imunes/proteassoma/OXPHOS/replicação/p53 → robustas à composição.

## 6. PPI (congelado)

- 177 genes → 148 nós (GCC) / 727 arestas / 8 comunidades.
- Hubs: FLNC (betweenness 0,27), ITGA2B (0,22), FN1 (0,16).
- FN1 e ITGA2 = hubs estáveis (top-2 por grau em todas as variantes de rede).

## 7. Ranking (congelado)

- Dominado por genes musculares (composicionais) no composite original.
- Após remoção de músculo: **FN1 = #1**, seguido por LAMA2, ITGA2B, ACTN2, ITGA2, CTSS, CCND1.

## 8. FN1 e ITGA2 (congelados como hipótese N5)

- FN1: TPM +4,32 / voom +4,96 / DESeq2 +6,10.
- ITGA2: TPM +2,47 / voom +2,86 / DESeq2 +3,36.
- Interpretação permitida: **hipótese molecular/translacional para investigação futura** (N5/N6).
- **NÃO** são alvos terapêuticos validados, biomarcadores clínicos validados ou alvos de drug delivery demonstrados.

## 9. Critérios de validação (pré-definidos)

Para cada gene candidato congelado, em cada dataset GEO:
- **INDEPENDENTLY REPLICATED:** mesma direção E FDR<0,05.
- **PARTIALLY REPLICATED:** mesma direção E FDR<0,05 em apenas um dataset (ou direção consistente com FDR≥0,05).
- **NOT REPLICATED:** direção inconsistente/ausente.
- **CONTRADICTED:** direção oposta significativa.
- **NOT TESTABLE:** gene não mapeado/ausente na plataforma.

Critérios aplicados **identicamente** a todos os candidatos (sem tratamento especial a FN1/ITGA2).

## 10. SHA256 (snapshot em `results/discovery_freeze/`)

Ver `results/discovery_freeze/SHA256SUMS.txt`.

---

*Após este freeze, nenhum candidato, threshold, probe, peso de ranking ou contraste foi
alterado com base nos resultados GEO.*
