# RECOVERY STATUS — thyroid-volcano-ppi (Fase 2 global)

**Data da auditoria:** 2026-09-06
**R:** 4.6.1 · **Ambiente:** Windows 11 x64
**Método de auditoria:** evidência de artefatos (timestamps, dimensões, checksums, logs), não apenas nomes de arquivos.

---

## 0. RESUMO DO PONTO DE PARADA

A execução anterior foi interrompida **logo após criar** os scripts exploratórios
`scripts/_recount3_list.R` (08:56) e `scripts/_recount3_explore.R` (09:00), que
investigavam o **recount3** como fonte alternativa de contagens — porque o arquivo
Xena `TcgaTargetGtex_rsem_gene_count.gz` retornou **HTTP 403 Forbidden**.

**Evidência de interrupção:** timestamps de arquivos — a atividade termina em
09:00 (criação de `_recount3_explore.R`); nenhum resultado de contagem existia
antes desta recuperação.

**O que já estava COMPLETO e VALIDADO antes da interrupção:** toda a Fase 2
global sobre `log2(TPM+0.001)` (aquisição → QC → DEG → GSEA → redundância → PPI
→ priorização → sensibilidade → composição), documentada em `RELATORIO_FASE2.md`.

---

## 1. CLASSIFICAÇÃO DE ETAPAS

| ETAPA | STATUS | EVIDÊNCIA |
|---|---|---|
| Aquisição matriz global TPM | COMPLETA E VALIDADA | `data/global/TCGA_GTEx_thyroid_tpm.tsv` (58581×783, 352MB) + `MANIFEST.tsv` + `EXECUCAO_FASE2_LOG.md` |
| QC global (TPM) | COMPLETA E VALIDADA | `02_qc/QC_summary.txt`, `QC_sample_metrics.tsv` (783 linhas), `QC_pca_variance.tsv` |
| DEG genoma-wide (limma, TPM) | COMPLETA E VALIDADA | `04_differential_expression/DEG_full_results.tsv` (20376 genes) + `DEG_summary.tsv` |
| GSEA global (fgsea KEGG) | COMPLETA E VALIDADA | `05_gsea/GSEA_KEGG_all.tsv` (335 vias; 22 padj<0.05, 50 padj<0.25) |
| Painel a priori (10 vias) | COMPLETA E VALIDADA | `06_kegg_panel/KEGG_a_priori_panel.tsv` (só p53 significativa) |
| Redundância de vias | COMPLETA E VALIDADA | `07_pathway_redundancy/*` (22 vias → 6 módulos → 7 representativas) |
| PPI (STRING) | COMPLETA E VALIDADA | `08_ppi/*` (177 genes → 148 nós/727 arestas, 8 comunidades) |
| Priorização multicritério | COMPLETA E VALIDADA | `09_target_prioritization/*` |
| Sensibilidade (DEG/ranking/STRING) | COMPLETA E VALIDADA | `scripts/08_sensitivity.R` + ranking_stability |
| Composição (marcadores) | COMPLETA E VALIDADA | `results/composition/*` (músculo: THCA −2.19 vs Normal +1.52) |
| Sensibilidade à composição | COMPLETA E VALIDADA | `results/sensitivity/*` (GSEA sem músculo) |
| **Count matrix original (Xena)** | **FALHOU (HTTP 403)** | `TcgaTargetGtex_rsem_gene_count.gz` → 403; `_count` e `_norm_count` → 403 |
| **Count matrix alternativa (recount3)** | **EXECUTADA NESTA RECUPERAÇÃO** | `data/global/TCGA_GTEx_thyroid_counts.tsv` (56937×782) + `MANIFEST_counts.tsv` |
| QC de counts | EXECUTADA NESTA RECUPERAÇÃO | `scripts/12_counts_deg.R` (lib size, filtro CPM>1 em ≥10%) |
| limma-voom | EXECUTADA NESTA RECUPERAÇÃO | `results/counts_deg/voom_full_results.tsv` (22118 genes, 7185 DEGs) |
| DESeq2 | EXECUTADA NESTA RECUPERAÇÃO | `results/counts_deg/deseq2_full_results.tsv` (22118 genes, 6976 DEGs) |
| Comparação TPM×voom×DESeq2 | EXECUTADA NESTA RECUPERAÇÃO | `results/counts_deg/method_concordance.tsv` |
| Contribuição da composição | EXECUTADA NESTA RECUPERAÇÃO | `results/composition/composition_contribution.tsv` |
| Deconvolução formal (CIBERSORTx etc.) | NÃO EXECUTADA (limitação registrada) | marker-based scoring como estimativa; ver §6 |
| Validação externa (GEO GSE33630/GSE60542) | EXECUTADA NESTA RECUPERAÇÃO | `10_validation/VALIDATION_REPORT.md` + DEGs GEO |
| Nanomedicina | AVALIADA (não sustentada como proposta específica) | ver `12_nanomedicine/` |
| Testes automatizados | EXECUTADA NESTA RECUPERAÇÃO | `scripts/15_validate_outputs.R` (0 falhas) |
| TRACEABILITY_MATRIX | EXECUTADA NESTA RECUPERAÇÃO | `documentation/TRACEABILITY_MATRIX.tsv` |

---

## 2. ESTADO DOS DADOS

| Matriz | Genes | Amostras | Tipo | Unidade | Fonte | Status |
|---|---|---|---|---|---|---|
| TPM global | 58.581 | 783 (504+279) | gene-level | log2(TPM+0.001) | UCSC Xena TOIL | disponível, validada |
| Count (Xena) | — | — | — | — | Xena (403) | INACESSÍVEL |
| **Count (recount3)** | **56.937** | **782 (504+278)** | gene-level | raw counts (STAR, G026) | recount3 | disponível, validada |

**Proveniência recount3 (registrada):**
- TCGA THCA (`project="THCA"`, `file_source="tcga"`) → 572 amostras; `raw_counts`.
- GTEx THYROID (`project="THYROID"`, `file_source="gtex"`) → 706 amostras; `raw_counts`.
- Anotação GENCODE G026 (63856 features); símbolos colapsados por SOMA (2699 duplicados).
- **Correspondência:** 504/504 TCGA + 278/279 GTEx. Única ausente: `GTEX-SUCS-0226-SM-5CHQG`.
- **Caveat:** contagens recount3 são *gene sums* (soma de cobertura por base) — library sizes
  ~5e8–1.4e10. São proporcionais a read counts e válidas para voom/DESeq2 (métodos que
  normalizam por library size), mas **não são "raw read counts" no sentido RSEM**. Diferente
  também da quantificação TOIL RSEM (GENCODE v23). Tratado como análise de sensibilidade.

---

## 3. RESULTADOS ESTATÍSTICOS (TPM × voom × DESeq2)

| Método | Genes testados | DEGs | up | down | Observações |
|---|---|---|---|---|---|
| limma (TPM) | 20.376 | 8.161 | 1.530 | 6.631 | `trend=TRUE`, contraste THCA−Normal |
| limma-voom | 22.118 | 7.185 | 2.710 | 4.475 | TMM + voom + eBayes |
| DESeq2 | 22.118 | 6.976 | 2.888 | 4.088 | raw counts, design ~condition |

**Concordância entre métodos (genes comuns = 18.459):**

| Comparação | logFC Pearson | logFC Spearman | direção | Jaccard DEG | ranking Spearman |
|---|---|---|---|---|---|
| TPM × voom | 0.895 | 0.922 | 0.813 | 0.602 | 0.891 |
| TPM × DESeq2 | 0.871 | 0.906 | 0.808 | 0.557 | 0.907 |
| voom × DESeq2 | 0.954 | 0.978 | 0.964 | 0.857 | 0.983 |

**Classificação de estabilidade (genes comuns):** ROBUSTO 4293 · PARCIALMENTE ROBUSTO 1218 ·
MÉTODO-DEPENDENTE 2471 · NÃO ROBUSTO 10477.

---

## 4. COMPOSIÇÃO TECIDUAL (FASE 9/10)

| Categoria | marcadores testados | DEGs | up | down |
|---|---|---|---|---|
| muscle_skeletal | 31 | 24 | 0 | **24** |
| thyroid_epithelial | 13 | 7 | 1 | 6 |
| fibroblast_ecm | 16 | 9 | 1 | 8 |
| endothelial | 9 | 9 | 0 | 9 |
| immune | 18 | 5 | 3 | 2 |

- **Score muscular:** THCA = −2,19 vs Normal = **+1,52** (log2 TPM; p=8,2e-114). → sinal muscular
  **concentrado nos GTEx normais**, compatível com músculo esquelético adjacente (strap muscle).
- **PC1 correlaciona −0,845** com score muscular (e −0,657 epitélio tireoidiano, −0,681 endotelial).
- Os 6 genes musculares levantados (MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM) são **ROBUSTOS** nos 3
  métodos (todos down em THCA), mas sua magnitude é **composicional**, não biologia tumoral
  intrínseca.
- **Conclusão:** o sinal muscular é quase certamente artefato de composição tecidual (diferença
  TCGA-tumor vs GTEx-normal), NÃO regulação oncogênica.

---

## 5. GENES-CHAVE (FN1/ITGA2 e musculares)

| Gene | logFC TPM | logFC voom | logFC DESeq2 | Robustez |
|---|---|---|---|---|
| FN1 | +4,32 | +4,96 | +6,10 | ROBUSTO |
| ITGA2 | +2,47 | +2,86 | +3,36 | ROBUSTO |
| MYH7 | −7,01 | −7,04 | −5,91 | ROBUSTO |
| MYL2 | −8,69 | −9,43 | −6,29 | ROBUSTO |
| ACTA1 | −7,17 | −5,41 | −5,53 | ROBUSTO |
| TNNT3 | −4,77 | −4,72 | −3,74 | ROBUSTO |
| CKM | −6,92 | −6,65 | −5,81 | ROBUSTO |

FN1/ITGA2 permanecem superexpressos no tumor de forma robusta; musculares permanecem down de
forma robusta mas são composicionais.

---

## 6. LIMITAÇÕES QUE PERMANECEM

1. **Count original (Xena) = 403** — a contagem usada é recount3 (STAR G026), não TOIL RSEM.
2. **1 amostra GTEx ausente** (`GTEX-SUCS-0226-SM-5CHQG`) na matriz recount3 (782 vs 783).
3. **Batch confundido com condição** (TCGA=tumor, GTEx=normal) — limitação estrutural.
4. **Deconvolução formal não executada** (sem CIBERSORTx/assinatura de referência); usado
   marker-based scoring (estimativa/modelo).
5. **Sem validação externa independente** (GEO GSE33630/GSE60542 pendentes).
6. Unidades TPM (log2 TPM+0.001) vs count (gene sums) — magnitudes não diretamente comparáveis;
   comparáveis em direção/ranking.

---

## 7. PRÓXIMAS AÇÕES (ordem)

1. ✅ Contagem recount3 + voom/DESeq2 + comparação (feito nesta recuperação).
2. ✅ Validação externa (GSE33630 primário, GSE60542 secundário) — FN1/ITGA2/CTSS/HLA-DPA1/CCND1 replicados; musculares NÃO replicados (confirma composição).
3. ⬜ Deconvolução formal (CIBERSORTx / assinatura de referência) — opcional, justificado.
4. ⬜ Atualizar GSEA/redundância/PPI/priorização incorporando a sensibilidade à composição (já
   parcialmente feito em `results/sensitivity/`).
