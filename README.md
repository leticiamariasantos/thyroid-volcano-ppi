# thyroid-volcano-ppi

**Análise transcriptômica do carcinoma papilífero da tireoide com painel de 30 vias, convergência molecular e priorização de ITGA2 como candidato translacional**

> **Versão:** 4.1.0 (Fase 2 — reboot) | **Data da execução:** 2026-09-06 | **Tipo de estudo:** Exploratório, gerador de hipóteses

[![R](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![renv](https://img.shields.io/badge/renv-locked-blueviolet)](renv.lock)
[![Docker](https://img.shields.io/badge/Docker-ready-2496ED?logo=docker)](Dockerfile)
[![CI](https://github.com/leticiamariasantos/thyroid-volcano-ppi/actions/workflows/ci.yml/badge.svg)](https://github.com/leticiamariasantos/thyroid-volcano-ppi/actions)

---

## 1. Sobre esta versão

> **Esta é a fase atual e única documentada neste README: uma reconstrução independente do
> pipeline com expansão do painel de vias de 10 para 30 vias.**

A análise anterior (Fase 1 direcionada à via KEGG hsa04919) e a Fase 2 original (painel de 10
vias) foram **removidas da camada analítica** e preservadas apenas como documentação histórica
em `documentation/` e `scripts/legacy/`. Todo o pipeline foi **reconstruído do zero** a partir
dos dados brutos, sem reutilizar resultados estatísticos, rankings, DEGs ou conclusões das
fases anteriores.

**Pergunta científica:** quais programas biológicos emergem do transcriptoma do carcinoma
papilífero da tireoide (PTC) em escala global e, após a expansão **pré-especificada** do espaço
de hipóteses (10 → 30 vias), quais genes emergem como componentes convergentes dos programas
tumorais — e quais podem ser priorizados como candidatos translacionais (não alvos validados)?

---

## 2. Justificativas metodológicas (pré-registradas)

Todas as decisões metodológicas foram tomadas **antes** da interpretação dos resultados desta
execução, para evitar circularidade:

| Decisão | Justificativa |
|---|---|
| Painel de 10 vias preservado | As vias originais (hsa05216, hsa04919, hsa04010, hsa04151, hsa04150, hsa04115, hsa04210, hsa04110, hsa04310, hsa04064) são hipótese **a priori** do projeto |
| 20 vias adicionais selecionadas **a priori** | Ampliar o espaço biológico (adesão/ECM, EMT/TGF-β, Hippo, JAK-STAT, TNF, interferon, imunidade, ferroptose, senescência, proteassoma, OXPHOS, reparo/replicação de DNA, hipóxia, autofagia) **sem** consultar quais vias seriam significativas |
| `source ≡ condition` documentado | TCGA ≡ tumor e GTEx ≡ normal estão perfeitamente confundidos; **nenhuma** correção de batch é alegada como eliminadora |
| Controle de composição muscular | A assinatura de músculo estriado (MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM e 48 marcadores) é **artefato de composição** do tecido normal GTEx, não regulação tumoral |
| ITGA2/FN1/CCND1 avaliados **após** o pipeline | Não se selecionam vias nem genes por esses candidatos; a convergência é verificada a posteriori |
| Nanotecnologia como consequência | A nanotecnologia **não** é premissa da seleção; é consequência translacional hipotética (biologia → candidato → localização → acessibilidade → hipótese de targeting → nanomedicina) |

A seleção das 20 vias adicionais está fundamentada em
`results/phase2/pathways/PANEL_SELECTION_RATIONALE.md` e na tabela
`results/phase2/pathways/PANEL_30_PATHWAYS.tsv`.

---

## 3. Dados

| Camada | Arquivo | Genes | Amostras | Unidade |
|---|---|---|---|---|
| Principal (TPM) | `data/global/TCGA_GTEx_thyroid_tpm.tsv` | 58.581 | 783 (504 TCGA-THCA, 279 GTEx) | log2(TPM+0.001) |
| Contagens (sensibilidade) | `data/global/TCGA_GTEx_thyroid_counts.tsv` | 56.937 | 782 (504 TCGA, 278 GTEx) | raw counts (recount3, GENCODE) |
| Bruto | `data/global/TcgaTargetGtex_rsem_gene_tpm.gz` | 60.498 | 19.131 (global TOIL) | superconjunto |

Proveniência em `data/manifests/` (MANIFEST_TPM.tsv, MANIFEST_counts.tsv). Auditoria de
integridade em `results/phase2/audit/DATA_AUDIT.txt` (sem duplicatas, sem NA/Inf, contagens
inteiras).

---

## 4. Principais resultados

| Etapa | Resultado |
|---|---|
| **DE (limma, principal)** | 26.011 genes testados → 12.200 DEGs (2.485 Up / 9.715 Down; \|log2FC\|>1, FDR<0,05) |
| **Concordância** | Spearman(log2FC) 0,88–0,96 entre limma/voom/DESeq2 |
| **Painel de 30 vias (GSEA)** | 6 significativas e robustas (todas Up): Proteasome, Antigen processing, OXPHOS, DNA replication + p53, Cell cycle |
| **GSEA global** | Ribosome/tradução, Proteasome, resposta imune, Cell cycle, OXPHOS, p53 (Up); Myogenesis (Down) |
| **Composição** | Score muscular tumor −2,91 vs normal +0,03; as 6 vias persistem após remoção dos 48 marcadores |
| **Convergência** | 18 genes no leading edge de ≥2 vias robustas: eixo ciclo celular/p53 (TP53, CCND1, CDK1, CDKN1A, CDKN2A, MDM2, CHEK1, MCM/PCNA) |
| **PPI** | Rede STRING (excluindo músculo): 70 nós / 65 arestas; hub FN1 (grau 9), com CCND1, CDT1, CDH2, SDC1, TK1, PCLAF, CENPM, MET, FGF17 |
| **ITGA2 / FN1 / CCND1** | Todos Up (logFC +2,47 / +4,32 / +2,26), consistentes nos 3 métodos, replicados em GSE33630/GSE60542/GSE224356 |

**Sobre a rede PPI:** na Fase 1 (análise direcionada à via hsa04919, 29 DEGs), **PRKCA** era o
hub. Na análise **genoma-wide** (12.200 DEGs), a rede tumor-relevante é dominada por genes de
**proliferação e adesão/EMT** — **FN1** é o hub principal, e PRKCA **não** é mais hub. A
centralidade é reportada como **propriedade topológica**, não como relevância funcional ou alvo
terapêutico.

**Conclusão oficial:** a expansão 10→30 **não** convergiu para ITGA2/FN1/CCND1 como componentes
dos programas mais enriquecidos (proteassoma/antígeno/OXPHOS/replicação). Eles emergem como
**DEGs individuais robustos e replicados** (CCND1 com core enrichment restrito a vias
**originais** Cell cycle/p53). **ITGA2 permanece candidato translacional para investigação de
direcionamento molecular — nunca alvo validado.** Relatórios completos em
`results/phase2/reports/`.

---

## 5. Estrutura do repositório

```
thyroid-volcano-ppi/
├── .github/workflows/ci.yml     # CI (structure-check, smoke-test, lint)
├── data/
│   ├── global/                  # matriz TPM, contagens e bruto (.gz)
│   ├── external/                # GEO (GSE33630, GSE60542, GSE224356, GSE232237, GSE182416) e Reactome GMT
│   ├── raw/                     # XENA_THCA.tsv (Fase 1)
│   └── manifests/               # manifestos de proveniência
├── docs/                        # protocolo de análise, dicionário de dados, literatura
├── documentation/               # rastreabilidade histórica (Fases 1 e 2 originais)
├── R/                           # funções R legadas (Fase 1)
├── renv/ + renv.lock            # ambiente reprodutível
├── results/phase2/              # resultados da NOVA fase (audit, preprocessing,
│   │                            #   differential_expression, pathways, gsea, ppi,
│   │                            #   validation, sensitivity, figures, reports)
├── scripts/
│   ├── phase2/                  # pipeline da nova fase (00–17)
│   ├── legacy/                  # scripts legados (Fases 1 e 2 originais)
│   ├── 15_validate_outputs.R    # validação técnica (delega para phase2/17_validate.R)
│   ├── download_data.R          # aquisição de dados
│   └── setup_renv.R             # instalação do ambiente
├── tests/                       # testes unitários
├── Dockerfile                   # container reprodutível
├── README.md
└── LICENSE
```

---

## 6. Como executar (cada etapa)

Pré-requisitos: R ≥ 4.1 (testado em R 4.6.1), `renv::restore()` para instalar os pacotes
(versões em `renv.lock`), e conexão de internet para STRING, KEGG, Reactome, MSigDB e cBioPortal.

Os scripts são executados **em ordem** (cada um depende do anterior via `results/phase2/data/*.rds`
e `*.tsv`):

| # | Script | O que faz | Saídas |
|---|---|---|---|
| 00 | `00_config.R` | Parâmetros, thresholds, seeds, caminhos, helpers | — |
| 01 | `01_audit_data.R` | Auditoria de integridade (dimensões, IDs, duplicatas, NA, distribuição) | `audit/DATA_AUDIT.txt`, `audit/sample_metadata.tsv`, `data/tpm_matrix.rds`, `data/counts_matrix.rds` |
| 02 | `02_qc.R` | QC: distribuição, PCA, clustering, outliers, `source≡condition` | `preprocessing/`, figuras QC |
| 03 | `03_de.R` | DE: limma (TPM), limma-voom e DESeq2 (contagens) + comparação de métodos | `differential_expression/` |
| 04 | `04_panel.R` | Define o painel de 30 vias (10 + 20 *a priori*) e obtém gene sets (KEGG/Reactome) | `pathways/PANEL_30_PATHWAYS*.tsv`, `data/panel30_genesets.rds` |
| 05 | `05_gsea_global.R` | GSEA global (KEGG, Reactome, Hallmark) | `gsea/GSEA_GLOBAL_*.tsv` |
| 06 | `06_gsea_panel.R` | GSEA do painel de 30 vias para os 3 métodos | `pathways/PANEL_30_GSEA_*.tsv` |
| 07 | `07_redundancy.R` | Redundância (Jaccard) e clustering das 30 vias | `pathways/PATHWAY_REDUNDANCY.tsv`, `PATHWAY_MODULES.tsv` |
| 08 | `08_composition.R` | Sensibilidade de composição muscular (FULL vs removido) | `sensitivity/` |
| 09 | `09_ppi.R` | Rede PPI (STRING ≥700) sobre DEGs tumor-relevantes (excluindo músculo) | `ppi/` |
| 10 | `10_validation.R` | Validação externa (GSE33630, GSE60542, GSE224356) | `validation/` |
| 11 | `11_singlecell.R` | Localização celular (GSE232237, marker-based) | `validation/singlecell_*` |
| 12 | `12_protein_mutation.R` | RPPA (CCND1) e mutação/CNV (cBioPortal THCA) | `validation/RPPA_summary.tsv`, `mutation_frequency.tsv` |
| 13 | `13_robustness.R` | Matriz de robustez das 30 vias | `pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv` |
| 14 | `14_convergence.R` | Convergência do painel + candidatos ITGA2/FN1/CCND1 | `validation/panel_convergence_genes.tsv`, `convergence_candidates.tsv` |
| 15 | `15_figures.R` | Figuras (volcano, heatmap, GSEA, dotplot, NES, PPI, single-cell) | `figures/` |
| 17 | `17_validate.R` | Validação técnica (0 failures) | — |

**Execução em lote** (na ordem acima):

```bash
cd thyroid-volcano-ppi
Rscript scripts/phase2/00_config.R && Rscript scripts/phase2/01_audit_data.R && \
Rscript scripts/phase2/02_qc.R      && Rscript scripts/phase2/03_de.R          && \
Rscript scripts/phase2/04_panel.R   && Rscript scripts/phase2/05_gsea_global.R && \
Rscript scripts/phase2/06_gsea_panel.R && Rscript scripts/phase2/07_redundancy.R && \
Rscript scripts/phase2/08_composition.R && Rscript scripts/phase2/09_ppi.R      && \
Rscript scripts/phase2/10_validation.R  && Rscript scripts/phase2/11_singlecell.R && \
Rscript scripts/phase2/12_protein_mutation.R && Rscript scripts/phase2/13_robustness.R && \
Rscript scripts/phase2/14_convergence.R && Rscript scripts/phase2/15_figures.R && \
Rscript scripts/phase2/17_validate.R
```

Validação técnica (também acessível por `Rscript scripts/15_validate_outputs.R`).

---

## 7. Pacotes utilizados

| Pacote | Versão | Função | Fonte |
|---|---|---|---|
| R | 4.6.1 | Runtime | CRAN |
| limma | 3.68.0 | DE (lmFit/eBayes), voom | Bioconductor |
| edgeR | 4.10.0 | DGEList, filterByExpr, normLibSizes | Bioconductor |
| DESeq2 | 1.51.7 | DE por contagens (Wald) | Bioconductor |
| fgsea | 1.37.4 | GSEA (fgseaMultilevel) | Bioconductor |
| clusterProfiler | 4.20.0 | Enriquecimento (apoio) | Bioconductor |
| msigdbr | 26.1.0 | Gene sets MSigDB (KEGG_LEGACY, Hallmark) | CRAN |
| KEGGREST | 1.52.0 | Gene sets KEGG (painel) | Bioconductor |
| ReactomePA / reactome.db | 1.56.0 / 1.96.0 | Gene sets Reactome | Bioconductor |
| org.Hs.eg.db | 3.23.1 | Anotação gênica humana | Bioconductor |
| AnnotationDbi | 1.74.0 | Interface de anotação | Bioconductor |
| hgu133plus2.db | 3.13.0 | Mapeamento probe→símbolo (GPL570) | Bioconductor |
| data.table | 1.18.4 | I/O e manipulação de dados | CRAN |
| dplyr / tidyr | 1.2.1 / 1.3.2 | Manipulação/organização | CRAN |
| ggplot2 / ggrepel | 4.0.3 / 0.9.8 | Visualização | CRAN |
| pheatmap / RColorBrewer | 1.0.13 / 1.1-3 | Heatmaps | CRAN |
| igraph | 2.3.0 | Rede PPI, centralidade, walktrap | CRAN |
| httr / jsonlite | 1.4.8 / 2.0.0 | APIs (STRING, cBioPortal) | CRAN |
| here | 1.0.2 | Caminhos portáveis | CRAN |
| readxl | 1.5.0 | Leitura de listas de DEGs (GSE224356) | CRAN |
| reshape2 / MASS | 1.4.5 / — | Manipulação / Mahalanobis | CRAN |

Versões completas em `renv.lock`.

---

## 8. Declaração de uso de Inteligência Artificial

Em conformidade com a Portaria CNPq nº 2.664/2026, declara-se que este projeto utilizou
ferramentas de IA como **suporte técnico e metodológico** (geração/depuração de código R e
Python, revisão de documentação e auditoria de qualidade científica), e não como autora do
conteúdo científico.

- **Nenhuma conclusão científica foi derivada exclusivamente por IA.** Hipóteses, interpretação
  e discussão foram formuladas pelos autores com base nos outputs do pipeline e na literatura.
- **As análises estatísticas foram executadas pelo pipeline em R**, com `set.seed(42)`
  garantindo reprodutibilidade determinística.
- **Nenhum texto científico final foi redigido por IA** — as ferramentas atuaram como
  assistentes de programação e revisão metodológica.
- **As ferramentas não substituem o julgamento científico**; os autores assumem responsabilidade
  integral pela acurácia, validade e adequação das conclusões.

---

## 9. Contribuições (CRediT)

| Autor | Contribuição |
|---|---|
| **Letícia Maria Dias Freitas** | Conceitualização (liderança); Metodologia; Software; Análise formal; Curadoria de dados; Validação; Visualização; Investigação; Redação — rascunho original |
| **Ryan de Paulo Santos** | Conceitualização (suporte); Metodologia; Software; Análise formal; Curadoria; Validação; Visualização; Administração do projeto; Revisão e edição |
| **Thais Faria Coutinho da Silva Pereira** | Supervisão (liderança); Revisão científica (liderança); Validação (suporte) |

## 10. Autores

| Autor | ORCID | Afiliação |
|---|---|---|
| Letícia Maria Dias Freitas | [0009-0009-9930-9588](https://orcid.org/0009-0009-9930-9588) | ETEJBM, Campos dos Goytacazes, RJ |
| Ryan de Paulo Santos | [0009-0005-6770-2001](https://orcid.org/0009-0005-6770-2001) | IFF, Campus Campos Guarus, RJ |
| Thais Faria Coutinho da Silva Pereira | [0009-0005-7091-2480](https://orcid.org/0009-0005-7091-2480) | ETEJBM, Campos dos Goytacazes, RJ |

---

## 11. Limitações

1. **`source ≡ condition`** — TCGA (tumor) e GTEx (normal) estão perfeitamente confundidos;
   nenhuma correção de batch elimina o problema sem remover o sinal biológico.
2. **Composição muscular do normal** — a assinatura de músculo estriado é artefato do GTEx.
3. **PPI esparso e in silico** — interações STRING são inferidas; centralidade é topológica.
4. **RPPA tumor-only** — ITGA2 e FN1 fora do painel; CCND1 com correlação RNA-proteína fraca (0,051).
5. **TERT promotor não captado** — exoma não detecta mutações de promotor.
6. **Single-cell marker-based** — classificação do GSE232237 sem anotação autoral.
7. **Ausência de alteração genômica** que explique a expressão de ITGA2/FN1/CCND1.
8. **Estudo exploratório, não confirmatório** — resultados restritos ao PTC no contexto TCGA/GTEx.

---

## 12. Licença e citação

MIT License — veja [LICENSE](LICENSE).

```bibtex
@software{freitas2026thyroid30,
  title  = {thyroid-volcano-ppi: painel de 30 vias e convergência molecular no
            carcinoma papilífero da tireoide},
  author = {Letícia Maria Dias Freitas and Ryan de Paulo Santos and
            Thais Faria Coutinho da Silva Pereira},
  year   = {2026},
  url    = {https://github.com/leticiamariasantos/thyroid-volcano-ppi},
  note   = {v4.1.0}
}
```
