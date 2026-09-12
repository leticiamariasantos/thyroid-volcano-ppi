# Análise transcriptômica robusta do carcinoma papilífero da tireoide com controle de batch e composição

> **Versão:** 5.0.0 (upgrade metodológico 2026) | **Protocolo congelado:** 2026-09-08 | **Tipo de estudo:** Exploratório, gerador de hipóteses

[![R](https://img.shields.io/badge/R-4.6.1-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![renv](https://img.shields.io/badge/renv-locked-blueviolet)](renv.lock)
[![Docker](https://img.shields.io/badge/Docker-ready-2496ED?logo=docker)](Dockerfile)
[![CI](https://github.com/leticiamariasantos/thyroid-volcano-ppi/actions/workflows/ci.yml/badge.svg)](https://github.com/leticiamariasantos/thyroid-volcano-ppi/actions)

---

## 1. Sobre esta versão

> **Esta é a fase atual e única documentada neste README: uma reconstrução independente do
> pipeline com um painel de 30 vias.**

A análise anterior (Fase 1 direcionada à via KEGG hsa04919) e a Fase 2 original (painel de 10
vias) foram **removidas da camada analítica** e preservadas apenas como documentação histórica
em `documentation/` e `scripts/legacy/`. Todo o pipeline foi **reconstruído do zero** a partir
dos dados brutos, sem reutilizar resultados estatísticos, rankings, DEGs ou conclusões das
fases anteriores.

**Questão norteadora:** quais programas biológicos emergem do transcriptoma do carcinoma
papilífero da tireoide (PTC) em escala global e, com um painel **pré-especificado** de 30 vias,
quais genes emergem como componentes convergentes dos programas
tumorais, e quais podem ser priorizados como candidatos translacionais (não alvos validados)?

---

## 2. Justificativas metodológicas

Todas as decisões metodológicas foram tomadas antes da interpretação dos resultados desta
execução, para evitar circularidade:

| Decisão | Justificativa |
|---|---|
| Painel de 30 vias definido **a priori** | As 30 vias foram pré-especificadas antes da interpretação dos resultados (sinalização, adesão/ECM, EMT/TGF-β, Hippo, JAK-STAT, TNF, interferon, imunidade, ferroptose, senescência, proteassoma, OXPHOS, reparo/replicação de DNA, hipóxia e autofagia), sem consultar quais seriam significativas |
| Contraste TCGA-matched primário | 505 tumores e 59 normais adjacentes recount3; os 59 participantes com ambos os tecidos também formam uma sensibilidade pareada |
| Batch identificável | Normais TCGA quebram o confundimento perfeito; TCGA×GTEx bruto continua sendo apenas sensibilidade, e correção ainda pode sobrecorrigir |
| Deconvolução e composição | xCell + EPIC formais, com consenso multiassinatura auditável para músculo, fibroblasto, endotélio, imune e epitélio |
| Cinco métodos de DE | limma, voom, voom com pesos de qualidade, edgeR-QL e DESeq2 em versões raw, batch-corrected, composition-adjusted e TCGA-matched |
| ITGA2/FN1/CCND1 avaliados **após** o pipeline | Não se selecionam vias nem genes por esses candidatos; a convergência é verificada a posteriori |
| Nanotecnologia como consequência | A nanotecnologia **não** é premissa da seleção; é consequência translacional hipotética (biologia → candidato → localização → acessibilidade → hipótese de targeting → nanomedicina) |

A seleção das 30 vias está fundamentada em
`results/phase2/pathways/PANEL_SELECTION_RATIONALE.md` e na tabela
`results/phase2/pathways/PANEL_30_PATHWAYS.tsv`.

---

## 3. Dados

| Camada | Arquivo | Genes | Amostras | Unidade |
|---|---|---|---|---|
| Principal (TPM) | `data/global/TCGA_GTEx_thyroid_tpm.tsv` | 58.581 | 783 (504 TCGA-THCA, 279 GTEx) | log2(TPM+0.001) |
| Contagens (sensibilidade) | `data/global/TCGA_GTEx_thyroid_counts.tsv` | 56.937 | 782 (504 TCGA, 278 GTEx) | raw counts (recount3, GENCODE) |
| Bruto | `data/global/TcgaTargetGtex_rsem_gene_tpm.gz` | 60.498 | 19.131 (global TOIL) | superconjunto |
| Primário TCGA-matched | cache recount3 + `upgrade_2026/matrices/` | GENCODE G026 filtrado | 505 tumores + 59 normais adjacentes | raw counts |

Proveniência em `data/manifests/` (MANIFEST_TPM.tsv, MANIFEST_counts.tsv). Auditoria de
integridade em `results/phase2/audit/DATA_AUDIT.txt` (sem duplicatas, sem NA/Inf, contagens
inteiras).

---

## 4. Principais resultados

| Etapa | Resultado |
|---|---|
| **DE (limma, principal)** | 26.011 genes testados → 12.200 DEGs (2.485 Up / 9.715 Down; \|log2FC\|>1, FDR<0,05) |
| **Concordância** | Spearman(log2FC) 0,88 a 0,96 entre limma/voom/DESeq2 |
| **Painel de 30 vias (GSEA)** | 6 significativas e robustas (todas Up): Proteasome, Antigen processing, OXPHOS, DNA replication + p53, Cell cycle |
| **GSEA global** | Ribosome/tradução, Proteasome, resposta imune, Cell cycle, OXPHOS, p53 (Up); Myogenesis (Down) |
| **Composição** | Score muscular tumor −2,91 vs normal +0,03; as 6 vias persistem após remoção dos 48 marcadores |
| **Convergência** | 18 genes no leading edge de ≥2 vias robustas: eixo ciclo celular/p53 (TP53, CCND1, CDK1, CDKN1A, CDKN2A, MDM2, CHEK1, MCM/PCNA) |
| **PPI** | Rede STRING (excluindo músculo): 70 nós / 65 arestas; hub FN1 (grau 9), com CCND1, CDT1, CDH2, SDC1, TK1, PCLAF, CENPM, MET, FGF17 |
| **ITGA2 / FN1 / CCND1** | Todos Up (logFC +2,47 / +4,32 / +2,26, respectivamente), consistentes nos 3 métodos, replicados em GSE33630/GSE60542/GSE224356 |

**Sobre a rede PPI:** Na análise **genoma-wide** (12.200 DEGs), a rede tumor-relevante é dominada por genes de
**proliferação e adesão/EMT**. **FN1** é o hub principal, e PRKCA **não** é mais hub. A
centralidade é reportada como **propriedade topológica**, não como relevância funcional ou alvo
terapêutico.

**Conclusão da versão 5.0.0:** o painel de 30 vias não convergiu para ITGA2/FN1/CCND1 como componentes
dos programas mais enriquecidos (proteassoma/antígeno/OXPHOS/replicação). Eles emergem como
**DEGs individuais robustos e replicados**. Na execução 5.0.0, **ITGA2 passou em 5 dos 6
critérios pré-especificados** (falhou apenas em `surface_and_protein`, porque a localização
subcelular de membrana no HPA não estava anotada para ITGA2) e foi classificado como
**“DEG robusto de interesse”** — nunca “alvo validado”. Relatórios completos em
`results/phase2/reports/`.

### Resultados 5.0.0 (upgrade metodológico 2026)

**Equivalência formal (voom quality-weights).** A reimplementação do passo `arrayWeights`
(genebygene) reproduziu a referência oficial do limma 3.68.0: diferenças máximas de ~1e-13
em logFC, estatística t, p-valor, FDR e pesos de qualidade, com ranking idêntico
(`pipeline/DE_acceleration_equivalence.tsv`, `passed=TRUE`).

**Expressão diferencial (5 variantes × 5 métodos).** Genes testados e DEGs do método
`voom_qw` por variante (limiar |log2FC|≥1, FDR<0,05):

| Variante | Genes testados | DEGs (voom_qw) | Up | Down |
|---|---:|---:|---:|---:|
| tcga_matched | 22.787 | 4.015 | 1.961 | 2.054 |
| raw | 22.103 | 3.677 | 1.688 | 1.989 |
| batch_corrected | 21.446 | 3.500 | 2.348 | 1.152 |
| composition_adjusted | 22.103 | 2.592 | 1.132 | 1.460 |
| tcga_paired | 20.984 | 3.080 | 1.572 | 1.508 |

**GSEA e robustez.** 14 vias consistentes após 100 bootstraps (fração de direção e de
significância ≥0,80/0,70): `hsa03050` proteassoma, `hsa04110` ciclo celular, `hsa04115`
p53, `hsa03030` replicação de DNA, `hsa04510` adesão focal, `hsa04512` interação
ECM-receptor, `hsa04668` TNF, `hsa05216` câncer de tireoide, `hsa04210` apoptose,
`hsa04218` senescência celular, `hsa04010` MAPK, `hsa04151` PI3K-Akt,
`R-HSA-1474244` organização da matriz extracelular e `R-HSA-913531` sinalização de
interferon.

**Candidatos (tcga_matched, voom_qw).** ITGA2 log2FC 2,76 (FDR 9,3e-40); FN1 5,35
(FDR 1,3e-108); CCND1 1,83 (FDR 3,0e-109); CDH2 2,70 (FDR 8,1e-16), todos superexpressos.

**Meta-análise e gate ITGA2.** Meta-análise REML: log2FC 2,42 (IC 95% 1,94–2,91),
FDR 5,6e-22, I² 66,7. O gate pré-especificado de acessibilidade de superfície resultou
em 5/6 critérios (`surface_and_protein` falhou por ausência de anotação de membrana
plasmática no HPA) → classificação **“DEG robusto de interesse”**.

**PPI (3 universos, STRING experimental + coexpressão |rho|≥0,30).**
`genome_wide_top500` (34 nós), `panel30` (1.333 nós), `robust_leading_edge` (1.109 nós);
centralidade dominada por componentes da fosforilação oxidativa e da tradução ribossomal,
além de hubs de proliferação e adesão.

**Deconvolução e composição.** xCell + EPIC + consenso multiassinatura; 51 amostras EPIC
não-convergentes divulgadas por amostra; R² parcial dos 100 top DEGs em
`deconvolution/top_deg_celltype_variance.tsv`.

**Validação técnica.** `31_validate_upgrade.R` → 42/42 checks `TRUE` (`0 failures`);
manifesto de checksums em `validation/OUTPUT_FILE_CHECKSUMS.tsv` (125 artefatos).

### Antes vs depois

| Métrica | 4.1.0 | 5.0.0 |
|---|---:|---:|
| Contraste primário | TCGA tumor × GTEx normal (confundido) | TCGA tumor × normal adjacente; sensibilidade pareada |
| Métodos de DE | 3 | 5 |
| Versões analíticas | 1 principal + sensibilidades musculares | raw, batch-corrected, composition-adjusted, TCGA-matched e TCGA-paired |
| DEGs / vias robustas / ranking ITGA2 | 12.200 / 6 / resultado 4.1.0 | 5 variantes × 5 métodos; 14 vias robustas; ITGA2 = “DEG robusto de interesse” (5/6 critérios) |
| Composição | score/remoção muscular | xCell, EPIC e consenso multiassinatura + R² parcial |
| PPI | uma rede | três redes; evidência experimental + coexpressão; Walktrap/Louvain/Leiden |

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
│   ├── phase2/                  # pipeline base (00 a 16) + upgrade (20 a 31)
│   ├── legacy/                  # scripts legados (Fases 1 e 2 originais)
│   ├── 15_validate_outputs.R    # validação técnica (delega para phase2/16_validate.R)
│   ├── download_data.R          # aquisição de dados
│   └── setup_renv.R             # instalação do ambiente
├── tests/                       # testes unitários
├── Dockerfile                   # container reprodutível
├── README.md
└── LICENSE
```

O protocolo a priori completo está em `docs/METHODOLOGICAL_UPGRADE_2026.md`. O painel de
30 vias não foi alterado. `documentation/` e `scripts/legacy/` permanecem históricos.

---

## 6. Como executar (passo a passo completo)

Esta seção descreve, em detalhe, como reproduzir integralmente a análise. Foi escrita para ser
seguida tanto por uma pessoa quanto por um agente de IA, com os comandos exatos e o que verificar
em cada etapa.

### 6.1. Pré-requisitos

| Item | Requisito |
|---|---|
| Sistema operacional | Windows, macOS ou Linux |
| R | 4.6.1 (versão registrada e testada em 2026-09-08) |
| Gerenciador de pacotes | renv (o projeto usa `renv.lock` com versões exatas) |
| Internet | Necessária para STRING, KEGG, Reactome (MSigDB), cBioPortal e GEO |
| Espaço em disco | ~3 GB (matrizes grandes versionadas via Git LFS) |
| Git LFS | Necessário para clonar as matrizes grandes |

### 6.2. Etapa 0: obter o código e os dados

```bash
git lfs install
git clone https://github.com/leticiamariasantos/thyroid-volcano-ppi.git
cd thyroid-volcano-ppi
git lfs pull
```

Os dados brutos ficam em `data/global/` (TPM, contagens, arquivo `.gz`), `data/external/`
(GEO, single-cell, Reactome GMT) e `data/raw/`. Os manifestos de proveniência estão em
`data/manifests/`.

### 6.3. Etapa 1: restaurar o ambiente R

```bash
Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore()'
```

Confirme que os pacotes foram instalados nas versões registradas em `renv.lock`.

### 6.4. Etapa 2: executar o pipeline (em ordem)

Os scripts estão em `scripts/phase2/` e devem rodar em ordem numérica (00 a 16), pois cada um
consome artefatos do anterior (`results/phase2/data/*.rds` e `*.tsv`).

**Forma recomendada (orquestrador):**

```bash
Rscript scripts/run_phase2.R
```

**Forma manual (um a um, para depuração):**

| # | Comando | O que faz e o que verificar |
|---|---|---|
| 00 | `Rscript scripts/phase2/00_config.R` | Carrega parâmetros, thresholds, seed e caminhos; verifique se a raiz foi detectada |
| 01 | `Rscript scripts/phase2/01_audit_data.R` | Auditoria de integridade; deve terminar com 0 falhas e gerar `audit/DATA_AUDIT.txt` |
| 02 | `Rscript scripts/phase2/02_qc.R` | QC, PCA, clustering e outliers; gera `preprocessing/QC_summary.txt` |
| 03 | `Rscript scripts/phase2/03_de.R` | DE com limma, limma-voom e DESeq2; gera `differential_expression/` (12.200 DEGs no limma) |
| 04 | `Rscript scripts/phase2/04_panel.R` | Define o painel de 30 vias e obtém gene sets; gera `pathways/PANEL_30_PATHWAYS.tsv` |
| 05 | `Rscript scripts/phase2/05_gsea_global.R` | GSEA global (KEGG, Reactome, Hallmark); gera `gsea/GSEA_GLOBAL_*.tsv` |
| 06 | `Rscript scripts/phase2/06_gsea_panel.R` | GSEA das 30 vias para os 3 métodos; gera `pathways/PANEL_30_GSEA_*.tsv` |
| 07 | `Rscript scripts/phase2/07_redundancy.R` | Redundância (Jaccard) e clustering; gera `pathways/PATHWAY_REDUNDANCY.tsv` |
| 08 | `Rscript scripts/phase2/08_composition.R` | Sensibilidade de composição muscular; gera `sensitivity/` |
| 09 | `Rscript scripts/phase2/09_ppi.R` | Rede PPI (STRING, escore >=700); gera `ppi/` |
| 10 | `Rscript scripts/phase2/10_validation.R` | Validação externa (GSE33630, GSE60542, GSE224356); gera `validation/` |
| 11 | `Rscript scripts/phase2/11_singlecell.R` | Localização single-cell (GSE232237); gera `validation/singlecell_*` |
| 12 | `Rscript scripts/phase2/12_protein_mutation.R` | RPPA e mutação/CNV (cBioPortal THCA); gera `validation/RPPA_summary.tsv` |
| 13 | `Rscript scripts/phase2/13_robustness.R` | Matriz de robustez; gera `pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv` |
| 14 | `Rscript scripts/phase2/14_convergence.R` | Convergência do painel e candidatos; gera `validation/` |
| 15 | `Rscript scripts/phase2/15_figures.R` | Figuras; gera `figures/` |
| 16 | `Rscript scripts/phase2/16_validate.R` | Validação técnica; deve terminar com 0 failures |

**Execução em lote alternativa (linha única, sem o orquestrador):**

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
Rscript scripts/phase2/16_validate.R
```

**Upgrade metodológico 2026 (20–31):**

```bash
Rscript scripts/run_phase2.R --from=20
```

| # | Função principal |
|---:|---|
| 20–22 | Congela configuração; integra normais TCGA; filtra CPM ≥1 em ≥25% do menor grupo e corrige batch |
| 23 | Executa xCell, EPIC, consenso multiassinatura e matriz ajustada por composição |
| 24–26 | Executa cinco métodos de DE, GSEA/GSVA/leading edge e R² parcial por compartimento |
| 27 | Meta-análise REML, OS/PFS/DSS estratificados e gate quantitativo de superfície |
| 28–29 | Três redes PPI filtradas e 100 bootstraps de GSEA + remoção de outliers |
| 30–31 | Renderiza relatório automatizado e executa validação técnica integral |

O orquestrador grava assinaturas de script + inputs, checksums dos outputs e tempos em
`results/phase2/upgrade_2026/pipeline/`. Use `--force` para uma rodada limpa e
`--adopt-existing` apenas para registrar outputs históricos já verificados.

### 6.5. Etapa 3: validação técnica

```bash
Rscript scripts/phase2/16_validate.R
# ou
Rscript scripts/15_validate_outputs.R
```

Saída esperada: `VALIDAÇÃO OK (0 failures)`.

### 6.6. Etapa 4: testes unitários (opcional)

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### 6.7. Dependências entre etapas (para depuração)

- `00_config.R` é carregado por todos os scripts via `source(...)`, portanto roda primeiro.
- `01_audit_data.R` gera `data/tpm_matrix.rds` e `data/counts_matrix.rds`, consumidos por `02_qc.R` e `03_de.R`.
- `03_de.R` gera rankings e resultados de DE, consumidos a partir de `04_panel.R`.
- `04_panel.R` gera `data/panel30_genesets.rds`, consumido por `06_gsea_panel.R`, `07_redundancy.R`, `08_composition.R` e `14_convergence.R`.
- `05_gsea_global.R` persiste os gene sets globais (RDS), consumidos por `13_robustness.R`.
- `09_ppi.R` consulta a API do STRING; `10_validation.R` e `12_protein_mutation.R` dependem de internet e são independentes entre si.

### 6.8. Tempo estimado

| Etapa | Tempo aproximado |
|---|---|
| 00 a 02 | 2 a 5 minutos |
| 03 (inclui DESeq2) | 10 a 15 minutos |
| 04 a 08 | 5 a 10 minutos |
| 09 a 12 | 3 a 8 minutos (depende da rede) |
| 13 a 16 | 2 a 5 minutos |
| 20 a 21 | 3 a 10 minutos na primeira execução; cache nas seguintes |
| 22 (ComBat-seq completo) | potencialmente 1–2 horas; retomado por checksum |
| 23 a 31 | dependente de deconvolução, bootstrap e APIs; cache local versionado |

### 6.9. Resolução de problemas comuns

| Sintoma | Causa provável e solução |
|---|---|
| Erro `pandoc not found` | Instale o Pandoc ou adicione-o ao PATH |
| Erro de pacote ausente | Rode `renv::restore()` novamente |
| Falha em STRING/cBioPortal/HPA/Surfaceome | O pipeline reutiliza cache válido; sem rede e sem cache válido, falha explicitamente |
| Arquivos grandes não vieram no clone | Rode `git lfs pull` |

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
| reshape2 / MASS | 1.4.5 / n/a | Manipulação / Mahalanobis | CRAN |
| sva | 3.60.0 | ComBat-seq e variáveis substitutas | Bioconductor |
| xCell / EPIC | 1.1.0 / 1.1.7 | Deconvolução formal | GitHub (versões fixadas) |
| GSVA | 2.6.6 | Scores de vias por amostra | Bioconductor |
| metafor / survival | 5.0-1 / 3.8-6 | Meta-análise REML e OS/PFS/DSS | CRAN |
| digest | 0.6.39 | Checksums de inputs e caches | CRAN |

Versões completas em `renv.lock`.

---

## 8. Declaração de uso de Inteligência Artificial

Em conformidade com a Portaria CNPq nº 2.664/2026, este projeto declara o uso de ferramentas de
inteligência artificial (IA) exclusivamente como suporte técnico, e nunca
como autora do conteúdo, das hipóteses ou das conclusões.

### 8.1. Natureza do uso de IA

| Atividade | Papel da IA | Papel dos pesquisadores |
|---|---|---|
| Geração, depuração e otimização de código R e Python | Assistência na implementação do pipeline | Especificação dos métodos, revisão e validação de cada script |
| Auditoria de qualidade científica e linguística | Busca automática de padrões (overclaiming, inconsistências) | Decisão sobre o que corrigir e redação das correções |
| Revisão de documentação | Sugestões de estrutura e clareza | Redação e aprovação final de todo o texto |
| Interpretação dos resultados | Apoio na checagem de consistência e identificação de ressalvas | Formulação de hipóteses, leitura dos outputs, decisão interpretativa e discussão |

### 8.2. Garantias de integridade

- **Nenhuma conclusão científica foi derivada exclusivamente por IA.** As hipóteses biológicas, a
  interpretação dos resultados e a discussão sobre relevância translacional foram formuladas
  pelos autores, com base nos outputs reproduzíveis do pipeline e na literatura.
- **As análises estatísticas foram executadas integralmente pelo pipeline em R**, com
  `set.seed(42)` e versões fixadas em `renv.lock`, garantindo reprodutibilidade determinística.
- **A IA auxiliou na redação técnica e revisão metodológica.** Todo texto e toda interpretação
  permanecem sujeitos à revisão, aprovação e responsabilidade dos autores; assistência de IA
  não constitui autoria científica.
- **As ferramentas de IA não substituem o julgamento científico.** Os autores assumem
  responsabilidade integral pela acurácia dos dados, pela validade das análises e pela adequação
  das conclusões apresentadas.

### 8.3. Rastreabilidade

O registro das tarefas assistidas por IA, a trilha de auditoria do pipeline e a matriz de
rastreabilidade estão documentados em `documentation/TRACEABILITY_MATRIX.tsv`,
`documentation/TRANSPARENCY_LOG.md` e `results/phase2/reports/PHASE2_REPRODUCIBILITY.md`.

---

## 9. CRediT

| Autor | ORCID | Afiliação | Contribuição (CRediT) |
|---|---|---|---|
| Letícia Maria Dias Freitas | [0009-0009-9930-9588](https://orcid.org/0009-0009-9930-9588) | ETEJBM, Campos dos Goytacazes, RJ | Conceitualização (liderança); Metodologia; Software; Análise formal; Curadoria de dados; Validação; Visualização; Investigação; Redação (rascunho original) |
| Ryan de Paulo Santos | [0009-0005-6770-2001](https://orcid.org/0009-0005-6770-2001) | IFF, Campus Campos Guarus, RJ | Conceitualização (suporte); Metodologia; Software; Análise formal; Curadoria; Validação; Visualização; Administração do projeto; Revisão e edição |
| Thais Faria Coutinho da Silva Pereira | [0009-0005-7091-2480](https://orcid.org/0009-0005-7091-2480) | ETEJBM, Campos dos Goytacazes, RJ | Supervisão (liderança); Revisão científica (liderança); Validação (suporte) |

---

## 10. Limitações da análise 4.1.0 (históricas)

1. **`source ≡ condition`**: TCGA (tumor) e GTEx (normal) estão perfeitamente confundidos;
   nenhuma correção de batch elimina o problema sem remover o sinal biológico.
2. **Composição muscular do normal**: a assinatura de músculo estriado é artefato do GTEx.
3. **PPI esparso e in silico**: interações STRING são inferidas; centralidade é topológica.
4. **RPPA tumor-only**: ITGA2 e FN1 fora do painel; CCND1 com correlação RNA-proteína fraca (0,051).
5. **TERT promotor não captado**: exoma não detecta mutações de promotor.
6. **Single-cell marker-based**: classificação do GSE232237 sem anotação autoral.
7. **Ausência de alteração genômica** que explique a expressão de ITGA2/FN1/CCND1.
8. **Estudo exploratório**: resultados restritos ao PTC no contexto TCGA/GTEx.

### Limitações residuais após as melhorias

1. Normais adjacentes reduzem o confundimento de fonte, mas podem conter *field effects* e
   são muito menos numerosos que os tumores.
2. O desenho ampliado torna `source` estimável, porém ComBat-seq, SVA e ajustes contínuos
   ainda podem remover sinal biológico verdadeiro; nenhuma matriz corrigida é “ground truth”.
3. Deconvolução depende de assinaturas de referência e não demonstra origem celular,
   pureza, estado funcional ou causalidade.
4. A composição residual, diferenças pré-analíticas e covariáveis clínicas incompletas podem
   permanecer após o ajuste.
5. Meta-análise combina plataformas heterogêneas; alguns erros-padrão são aproximados a
   partir de logFC e p-value.
6. THCA tem poucos eventos de OS/PFS/DSS, sobretudo após estratificação por BRAF/RAS e
   histologia, limitando precisão e poder.
7. Human Protein Atlas/Surfaceome indicam localização/acessibilidade potencial, não
   internalização, especificidade tumoral, segurança ou janela terapêutica.
8. Mesmo após batch correction e deconvolução, o estudo permanece **exploratório e gerador
   de hipóteses**. Nenhum resultado valida alvo terapêutico ou estratégia nanomédica.

---

## 11. Licença

MIT License. Veja [LICENSE](LICENSE).

```bibtex
@software{freitas2026thyroid30,
  title  = {thyroid-volcano-ppi: painel de 30 vias e convergência molecular no
            carcinoma papilífero da tireoide},
  author = {Letícia Maria Dias Freitas and Ryan de Paulo Santos and
            Thais Faria Coutinho da Silva Pereira},
  year   = {2026},
  url    = {https://github.com/leticiamariasantos/thyroid-volcano-ppi},
  note   = {v5.0.0}
}
```
