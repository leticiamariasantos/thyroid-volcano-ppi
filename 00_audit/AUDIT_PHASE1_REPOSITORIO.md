# FASE 1 — AUDITORIA INTEGRAL DO REPOSITÓRIO
## thyroid-volcano-ppi (v3.1.0)

**Auditor:** revisor científico independente (papel: metodologista/cético)
**Data da auditoria:** 2026-09-06
**Origem auditada:** `https://github.com/leticiamariasantos/thyroid-volcano-ppi.git` (branch `main`, HEAD `f4c0e1f`)
**Método:** inspeção estática integral (código, dados, resultados, logs, docs, lockfiles, git). **Não** houve re-execução do pipeline neste ambiente (R ≥ 4.6 não disponível localmente — ver §8).

> **Regra aplicada nesta auditoria:** separar rigorosamente (i) *código existe* de *análise executada*; (ii) *resultado produzido* de *resultado esperado*; (iii) *afirmação documentada* de *afirmação verificada*.

---

## 1. INVENTÁRIO DO REPOSITÓRIO

| Item | Presente | Observações |
|---|---|---|
| `README.md` | Sim | v3.1.0, extenso, com resultados auto-descritos |
| `run_pipeline.R` | Sim | Orquestrador mestre |
| `R/00_setup.R` | Sim | Parâmetros, pacotes, seed 42 |
| `R/01_functions.R` | Sim | Funções (KEGG, STRING, centralidade, export) |
| `R/02_import.R` | Sim | Importação/validação |
| `R/03_deg.R` | Sim | Expressão diferencial (limma) |
| `R/03b_pca.R` | Sim | PCA/MDS (QC) — chamado pelo pipeline |
| `R/03c_heatmap.R` | Sim | **NÃO chamado pelo pipeline** |
| `R/03d_qc_outliers.R` | Sim | **Chamado, mas FALHA em runtime** |
| `R/03e_umap_qc.R` | Sim | UMAP (QC) — chamado pelo pipeline |
| `R/04_volcano.R` | Sim | Volcano plot — chamado |
| `R/05_ppi.R` | Sim | Rede PPI — chamado |
| `R/06_supplementary.R` | Sim | **NÃO chamado pelo pipeline** (tabelas S1–S4 geradas à parte) |
| `scripts/download_data.R` | Sim | Baixa o arquivo **completo** (`TcgaTargetGtex_rsem_gene_tpm.gz`) |
| `scripts/setup_renv.R` | Sim | Inicialização renv |
| `scripts/audit_docx_*.py`, `lock_docx_final.py`, `linguistic_audit.py` | Sim | Auditorias de manuscrito .docx (fora do pipeline R) |
| `data/raw/XENA_THCA.tsv` | Sim | **640 KB — matriz RESTRITA a 121 genes** |
| `data/README.md` | Sim | Documentação do diretório de dados |
| `docs/analysis_protocol.md` | Sim | Protocolo |
| `docs/data_dictionary.md` | Sim | Dicionário |
| `docs/figure_specs.md` | Sim | Especificação de figuras |
| `results/figures/`, `results/tables/`, `results/network/` | Sim | Saídas commitadas |
| `logs/` (pipeline + session_info) | Sim | 22 logs de pipeline + 22 logs de mensagens + session_info |
| `results/CHECKSUMS.md` | Sim | MD5 dos outputs |
| `renv.lock`, `renv/` | Sim | Lockfile |
| `Dockerfile` | Sim | **Usa `rocker/r-ver:4.4.0`** |
| `tests/testthat/` | Sim | Testes com dados simulados |
| `CITATION.cff`, `LICENSE` | Sim | MIT |
| `AUDIT_REPORT.md` | **AUSENTE** | Referenciado em `R/02_import.R` (linhas 7 e 78) e no `run_pipeline.R` — **referência pendurada** |
| `data/processed/` | **AUSENTE** | Documentado como "dados intermediários", não versionado |

---

## 2. ACHADO CENTRAL (MUDANÇA DE PARADIGMA PARA O PROJETO)

**O dataset commitado NÃO é transcriptoma global. É uma matriz de UMA ÚNICA VIA.**

- `data/raw/XENA_THCA.tsv` tem **783 linhas (amostras)** e **127 colunas**:
  - 6 colunas de metadados: `sample`, `samples`, `_study`, `TCGA_GTEX_main_category`, `_sample_type`, `_primary_site`;
  - **121 colunas de genes** (todas únicas).
- As 121 colunas correspondem, uma a uma, aos genes da via **KEGG hsa04919 (sinalização do hormônio tireoidiano)**. Confirmado por: (a) inspeção da lista de símbolos; (b) `T05_kegg_missing_genes.tsv` (apenas PIK3R3 e PRKACG ausentes após filtro); (c) log de execução `KEGG: 121 genes | 119 in data`.
- Portanto, **todas as análises do repositório (DEG, volcano, PPI) foram realizadas sobre 119 genes pré-selecionados de uma via candidata**, não sobre o transcriptoma (~20.000 genes).

**Consequências diretas:**
1. O título/descrição "Análise Transcriptômica" é **impreciso**: trata-se de *análise de expressão diferencial restrita a uma via candidata*.
2. **GSEA não pode ser executado com o dado atual.** GSEA exige um *ranking global* de todos os genes testados (idealmente genoma inteiro). Uma lista de 119 genes de uma única via é estruturalmente incompatível com GSEA.
3. **A nova etapa exige re-obter a matriz completa.** `scripts/download_data.R` já aponta para o arquivo completo `TcgaTargetGtex_rsem_gene_tpm.gz` (TOIL RSEM, TCGA + GTEx). O arquivo commitado é um subconjunto derivado, **sem script de derivação versionado** (gap de reprodutibilidade).

> **Classificação:** LACUNA CRÍTICA. IMPACTO: bloqueia a Fase 6 (GSEA) do novo projeto até que a matriz completa seja obtida e validada.

---

## 3. TABELA DE AUDITORIA POR ETAPA (EXIGIDA)

Formato: **ETAPA | EXISTE | EXECUTADA | DOCUMENTADA | EVIDÊNCIA | LACUNA | RISCO**

| ETAPA | EXISTE (código) | EXECUTADA (evidência) | DOCUMENTADA | EVIDÊNCIA | LACUNA | RISCO |
|---|---|---|---|---|---|---|
| Obtenção de dados | Sim (`download_data.R`) | Parcial (arquivo completo NÃO commitado; commitado só subconjunto 121 genes) | Sim | `data/raw/XENA_THCA.tsv` (640 KB) | Sem script versionado que gera o subconjunto de 121 genes | ALTO — não reprodutível do zero |
| Importação/validação | Sim (`02_import.R`) | Sim | Sim | `logs/pipeline_*.log`, T01 | — | BAIXO |
| QC — escala log2 | Sim (`validate_expression_scale`) | Sim | Sim | log "QC: max=… non-int=…" | — | BAIXO |
| QC — PCA/MDS | Sim (`03b_pca.R`) | Sim | Sim | `Fig_QC_PCA_*.png/pdf`, log "PC1 35.0%" | — | MÉDIO (batch confundido) |
| QC — UMAP | Sim (`03e_umap_qc.R`) | Sim | Sim | `Fig_QC_UMAP_*.png/pdf` | `uwot` fora do renv.lock | MÉDIO |
| QC — outliers | Sim (`03d_qc_outliers.R`) | **FALHOU (skipped)** | **Não declarado como falha** | log: "QC skipped … argumentos implicam em número de linhas distintos: 0, 783" | `T_QC_outlier_candidates.tsv` e `T_QC_sample_composition.tsv` AUSENTES | ALTO — README lista outlier como QC executado |
| Heatmap DEGs | Sim (`03c_heatmap.R`) | **NÃO executado** (não chamado) | **Não declarado** | `FigS2_Heatmap_*` AUSENTE | `pheatmap` fora do renv.lock | MÉDIO |
| Expressão diferencial (limma) | Sim (`03_deg.R`) | Sim | Sim | `T03_deg_full_results.tsv` (119 linhas) | Escopo = 119 genes | ALTO — escopo restrito |
| Volcano plot | Sim (`04_volcano.R`) | Sim | Sim | `Fig1_Volcano_*.png/pdf` | — | BAIXO |
| Anotação KEGG | Sim (`fetch_kegg_genes`) | Sim | Sim | T02/T03/T05 | — | BAIXO |
| Rede PPI (STRING) | Sim (`05_ppi.R`) | Sim | Sim | `N01–N04`, `T06/T07`, `Fig2_PPI_*` | Rede mínima (8 nós/12 arestas) | ALTO — rede subdimensionada |
| Centralidade | Sim (`compute_centrality`) | Sim | Sim | `N03_centrality_metrics.tsv` | — | MÉDIO |
| Tabelas suplementares S1–S4 | Sim (`06_supplementary.R`) | Sim (geradas à parte) | Sim | `S1–S4` presentes | Não orquestrado no pipeline | BAIXO |
| Testes unitários | Sim (`tests/`) | **Não verificado neste ambiente** | Sim | código presente | — | BAIXO |
| Ambiente reprodutível (renv/Docker) | Sim | — | Parcial | renv.lock; Dockerfile 4.4.0 vs README 4.6.0 | **inconsistência de versão R** | ALTO |
| Session info | Sim | Sim | Sim | `logs/session_info.txt` | — | BAIXO |
| Checksums | Sim | Sim | Sim | `results/CHECKSUMS.md` | cobre só parte dos outputs | BAIXO |

---

## 4. FASE 2 — CARACTERIZAÇÃO DOS DADOS (VERIFICADO vs A CONFIRMAR)

| Atributo | Valor | Status |
|---|---|---|
| Fonte | UCSC Xena Browser (TOIL recompute) | Documentado |
| Dataset (arquivo completo) | `TcgaTargetGtex_rsem_gene_tpm.gz` | A CONFIRMAR (não commitado) |
| Dataset (commitado) | `XENA_THCA.tsv` — **121 genes da via hsa04919** | VERIFICADO |
| Identificador/bookmark | `c486b845ee2e750c3a9d2fc5145c8426` | Documentado |
| Tipo de dado | RNA-seq (expressão log₂, provável TPM/RSEM normalizado) | A CONFIRMAR (TPM vs expected_count — doc contraditória) |
| Plataforma | Illumina (TCGA/GTEx) | Documentado |
| Tecido | Tireoide | VERIFICADO |
| Tipo tumoral | Carcinoma papilífero de tireoide (THCA) | Documentado |
| Nº indivíduos | **Não discriminado** (só nº amostras) | INFORMAÇÃO A CONFIRMAR |
| Nº amostras | **783** (504 THCA + 279 Normal) | VERIFICADO no arquivo bruto |
| Nº por grupo | THCA=504; Normal=279 | VERIFICADO |
| Critérios de inclusão | categoria == "TCGA Thyroid Carcinoma" ou "GTEX Thyroid" | VERIFICADO no código |
| Critérios de exclusão | condição indefinida (NA) | VERIFICADO no código (0 excluídos aqui) |
| Amostras pareadas | **Não** (TCGA ≠ GTEx) | VERIFICADO |
| Amostras repetidas/duplicadas | 0 duplicatas de ID | VERIFICADO (783 IDs únicos) |
| Variáveis clínicas | Apenas `study`, `sample_type`, `primary_site`, `main_category` | VERIFICADO (sem idade/sexo/estadiamento) |
| Fatores de confusão | **Batch (coorte) perfeitamente confundido com condição** | VERIFICADO — limitação reconhecida no README |

> **Conclusão Fase 2:** tamanho amostral disponível e confirmado = **783 amostras**. Indivíduos, idade, sexo, subtipo histológico e estadiamento = **INFORMAÇÃO A CONFIRMAR** (não presentes no arquivo commitado).

---

## 5. FASE 3–5 — QC, PRÉ-PROCESSAMENTO E DEG (O QUE EXISTE E O QUE FALTA)

**Executado (verificado nos logs/resultados):**
- Validação de escala log₂ (`validate_expression_scale`).
- Filtro de baixa expressão: `rowMeans(E > 0.5) >= 0.1` → reteve **119/121** genes (removeu PIK3R3, PRKACG).
- DEG com `limma`: design `~0 + condition`, contraste `THCA − Normal`, `eBayes`, BH.
- Resultado DEG: **29 genes (9↑, 20↓)** — consistente entre T02, T03, T04, T07 e log.

**Não executado / lacunas de QC:**
- **Outliers:** script falhou e foi pulado (bug de dimensões no `03d_qc_outliers.R`).
- **Heatmap:** não orquestrado.
- **Não há análise formal de batch effect além de PCA/UMAP descritivo.** Batch está confundido com condição (TCGA=tumor, GTEx=normal), logo qualquer "correção" removeria o sinal biológico — a decisão de *não corrigir* é defensável, mas torna o contraste **não causalmente interpretável**.
- **Não há deconvolução de composição celular** (reconhecido como limitação).

**Observações metodológicas sobre o DEG atual:**
- `limma` sobre log₂-TPM/RSEM **sem `voom`, sem `trend=TRUE` e sem quality weights** é aceitável como aproximação, mas não é a prática mais robusta para RNA-seq. Método defensável para dado já em escala log₂, porém deve ser explicitado como escolha.
- O cutoff é |log₂FC|>1 e FDR<0.05. **Significância estatística (FDR) e tamanho de efeito (logFC) estão corretamente separados** na classificação.

---

## 6. DIFERENÇAS ENTRE "DOCUMENTADO" E "REAL" (INCONSISTÊNCIAS)

1. **Escopo do dado.** README descreve "RNA-seq do TCGA e GTEx (783 amostras)" sem deixar explícito que a matriz commitada contém **apenas 121 genes de uma via**. → Redação enganosa por omissão.
2. **Outlier detection** listada no README como etapa do pipeline, mas **falhou em runtime** e não foi registrada como falha.
3. **Heatmap (`03c_heatmap.R`) e `06_supplementary.R`** não são orquestrados por `run_pipeline.R`.
4. **`AUDIT_REPORT.md`** é referenciado em código e no disclaimer final, mas **não existe**.
5. **Versão R:** README/session_info = R 4.6.0; **Dockerfile = `rocker/r-ver:4.4.0`**.
6. **`uwot` e `pheatmap`** são usados por scripts mas **ausentes do `renv.lock`**.
7. **Git:** 6 commits "Add files via upload" — **sem histórico semântico** (impede auditoria de evolução).
8. **Contas GitHub divergentes:** remote = `leticiamariasantos`, README/CITATION = `santosry`.
9. **Transformação do dado:** `download_data.R` diz "RSEM expected_count", `docs` diz "RSEM expected_count, upper quartile normalized", URL diz `gene_tpm`. → **TPM vs count vs upper-quartile: A CONFIRMAR**.
10. **Versões de pacotes** (ggplot2 4.0.3, dplyr 1.2.1, igraph 2.3.2, etc.) internamente consistentes entre S1 e renv.lock, mas **não verificáveis** contra CRAN/Bioc neste ambiente → **INFORMAÇÃO A CONFIRMAR**.
11. **FigS1_PCA_THCA_vs_Normal** presente em `results/figures/`, mas **nenhum script atual o gera** (provável resíduo de versão anterior).

---

## 7. LACUNAS PRIORIZADAS (PARA O NOVO PROJETO)

| # | Lacuna | Severidade | Bloqueia |
|---|---|---|---|
| L1 | Matriz de expressão **global** ausente (só 121 genes) | CRÍTICA | Fase 6 (GSEA) |
| L2 | Sem GSEA / enriquecimento de vias | CRÍTICA | Objetivo geral novo |
| L3 | Sem análise de redundância/seleção racional de vias | CRÍTICA | Fase 7 |
| L4 | PPI restrito a 29 DEGs de 1 via (rede de 8 nós) | ALTA | Fase 8–9 |
| L5 | Sem priorização multicritério de alvos | ALTA | Fase 9 |
| L6 | Sem validação independente / sensibilidade | ALTA | Fase 14–15 |
| L7 | Bug em `03d_qc_outliers.R` (dimensões) | MÉDIA | Fase 3 |
| L8 | Orquestração incompleta (heatmap/suplementares) | MÉDIA | Reprodução |
| L9 | Inconsistência renv/Docker/versões | MÉDIA | Reprodução |
| L10 | Batch confundido com condição | INERENTE | Interpretação causal |

---

## 8. LIMITES DESTA AUDITORIA

- **Não re-executei o pipeline R** (R não disponível neste ambiente; a execução original usou R 4.6.0 + internet). A coluna "EXECUTADA" baseia-se em **evidência de artefatos** (logs, tabelas, figuras, checksums), não em re-execução.
- **Não validei as chamadas STRING/KEGG em tempo real** (requer internet + R). Os valores numéricos de DEG/PPI foram conferidos **por consistência interna** entre T02–T07, N01–N04 e logs.
- Re-execução de ponta a ponta é **recomendada** assim que houver ambiente R ≥ 4.1 + renv restaurado.

---

## 9. VEREDITO DA FASE 1

**O pipeline existente está funcional para o que faz (DEG + volcano + PPI de uma única via candidata), mas NÃO é uma análise transcriptômica global e NÃO suporta GSEA.**

A nova etapa **não pode ser implementada sobre o dado atual**. É obrigatório, antes de qualquer GSEA:

1. **Obter e validar a matriz de expressão global** (o próprio `download_data.R` já indica a fonte correta).
2. **Re-executar QC completo** (incluindo o outlier QC corrigido) sobre a matriz global.
3. **Re-executar DEG genoma-wide** (limma com design adequado) para produzir o **ranking global** necessário ao GSEA.
4. **Registrar cada decisão** conforme Fase 17/18.

Sem isso, qualquer "GSEA" seria tecnicamente inválido (lista de genes de via única ≠ ranking genômico).

---

*Fim da Fase 1. A implementação da nova análise fica condicionada à obtenção da matriz global e à re-execução do DEG genoma-wide.*
