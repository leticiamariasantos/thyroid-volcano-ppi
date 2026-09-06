# PLANO OPERACIONAL — FASE 2 (TRANSCRIPTOMA GLOBAL)
## thyroid-volcano-ppi → extensão para caracterização transcriptômica global do carcinoma de tireoide

**Status:** PLANO PRÉ-EXECUÇÃO (nada de GSEA/DEG/PPI novo foi executado ainda).
**Data:** 2026-09-06
**Precedência decisória (imutável):** validade > qualidade/compatibilidade dos dados > reprodutibilidade > rigor estatístico > coerência biológica > robustez/sensibilidade > inovação > aplicabilidade translacional > adequação ao evento.

---

## 0. ESTADO DO AMBIENTE (VERIFICADO NESTA SESSÃO)

| Recurso | Estado | Implicação |
|---|---|---|
| R / Rscript | **DISPONÍVEL — R 4.6.1** (`C:\Program Files\R\R-4.6.1\bin\R.exe`, "Happy Hop") | Fase 2 pode ser conduzida em R, mantendo continuidade com a Fase 1 |
| Pacotes R para GSEA/DE/PPI | limma 3.68.0, edgeR 4.10.0, DESeq2 1.51.7, KEGGREST 1.52.0, org.Hs.eg.db 3.23.1, igraph 2.3.0, ggraph 2.2.2, fgsea 1.37.4, clusterProfiler 4.20.0, msigdbr 26.1.0, data.table 1.18.4 | OK — stack completo |
| Python | 3.14.5 | Alternativa (não necessária para a Fase 2) |
| pandas / numpy / scipy | 2.3.3 / 2.4.6 / 1.16.3 | OK |
| statsmodels | 0.14.6 | Modelos lineares / moderados |
| gseapy | 1.3.1 (instalado) | GSEA |
| networkx | 3.6.1 | PPI / centralidade |
| scikit-learn | 1.9.0 | PCA/clustering (se necessário) |
| matplotlib / seaborn | 3.11.0 / 0.13.2 | Figuras |
| Internet | OK | Download + KEGG REST |
| KEGG REST (HTTPS) | OK | Gene sets de vias |
| Disco livre | 265 GB | Suficiente |
| RAM | 8,3 GB total (~1,4 GB livre) | **Matriz global requer extração streaming** |

**Decisão de software (registrada, não assumida):** com R 4.6.1 disponível e com o stack completo (limma/edgeR/DESeq2 + fgsea/clusterProfiler/msigdbr + igraph), **a Fase 2 será conduzida em R**, mantendo continuidade metodológica com a Fase 1. **Divergência de versão registrada:** o projeto documenta R 4.6.0 (renv.lock/README); a instalação disponível é **R 4.6.1**. Pacotes-chave em versões de patch ligeiramente diferentes do `renv.lock` (ex.: limma 3.68.0 vs 3.68.4). Isto será registrado em `sessionInfo()` e na trilha de execução — não é bloqueante, mas é uma diferença de ambiente a declarar. A reprodutibilidade da Fase 2 será garantida por `sessionInfo()` + seeds + parâmetros + checksums (espelhando o padrão da Fase 1).

---

## 1. MATRIZ GLOBAL — AQUISIÇÃO E VALIDAÇÃO (BLOQUEIO DE ENTRADA)

**Origem:** UCSC Xena — `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz`
**Tamanho:** 1.323.254.426 bytes (gz) · **Last-Modified:** 2021-04-09

| Campo do protocolo | Conteúdo |
|---|---|
| PERGUNTA | Qual é a matriz transcriptômica global TCGA-THCA + GTEx-tireoide e qual a sua natureza (unidade, escala, dimensões)? |
| HIPÓTESE | A matriz completa contém a expressão genoma-wide (genes × amostras) em TPM/RSEM normalizado, da qual se extraem as mesmas 783 amostras já usadas na Fase 1. |
| DADOS DE ENTRADA | URL acima (gz); `data/raw/XENA_THCA.tsv` (Fase 1) como referência de IDs das 783 amostras |
| MÉTODO | Download com verificação de integridade (hash/byte-count); leitura **streaming** (linha a linha) extraindo apenas as colunas das 783 amostras; conversão para float32; checagem de dimensões, símbolos gênicos, duplicatas, NA |
| PARÂMETROS | IDs das amostras = os 783 presentes na Fase 1 (504 THCA + 279 GTEx Normal); genes = todas as linhas da matriz |
| OUTPUT ESPERADO | `data/global/TCGA_GTEx_thyroid_tpm.tsv` (~60k genes × 783 amostras) + manifesto (dimensões, checksum, data de acesso, unidade) |
| CRITÉRIO DE SUCESSO | (a) arquivo íntegro; (b) nº de genes ≈ genoma (dezenas de milhares); (c) 783 colunas com IDs exatamente iguais aos da Fase 1; (d) valores ≥ 0 e em escala plausível |
| CRITÉRIO DE PARADA | Download corrompido; IDs não coincidirem com a Fase 1; matriz não for genoma-wide; unidade não puder ser estabelecida |
| RISCO DE VIÉS | Extração de colunas errada (offset); duplicação de IDs; genes com sufixo `|` (parálogos) mal tratados |
| SENSIBILIDADE | Comparar IDs e contagens com a Fase 1; conferir com a descrição oficial do TOIL (TPM, upper-quartile vs não) |
| INTERPRETAÇÃO PERMITIDA | "A matriz global contém N genes e 783 amostras, com unidade X, verificada por …" |
| INTERPRETAÇÃO NÃO PERMITIDA | Assumir TPM vs counts sem verificar; tratar amostra como paciente |

> **Verificação de unidade (obrigatória, antes de qualquer DEG):** inspecionar distribuição (min/max/quantis), presença de valores inteiros vs contínuos, e cruzar com a documentação TOIL. Registrar explicitamente **TPM | expected_count | FPKM | upper-quartile**. Se não for determinável → **INFORMAÇÃO A CONFIRMAR** e **parar** antes do DEG.

---

## 2. VALIDAÇÃO DO PROCESSAMENTO (PRÉ-QC)

- Comparabilidade TCGA × GTEx: registrar que **origem = condição** (TCGA=tumor, GTEx=normal). Não é batch convencional.
- Harmonização: **não aplicar** correção de batch automaticamente. Se algo for aplicado, justificar problema/pressupostos/perda biológica/método/sensibilidade.
- Anotação gênica: usar símbolos HGNC; tratar sufixos `|` do TOIL (parálogos) com regra explícita e registrada.
- Duplicatas de genes/amostras: checar e registrar regra de colapso (ex.: média, ou maior variância — **decidir antes**).

---

## 3. CONTROLE DE QUALIDADE (QC)

| Campo | Conteúdo |
|---|---|
| PERGUNTA | As amostras e genes atendem a critérios mínimos de qualidade para DEG genoma-wide? |
| MÉTODO | distribuição por amostra; library-size/equivalente (soma TPM); nº genes expressos; PCA (top genes variáveis); clustering hierárquico + correlação; outliers (mediana/IQR/% zeros); separação tumor×normal; estruturas técnicas (coorte) |
| PARÂMETROS | limiares de outliers **definidos antes** (ex.: >3 IQR em ≥2 métricas OU %genes detectados < limiar) |
| OUTPUT | `02_qc/` figuras + tabela de amostras flaggeadas + justificativa de qualquer exclusão |
| CRITÉRIO DE SUCESSO | QC executado sem erro; cada exclusão documentada; separação biológica visível sem dominância técnica total |
| CRITÉRIO DE PARADA | dados inadequados; batch técnico dominante; script QC falhar |
| RISCO DE VIÉS | excluir amostra para "melhorar" resultado; ignorar estrutura de coorte |
| SENSIBILIDADE | PCA com/sem top genes; limiares de outliers distintos |
| INTERPRETAÇÃO PERMITIDA | descritiva: separação, outliers, efeito de coorte |
| INTERPRETAÇÃO NÃO PERMITIDA | "QC passou" se o script falhar; afirmar ausência de batch sem teste |

> **Correção obrigatória:** o script de outlier da Fase 1 (`03d_qc_outliers.R`) falhou por bug de dimensões. A Fase 2 reimplementa o QC de outliers em Python **do zero**, com log explícito.

---

## 4. EXPRESSÃO DIFERENCIAL GENOMA-WIDE

| Campo | Conteúdo |
|---|---|
| PERGUNTA | Quais genes estão diferencialmente expressos (THCA vs Normal) em escala global? |
| HIPÓTESE | Existe um conjunto genoma-wide de genes diferencialmente expressos além da via hsa04919. |
| MÉTODO | Dado log₂-TPM/RSEM → modelo linear por gene (OLS moderado: `statsmodels` com regularização empírica aproximada, OU `limma`-like shrinkage manual). **Decisão registrada após validar a unidade.** |
| PARÂMETROS | contraste THCA − Normal; FDR (BH) < 0.05; |log₂FC| ≥ 1 para classificação; filtro de baixa expressão pré-definido |
| OUTPUT | tabela completa ranqueada (gene, logFC, estatística, p, FDR); nº avaliados; nº significativos; up/down; volcano; MA; heatmap |
| CRITÉRIO DE SUCESSO | ranking global gerado; nº de genes avaliados = genoma (dezenas de milhares); contagens up/down reportadas |
| CRITÉRIO DE PARADA | método incompatível com a unidade dos dados; nº de genes avaliados ≈ 119 (indicaria matriz errada) |
| RISCO DE VIÉS | confundimento TCGA×GTEx (estrutural); escolher método por popularidade |
| SENSIBILIDADE | thresholds de expressão; FDR 0.01/0.05/0.10; método (OLS vs moderado) |
| INTERPRETAÇÃO PERMITIDA | "associação estatística de expressão" (Nível 2) |
| INTERPRETAÇÃO NÃO PERMITIDA | causalidade; "driver"; "biomarcador clínico" |

---

## 5. RANKING PARA GSEA

- **Definir ANTES de ver vias.** Ranking = estatística que preserva magnitude e direção: `rank_stat = sign(logFC) × (-log10(p))` (ou estatística t). Não usar lista de DEGs.
- Documentar a métrica e sua motivação.

---

## 6. GSEA GLOBAL (EXPLORATÓRIO)

| Campo | Conteúdo |
|---|---|
| PERGUNTA | Quais conjuntos gênicos/funcionais estão associados à condição no transcriptoma global? |
| MÉTODO | GSEA pré-ranqueado (gseapy `prerank`); universo = **todos os gene sets KEGG disponíveis** (via KEGG REST) + opcional MSigDB Hallmark |
| PARÂMETROS | método `fgsea`; permutações ≥ 10.000; min=15, max=500; ranking da §5; FDR<0.25 (triagem) e <0.05 |
| OUTPUT | tabela GSEA completa (term, NES, p nominal, FDR, leading edge) + figuras de enriquecimento |
| CRITÉRIO DE SUCESSO | GSEA executado sobre ranking global; todos os termos reportados (não só os significativos) |
| CRITÉRIO DE PARADA | sem resultados robustos; ranking inválido |
| RISCO DE VIÉS | restringir universo a vias "desejadas"; ler NES como causalidade |
| SENSIBILIDADE | database (KEGG vs Hallmark); min/max gene set; métrica de ranking |
| INTERPRETAÇÃO PERMITIDA | "enriquecimento associado à condição" (Nível 2) |
| INTERPRETAÇÃO NÃO PERMITIDA | "via ativada/inibida" sem evidência proteica |

---

## 7. PAINEL A PRIORI — 10 VIAS KEGG (PRÉ-ESPECIFICADO, IMOBILIZADO)

Vias (fixadas antes de ver resultados, **não alteráveis**):
hsa05216, hsa04919, hsa04010, hsa04151, hsa04150, hsa04115, hsa04210, hsa04110, hsa04310, hsa04064.

| Campo | Conteúdo |
|---|---|
| PERGUNTA | O que o transcriptoma global mostra **especificamente** para cada uma das 10 vias pré-especificadas? |
| MÉTODO | Para cada via: nº genes na matriz; nº DEGs; genes alterados; direção; estatística; FDR; NES (do GSEA quando aplicável); cobertura; leading-edge; sobreposição com as demais |
| OUTPUT | tabela do painel (1 linha por via) + interpretação por via |
| CRITÉRIO DE SUCESSO | todas as 10 vias reportadas, inclusive as **não enriquecidas** |
| CRITÉRIO DE PARADA | via sem genes na matriz (reportar como NA) |
| RISCO DE VIÉS | apresentar o painel como "top 10"; adicionar/remover vias após ver resultados |
| SENSIBILIDADE | cobertura com min/max gene set distintos |
| INTERPRETAÇÃO PERMITIDA | "enriquecimento associado à condição"; "alteração transcriptômica compatível com"; "sinal de enriquecimento"; "perturbação do programa molecular" |
| INTERPRETAÇÃO NÃO PERMITIDA | "via ativada"; inferência de atividade proteica a partir de RNA |

---

## 8. ANÁLISE DE REDUNDÂNCIA ENTRE VIAS

- Matriz de sobreposição (Jaccard) entre gene sets significativos + painel.
- Rede de similaridade de vias; identificação de **módulos de vias relacionadas**.
- Não tratar 10 vias sobrepostas como 10 evidências independentes; **não inflar o nº de achados**.

---

## 9. SELEÇÃO DE VIAS PARA PPI (CRITÉRIOS PRÉ-DEFINIDOS)

Selecionar para PPI apenas vias com: FDR relevante **e** NES coerente **e** baixa redundância **e** coerência com a biologia do carcinoma de tireoide **e** nº de genes compatível com rede. Número **não imposto** (pode ser 3, 5, 7…). Critérios congelados **antes** de ver a PPI.

---

## 10. PPI

| Campo | Conteúdo |
|---|---|
| PERGUNTA | Quais módulos de interação emergem dos genes das vias selecionadas? |
| MÉTODO | STRING (REST API), Homo sapiens; genes das vias selecionadas (conjunto de entrada justificado); centralidade (degree, betweenness, closeness); comunidades |
| PARÂMETROS | score mínimo ≥ 0.700 (alta confiança), com sensibilidade 0.400/0.900; registrar versão STRING |
| OUTPUT | rede(s) + métricas (nós, arestas, componentes, modularidade, módulos) |
| CRITÉRIO DE SUCESSO | rede construída com conjunto justificado; métricas reportadas |
| CRITÉRIO DE PARADA | rede insuficiente (nós/arestas abaixo de limiar pré-definido) |
| RISCO DE VIÉS | escolher genes "interessantes" após ver resultados; chamar hub de alvo |
| SENSIBILIDADE | score STRING; métrica de centralidade |
| INTERPRETAÇÃO PERMITIDA | "associação em rede" (Nível 3); hub = topologia |
| INTERPRETAÇÃO NÃO PERMITIDA | hub = biomarcador/driver/alvo/causal |

---

## 11. MÓDULOS

- Módulos = comunidades da rede; coerência funcional descritiva; **não** inferir função nova sem anotação.

---

## 12. PRIORIZAÇÃO MULTICRITÉRIO DE ALVOS

| Campo | Conteúdo |
|---|---|
| PERGUNTA | Quais genes/proteínas apresentam convergência de evidências (transcriptoma + via + rede + anotação)? |
| MÉTODO | ranking multicritério **definido antes**; pesos justificados e submetidos a sensibilidade |
| CRITÉRIOS | |logFC|, FDR, contribuição GSEA, presença em vias, centralidade, módulo, localização celular, membrana/extracelular, função, evidência externa, consistência entre datasets, plausibilidade de mensuração/intervenção |
| OUTPUT | tabela final (gene, logFC, FDR, via(s), NES, grau, betweenness, módulo, localização, evidência, estabilidade, plausibilidade diagnóstica/terapêutica/nanomédica, nível de evidência) |
| CRITÉRIO DE SUCESSO | ranking estável em sensibilidade de pesos |
| CRITÉRIO DE PARADA | candidatos instáveis; dependência de escolha arbitrária |
| RISCO DE VIÉS | soma arbitrária; pesos pós-hoc |
| INTERPRETAÇÃO PERMITIDA | "candidato priorizado por convergência computacional" |
| INTERPRETAÇÃO NÃO PERMITIDA | "alvo terapêutico"; "biomarcador clínico" |

> **Campos não disponíveis = `NA` ou `INFORMAÇÃO A CONFIRMAR`.** Proibido inventar.

---

## 13. VALIDAÇÃO EXTERNA

- Buscar dataset independente (outra coorte GEO/proteômica) **com critério definido antes** de observar o resultado.
- Avaliar direção/magnitude/consistência de genes, vias e candidatos.
- Se não houver dataset metodologicamente comparável → declarar **"Não foi realizada validação externa independente."**

---

## 14. ANÁLISE DE SENSIBILIDADE

- Thresholds de expressão, critérios de DEG, FDR, parâmetros GSEA, gene sets, score PPI, métricas de centralidade, influência de hubs, estabilidade dos candidatos e dos pesos.
- Resultado que desaparece sob pequena variação = **instável**.

---

## 15. PLAUSIBILIDADE NANOMÉDICA (OBRIGATORIAMENTE POSTERIOR)

- Só após: transcriptoma global → DEG → GSEA → vias → redundância → PPI → módulos → priorização → validação/sensibilidade.
- Pergunta-guia: *"Quais características moleculares observadas poderiam, de forma racional e explicitamente hipotética, justificar investigação futura por nanomedicina?"*
- **Não** escolher nanopartícula/magnética/hipertermia/biossensor a priori.
- Para cada possibilidade, separar: **EVIDÊNCIA DIRETA × INFERÊNCIA MOLECULAR × HIPÓTESE NANOMÉDICA**.
- Se não houver fundamento suficiente → declarar explicitamente e classificar como **hipótese futura**.

---

## 16. ESTRUTURA DE DIRETÓRIOS DA FASE 2

```
00_audit/ 01_data/ 02_qc/ 03_preprocessing/ 04_differential_expression/
05_gsea/ 06_kegg_panel/ 07_pathway_redundancy/ 08_ppi/
09_target_prioritization/ 10_validation/ 11_sensitivity/
12_nanomedicine/ 13_reports/ 14_figures/
```

A Fase 1 (R) permanece intacta. Fase 2 em Python, com rastreabilidade explícita entre as duas fases.

---

## 17. PRÓXIMA AÇÃO CONCRETA (GATE)

1. **Confirmar a decisão Python-vs-R** (R ausente ⇒ Python é a única via executável aqui; registrada acima).
2. **Baixar e validar a matriz global** (1,32 GB) com extração streaming das 783 amostras.
3. Só então: QC → DEG → ranking → GSEA → painel → redundância → PPI → priorização → validação → sensibilidade → (nanomedicina por último).

> **Este plano não executa nenhuma análise nova.** A execução começa pela aquisição/validação da matriz global, após este registro.

*Plano operacional fechado em 2026-09-06. Nenhuma via do GSEA global, alvo ou nanotecnologia foi pré-escolhido além do painel a priori de 10 vias KEGG explicitamente especificado.*
