# thyroid-volcano-ppi

**Análise Transcriptômica da Via de Sinalização do Hormônio Tireoidiano no Carcinoma de Tireoide e Potenciais Implicações para a Enfermagem de Precisão**

> **Versão:** 3.1.0 | **Data:** 2026-06-24 | **Tipo de estudo:** Exploratório, gerador de hipóteses

[![R >= 4.1](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![renv](https://img.shields.io/badge/renv-locked-blueviolet)](renv.lock)
[![Docker](https://img.shields.io/badge/Docker-ready-2496ED?logo=docker)](Dockerfile)

---

## Sobre este estudo

Estudo exploratório de bioinformática aplicada à oncologia tireoidiana. Foram analisados dados públicos de RNA-seq do TCGA e GTEx (783 amostras: 504 de carcinoma papilífero de tireoide, 279 de tecido tireoidiano normal) com foco na via de sinalização do hormônio tireoidiano (KEGG hsa04919, 119 genes).

> **Nota sobre a contagem de genes da via KEGG:** O website do KEGG exibe 78 genes no diagrama da via hsa04919 (genes manualmente curados com posição na ilustração). A função `keggGet("hsa04919")` do pacote KEGGREST, utilizada neste estudo, retorna 121 genes — incluindo todos os parálogos humanos e membros de famílias gênicas que compartilham KEGG Orthology (KO) com os genes do diagrama. Por exemplo, se o diagrama utiliza "PRKC", a API retorna PRKCA, PRKCB, PRKCG e demais membros da família. Dos 121 genes retornados, 119 intersectam com a matriz de expressão TCGA-GTEx (2 são excluídos pelo filtro de baixa expressão). Esta abordagem é conservadora e padrão em estudos transcriptômicos baseados em vias KEGG, pois evita a perda de genes funcionalmente relacionados.

### Hipótese científica

O carcinoma de tireoide apresenta alterações transcriptômicas em genes da via de sinalização do hormônio tireoidiano. A identificação desses genes e de suas interações proteicas pode gerar hipóteses biológicas relevantes para futuras investigações em oncologia molecular e enfermagem de precisão.

### Resultados

**Expressão diferencial.** Dos 119 genes da via de sinalização do hormônio tireoidiano (KEGG hsa04919) analisados, 29 (24,4%) foram diferencialmente expressos entre carcinoma de tireoide (THCA, n = 504) e tecido tireoidiano normal (GTEx, n = 279), considerando |log2FC| > 1 e FDR < 0,05. Destes, 9 genes encontram-se superexpressos e 20 subexpressos no tumor.

| Gene | log2FC | FDR | Regulação | Função biológica (STRING) |
|------|--------|-----|-----------|---------------------------|
| MYH7 | −5,59 | 7,3×10⁻²²⁴ | Subexpresso | Cadeia pesada da miosina-7; contração muscular, atividade ATPase dependente de actina |
| RXRG | +5,28 | 3,3×10⁻¹⁴⁹ | Superexpresso | Receptor nuclear; mecanismos de regulação transcricional |
| DIO3 | −4,97 | 4,4×10⁻¹⁶⁸ | Subexpresso | Metabolismo de hormônios tireoidianos; inativação de iodotironinas |
| MYH6 | −4,42 | 1,5×10⁻¹⁵⁷ | Subexpresso | Cadeia pesada da miosina-6 (alfa); isoforma predominante no átrio cardíaco |
| ATP2A1 | −2,83 | 1,3×10⁻¹⁶⁹ | Subexpresso | Transporte de íons cálcio; homeostase intracelular de cálcio |
| CCND1 | +2,65 | 1,0×10⁻²¹⁴ | Superexpresso | Regulação do ciclo celular (ciclina D1) |

**Genes adicionais com elevada significância** incluem DIO1 (log2FC = −2,60; FDR = 2,4×10⁻²⁷), PLCD3 (log2FC = +2,38; FDR = 1,2×10⁻⁹⁵), PLN (log2FC = −2,29; FDR = 5,9×10⁻⁹⁷), ATP1A2 (log2FC = −2,30; FDR = 1,2×10⁻⁸⁵), KAT2A (log2FC = −1,80; FDR = 6,1×10⁻¹³⁶), FXYD2 (log2FC = −1,74; FDR = 1,9×10⁻⁷⁵), TBC1D4 (log2FC = −1,64; FDR = 2,1×10⁻⁹⁴), STAT1 (log2FC = +1,57; FDR = 5,8×10⁻⁷³) e PLCD4 (log2FC = −1,53; FDR = 2,9×10⁻⁷¹). A lista completa dos 29 DEGs está disponível em `results/tables/T03_deg_full_results.tsv`.

**Rede PPI.** A consulta à base STRING v12.0 (Homo sapiens, taxon 9606) mapeou todos os 29 DEGs, dos quais 24 apresentaram pelo menos uma interação proteína-proteína com escore combinado ≥ 700 (alta confiança). Cinco genes não formaram interações acima desse limiar: RXRG, DIO1, MED13, TBC1D4 e PFKFB2.

A rede conectada final reteve 8 proteínas organizadas em 2 módulos funcionais detectados pelo algoritmo walktrap, com 12 arestas, modularidade 0,109, coeficiente de agrupamento 0,656 e densidade 0,429. Os 21 DEGs restantes não integraram a rede conectada, seja por ausência de interações de alta confiança ou por formarem componentes desconexos de pequeno tamanho.

| Módulo | Proteínas | Função predominante |
|--------|-----------|---------------------|
| **Módulo 1** (6 proteínas) | PRKCA, PRKCG, PLCG1, PLCD3, PLCD4, PIK3R2 | Transdução de sinais intracelulares mediada por fosfolipases C e proteínas quinase C |
| **Módulo 2** (2 proteínas) | ACTG1, ITGAV | Organização estrutural do citoesqueleto de actina e interações célula-matriz |

**Centralidade.** PRKCA (Proteína Quinase C Alfa) destacou-se como gene de maior centralidade na rede (grau = 5, betweenness = 0,476, closeness = 8,52×10⁻⁴, hub score = 1,000), atuando como único elo entre os dois módulos por meio da interação PRKCA–ACTG1 (escore STRING = 918). Os demais hubs identificados foram PLCG1 (grau = 4, betweenness = 0,286), PLCD4 (grau = 4, betweenness = 0,095) e PRKCG (grau = 4, betweenness = 0), todos restritos ao Módulo 1. PLCD3 (grau = 3) não atingiu os critérios de hub. As métricas completas de centralidade estão em `results/network/N03_centrality_metrics.tsv`.

---

## Como interpretar os resultados

### Volcano Plot

| Elemento | Significado |
|----------|-------------|
| **Eixo X** | log2(fold change): direção e magnitude da alteração. Positivo = mais expresso no tumor. Negativo = menos expresso |
| **Eixo Y** | -log10(valor-p ajustado): significância estatística. Mais alto = mais confiável |
| **Pontos azuis** | Genes superexpressos no tumor (9 genes) |
| **Pontos magenta** | Genes subexpressos no tumor (20 genes) |
| **Pontos cinza** | Genes sem diferença significativa (90 genes) |
| **Linhas tracejadas** | Limiares estatísticos: |log2FC| = 1 (vertical), FDR = 0,05 (horizontal) |
| **Rótulos** | 15 genes selecionados por significância, magnitude de fold-change e relevância na rede |

### Rede PPI

| Elemento | Significado |
|----------|-------------|
| **Nós (círculos)** | Proteínas codificadas pelos genes diferencialmente expressos |
| **Cores** | Módulos funcionais detectados por walktrap |
| **Tamanho do nó** | Proporcional ao grau de conectividade |
| **Borda escura** | Gene com elevada centralidade na rede |
| **Linhas** | Interações proteína-proteína (STRING, escore >= 700) |

**Importante:** As métricas de centralidade são exploratórias. A rede PPI reflete conhecimento acumulado na literatura e não demonstra causalidade, ativação ou inibição direta. Os resultados geram hipóteses, não validação clínica.

---

## Como este estudo foi conduzido

### Ambiente computacional

Todas as análises foram executadas em **R versão 4.6.0** (2026-04-26, "Jitterbug Madness") utilizando **Visual Studio Code versão 1.125** como ambiente de desenvolvimento integrado (IDE), com suporte das extensões R e Git. O `sessionInfo()` completo de cada execução é registrado em `logs/session_info.txt`.

### Desenho do estudo

| Etapa | Método | Ferramenta |
|--------|--------|-----------|
| Obtenção dos dados | Download do conjunto integrado TCGA-GTEx | UCSC Xena Browser |
| Pré-processamento | Filtro de baixa expressão, validação de escala log2 | R base, dplyr |
| Controle de qualidade | PCA e UMAP por condição e por fonte | limma (plotMDS), uwot |
| Expressão diferencial | Modelo linear com moderação empírica de Bayes | limma (lmFit, contrasts.fit, eBayes) |
| Correção para múltiplos testes | Benjamini-Hochberg (FDR < 0,05) | limma (topTable) |
| Anotação de via | KEGG hsa04919 | KEGGREST |
| Rede PPI | Interações de alta confiança (escore >= 700) | STRING v12.0 via REST API |
| Centralidade | Betweenness, degree, closeness, hub score | igraph |
| Comunidades | Walktrap | igraph |
| Visualização | Volcano plot + rede PPI | ggplot2, ggrepel, ggraph |

### Fontes de dados

| Fonte | Conteúdo | Acesso |
|--------|----------|--------|
| TCGA THCA | 504 amostras de carcinoma papilífero de tireoide | UCSC Xena Browser |
| GTEx Thyroid | 279 amostras de tecido tireoidiano normal | UCSC Xena Browser |
| KEGG hsa04919 | Via de sinalização do hormônio tireoidiano (121 genes) | KEGG REST API |
| STRING v12.0 | Interações proteína-proteína (Homo sapiens, taxon 9606) | STRING REST API |

### Parâmetros

| Parâmetro | Valor |
|-----------|-------|
| |log2FC| mínimo | 1,0 |
| FDR máximo | 0,05 |
| Escore STRING mínimo | 700 |
| Semente de reprodutibilidade | 42 |

---

## Como executar

### Requisitos

- R >= 4.1 (este estudo utilizou R v4.6.0)
- Conexão com internet (para STRING e KEGG)
- VS Code com extensão R (este estudo utilizou VS Code v1.125; RStudio também compatível)

### Passo a passo

**1. Obter o repositório**

```bash
git clone https://github.com/santosry/thyroid-volcano-ppi.git
cd thyroid-volcano-ppi
```

**2. Obter os dados**

```bash
Rscript scripts/download_data.R
```

O arquivo `XENA_THCA.tsv` será baixado para `data/raw/`. Alternativamente, acesse https://xenabrowser.net/?bookmark=c486b845ee2e750c3a9d2fc5145c8426 e faça o download manual.

**3. Instalar os pacotes**

```r
install.packages("renv")
renv::restore()
```

O `renv.lock` garante as versões exatas utilizadas pelos autores. A instalação leva de 5 a 15 minutos na primeira execução.

**4. Executar o pipeline**

```r
source("run_pipeline.R")
```

Tempo estimado: 3 a 5 minutos. O pipeline gera:

- `results/figures/` — Volcano plot e rede PPI (PNG 600 dpi + PDF vetorial)
- `results/tables/` — Sete tabelas TSV com resultados completos
- `results/network/` — Metadados da rede PPI
- `logs/` — Logs de execução e sessionInfo()

**5. Executar os testes**

```r
testthat::test_dir("tests/testthat")
```

---

## Pacotes utilizados

| Pacote | Versão | Função | Fonte |
|--------|--------|--------|-------|
| limma | 3.68.4 | Expressão diferencial, eBayes, BH, PCA (plotMDS) | Bioconductor |
| ggplot2 | 4.0.3 | Visualização (volcano plot) | CRAN |
| ggrepel | 0.9.8 | Rótulos com repulsão | CRAN |
| igraph | 2.3.2 | Construção e análise de redes, centralidade, walktrap | CRAN |
| ggraph | 2.2.2 | Visualização de redes | CRAN |
| dplyr | 1.2.1 | Manipulação de dados | CRAN |
| tidyr | 1.3.2 | Organização de dados | CRAN |
| readr | 2.2.0 | Leitura e exportação de TSV | CRAN |
| tibble | 3.3.1 | Estruturas de dados | CRAN |
| stringr | 1.6.0 | Processamento de texto | CRAN |
| purrr | 1.2.2 | Programação funcional | CRAN |
| here | 1.0.2 | Portabilidade de caminhos | CRAN |
| KEGGREST | 1.52.2 | Consulta à via KEGG hsa04919 | Bioconductor |
| org.Hs.eg.db | 3.23.1 | Anotação gênica humana | Bioconductor |
| AnnotationDbi | 1.74.0 | Interface de anotação | Bioconductor |
| httr | 1.4.8 | Requisições HTTP à STRING API | CRAN |
| jsonlite | 2.0.0 | Processamento de JSON da STRING API | CRAN |
| uwot | — | UMAP para controle de qualidade | CRAN |

Versões completas em `renv.lock` e `logs/session_info.txt`.

---

## Estrutura do repositório

```
thyroid-volcano-ppi/
├── run_pipeline.R              Script principal
├── R/
│   ├── 00_setup.R              Parâmetros, pacotes, cores, verificação de internet
│   ├── 01_functions.R          Funções (KEGG, STRING, exportação)
│   ├── 02_import.R             Importação e validação
│   ├── 03_deg.R                Expressão diferencial (limma)
│   ├── 03b_pca.R               PCA (controle de qualidade)
│   ├── 03c_heatmap.R           Heatmap DEGs
│   ├── 03d_qc_outliers.R       Detecção de outliers
│   ├── 03e_umap_qc.R           UMAP (controle de qualidade)
│   ├── 04_volcano.R            Volcano plot
│   ├── 05_ppi.R                Rede PPI
│   └── 06_supplementary.R      Tabelas suplementares
├── data/
│   ├── raw/                    XENA_THCA.tsv
│   └── processed/              Dados intermediários
├── scripts/
│   ├── download_data.R         Download automático dos dados
│   └── setup_renv.R            Inicialização do renv
├── results/
│   ├── figures/                PNGs 600 dpi + PDFs vetoriais
│   ├── tables/                 Tabelas TSV
│   └── network/                Metadados da rede PPI
├── tests/
│   └── testthat/               Testes unitários
├── docs/                       Documentação complementar
├── logs/                       Logs e sessionInfo()
├── Dockerfile                  Container reprodutível
├── renv.lock                   Versões exatas dos pacotes
├── LICENSE                     MIT
└── CITATION.cff                Metadados de citação
```

---

## Reprodutibilidade, interoperabilidade e portabilidade

### Reprodutibilidade

- `set.seed(42)` em todos os scripts que utilizam aleatoriedade
- Parâmetros centralizados em `R/00_setup.R`
- Caminhos via `here::here()` — nenhum caminho absoluto
- `renv.lock` com versões exatas de todos os pacotes R
- `sessionInfo()` capturado a cada execução em `logs/`
- `CHECKSUMS.md` com hashes MD5 de todos os outputs
- `Dockerfile` para ambiente Linux totalmente reprodutível
- Testes unitários em `tests/testthat/` com dados simulados (sem internet)

### Interoperabilidade

- Todos os outputs em TSV (valores separados por tabulação), legíveis por R, Python, Excel e qualquer linguagem
- Figuras em PNG (600 dpi, raster) e PDF (vetorial) para qualquer software de editoração
- Dados de entrada em TSV padronizado conforme exportação do UCSC Xena Browser
- Metadados das figuras documentados em `docs/figure_specs.md`
- Dicionário de dados em `docs/data_dictionary.md`

### Portabilidade

- `here::here()` resolve caminhos independentemente do sistema operacional
- `renv.lock` garante ambiente R idêntico em Windows, macOS e Linux
- `Dockerfile` permite execução em qualquer sistema com Docker
- Sem dependências de interfaces gráficas (totalmente executável em linha de comando)
- Verificação de conectividade (`check_internet()`) antes de chamadas a APIs externas

---

## Declaração de uso de Inteligência Artificial

Em conformidade com a Portaria CNPq nº 2.664/2026, declaramos:

Este projeto utilizou ferramentas de inteligência artificial como suporte técnico e metodológico durante as etapas de desenvolvimento de código, depuração, revisão de documentação e auditoria de qualidade científica.

### Ferramentas empregadas

| Ferramenta | Desenvolvedor | Etapa |
|------------|---------------|-------|
| DeepSeek-v4-pro | DeepSeek AI | Otimização de código R, auditoria de namespaces, revisão de funções estatísticas |
| Codex | OpenAI | Geração e depuração de scripts R, documentação técnica |
| ChatGPT 5.5 | OpenAI | Revisão textual, estruturação de documentação |
| Grok | xAI | Análise exploratória, prototipagem de visualizações |

### Natureza da participação humana

Em todas as etapas, a participação humana foi integral e soberana:

- **Nenhuma conclusão científica foi derivada exclusivamente por IA.** As hipóteses biológicas, a interpretação dos resultados e as discussões sobre relevância para a enfermagem de precisão foram formuladas pelos autores com base nos outputs do pipeline, na literatura científica e na experiência clínica da equipe.
- **As análises estatísticas foram executadas pelo pipeline em R**, com `set.seed(42)` garantindo reprodutibilidade determinística. Todos os resultados numéricos foram verificados manualmente pelos autores contra os arquivos de output em `results/tables/`.
- **Nenhum texto científico final foi gerado por IA.** As ferramentas de IA atuaram como assistentes de programação (geração, depuração e otimização de código R/Python) e revisão metodológica, não como redatores do conteúdo científico. Todo o conteúdo textual do resumo expandido e da documentação foi redigido, revisado e aprovado pelos autores.
- **As ferramentas de IA não substituem o julgamento científico.** Os autores assumem responsabilidade integral pela acurácia dos dados, pela validade das análises e pela adequação das conclusões apresentadas.

### Rastreabilidade

O registro completo das tarefas assistidas por IA está disponível em `results/tables/S4_ai_assisted_tasks.tsv`, contendo para cada tarefa: ferramenta utilizada, natureza da participação humana, método de validação e status de conformidade. A trilha de auditoria do pipeline (etapas, validações, dependências) está documentada em `results/tables/S3_pipeline_audit_trail.tsv`, permitindo verificação independente de cada etapa da análise.

---

## Limitações

1. **Comparação TCGA vs GTEx sem correção de batch effect.** Amostras tumorais (TCGA) e normais (GTEx) provêm de coortes distintas com protocolos de sequenciamento, processamento e perfis demográficos diferentes. Como TCGA corresponde a tumor e GTEx a normal, batch e condição estão perfeitamente confundidos, impedindo correção sem remover o sinal biológico. Os resultados devem ser interpretados como exploratórios.

2. **Estudo exploratório, não confirmatório.** A análise de expressão diferencial indica associação estatística, não relação causal. Os genes diferencialmente expressos podem ser consequência, e não causa, do processo neoplásico.

3. **Rede PPI in silico.** As interações proteína-proteína são preditas ou inferidas pelo STRING a partir de evidência combinada (mineração de texto, experimentos, coexpressão). Não há validação experimental direta.

4. **Métricas de centralidade exploratórias.** Betweenness, degree, closeness e hub score são métricas de rede. Não constituem validação de relevância funcional, não identificam alvos terapêuticos e não substituem ensaios experimentais.

5. **Foco em uma única via.** Apenas a via KEGG hsa04919 foi analisada. Genes fora dessa via, potencialmente relevantes, não foram considerados.

6. **Generalização limitada.** Os resultados aplicam-se ao carcinoma papilífero de tireoide (THCA) no contexto dos dados TCGA-GTEx. A extrapolação para outros subtipos histológicos requer validação independente.

---

## Autores

| Autor | ORCID | Afiliação |
|--------|-------|-----------|
| **Leticia Maria Dias Freitas** (autora correspondente) | [0009-0009-9930-9588](https://orcid.org/0009-0009-9930-9588) | Escola Técnica Estadual João Barcelos Martins (FAETEC), Campos dos Goytacazes, RJ |
| Ryan de Paulo Santos | [0009-0005-6770-2001](https://orcid.org/0009-0005-6770-2001) | Instituto Federal Fluminense (IFF), Campus Campos Guarus, Campos dos Goytacazes, RJ |
| Thais Faria Coutinho da Silva Pereira | [0009-0005-7091-2480](https://orcid.org/0009-0005-7091-2480) | Escola Técnica Estadual João Barcelos Martins (FAETEC), Campos dos Goytacazes, RJ |

**Contato:** leticiamariadiasfreitas@gmail.com

---

## Contribuições (CRediT)

| Autor | Contribuição |
|--------|-------------|
| Leticia Maria Dias Freitas | Conceitualização (Liderança); Metodologia (Igual); Software (Igual); Análise Formal (Igual); Curadoria de Dados (Igual); Validação (Igual); Visualização (Igual); Investigação (Igual); Redação — Rascunho Original (Liderança); Administração do Projeto (Suporte) |
| Ryan de Paulo Santos | Conceitualização (Suporte); Metodologia (Igual); Software (Igual); Análise Formal (Igual); Curadoria de Dados (Igual); Validação (Igual); Visualização (Igual); Investigação (Igual); Redação — Rascunho Original (Igual); Administração do Projeto (Liderança); Redação — Revisão e Edição (Suporte) |
| Thais Faria Coutinho da Silva Pereira | Supervisão (Liderança); Revisão Científica (Liderança); Validação (Suporte) |

---

## Licença

MIT License. Veja [LICENSE](LICENSE).

---

## Como citar

```bibtex
@software{freitas2026thyroid,
  title        = {thyroid-volcano-ppi: Análise Transcriptômica da Via de
                  Sinalização do Hormônio Tireoidiano no Carcinoma de
                  Tireoide e Potenciais Implicações para a Enfermagem
                  de Precisão},
  author       = {Leticia Maria Dias Freitas and Ryan de Paulo Santos
                  and Thais Faria Coutinho da Silva Pereira},
  year         = {2026},
  url          = {https://github.com/santosry/thyroid-volcano-ppi},
  note         = {v3.1.0}
}
```

---

## Referências

1. LEE, K.; ANASTASOPOULOU, C.; CHANDRAN, C.; CASSARO, S. Thyroid Cancer. In: STATPEARLS [Internet]. Treasure Island (FL): StatPearls Publishing, 2023. Disponível em: https://www.ncbi.nlm.nih.gov/books/NBK45929/. Acesso em: 24 jun. 2026.

2. TSAI, W. H. et al. Association between thyroid cancer and cardiovascular disease: a meta-analysis. *Frontiers in Cardiovascular Medicine*, v. 10, 2023. DOI: 10.3389/fcvm.2023.1075842.

3. ZHANG, B. et al. Integrated bioinformatics analysis for the identification of key genes and signaling pathways in thyroid carcinoma. *Experimental and Therapeutic Medicine*, v. 21, n. 3, 2021. DOI: 10.3892/etm.2021.9664.

4. FU, M. R. et al. Precision health: a nursing perspective. *International Journal of Nursing Sciences*, v. 7, n. 1, p. 5-12, 2020. DOI: 10.1016/j.ijnss.2019.12.008.

5. GOLDMAN, M. J. et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology*, v. 38, n. 6, p. 675-678, 2020. DOI: 10.1038/s41587-020-0546-8.

6. RITCHIE, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, v. 43, n. 7, e47, 2015. DOI: 10.1093/nar/gkv007.

7. SZKLARCZYK, D. et al. The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research*, v. 51, n. D1, p. D638-D646, 2023. DOI: 10.1093/nar/gkac1000.

8. CANCER GENOME ATLAS RESEARCH NETWORK. Integrated genomic characterization of papillary thyroid carcinoma. *Cell*, v. 159, n. 3, p. 676-690, 2014. DOI: 10.1016/j.cell.2014.09.050.

9. KANEHISA, M. et al. KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Research*, v. 51, n. D1, p. D587-D592, 2023. DOI: 10.1093/nar/gkac963.

10. CSARDI, G. et al. *igraph*: network analysis and visualization. Version 2.0.3. [Software]. CRAN, 2024. Disponível em: https://CRAN.R-project.org/package=igraph. Acesso em: 24 jun. 2026.

11. PEDERSEN, T. L. *ggraph*: an implementation of grammar of graphics for graphs and networks. Version 2.2.1. [Software]. CRAN, 2024. Disponível em: https://CRAN.R-project.org/package=ggraph. Acesso em: 24 jun. 2026.

---

*Pipeline mantido por [Ryan de Paulo Santos](https://github.com/santosry), ORCID: [0009-0005-6770-2001](https://orcid.org/0009-0005-6770-2001)*
