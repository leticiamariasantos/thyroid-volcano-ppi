# PROTOCOLO DAS NOVAS ANÁLISES (PRÉ-REGISTRO METODOLÓGICO)
## GSEA → seleção de vias → PPI → priorização molecular → (só depois) avaliação de plausibilidade nanotecnológica

**Status:** PROTOCOLO. Nenhuma via, gene, alvo ou tecnologia foi pré-selecionado.
**Princípio:** validade científica > reprodutibilidade > coerência biológica > inovação > adequação ao evento.
**Condição de entrada (gate):** obter a matriz de expressão **global** (ver `AUDIT_PHASE1_REPOSITORIO.md`, lacuna L1). Sem matriz global, **nenhuma etapa abaixo é executável de forma válida**.

---

## 0. PERGUNTA CIENTÍFICA VINCULANTE

- **Principal:** Quais alterações transcriptômicas e vias biológicas caracterizam o carcinoma de tireoide e quais alvos moleculares podem ser priorizados computacionalmente para futuras aplicações em nanomedicina?
- **Secundária:** Os alvos identificados apresentam características que permitem uma conexão **biologicamente plausível** com estratégias nanotecnológicas?

> A pergunta secundária só é respondida **após** a obtenção dos resultados moleculares. A hipótese nanotecnológica é **exploratória** e pode ser refutada (resultado válido).

---

## 1. FASE 2 — CARACTERIZAÇÃO DOS DADOS (MATRIZ GLOBAL)

| Item a registrar | Obrigatório |
|---|---|
| Fonte / dataset / identificador | UCSC Xena `TcgaTargetGtex_rsem_gene_tpm.gz` (bookmark a confirmar) |
| Tipo de dado / plataforma | RNA-seq TOIL RSEM (TPM vs count — **verificar antes**) |
| Tecido / tipo tumoral | Tireoide / THCA (e normal GTEx) |
| Nº de amostras / por grupo | A CONFIRMAR após download (esperado 504 + 279, **não assumir**) |
| Nº de indivíduos | A CONFIRMAR (não confundir amostra com paciente) |
| Critérios inclusão/exclusão | Documentar explicitamente |
| Pareamento / duplicatas | Verificar IDs; tratar duplicatas |
| Variáveis clínicas / confundidores | Inventariar o que existir; registrar ausências |

**Regra:** qualquer número não verificável = **INFORMAÇÃO A CONFIRMAR**. Proibido estimar.

---

## 2. FASE 3 — CONTROLE DE QUALIDADE (matriz global)

Procedimentos obrigatórios (adequados a RNA-seq):
- distribuição da expressão por amostra (densidade / boxplot);
- fração de genes pouco expressos; genes não detectados;
- outliers por amostra (métricas: mediana, IQR, % zeros; **corrigir o bug do `03d_qc_outliers.R`** antes de reusar);
- PCA/MDS por **condição** e por **coorte** (avaliação de batch);
- clustering hierárquico + correlação amostral;
- duplicatas de ID; valores ausentes (NA);
- batch effect: quantificar (não apenas descrever).

**Critério objetivo de exclusão (definido ANTES de olhar o resultado):**
- ex.: amostra com fração de genes detectados < limiar pré-definido OU desvio > 3 IQR em ≥2 métricas de QC.
- **Nunca excluir amostra apenas porque altera o resultado.** Registrar cada exclusão + justificativa.

---

## 3. FASE 4 — PRÉ-PROCESSAMENTO

- Normalização/transformação: **decidir após inspeção** (TPM log₂ vs counts → voom/DESeq2).
- Filtragem: genes de baixa expressão (critério pré-definido).
- Anotação: símbolo ↔ Entrez (org.Hs.eg.db) para GSEA; tratar genes duplicados/1-para-muitos.
- Valores ausentes: registrar e tratar (ou excluir) com critério.
- Batch: **não corrigir cegamente**. TCGA=tumor e GTEx=normal confundem batch com condição. Se corrigir, o sinal biológico é destruído; documentar a decisão e suas consequências.

---

## 4. FASE 5 — EXPRESSÃO DIFERENCIAL (GENOMA-WIDE)

- **Método:** decidir conforme o dado:
  - counts → DESeq2 ou edgeR (ou limma-voom);
  - log₂ normalizado (TPM/RSEM) → limma (`lmFit + eBayes`), com `trend=TRUE`/quality weights se aplicável.
  - **Não escolher ferramenta por popularidade.**
- **Design/contraste:** `~ 0 + condition`; contraste `THCA − Normal`.
- Covariáveis: incluir se houver variáveis clínicas disponíveis (senão, declarar ausência).
- **FDR:** Benjamini–Hochberg; **cutoff** pré-definido (FDR < 0.05; |log₂FC| ≥ 1 para classificação de DEG).
- Saídas: tabela completa (ranking), genes significativos, volcano, MA, heatmap, diagnósticos.
- **Separar:** significância (FDR) × tamanho de efeito (logFC) × relevância biológica.
- **Produto crítico:** uma **tabela ranqueada de TODOS os genes testados** (logFC, estatística, p, FDR) — insumo do GSEA.

---

## 5. FASE 6 — GSEA (PROTOCOLO)

**Pré-condição:** ranking global de todos os genes testados (não apenas DEGs).

| Parâmetro | Definição |
|---|---|
| Métrica de ranking | `-log10(p) × sign(logFC)` (ou estatística t moderada). **Não usar apenas lista de DEGs.** |
| Método | GSEA clássico (pre-ranked) — `fgsea` (recomendado) ou `clusterProfiler::GSEA` |
| Database | MSigDB (Hallmark, C2 CP, C5 GO BP, C6 oncogenic) **e/ou** KEGG/Reactome via `msigdbr`/`clusterProfiler` |
| Versão do database | Registrar exatamente (ex.: MSigDB v2024.1.Hs; KEGG release) |
| Nº permutações | ≥ 10.000 (fgsea: default/`nPermSimple`) |
| Tamanho min/max do gene set | ex.: 15–500 (justificar) |
| FDR | < 0.25 (convenção GSEA) para triagem; reportar também < 0.05 |
| NES | reportar; usar na seleção de vias |

**Proibições:**
- Não usar como ranking apenas os genes que passaram cutoff de DE.
- Não re-executar DE depois do GSEA.
- Não pré-selecionar MAPK, PI3K-AKT, p53, JAK-STAT, etc. Essas só entram **se** emergirem dos dados.

---

## 6. FASE 7 — SELEÇÃO RACIONAL DAS VIAS

**Critérios de seleção (aplicados em ordem, definidos ANTES de ver resultados):**
1. FDR do GSEA;
2. NES (direção e magnitude);
3. tamanho do gene set (nem minúsculo nem gigante);
4. coerência biológica com carcinoma de tireoide (avaliada APÓS a lista surgir);
5. **redundância** (overlap de genes entre vias; via similarity/leading-edge analysis);
6. independência entre vias;
7. potencial de gerar hipóteses moleculares.

**Regras anti-viés:**
- Identificar **todas** as vias significativas primeiro; depois filtrar.
- **Não** escolher "as 5 menores p-values" por padrão.
- Se houver 5 vias excessivamente redundantes, reduzir. Se houver só 3 robustas, trabalhar com 3.
- O número de vias **não é imposto artificialmente**.

---

## 7. FASE 8 — PPI (ESTRATÉGIA)

- Ferramenta: STRING (REST API), versão a registrar; espécie 9606.
- Score mínimo: **pré-definido** (ex.: ≥ 0.700 alta confiança), com **análise de sensibilidade** (0.400 / 0.700 / 0.900).
- Registrar: nº nós, nº arestas, % mapeamento, componentes.
- Métricas de centralidade: degree, betweenness, closeness (+ eigenvector/hub se justificado).
- **Semântica obrigatória:** hub **não é** biomarcador **nem** alvo terapêutico. Apenas topologia.
- Se a rede for insuficiente (< limiar pré-definido de nós/arestas), **interromper** e registrar.

---

## 8. FASE 9 — PRIORIZAÇÃO MOLECULAR (CRITÉRIOS OBJETIVOS)

Pontuação multicritério (pesos definidos **antes** de aplicar; registrar):
1. magnitude DE (|logFC|) e FDR;
2. NES da via à qual pertence;
3. centralidade na PPI (degree/betweenness);
4. localização subcelular (membrana/extracelular > citoplasma > núcleo) — se disponível;
5. estrutura experimental disponível (PDB) / predita (AlphaFold);
6. domínios funcionais;
7. evidência experimental prévia e em carcinoma de tireoide;
8. consistência entre datasets (se houver validação).

**Vocabulário (não são sinônimos):** gene diferencialmente expresso ≠ hub ≠ biomarcador ≠ alvo molecular ≠ alvo terapêutico ≠ alvo de direcionamento nanotecnológico.

**Critério de interrupção:** se nenhum candidato atingir os limiares mínimos pré-definidos, declarar "sem candidato prioritário suficiente" — resultado válido.

---

## 9. FASE 10 — CARACTERIZAÇÃO ESTRUTURAL (CONDICIONAL)

Só executar **se** responder uma pergunta científica definida:
- Estrutura experimental (PDB) / predita (AlphaFold), qualidade, domínios, sítios, acessibilidade.
- **Docking só se houver hipótese molecular clara + justificativa de ligante.** Não adicionar docking "para complexidade".
- Caso contrário: **não executar** e registrar o motivo.

---

## 10. FASE 11–12 — NANOTECNOLOGIA COMO RESULTADO EMERGENTE (ADIADA)

**Gate:** só prosseguir se a resposta à pergunta for **SIM**:
> "Os resultados obtidos permitem estabelecer uma conexão científica defensável com nanotecnologia?"

Se **NÃO**: registrar "não há evidência suficiente para conexão nanotecnológica" e encerrar essa camada (resultado científico válido, não fracasso).

Se **SIM**, para cada conexão proposta, preencher obrigatoriamente:
| ALVO | EVIDÊNCIA MOLECULAR | CARACTERÍSTICA RELEVANTE | APLICAÇÃO POSSÍVEL | EVIDÊNCIA NA LITERATURA | NÍVEL DE EVIDÊNCIA | LIMITAÇÃO | PRÓXIMO EXPERIMENTO |

- A tecnologia é escolhida **após** a caracterização dos alvos, nunca antes.
- "Nanomedicina personalizada" só com dados individuais de pacientes; caso contrário usar "estratificação molecular / priorização computacional / nanomedicina de precisão potencial".

---

## 11. FASE 13–15 — ESTRATIFICAÇÃO, VALIDAÇÃO, SENSIBILIDADE

- **Estratificação (clustering/subtipos):** só se os dados suportarem (PCA/clustering/estabilidade/coerência). Evitar data leakage. Não criar subtipos artificiais.
- **Validação:** buscar dataset independente (outra coorte GEO, proteômica, evidência publicada). Separar descoberta × validação × interpretação. Sem validação → declarar explicitamente.
- **Sensibilidade:** variar FDR, cutoff, método DE, score STRING, métrica de centralidade, database, parâmetros GSEA. Resultado que desaparece com pequena mudança = **instável**.

---

## 12. FASE 16 — CRITÉRIOS DE INTERRUPÇÃO (APLICÁVEIS A CADA ETAPA)

Interromper/reformular quando: dados insuficientes; QC inadequado; batch não controlável; método incompatível; sem significância após correção; vias excessivamente redundantes; rede insuficiente; resultado dependente de decisão arbitrária única; circularidade; data leakage; sem evidência para conclusão; análise não responde à pergunta; etapa adicionada só por aparência de complexidade; conclusão exige extrapolação não sustentada.

Registrar sempre: **MOTIVO | IMPACTO | ALTERNATIVA | EFEITO SOBRE O OBJETIVO**.

---

## 13. FASE 17 — PROTOCOLO OPERACIONAL (TEMPLATE POR ANÁLISE)

Para **cada** análise, preencher antes de executar:
- PERGUNTA · HIPÓTESE · DADOS DE ENTRADA · MÉTODO · SOFTWARE · VERSÃO · PACOTES · VERSÕES · PARÂMETROS · SAÍDA ESPERADA · CRITÉRIO DE SUCESSO · CRITÉRIO DE INTERRUPÇÃO · RISCO DE VIÉS · ANÁLISE DE SENSIBILIDADE · FORMA DE INTERPRETAÇÃO.

---

## 14. FASE 18 — REGISTRO DE EXECUÇÃO

Registrar: data, dataset, versão, script, commit, software, pacotes, parâmetros, arquivos de entrada/saída, logs, erros, alterações, decisões. **Nunca sobrescrever resultados anteriores.** Estrutura: `data/ metadata/ scripts/ results/ figures/ tables/ logs/ documentation/`.

---

## 15. FASE 19 — HIERARQUIA DA EVIDÊNCIA

| Nível | Significado |
|---|---|
| N1 | observado diretamente nos dados |
| N2 | associação estatística |
| N3 | associação de rede |
| N4 | inferência funcional |
| N5 | hipótese nanotecnológica |
| N6 | hipótese translacional |

**N5/N6 nunca apresentados como evidência observacional.**

---

## 16. FASE 20–21 — CONTROLES ANTI-CHERRY-PICKING E ANTI-"NANOTECNOLOGIZAÇÃO"

- Critérios definidos **antes** de ver resultados; se definidos depois, declarar "análise exploratória".
- Preservar resultados negativos.
- Para cada elemento nanotecnológico, responder: (1) qual resultado o justifica? (2) qual característica molecular sustenta? (3) há literatura? (4) conexão direta ou especulativa? (5) necessária ou decorativa? (6) o trabalho permanece válido sem ela? — Se a resposta a (6) for SIM e não houver justificativa forte, **remover**.

---

## 17. ORDEM DE EXECUÇÃO PROPOSTA (NENHUMA ETAPA PULADA)

1. Obter matriz global → 2. QC completo → 3. pré-processamento → 4. DEG genoma-wide (ranking global) → 5. GSEA → 6. análise de redundância e seleção de vias → 7. PPI das vias selecionadas → 8. centralidade → 9. priorização multicritério → 10. caracterização funcional/estrutural (se justificada) → 11. sensibilidade/validação → 12. **só então** avaliação de plausibilidade nanotecnológica.

> **Bloqueio atual:** a ordem acima não pode começar sem a matriz global (Fase 1, lacuna L1).

---

*Protocolo fechado em 2026-09-06. Nenhuma via, alvo ou nanotecnologia foi pré-escolhida.*
