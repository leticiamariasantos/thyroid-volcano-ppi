# Painel de 30 vias — Fundamentação da Seleção (a priori)

**Data da definição do painel:** 2026-09-06 (a mesma da nova execução analítica)
**Natureza da seleção:** pré-especificada ANTES da interpretação dos resultados desta execução.

## Regra de reprodutibilidade

As 20 vias adicionais foram selecionadas **antes** de qualquer consulta aos resultados de
expressão diferencial, GSEA ou enriquecimento desta nova fase. Nenhuma via foi incluída por
ter se mostrado significativa nesta execução. A seleção permaneceria defensável mesmo se os
genes ITGA2, FN1 ou CCND1 não fossem diferencialmente expressos.

## Critérios objetivos de seleção

1. **Relevância para carcinoma de tireoide** — vias com papel documentado em biologia
   tireoidiana ou oncogênese tireoidiana.
2. **Relevância para biologia tumoral** — vias canônicas de progressão, microambiente,
   metabolismo e resposta a estresse.
3. **Complementaridade às 10 vias originais** — ampliar o espaço biológico, não duplicar
   (MAPK/PI3K/mTOR/p53/apoptose/ciclo/Wnt/NF-kB já cobertos).
4. **Mecanismos celulares distintos** — cada via deve representar um processo celular
   distinguível (adesão, imunidade, metabolismo, reparo, etc.).
5. **Qualidade e estabilidade da definição** — vias bem anotadas e estáveis.
6. **Disponibilidade confiável** — presentes no KEGG ou Reactome com IDs resolvíveis.
7. **Baixa redundância excessiva** — vias relacionadas são mantidas apenas quando
   representam facetas funcionalmente distintas (sinalização vs. estrutura vs. remodelamento).

## As 10 vias originais (INALTERADAS)

As 10 vias originais do projeto foram preservadas sem alteração de composição:

`hsa05216` Thyroid cancer · `hsa04919` Thyroid hormone signaling · `hsa04010` MAPK ·
`hsa04151` PI3K-Akt · `hsa04150` mTOR · `hsa04115` p53 · `hsa04210` Apoptosis ·
`hsa04110` Cell cycle · `hsa04310` Wnt · `hsa04064` NF-kappa B.

## Justificativa individual das 20 vias adicionais

### Bloco A — Interação célula-matriz / adesão / ECM

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04510` | Focal adhesion | Plataforma de sinalização por integrinas que converte contato com a matriz em sinais de sobrevivência, migração e invasão. Complementa as vias de crescimento (MAPK/PI3K) cobrindo o *input* mecânico/adesivo ausente do painel original. |
| `hsa04512` | ECM-receptor interaction | Define os receptores estruturais (integrinas, colágeno, laminina) que medeiam a interação tumor–estroma. Distinta de Focal adhesion por representar a face de reconhecimento de ligantes. |
| `R-HSA-1474244` | Extracellular matrix organization | Cobre a biossíntese, montagem e degradação de colágeno e matriz, o componente *estrutural* que as duas vias acima apenas sinalizam. |

### Bloco B — EMT / TGF-beta / Hippo

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04350` | TGF-beta signaling | Indutor canônico de transição epitélio-mesenquimal e imunossupressão no microambiente tumoral; mecanismo de plasticidade ausente do painel original. |
| `hsa04390` | Hippo signaling | Eixo YAP/TAZ de controle de tamanho de órgão e supressão tumoral; complementa Wnt com um regulador de crescimento mecanicamente distinto. |

### Bloco C — JAK-STAT / TNF / inflamação

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04630` | JAK-STAT signaling | Transdução de citocinas com papel em inflamação e imunidade antitumoral; via de sinalização própria, distinta de MAPK/PI3K. |
| `hsa04668` | TNF signaling | Inflamação crônica, sobrevivência e morte mediadas por TNF; relacionada mas distinta de NF-kB original (foco no ligante/receptor). |
| `hsa04657` | IL-17 signaling | Eixo pró-inflamatório Th17; amplia o espaço inflamatório além de TNF/NF-kB. |

### Bloco D — Imunidade / interferon / apresentação de antígeno

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04612` | Antigen processing and presentation | Mecanismo MHC/processamento que conecta imunidade adaptativa ao reconhecimento tumoral; não representado no painel original. |
| `R-HSA-913531` | Interferon Signaling | Programa imunológico de interferon tipo I, modulador central do microambiente imune tumoral. |
| `hsa04620` | Toll-like receptor signaling | Sensoriamento imune inato de PAMPs; complementa interferon com o braço de detecção inata. |

### Bloco E — Morte celular regulada / senescência / degradação

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04216` | Ferroptosis | Morte celular dependente de ferro, distinta de apoptosis; vulnerabilidade metabólica emergente em câncer. |
| `hsa04218` | Cellular senescence | Parada proliferativa irreversível ligada a senescência tumoral e resposta a terapia. |
| `hsa03050` | Proteasome | Degradação proteica ubiquitina-proteassoma; regula turnover de oncoproteínas (mecanismo distinto de morte). |

### Bloco F — Metabolismo / OXPHOS / estresse oxidativo

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa00190` | Oxidative phosphorylation | Metabolismo mitocondrial e reprogramação energética tumoral (ausente do painel original, dominado por sinalização). |
| `R-HSA-3299685` | Detoxification of Reactive Oxygen Species | Homeostase redox e defesa antioxidante; distinta de OXPHOS por representar a resposta ao estresse oxidativo. |

### Bloco G — Replicação / reparo de DNA

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa03440` | Homologous recombination | Reparo de quebras de fita dupla e instabilidade genômica. |
| `hsa03030` | DNA replication | Maquinaria replicativa e estresse replicativo (proliferação). |

### Bloco H — Hipóxia / autofagia

| ID | Via | Justificativa |
|----|-----|---------------|
| `hsa04066` | HIF-1 signaling | Resposta a hipóxia e angiogênese; adaptação tumoral ao microambiente hipóxico (distinta de OXPHOS). |
| `hsa04140` | Autophagy - animal | Reciclagem celular e sobrevivência sob estresse metabólico; mecanismo distinto de morte celular. |

## Nota sobre não-circularidade

As vias não foram escolhidas por conterem ITGA2, FN1 ou CCND1. Ainda que Focal adhesion,
ECM-receptor e ECM organization contenham esses genes, sua inclusão se justifica por
representarem programas biológicos de adesão/matriz com relevância independente para a
biologia tumoral. A avaliação de convergência desses candidatos (seção dedicada) é feita
**após** o pipeline, sem usar a composição do painel como evidência.
