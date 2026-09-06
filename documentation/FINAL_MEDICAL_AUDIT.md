# AUDITORIA MÉDICA FINAL — Revisão multidisciplinar

**Data:** 2026-09-06

Revisão sob as perspectivas de endocrinologista, oncologista, patologista, bioinformata e
pesquisador translacional.

| # | Pergunta | Resposta | Evidência |
|---|---|---|---|
| 1 | Demonstra alteração molecular? | **SIM** | 8.161 DEGs; FN1/ITGA2 robustos em 3 métodos |
| 2 | O sinal pode ser composição? | **SIM, parcialmente** (músculo) | muscle score −3,71; genes musculares colapsam e não replicam em coortes pareadas |
| 3 | O resultado replica? | **SIM** (FN1/ITGA2/CCND1/CTSS/HLA-DPA1) | GSE33630 + GSE60542 + GSE224356 (RNA-seq 3/3) |
| 4 | Existe proteína? | **SIM** (FN1, ITGA2, CCND1) | TCGA RPPA; RNA-proteína r=0.53/0.41 |
| 5 | Onde está a proteína? | **Parcialmente definido** | FN1=ECM/estroma (provável); ITGA2=membrana (provável); single-cell pendente |
| 6 | Qual célula expressa? | **NÃO resolvido** | single-cell não extraído |
| 7 | Coerência multiômica? | **SIM** | RNA↑ + proteína↑ (correlacionada); sem CNV/mutação dirigindo → transcricional |
| 8 | Plausibilidade funcional? | **Plausível** | FN1/ITGA2 = ECM/adesão, coerente com remodelamento de microambiente |
| 9 | Acessibilidade? | **Parcial** | ITGA2 membrana (acessível); FN1 extracelular (acessível ao compartimento); CCND1 intracelular (não) |
| 10 | Plausibilidade para targeting? | **Hipotética** | sem binding/especificidade/internalização demonstrados |
| 11 | Evidência suficiente para nanomedicina? | **NÃO** | NANO-TIER B (plausível, não validada); nenhum NANO-TIER A |

## Veredito médico

O estudo demonstra **de forma robusta e replicada** uma assinatura transcriptômica de
THCA com componentes tumor-intrínsecos (FN1/ITGA2/CCND1) e um artefato de composição
(músculo) corretamente identificado e excluído.

**Nenhuma afirmação terapêutica, prognóstica ou nanomédica é suportada.** Os candidatos
FN1 e ITGA2 são **hipóteses moleculares prioritárias para investigação experimental futura**,
não biomarcadores nem alvos validados.

## Teste de honestidade científica (revisor hostil)

**Principal crítica potencial:** "O estudo compara TCGA (tumor) com GTEx (normal) — batch
totalmente confundido com condição; qualquer DEG pode ser artefato técnico; e FN1/ITGA2 podem
ser apenas sinal de estroma."

**A crítica foi resolvida?** **Parcialmente, e de forma explícita:**
- O confundimento foi **quantificado formalmente** (source≡condition, VIF=∞, R²=0,24) — não negado.
- A triangulação em **3 coortes independentes com desenho controlado** (GSE33630, GSE60542
  microarray pareado; GSE224356 RNA-seq pareado) **replica** FN1/ITGA2/CCND1, escapando do
  confundimento.
- O sinal muscular (que era o principal artefato de composição) **não replica** nas coortes
  pareadas — confirmando que a composição foi corretamente identificada.
- **Limitação remanescente:** a origem celular (tumor vs estroma) de FN1/ITGA2 permanece não
  resolvida por single-cell; esta é a principal razão para **não** elevar FN1/ITGA2 a alvo.
