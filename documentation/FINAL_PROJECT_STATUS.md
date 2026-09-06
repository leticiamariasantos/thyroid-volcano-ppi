# FINAL PROJECT STATUS — thyroid-volcano-ppi

**FINAL SCIENTIFIC STATUS** — 2026-09-06

---

## PENDING → ACTION → RESULT → STATUS

| Pendência | Ação | Resultado | Status |
|---|---|---|---|
| Deconvolução formal | 4 refs scRNA-seq investigadas (GSE182416 tem smooth muscle/pericyte) | marker-based validado; ferramenta formal pendente | PARTIALLY_RESOLVED |
| Bulk RNA-seq independente | GSE224356 (3 pares) baixado e analisado | FN1/ITGA2/CCND1 replicados 3/3 | RESOLVED |
| Proteína | RPPA TCGA + correlação RNA-proteína | FN1 (r=0.53), ITGA2 (r=0.41) detectados e correlacionados | RESOLVED |
| Mutação | MAF cBioPortal | BRAF 57.4%, RAS 11.8%, TP53 0.4%; FN1/ITGA2 não mutados | RESOLVED |
| CNV | GISTIC cBioPortal | candidatos diploides (expressão não dirigida por CNV) | RESOLVED |
| IHC/metilação | avaliado, não extraído | — | PARTIALLY_RESOLVED |
| renv/Dockerfile | snapshot + Dockerfile v2 | R 4.6.1 congelado | RESOLVED |
| Batch TCGA×GTEx | VIF/rank/R² formal | confundimento estrutural total; mitigado por triangulação | RESOLVED (caracterizado) |

## TOP MOLECULAR CANDIDATES (Rank A — evidência molecular)

FN1 (0.689) · ITGA2 (0.647) · CCND1 (0.625) · LAMA2 (0.601) · RRM2 (0.598) · ITGA2B (0.568) · DCN (0.565) · FBLN1 (0.558) · COL27A1 (0.558) · CTSS (0.555)

## TOP TRANSLATIONAL CANDIDATES (Rank B — prioridade translacional)

ITGA2 (0.711) · FN1 (0.698) · CCND1 (0.686) · LAMA2 (0.650) · RRM2 (0.640) · ITGA2B (0.632) · CTSS (0.617) · FBLN1 (0.612) · COL27A1 (0.610) · SDC4 (0.597)

## TOP NANOMEDICINE HYPOTHESES (Rank C — hipótese nanomédica)

ITGA2 (0.788) · FN1 (0.773) · LAMA2 (0.719) · ITGA2B (0.709) · FBLN1 (0.684) · SDC4 (0.672) · CTSS (0.662) · CCND1 (0.655) · DCN (0.647) · ADH1B (0.641)

**Classificação nanomédica:** todos os top estão em **NANO-TIER B** (hipótese plausível, **não** validada).
**Nenhum candidato atinge NANO-TIER A** (evidência para investigação experimental direcionada) — porque
faltam binding/especificidade/internalização/eficácia e a origem celular (tumor vs estroma) não foi
resolvida por single-cell.

## Genes musculares (controle negativo validado)

MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM → **TIER 5 / NANO-TIER E** (artefato de composição). O modelo
corretamente os desprioriza, demonstrando que robustez estatística ≠ relevância tumoral.

## REMAINING LIMITATIONS

1. **Origem celular de FN1/ITGA2** (tumor vs estroma) — single-cell não extraído.
2. **Deconvolução formal** (MuSiC/BisqueRNA) não executada (referência GSE182416 identificada).
3. **RPPA tumor-only** — sem comparação proteica tumor vs normal.
4. **GSE224356 n=3** — replicação independente com poder limitado.
5. **IHC/metilação** não extraídos de dados primários.
6. **Batch TCGA×GTEx** estruturalmente confundido (irremediável no desenho discovery).
7. **Associação candidato × status mutacional** não testada (poder amostral).

## REQUIRED FUTURE EXPERIMENTS

1. **Single-cell** (GSE182416 + GSE232237) para definir a célula que expressa FN1/ITGA2.
2. **IHC** (HPA/coorte própria) para localização proteica tecidual (tumor vs estroma vs ECM).
3. **Proteômica tumor vs normal** (CPTAC/independente) para confirmação proteica diferencial.
4. **Validação funcional** (knockdown/CRISPR de ITGA2; modelo de ECM para FN1) antes de qualquer
   claim translacional.
5. **Associação a BRAF/RAS/TERT** em coorte maior.
6. **Deconvolução formal** com referência single-cell em ambiente com mais memória.

---

## CONCLUSÃO CIENTÍFICA

**O QUE FOI DEMONSTRADO:** forte diferença transcriptômica THCA×normal; 8.161 DEGs (limma TPM).

**O QUE FOI REPLICADO:** FN1, ITGA2, CCND1, CTSS, HLA-DPA1 — direção e significância em microarray
(GSE33630/GSE60542) **e** RNA-seq independente (GSE224356, 3/3 para FN1/ITGA2/CCND1).

**O QUE É ROBUSTO:** FN1/ITGA2/CCND1 (3 métodos + composição + replicação); vias p53/nucleotídeos/
replicação/imune/proteassoma.

**O QUE É HIPÓTESE:** a interpretação biológica do aumento de FN1/ITGA2 (microambiente vs tumor).

**O QUE É HIPÓTESE TRANSLACIONAL:** FN1 e ITGA2 como candidatos prioritários para investigação
experimental (N5/N6) — **não** alvos terapêuticos validados.

**O QUE É HIPÓTESE NANOMÉDICA:** ITGA2 (membrana, acessível) e FN1 (ECM) como hipóteses de
reconhecimento/targeting (NANO-TIER B) — **nenhuma** estratégia nanomédica validada.

**O QUE CONTINUA DESCONHECIDO:** célula de origem, localização proteica tecidual definitiva,
especificidade tumoral, e validade funcional.

---

*ITGA2 emerges as a prioritized candidate for experimental investigation of targeted delivery —*
*isto é uma priorização de hipótese, não uma declaração de eficácia terapêutica.*

---

## ATUALIZAÇÃO FINAL — Origem celular (single-cell, GSE182416)

**Resolução das duas principais incertezas:**

| Gene | Origem celular (tireoide normal) | Conclusão |
|---|---|---|
| FN1 | baixa e ampla (endotélio 0.94, mieloide 0.90, fibroblasto 0.51, epitelial 0.40) | **stromal/ECM, não tumor-intrínseco** → hipótese de microambiente |
| ITGA2 | essencialmente ausente (<4% células, mean <0.05 em todas) | **linha de base normal muito baixa** → suporte moderado p/ targeting; célula de origem tumoral não resolvida |
| MYH7/MYL1/MYL2/ACTA1/TNNT3/CKM | **completamente ausentes** (<0.5% células) | **confirma definitivamente artefato de composição (strap muscle no GTEx)** |

**Mudança conceitual (ranking recalibrado):**
- **FN1** deixa de ser tratado como candidato a targeting de célula tumoral → passa a
  **hipótese de reconhecimento/modulação do microambiente/ECM**.
- **ITGA2** permanece como candidato de **targeting de membrana** (linha de base normal baixa +
  superexpressão replicada + proteína), porém **especificidade tumoral não demonstrada**.
- Musculares: **descartados definitivamente** (ausentes na tireoide em nível single-cell).

## Resposta final (uma frase — conservadora, pós-auditoria de falsificação)

**ITGA2 permanece como candidato translacional prioritário (não alvo validado): superexpressão
tumoral consistente em 3 coortes independentes e 3 métodos, com detecção proteica por RPPA e baixa
expressão na tireoide normal em single-cell; localização celular tumoral, especificidade e evidência
funcional permanecem não demonstradas.**

**O primeiro experimento para falsificar essa hipótese é o perfil single-cell (GSE232237, PTC) ou
IHC multiplexada para determinar se ITGA2 está nas células tumorais ou apenas no estroma/vasos do
tumor.**

> Ver `documentation/FINAL_FALSIFICATION_AUDIT.md` para a tentativa formal de falsificação, as
> incertezas (localização tumoral, especificidade, evidência funcional) e a correção do viés de
> acessibilidade no score nanomédico.

---

## ATUALIZAÇÃO DEFINITIVA — Single-cell tumoral (GSE232237)

**Descoberta que altera a interpretação final:** análise de scRNA-seq tumoral (7 PTC + 5 ATC,
marker-based) mostrou que **ITGA2 é predominantemente expresso nas células epiteliais/tumorais**
do PTC (52,7% das células; 5,5× endotélio, 22× fibroblasto), e **FN1 também é fortemente tumoral**
(89,3% das células epiteliais), corrigindo a visão anterior de "FN1 exclusivamente estromal".

### Classificação final (atualizada)

| Candidato | Classificação |
|---|---|
| ITGA2 | **Candidato translacional prioritário / hipótese de targeting** — expresso em célula tumoral (single-cell), superexpresso e replicado, proteína presente; **especificidade tumoral e evidência funcional ainda não demonstradas** |
| FN1 | **Componente tumoral + estromal (EMT/microambiente)** — não exclusivamente ECM |
| CCND1 | Candidato molecular de proliferação (intracelular, não acessível a targeting) |

### Frase final conservadora (atualizada)

> "ITGA2 apresentou superexpressão tumoral consistente em três coortes independentes e por três
> métodos, com detecção proteica por RPPA, baixa expressão na tireoide normal e expressão
> predominantemente epitelial/tumoral em single-cell de PTC, sustentando sua priorização como
> candidato translacional para investigação de targeted delivery; a especificidade tumoral, a
> acessibilidade funcional e a evidência funcional permanecem não demonstradas, impedindo sua
> classificação como alvo terapêutico validado."
