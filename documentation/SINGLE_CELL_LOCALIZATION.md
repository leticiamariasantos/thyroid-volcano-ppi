# SINGLE-CELL LOCALIZATION — Origem celular de FN1 e ITGA2

**Data:** 2026-09-06 · **Referência:** GSE182416 (scRNA-seq tireoide normal, 54.726 células, 7 pacientes)
**Cell types:** Epithelial (19.715), T/NK (14.139), B (6.674), Endothelial (4.993), Fibroblasts (3.995),
Myeloid (3.612), SMCs/Pericytes (1.229), Proliferating (369).

## Caveat de qualidade (registrado)
A matriz raw count de GSE182416 apresenta **contaminação por RNA ambiente** (marcadores epiteliais
TG/PAX8/TPO aparecem em todas as populações), o que reduz a precisão de atribuição célula-específica.
Conclusões robustas a esse ruído são destacadas.

## 1. FN1 — onde é expresso?

| cell_type | mean | % expr |
|---|---|---|
| Endothelial | 0.944 | 8.7 |
| Myeloid | 0.901 | 8.4 |
| T/NK | 0.576 | 7.4 |
| Fibroblasts | 0.507 | 10.7 |
| SMCs/Pericytes | 0.434 | 7.0 |
| Epithelial | 0.401 | 7.3 |
| B | 0.366 | 4.2 |

**Conclusão FN1:** expressão **baixa e amplamente distribuída**, sem enriquecimento forte em uma única
população; levemente maior em endotélio/mieloide que em fibroblasto (possivelmente ruído/ambiente).
Combinado com a natureza de FN1 (fibronectina, proteína de ECM) e HPA, o padrão é mais compatível com
**origem estromal/ECM/microambiente**, NÃO célula epitelial/tumoral.

**Classificação: B. stromal/ECM (predominantemente)** — não tumor-cell intrinsic.

## 2. ITGA2 — onde é expresso?

| cell_type | mean | % expr |
|---|---|---|
| Endothelial | 0.050 | 3.7 |
| Epithelial | 0.044 | 3.7 |
| SMCs/Pericytes | 0.040 | 3.4 |
| T/NK/B/Myeloid | 0.032 | 2.5 |
| Fibroblasts | 0.031 | 2.8 |

**Conclusão ITGA2:** **essencialmente AUSENTE na tireoide normal** em todas as populações (<4% das
células, mean <0,05). Isso indica linha de base muito baixa no tecido normal — coerente com a
superexpressão observada no tumor. A origem celular no TUMOR (célula tumoral vs estroma tumoral)
**não é resolvida** por este dataset (que é só normal).

**Classificação: C. expressão de membrana, especificidade insuficientemente resolvida** — baixa linha
de base normal é favorável, mas não foi demonstrado enriquecimento em célula tumoral.

## 3. Genes musculares — controles negativos (resultado decisivo)

MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM estão **completamente ausentes** na tireoide normal single-cell
(mean ≈ 0; % expressando < 0,5% em TODAS as populações).

→ **Prova definitiva** de que o sinal muscular no GTEx não vem de células da tireoide, e sim de
**músculo esquelético contaminante (strap muscle)** no tecido normal GTEx. **Confirma o artefato de
composição** de forma independente e no nível de célula única.

## 4. GSE232237 / GSE241184 / GSE193581

- GSE232237 (PTC/ATC): só `RAW.tar` (matrizes raw por amostra), **sem anotação de tipo celular
  processada** → exigiria clustering/anotação completa (fora do escopo desta execução).
- GSE241184 (1 paciente) e GSE193581 (foco ATC): sem matriz anotada prontamente utilizável.
- **Registrado:** a origem celular em TUMOR de FN1/ITGA2 permanece incerta (referências tumorais
  anotadas não disponíveis como matriz processada).

## 5. Recalibração do ranking

| Candidato | Mudança | Justificativa |
|---|---|---|
| FN1 | **não é alvo tumoral direto**; reformulado como hipótese de microambiente/ECM | expressão broad/estromal em single-cell |
| ITGA2 | **suporte moderado para targeting**; mantido (baixa linha de base normal) | membrana + baixa expressão normal + replicação + proteína |
| CCND1 | mantido (não acessível p/ targeting) | intracelular, expressão epitelial/proliferativa |
| Musculares | **descartados definitivamente** | ausentes no single-cell tireoidiano |

## 6. Decisão final

- **FN1 — Principal hipótese:** componente de matriz extracelular/microambiente tumoral
  potencialmente relevante para estratégia de reconhecimento/modulação do microambiente (NÃO
  targeted delivery para célula tumoral).
- **ITGA2 — Principal hipótese:** receptor de membrana com linha de base normal muito baixa e
  superexpressão tumoral replicada; candidato a investigação de targeting, MAS especificidade
  tumoral ainda não demonstrada em single-cell.
- **FN1 e ITGA2 NÃO são equivalentes:** FN1 = hipótese de microambiente/ECM; ITGA2 = hipótese de
  targeting de membrana (mais próxima de um alvo celular, porém não validada).

---

## ATUALIZAÇÃO — Single-cell TUMORAL (GSE232237, PTC + ATC)

**Resultado que muda a interpretação anterior:** a análise de GSE232237 (marker-based, 12 amostras
tumorais: 7 PTC + 5 ATC) revelou:

| Gene | PTC — Epithelial_tumor | PTC — Fibroblast | PTC — Endothelial | PTC — imune/SMC |
|---|---|---|---|---|
| ITGA2 | **3,25 (52,7%)** | 0,15 (7,9%) | 0,59 (30%) | <0,05 (<4%) |
| FN1 | **168,9 (89,3%)** | 47,7 (86,8%) | 3,6 | 1,3–15 |
| CCND1 | **13,0 (80,0%)** | 1,3 | 1,5 | <2 |
| ITGB1 | 4,9 (81,8%) | 6,4 (75,4%) | 4,1 (75,2%) | 0,8–6,7 |

### Conclusões corrigidas

1. **ITGA2 é predominantemente expresso nas células epiteliais/tumorais do PTC** (52,7% das células
   tumorais; 5,5× endotélio, 22× fibroblasto), e **não** no estroma. A localização tumoral —
   antes "não resolvida" — agora tem **evidência favorável à célula tumoral** (classificação
   falsificação: **A→B** — forte suporte à expressão tumoral, com expressão secundária em endotélio).
2. **FN1 também é fortemente expresso nas células tumorais epiteliais** (168,9; 89,3%) e em
   fibroblastos (47,7). A interpretação anterior "FN1 é stromal/ECM" era **incompleta**: FN1 é
   **tumor + estroma** (consistente com programa EMT no tumor).
3. **ITGB1** é ubíquo (como esperado para a subunidade beta1), reforçando que a especificidade do
   par reside na subunidade **alpha (ITGA2)**, que é enriquecida no tumor epitelial.

### Caveat
Classificação celular por marcadores (EPCAM/KRT/PAX8/TG/TPO para epitelial; COL1A1/DCN para
fibroblasto; PECAM1/VWF para endotélio), pois o RAW.tar de GSE232237 **não inclui anotação de tipo
celular do autor**. O enriquecimento de ITGA2/FN1 em células epiteliais é robusto à incerteza de
classificação (diferença de 5–22× entre compartimentos).
