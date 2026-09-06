# EXTERNAL VALIDATION REPORT — Validação externa independente (GEO)

**Data:** 2026-09-06 · **R:** 4.6.1
**Discovery:** TCGA-THCA (n=504) vs GTEx normal thyroid (n=279), `log2(TPM+0.001)`; complemento por counts (recount3, voom/DESeq2).

---

## 1. Objetivo

Avaliar, de forma independente e sem data leakage, se os candidatos e vias do discovery
sobrevivem à validação em coortes externas de microarray (PTC vs normal), e testar a hipótese
de que o sinal muscular era artefato de composição tecidual.

## 2. Discovery congelado

Ver `documentation/DISCOVERY_FREEZE.md`. Candidatos congelados: **FN1, ITGA2, CTSS, HLA-DPA1, CCND1**
(tumor-intrínsecos) e **MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM** (classificados composicionais).

## 3. Deconvolução

- **Deconvolução formal (CIBERSORTx/EPIC/MuSiC/BisqueRNA/MCPcounter/xCell) NÃO aplicável:**
  exigem referência single-cell (MuSiC/BisqueRNA) ou Docker (CIBERSORTx) ou não possuem assinatura
  de músculo esquelético/epitélio tireoidiano (EPIC/MCPcounter/xCell) — a hipótese central é
  justamente contaminação muscular. **Registrada a impossibilidade, sem simulação.**
- **Método usado:** marker-based scoring (MCPcounter-like) — média de expressão log2(TPM+0.001)
  de marcadores canônicos. Tratado como **estimativa/modelo**, não observação experimental.
- **Resultado (discovery):** muscle_skeletal enriquecido no **normal GTEx** (Cohen's d = −3,75,
  FDR = 4e-113); thyroid_epithelial, fibroblast_ecm e endothelial também enriquecidos no normal;
  immune sem diferença (p = 0,34). Ver `results/composition/composition_summary.tsv`.

## 4. GSE33630 (validação primária)

- **Proveniência:** GEO, Affymetrix HG-U133 Plus 2.0 (GPL570); 105 amostras; SHA256 em
  `results/validation/provenance.tsv`.
- **Amostras:** 49 PTC, 45 normal (patient-matched), 11 ATC (**excluídos** do contraste primário).
- **Contraste:** PTC − Normal (limma, trend=TRUE).
- **Normalização:** valores do series matrix já são RMA-normalizados (log2); não foi re-aplicada RMA.
- **Mapeamento de probes:** probe de **maior expressão média** por gene (regra pré-registrada e
  uniforme, via `hgu133plus2.db`).

## 5. GSE60542 (validação secundária/robustez)

- **Amostras:** 33 PTC, 30 normal thyroid (e metástases/linfonodos — não misturados).
- **Contraste primário:** PTC − Normal (mesmo pipeline/mesma regra de mapeamento do GSE33630).
- Usado também como **stress test de composição** (normal pareado, sem strap muscle).

## 6. QC

- Distribuições/padrão de microarray confirmados (valores log2). Sem exclusão de amostras no
  contraste primário (PTC/Normal). ATC e amostras não-tireoidianas foram apenas **excluídas do
  contraste**, não removidas do dataset.

## 7. Normalização

- Affymetrix: series matrix já RMA-normalizado. **Não** foi usado TPM nem tratado microarray como RNA-seq.

## 8. Gene mapping

- `hgu133plus2.db`; regra **pré-registrada**: probe com maior expressão média por gene, aplicada
  a todos os genes (sem seleção pós-hoc). `select()` gerou aviso de mapeamento 1:many (esperado,
  resolvido pela regra de colapso).

## 9. Validação de genes (candidatos congelados)

| Gene | Discovery | GSE33630 (FDR) | GSE60542 (FDR) | Status |
|---|---|---|---|---|
| FN1 | +4,32 | +2,91 (2e-25) | +3,06 (1e-14) | INDEPENDENTLY REPLICATED |
| ITGA2 | +2,47 | +2,00 (2e-16) | +2,50 (2e-11) | INDEPENDENTLY REPLICATED |
| CTSS | +2,42 | +1,38 (5e-7) | +1,12 (3e-4) | INDEPENDENTLY REPLICATED |
| HLA-DPA1 | +1,22 | +1,06 (1e-4) | +0,81 (2e-3) | INDEPENDENTLY REPLICATED |
| CCND1 | +2,26 | +1,42 (1e-21) | +1,07 (7e-11) | INDEPENDENTLY REPLICATED |
| MYH7 | −7,01 | −0,01 (0,86) | +0,12 (0,08) | NOT REPLICATED (composicional) |
| MYL1 | −8,30 | +0,05 (0,79) | +0,11 (0,46) | NOT REPLICATED |
| MYL2 | −8,69 | −0,11 (0,51) | +0,05 (0,79) | NOT REPLICATED |
| ACTA1 | −7,17 | −0,12 (0,64) | +0,11 (0,49) | NOT REPLICATED |
| TNNT3 | −4,77 | −0,07 (0,27) | −0,07 (0,38) | NOT REPLICATED |
| CKM | −6,92 | +0,07 (0,70) | +0,15 (0,13) | NOT REPLICATED |

Concordância global de direção (4.584 genes discovery-DEG): GSE33630 0,67; GSE60542 0,72.

## 10. Validação de vias

| Pathway | Discovery NES | GSE33630 NES | GSE60542 NES | Conclusão |
|---|---|---|---|---|
| hsa04115 p53 | 1,87 | 2,24** | 2,16** | replicada (ambas) |
| hsa04612 Antigen presentation | 2,58 | 1,18 | 0,99 | direção consistente (mais fraca) |
| hsa03030 DNA replication | 2,00 | 0,65 | 1,89** | replicada (GSE60542) |
| hsa01232 Nucleotide metabolism | 1,44 | 1,68** | 1,69** | replicada (ambas) |
| hsa04820 Cytoskeleton in muscle | −1,47 | **+1,73** | **+1,44** | **direção invertida → composicional** |

`**` = padj < 0,05.

## 11. Composição (stress test)

| Dataset | muscle PTC | muscle Normal | diff |
|---|---|---|---|
| Discovery (TCGA/GTEx) | −2,19 | +1,52 | **−3,71** |
| GSE33630 | 5,43 | 5,52 | −0,10 |
| GSE60542 | 5,13 | 5,24 | −0,11 |

O sinal muscular de −3,71 no discovery **desaparece completamente** nos GEO (diff ≈ −0,1),
porque o normal pareado não contém o músculo esquelético adjacente do GTEx. **Confirmação
definitiva do artefato composicional.**

## 12. Concordâncias

FN1, ITGA2, CTSS, HLA-DPA1, CCND1; vias p53, metabolismo de nucleotídeos, replicação de DNA.

## 13. Discordâncias

Genes musculares (MYH7/MYL1/MYL2/ACTA1/TNNT3/CKM) e via "Cytoskeleton in muscle cells" —
não replicam e a via inverte a direção → composicionais.

## 14. FN1

Replicado em 2 coortes independentes (direção e significância). Permanece **hipótese molecular
(N5/N6)** — superexpresso no tumor e no estroma, hub de rede; **não** é alvo terapêutico validado.

## 15. ITGA2

Replicado em 2 coortes independentes. Permanece **hipótese molecular (N5/N6)** — receptor de
superfície/adesão; **não** é alvo terapêutico validado.

## 16. PPI

Núcleo estável: FN1 e ITGA2 são os **top-2 hubs por grau em todas as 4 variantes de rede**
(original, count-based, no-composition, robust candidates). Correlação de grau 0,88–0,97;
Jaccard de nós 0,63–0,87. Hubs musculares (MYL1/MYH7/TPM1) **desaparecem** ao remover composição.

## 17. Limitações

1. Cross-plataforma (RNA-seq discovery × microarray GEO) — magnitudes não comparáveis diretamente.
2. GSE33630/GSE60542 têm amostras menores (n≈45–49 e 30–33 por grupo).
3. Deconvolução formal não aplicável (sem referência de músculo); usado marker-based scoring.
4. Mapeamento de probes 1:many resolvido por regra pré-registrada (maior média).
5. MYL1 ausente na anotação G026 dos counts (presente apenas no TPM).
6. 1 amostra GTEx sem correspondência no recount3.

## 18. Implicações translacionais

FN1 e ITGA2 (superexpressos, replicados, de superfície/ECM) sustentam **hipótese exploratória**
de reconhecimento/targeting para investigação futura. Genes musculares **descartados**.

## 19. Nanotecnologia

**Nenhuma evidência suficiente** para propor estratégia nanomédica específica. A replicação de
FN1/ITGA2 **não autoriza** chamá-los de alvos terapêuticos, biomarcadores clínicos ou alvos de
drug delivery demonstrados. Máximo permitido: **hipótese molecular/translacional para investigação
futura** (N5/N6).

## 20. Conclusão

A validação externa independente **corrobora** os achados tumor-intrínsecos (FN1, ITGA2, CTSS,
HLA-DPA1, CCND1; vias p53/nucleotídeos/replicação) e **confirma** que o sinal muscular era
artefato de composição tecidual — um resultado negativo válido e importante, que reorienta a
priorização para candidatos genuínos.

---

*Arquivos: `results/validation/*`, `10_validation/*`, `results/discovery_freeze/*`.*
