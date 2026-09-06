# PHASE 2 (reboot) — Resultados

**Data:** 2026-09-06

## 1. Expressão diferencial

| Método | Testados | DEGs | Up | Down |
|---|---|---|---|---|
| limma (log2 TPM) | 26.011 | 12.200 | 2.485 | 9.715 |
| limma-voom (counts) | 39.702 | 19.485 | 7.235 | 12.250 |
| DESeq2 (counts) | 39.702 | 17.737 | 5.572 | 12.165 |

**Concordância entre métodos:** Spearman(logFC) 0,88–0,96; concordância de direção
0,84–0,98; Jaccard de DEGs 0,59–0,80 (voom↔DESeq2 mais concordantes, como esperado por
compartilharem a mesma matriz de contagens).

## 2. Painel de 30 vias (GSEA, limma — método principal)

Vias significativas (FDR<0,05, |NES| decrescente), todas **Up no tumor**:

| Via | Grupo | NES | padj |
|---|---|---|---|
| Proteasome (hsa03050) | ADDITIONAL_20 | +2,82 | 1,8e-7 |
| Antigen processing and presentation (hsa04612) | ADDITIONAL_20 | +2,44 | 3,2e-6 |
| Oxidative phosphorylation (hsa00190) | ADDITIONAL_20 | +2,18 | 3,3e-7 |
| DNA replication (hsa03030) | ADDITIONAL_20 | +2,01 | 4,8e-3 |
| p53 signaling (hsa04115) | ORIGINAL_10 | +1,89 | 8,4e-5 |
| Cell cycle (hsa04110) | ORIGINAL_10 | +1,75 | 3,2e-6 |

No painel original de 10 vias, **p53 e Cell cycle** foram as vias a priori mais robustas.
As demais vias do painel (incluindo MAPK, PI3K-Akt, mTOR, Wnt, NF-kB, Thyroid cancer,
Thyroid hormone signaling, Apoptosis) não atingiram FDR<0,05 no limma.

**Cross-método:** várias vias adicionais tornam-se significativas em voom/DESeq2
(Interferon Signaling, Focal adhesion, Hippo, Cellular senescence, Ferroptosis, Autophagy,
Detoxification of ROS, Apoptosis, mTOR, Thyroid cancer) — reportado como consistência
parcial, não como confirmação.

## 3. GSEA global (sem restrição ao painel)

Top programas emergentes (NES, padj):

- **KEGG:** Ribosome (+3,94), Proteasome (+2,77), Graft-versus-host (+2,81),
  Antigen processing (+2,40), Cell cycle (+1,92), DNA replication (+2,01),
  OXPHOS (+1,87), p53 (+1,97); **Dilated cardiomyopathy (−1,63)** (músculo, normal).
- **Reactome:** tradução/ribossomo (topo, +4,05), seguidos de resposta imune.
- **Hallmark:** Myogenesis (−1,65, Down no tumor), Protein secretion (+1,62),
  Coagulation (+1,55).

Interpretação: emergem **proliferação/tradução**, **imunidade/inflamação**, **proteostase**
e **metabolismo mitocondrial**; a assinatura muscular emerge **Down no tumor** (tecido normal).

## 4. Redundância entre as 30 vias

Pares mais redundantes (Jaccard): Focal adhesion↔ECM-receptor (0,29), TNF↔IL-17 (0,26),
PI3K-Akt↔Focal adhesion (0,26), TNF↔TLR (0,25), MAPK↔PI3K-Akt (0,19). Módulos coerentes:
(i) adesão/ECM; (ii) inflamação/imunidade; (iii) crescimento/sinalização.

## 5. Composição celular (músculo)

Score muscular: **tumor −2,91 vs normal +0,03** — assinatura muscular concentrada no tecido
normal GTEx (músculo esquelético adjacente), não é regulação tumoral intrínseca.

**Sensibilidade (remoção de 48 marcadores musculares):** as 6 vias significativas do painel
**permanecem robustas** (delta_NES < 0,06 em módulo). Nenhuma via significativa desapareceu.
MYH7 não replica como Down nas coortes pareadas externas (GSE33630/GSE60542), corroborando
que o sinal muscular é artefato de composição do GTEx.

## 6. PPI

Rede STRING (escore≥700) reconstruída excluindo os 48 marcadores de músculo estriado
(artefato de composição do tecido normal) e equilibrando 250 genes Up + 250 Down: **70 nós / 65
arestas**, 18 comunidades (modularidade 0,794). **FN1** é o hub de maior grau (9), seguido por
CDT1, CDH2, SDC1 (grau 5) e CCND1, TK1, PCLAF, CENPM, MET, FGF17 (grau 4) — programa de
proliferação (ciclo celular) e adesão/EMT. **A centralidade é reportada como propriedade
topológica, não como alvo terapêutico.** PRKCA (hub da Fase 1 direcionada à via hsa04919) não
é hub na análise genoma-wide.

## 7. Candidatos ITGA2 / FN1 / CCND1 (convergência independente)

| Gene | logFC (limma) | FDR | rank | direção | consistente (3 métodos) | vias (membro) | core enrichment | compartimento (sc) |
|---|---|---|---|---|---|---|---|---|
| CCND1 | +2,26 | 2,7e-144 | 949 | Up | sim | 10 | **4 vias** (Cell cycle, p53, Senescence, Thyroid cancer) | Epithelial_tumor |
| FN1 | +4,32 | 1,8e-95 | 2926 | Up | sim | 4 (ECM/adesão) | 0 | Epithelial_tumor (+ Fibroblast) |
| ITGA2 | +2,47 | 4,1e-89 | 3325 | Up | sim | 4 (ECM/adesão) | 0 | Epithelial_tumor |

**Convergência do painel (leading edge recorrente):** 18 genes figuram no core enrichment
(leading edge) de ≥2 vias robustas. O painel de 30 vias converge para o **eixo ciclo
celular/p53 (vias ORIGINAIS)** — TP53 (7 vias), CCND1, CDK1, CDKN1A, CCNB1, CCNB2, CCND2,
CDK4, CDKN2A, CHEK1, MDM2, SFN — e para a **maquinaria de replicação** (MCM2–MCM5, PCNA) e
**PSME1** (elo proteassoma–apresentação de antígeno). As 20 vias adicionais convergem para
proteassoma/apresentação de antígeno/OXPHOS/replicação de DNA. `validation/panel_convergence_genes.tsv`.

**Candidatos pré-especificados:** CCND1 figura no core enrichment apenas de vias ORIGINAIS
(Cell cycle, p53) e de Cellular senescence/Thyroid cancer; não figura no core enrichment das
20 vias adicionais. ITGA2 e FN1 **não** figuram no core enrichment — emergem como DEGs
individuais robustos e replicados dentro dos programas de adesão/ECM (Focal adhesion,
ECM-receptor, ECM organization), não como componentes centrais dos programas mais enriquecidos.

**Validação externa (direção concordante):** ITGA2, FN1 e CCND1 **Up** em GSE33630
(logFC +2,00 / +2,88 / +1,42), GSE60542 (+2,50 / +3,06 / +1,07) e GSE224356 (T1N1, T2N2, T3N3).

**Proteína (RPPA):** CCND1 (Cyclin D1) no painel — correlação RNA-proteína fraca
(Spearman 0,051, n=366), coerente com regulação pós-transcricional da ciclina D1.
ITGA2 e FN1 fora do painel RPPA (limitação).

**Mutação/CNV (THCA):** BRAF 57% (V600E), NRAS 7,8%, HRAS 3,2%, KRAS 0,8%, TP53 0,4%,
TERT 0,6% (promotor não captado por exoma). **ITGA2 e CCND1 sem mutações**; FN1 0,4%.
"Não foi identificada alteração genômica suficiente para explicar a expressão observada"
de ITGA2/FN1/CCND1 — a superexpressão é provavelmente regulada em nível transcricional/epigenético.

## 8. Matriz de robustez

**6 vias ROBUSTAS (escore 4/4):** Oxidative phosphorylation, DNA replication, Proteasome,
Cell cycle, p53 signaling, Antigen processing and presentation — significativas no painel
e no global, consistentes entre métodos e persistentes após remoção muscular.
`pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv`.
