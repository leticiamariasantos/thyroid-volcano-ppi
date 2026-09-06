# COMPOSITE RANKING METHODOLOGY — Ranking multidimensional de hipóteses

**Data:** 2026-09-06 · **Script:** `scripts/28_composite_ranking.R`

## 1. Objetivo

Priorizar **hipóteses para investigação experimental**, separando:
- **Rank A** — evidência molecular integrada;
- **Rank B** — prioridade translacional;
- **Rank C** — hipótese nanomédica (explicitamente não-validação).

## 2. Universo de candidatos (congelado, não-circular)

Union de: genes do PPI (contexto de rede), top 60 DEGs por |logFC|, candidatos prévios
(FN1/ITGA2/CTSS/HLA-DPA1/CCND1) e **controles negativos** (genes musculares). Nenhum candidato
foi selecionado usando validação externa ou literatura.

## 3. Dimensões (todas normalizadas 0–1; ausência ≠ negativo)

| Dimensão | Componentes | Fonte |
|---|---|---|
| A. Robustez transcriptômica | |logFC| TPM, −log10 FDR, concordância TPM/voom/DESeq2, consistência de direção | DEG TPM/voom/DESeq2 |
| B. Replicação independente | direção + FDR em GSE33630 e GSE60542 (1=replicado, 0.5=direção, 0=não) | GSE33630_deg, GSE60542_deg |
| C. Proteína | presença no RPPA TCGA (1=medido; NA=não medido → neutro) | cBioPortal RPPA |
| D. Composição | flag muscular/epitelial/endotelial/fibroblasto/imune | marcadores canônicos |
| E. Pathway context | nº de vias selecionadas (0–2) | pathways_selected |
| F. PPI | degree + betweenness normalizados | PPI_centrality |
| G. Especificidade | penaliza muscular/epitelial/endotelial/fibroblasto | marcadores |
| H. Acessibilidade | ECM/superfície=1, secretado=0.8, intracelular=0.3, desconhecido=0.5 | classificação funcional |

## 4. Scores

- **MOLECULAR** = 0.30 A + 0.20 B + 0.15 C + 0.15 E + 0.10 F + 0.10 |logFC|.
- **TRANSLATIONAL** = 0.25 A + 0.20 B + 0.15 C + 0.15 G + 0.10 E + 0.10 F + 0.05 |logFC|.
- **NANOMEDICINE** = 0.20 A + 0.20 B + 0.15 C + 0.15 H + 0.15 G + 0.10 F + 0.05 |logFC|.

## 5. Penalizações explícitas (mostradas, nunca ocultas)

- Muscular: −0.35.
- Replicação < 0.5: −0.15.
- Especificidade < 0.5: −0.10.

`FINAL = max(0, SCORE − PENALTY)`.

## 6. Sensibilidade de pesos

1000 combinações aleatórias de pesos (6 dimensões); `RANKING_STABILITY` = frequência no top-10.

## 7. Tiers

- **TIER 1–4** (evidência molecular decrescente); **TIER 5** = provável composicional.
- **NANO-TIER A–E** (A=não atribuído nesta execução; B=plausível; C=exploratória; D=insuficiente; E=não recomendável).

## 8. Limitações reconhecidas

1. Lista de marcadores musculares canônica **não é exaustiva** (MYH3, MYH11, MYL5/6/9 etc. não
   flaggados) → alguns genes musculares podem escapar da penalização.
2. Proteína: RPPA cobre ~200 proteínas; ausência no painel = **não medido**, não "negativo".
3. Acessibilidade/localização: estimada por categoria funcional, **não** por dado de localização
   subcelular primário (pendente).
4. Single-cell, IHC, metilação: **não disponíveis** (dimensões ausentes marcadas como NA).

## 9. Anti-circularidade

Nenhum componente usa informação derivada da própria seleção do candidato. O ranking é
reproduzível a partir dos dados congelados (`results/discovery_freeze/`).
