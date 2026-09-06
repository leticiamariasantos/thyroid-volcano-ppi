# DECONVOLUTION FINAL ASSESSMENT — Avaliação final da deconvolução

**Data:** 2026-09-06 · **Status:** PARTIALLY_RESOLVED

## 1. Referências single-cell investigadas

| Accession | Tipo | Células | Populações | Adequação |
|---|---|---|---|---|
| GSE182416 | scRNA-seq tireoide normal | 54.762 (7 pacientes) | epitelial, fibroblasto, endotelial, **smooth muscle/pericyte**, imune | **apropriada (tem músculo liso/pericito)** |
| GSE232237 | scRNA-seq PTC/ATC/normal | 3 normal + 7 PTC + 5 ATC | epitelial, imune, estromal | apropriada p/ tumor |
| GSE241184 | scRNA-seq PTC + metástase | 3 amostras (1 paciente) | epitelial, imune, estromal | limitada (n=1) |
| GSE193581 | scRNA-seq ATC/DTC | múltiplos | tumor, imune, estromal | foco em ATC |

**Conclusão:** GSE182416 é a referência mais adequada (inclui smooth muscle/pericytes, permitindo
modelar o sinal muscular). GSE232237 complementa com células tumorais.

## 2. Ferramentas de deconvolução

| Ferramenta | Disponibilidade | Viabilidade | Veredito |
|---|---|---|---|
| MuSiC | GitHub (não instalado) | exige referência scRNA-seq + bulk counts | não executado (dependência) |
| BisqueRNA | CRAN (não instalado) | exige ExpressionSet scRNA-seq | não executado |
| CIBERSORTx | Docker | infraestrutura ausente | não executado |
| Marker-based (MCPcounter-like) | implementado | executa sobre bulk TPM | **executado** |

## 3. Resultado alcançado (marker-based, já validado)

- muscle_skeletal: efeito −3,71 (Cohen's d −3,75, FDR 4e-113) → **composicional**.
- thyroid_epithelial/fibroblast/endothelial: enriquecidos no normal.
- immune: sem diferença.
- **Validação independente da hipótese composicional:** em GSE224356 (pares pareados, sem strap
  muscle) e GSE33630/GSE60542, os genes musculares NÃO são DE → confirma composição.

## 4. Pergunta principal respondida

"Até que ponto a diferença TCGA×GTEx é explicada por composição?"
→ O **componente muscular** (MYH7/MYL1/MYL2/ACTA1/TNNT3/CKM) é quase inteiramente composicional
(colapsa/inverte ao ajustar por score muscular; ausente em coortes pareadas). O **sinal
tumor-intrínseco** (FN1/ITGA2/CCND1/CTSS/HLA-DPA1) **não** é explicado por composição (robusto ao
ajuste).

## 5. Estado final

- **PARTIALLY_RESOLVED** — a hipótese composicional foi respondida empiricamente por
  marker-based scoring + covariável + triangulação, mas a deconvolução formal (MuSiC/BisqueRNA)
  com a referência single-cell GSE182416 **não foi executada** (limitação de dependências/memória).
- **Próximo experimento (recomendado):** rodar MuSiC/BisqueRNA com GSE182416 + GSE232237 como
  referência em ambiente com mais memória.

*Arquivos: `results/composition/composition_summary.tsv`, `formal_deconvolution.tsv`,
`composition_covariate_*.tsv`, `results/rnaseq_validation/RNA_SEQ_VALIDATION.tsv`.*
