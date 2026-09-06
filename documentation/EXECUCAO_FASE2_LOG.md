# REGISTRO DE EXECUÇÃO — FASE 2 (AQUISIÇÃO DA MATRIZ GLOBAL)

**Data:** 2026-09-06
**Script:** `scripts/01_acquire_global_data.R`
**R:** 4.6.1 · **data.table:** 1.18.4

## 1. O que foi feito (verificado, não assumido)

| Item | Resultado |
|---|---|
| Matriz global | `TcgaTargetGtex_rsem_gene_tpm.gz` (TOIL recompute) — 1.323.254.426 bytes, download íntegro |
| Dimensões brutas | 60.498 genes × 19.131 amostras |
| Amostras Fase 1 | 783/783 presentes na matriz global (nenhuma ausente) |
| Mapeamento Ensembl→símbolo | 100% (60.498/60.498) via `gencode.v23.annotation.gene.probemap` |
| Genes após colapso de símbolos duplicados | **58.581** (1.917 removidos; regra: maior média de expressão) |
| Matriz final | `data/global/TCGA_GTEx_thyroid_tpm.tsv` — 58.581 genes × 783 amostras (352 MB) |

## 2. UNIDADE/ESCALA ESTABELECIDA

- min = **−9,966** (= log2(0,001)), max = 17,543, Q25/mediana = −9,966, Q99 = 7,547.
- 99,1% dos valores são não-inteiros; **0% zeros**.
- **Conclusão: os valores já estão em escala log2, com pseudocontagem 0,001 → `log2(TPM + 0,001)`.** Não é count bruto (não é inteiro) e não é log2(count+1) (permite negativos).

## 3. ACHADO CRÍTICO — ESCALA DIFERE DA FASE 1

Cross-check gene-a-gene (amostra `GTEX-OHPL-2626-SM-2HMJA`):

| Gene | Fase 1 (`XENA_THCA.tsv`) | Global TPM | Δ log2 |
|---|---|---|---|
| ACTB | 17,28 | 11,56 | +5,72 |
| CCND1 | 9,69 | 2,71 | +6,98 |
| MYH7 | 9,48 | 1,98 | +7,50 |
| PRKCA | 10,34 | 2,34 | +8,00 |
| TP53 | 11,43 | 5,23 | +6,20 |

**Interpretação:** a Fase 1 usou `log2(norm_count + 1)` (contagem RSEM normalizada por quartil superior), enquanto a matriz global acessível é `log2(TPM + 0,001)`. As diferenças **não são constantes** (variam por gene), o que é compatível com a diferença de normalização count-vs-TPM (TPM divide pelo comprimento do gene). Portanto, **a Fase 1 NÃO é um subconjunto direto do arquivo TPM**; foi derivada de outra quantificação.

## 4. DECISÃO REGISTRADA

- Os arquivos `TcgaTargetGtex_rsem_gene_expected_count.gz` e `..._norm_count.gz` retornam **HTTP 403 AccessDenied** neste hub S3 (testado em 2026-09-06). A listagem do bucket também é negada.
- **Decisão:** a Fase 2 (global) será conduzida sobre `log2(TPM + 0,001)`, unidade **estabelecida e verificada**. A Fase 1 permanece como análise direcionada em `log2(norm_count+1)`.
- **Consequência documentada:** valores de DE e ranking da Fase 2 **não são diretamente comparáveis em magnitude absoluta** com a Fase 1; são comparáveis em direção/ranking. Esta limitação irá para a discussão.
- **Não foi possível reproduzir a quantificação exata da Fase 1** por indisponibilidade do arquivo de expected/norm count → **INFORMAÇÃO A CONFIRMAR** quanto à normalização exata da Fase 1 (consistente com a descrição `log2(norm_count+1)` do README, mas não verificável contra a fonte).

## 5. PRÓXIMO PASSO

QC formal da matriz global (58.581 × 783) → DEG genoma-wide → ranking → GSEA.
