# VALIDATION LEAKAGE AUDIT — Auditoria de vazamento (data leakage)

**Data:** 2026-09-06

Verificação obrigatória de que os datasets GEO (GSE33630, GSE60542) foram usados
**exclusivamente para validação**, sem influenciar a seleção de candidatos, probes,
thresholds, pesos ou contrastes do discovery.

| Item | Status | Evidência |
|---|---|---|
| Nenhum candidato foi escolhido usando GEO | OK | Candidatos (FN1, ITGA2, CTSS, HLA-DPA1, CCND1, musculares) definidos na Fase 2 discovery, documentados em `RELATORIO_FASE2.md` (timestamp 08:40) **antes** de qualquer análise GEO (início 09:38). |
| Nenhuma probe selecionada por desempenho em GEO | OK | Regra de mapeamento **pré-registrada e uniforme**: para cada gene, o probe com **maior expressão média** (aplicado a todos os genes, não apenas FN1/ITGA2). |
| Nenhum threshold ajustado após observar GEO | OK | Cutoffs congelados no discovery: FDR<0,05, |log2FC|>1, STRING≥700, GSEA padj<0,25. |
| Nenhum peso de ranking alterado após observar GEO | OK | Pesos (60% DE + 40% topologia) definidos no discovery (`scripts/07_prioritization.R`) e congelados em `DISCOVERY_FREEZE.md`. |
| Nenhum contraste escolhido com base no resultado | OK | Contraste principal **pré-especificado** (PTC vs normal), sem ATC misturado; secundários (N0/N+, metástase) declarados exploratórios. |

## Cronologia (timestamps dos arquivos)

1. `documentation/RELATORIO_FASE2.md` (08:40) — discovery completo, candidatos já definidos.
2. `scripts/16_geo_download.R` (09:38) — primeiro acesso a GEO (após o discovery).
3. `scripts/17_geo_validate.R` (09:49) — validação (candidatos já congelados).
4. `documentation/DISCOVERY_FREEZE.md` (09:55) — formalização do congelamento.

## Veredito

**NENHUM vazamento detectado.** Os resultados GEO não foram usados para selecionar,
ajustar ou favorecer qualquer candidato, probe, threshold, peso ou contraste do discovery.

Os genes que replicaram (FN1, ITGA2, CTSS, HLA-DPA1, CCND1) e os que não replicaram
(musculares) emergiram **naturalmente** da análise independente, sem pós-seleção.
