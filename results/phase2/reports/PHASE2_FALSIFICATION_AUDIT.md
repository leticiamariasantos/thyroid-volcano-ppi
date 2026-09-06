# PHASE 2 (reboot) — Auditoria de Falsificação

Tentativa explícita de refutar os principais resultados antes de concluir.

## Para ITGA2

| Pergunta de falsificação | Resposta | Conclusão |
|---|---|---|
| O sinal desaparece com controle de composição? | Não — ITGA2 não é muscular; estável após remoção dos 48 marcadores | passa |
| É específico de um único método? | Não — Up em limma (+2,47), voom (+2,98), DESeq2 (+3,37) | passa |
| Depende de uma única coorte? | Não — replicado em GSE33630 (+2,00) e GSE60542 (+2,50) | passa |
| Depende apenas de TCGA/GTEx? | Não — GSE33630/GSE60542 são coortes pareadas independentes | passa |
| Desaparece no single-cell? | Não — Epithelial_tumor é o compartimento predominante (2,53 vs 0,12 fibroblasto) | passa |
| É predominantemente estromal? | Não no single-cell analisado (epitelial predominante) | passa (com ressalva) |
| Há expressão normal relevante? | Sim — ITGA2 tem expressão basal em alguns tecidos normais; não há especificidade absoluta | **limitação mantida** |
| Literatura contradiz? | ITGA2 é conhecido em adesão/migração; superexpressão em PTC é plausível, não contradita | passa |
| Proteína confirma RNA? | Indeterminado — ITGA2 fora do painel RPPA | **não testável** |
| Risco on-target/off-tumor? | Sim — integrina α2 expressa em plaquetas/endotélio → risco real | **risco mantido** |

**Veredicto sobre ITGA2:** o sinal de superexpressão tumoral sobrevive a todos os controles
testáveis, mas **não** há evidência de especificidade tumoral absoluta, acessibilidade
funcional, internalização ou eficácia. Permanece **candidato translacional para investigação
de direcionamento molecular**, não alvo validado.

## Para as vias

| Pergunta | Resposta |
|---|---|
| Alguma via significativa é consequência de genes musculares? | Não — as 6 vias robustas permanecem após remoção dos 48 marcadores |
| Há redundância que faz vias parecerem independentes? | Parcialmente — módulos de adesão/ECM e inflamação/imune são correlacionados (documentado) |
| O resultado depende de uma coleção específica? | Não — Proteasome/OXPHOS/Cell cycle/p53/Antigen processing são consistentes entre KEGG, Reactome e Hallmark |
| O resultado desaparece em outro método? | As 6 vias robustas são consistentes em limma/voom/DESeq2; outras vias (Interferon, Focal adhesion) são método-dependentes (exploratórias) |

## Itens que NÃO sobreviveram à falsificação

1. **ITGA2/FN1 como componentes do core enrichment** — não figuram no leading edge de
   nenhuma via; são DEGs individuais dentro de programas de adesão/ECM. (CCND1, em contraste,
   figura no core enrichment de vias ORIGINAIS — Cell cycle, p53 — e de Cellular
   senescence/Thyroid cancer, não das 20 vias adicionais.)
2. **ITGA2/FN1 como alvos RPPA-validados** — fora do painel RPPA.
3. **Assinatura muscular como sinal tumoral** — refutada (composicional; não replica em
   coortes pareadas).
4. **PPI como rede tumoral coerente** — refutada (esparsa e dominada por composição).
