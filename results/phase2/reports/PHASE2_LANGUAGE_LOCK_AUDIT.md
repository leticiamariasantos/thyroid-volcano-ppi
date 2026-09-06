# PHASE 2 (reboot) — Auditoria Linguística (Language Lock)

**Data da auditoria:** 2026-09-06
**Escopo:** os 6 relatórios em `results/phase2/reports/` (FINAL_REPORT, METHODS, RESULTS, LIMITATIONS, FALSIFICATION_AUDIT, REPRODUCIBILITY).
**Critério:** busca exaustiva e automática por expressões proibidas ou de alto risco de *overclaiming*, com reescrita conservadora alinhada à conclusão oficial da Fase 2.

## 1. Expressões proibidas — busca exaustiva (violações duras)

| Expressão | Ocorrências |
|---|---|
| "descobrimos" | 0 |
| "novo alvo" | 0 |
| "primeiro" | 0 |
| "comprovamos" | 0 |
| "demonstramos que é alvo terapêutico" | 0 |
| "específico do tumor" | 0 |
| "seguro" | 0 |
| "eficaz" | 0 |
| "validação definitiva" | 0 |
| "causa" / "causal" | 0 |
| "alvo terapêutico validado" / "alvo validado" | 0 (apenas em negação: "não alvo validado") |
| centralidade PPI → relevância funcional/causalidade | 0 (todas qualificadas como "propriedade topológica") |

**Resultado:** nenhuma violação dura. As menções a "alvo" ocorrem exclusivamente em
negação ("não alvo validado", "não identifica alvo terapêutico"), e a centralidade de PPI
está sempre qualificada como "topológica".

## 2. Ocorrências de risco (overclaiming leve — framing do CCND1)

| # | Arquivo | Seção | Frase original (resumo) | Risco | Reescrita proposta |
|---|---|---|---|---|---|
| 1 | PHASE2_FINAL_REPORT.md | item 9 (CCND1) | "É gene de convergência do painel: core enrichment de Cell cycle, p53, Cellular senescence e Thyroid cancer." | Apresenta CCND1 como gene de convergência **do painel**, sem esclarecer que o core enrichment é via vias **ORIGINAIS** (não das 20 adicionais), em tensão com a conclusão oficial | "Entre os três candidatos, é o único que figura no core enrichment — restrito às vias ORIGINAIS (Cell cycle, p53) e a Cellular senescence/Thyroid cancer; não figura no core enrichment das 20 vias adicionais." |
| 2 | PHASE2_FINAL_REPORT.md | Conclusão | "…(com CCND1, TP53, CDK1… como nós centrais de convergência). CCND1 emerge como gene de convergência…" | "nós centrais de convergência" e "gene de convergência" superdimensionam o papel do CCND1 e omitem que se trata de vias ORIGINAIS | "…(TP53, CDK1, CDKN1A, CDKN2A, MDM2, CHEK1 e a maquinaria MCM/PCNA). CCND1 é o único dos três candidatos que figura no core enrichment — restrito às vias ORIGINAIS (Cell cycle, p53) e a Cellular senescence/Thyroid cancer; não figura no core enrichment das 20 vias adicionais." |
| 3 | PHASE2_FINAL_REPORT.md | item 15 | "Convergência molecular do painel: eixo ciclo celular/p53 (TP53, CCND1, CDK1…)" | Lista CCND1 como parte da convergência sem distinguir vias originais vs. adicionais | "Convergência molecular do painel (vias ORIGINAIS): eixo ciclo celular/p53 (TP53, CCND1, CDK1…)" |
| 4 | PHASE2_RESULTS.md | seção 7 (Convergência) | "O painel de 30 vias converge para o eixo ciclo celular/p53 — TP53 (7 vias), CCND1…" | Idem — não distingue que a convergência de CCND1 é via vias ORIGINAIS | "O painel de 30 vias converge para o eixo ciclo celular/p53 (vias ORIGINAIS) — TP53 (7 vias), CCND1…; as 20 vias adicionais convergem para proteassoma/apresentação de antígeno/OXPHOS/replicação de DNA." |
| 5 | PHASE2_RESULTS.md | seção 7 (Candidatos) | "CCND1 **converge** (core enrichment de Cell cycle, p53, Cellular senescence e Thyroid cancer)." | "converge" sem qualificação de que é via vias ORIGINAIS | "CCND1 figura no core enrichment apenas de vias ORIGINAIS (Cell cycle, p53) e de Cellular senescence/Thyroid cancer; não figura no core enrichment das 20 vias adicionais." |
| 6 | PHASE2_FALSIFICATION_AUDIT.md | Itens não sobreviveram | "(CCND1, em contraste, É core enrichment de Cell cycle, p53, Cellular senescence e Thyroid cancer.)" | Afirmação absoluta sem distinguir vias originais | "(CCND1, em contraste, figura no core enrichment de vias ORIGINAIS — Cell cycle, p53 — e de Cellular senescence/Thyroid cancer.)" |

## 3. Contagem total

- **Violações duras (expressões proibidas):** 0
- **Ocorrências de risco leve (framing):** 6
- **Total de problemas encontrados:** 6

## 4. Declaração final

Após correção, o texto está livre de *overclaiming* e alinhado à conclusão oficial da
Fase 2: a expansão das 20 vias adicionais **não** convergiu para ITGA2/FN1/CCND1 como
componentes dos programas mais enriquecidos; os três emergem como DEGs individuais robustos
e replicados (CCND1 com core enrichment restrito a vias **originais**), e ITGA2 permanece
"candidato translacional para investigação de direcionamento molecular" — nunca alvo validado.
