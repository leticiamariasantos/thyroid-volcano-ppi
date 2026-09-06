# FINAL FALSIFICATION AUDIT — Tentativa de falsificar a conclusão

**Data:** 2026-09-06 · **Hipótese testada:** "ITGA2 é a hipótese molecular mais defensável."

## 1. Auditoria elo a elo da cadeia ITGA2

| Elo | Classificação | Evidência / ressalva |
|---|---|---|
| Diferencialmente expresso | SUPPORTED | TPM +4,32 / voom +2,86 / DESeq2 +3,36 (FDR <1e-86) |
| Replicação por múltiplos métodos | SUPPORTED | 3 métodos concordam em direção e significância |
| Replicação em microarray | SUPPORTED | GSE33630 +2,00 (2e-16); GSE60542 +2,50 (2e-11) |
| Replicação em RNA-seq independente | PARTIALLY_SUPPORTED | GSE224356 3/3 pares Up (+2,47/+2,60/+1,47), mas **n=3** |
| Proteína detectada | SUPPORTED | RPPA n=369 (mas **tumor-only**, inclui estroma) |
| Baixa expressão no normal | SUPPORTED | single-cell 3,1% positivas, consistente nos 7 pacientes (1,2–5,8%) |
| Proteína de membrana | SUPPORTED | integrina α2 é receptor transmembrana (biologia estabelecida) |
| Plausibilidade de targeting | PARTIALLY_SUPPORTED | membrana + baixa linha de base, mas sem binding/especificidade |
| Localização tumoral | **NOT_TESTABLE** | GSE182416 é só normal; GSE232237 sem anotação processada |
| Especificidade tumoral | **NOT_SUPPORTED** | não demonstrada; pode ser estroma/endotélio/plaqueta |

## 2. Perguntas adversariais (18) — respostas

1. **Produzido pelo estroma?** POSSÍVEL — não resolvido (bulk inclui estroma; single-cell normal não mostra origem tumoral).
2. **Associado a endotélio?** POSSÍVEL — ITGA2 (CD49b) é expresso em endotélio/plaquetas.
3. **Fibroblastos?** POSSÍVEL.
4. **Células imunes?** POSSÍVEL (leucócitos/plaquetas expressam integrinas).
5. **Marcador de composição?** PARCIALMENTE REJEITADO — ITGA2 é ROBUSTO após ajuste por músculo (logFC +3,46), não é marcador muscular; porém "composição estromal" genérica não foi formalmente testada.
6. **Dropout do single-cell?** **REJEITADO** — expressão baixa consistente nos 7 pacientes (1,2–5,8%); gene é detectável (max=28 counts) → não é artefato de sensibilidade.
7. **RPPA de população não tumoral?** POSSÍVEL — RPPA é bulk tumoral (inclui estroma); não distingue célula.
8. **Aumento bulk por composição?** POSSÍVEL — tumor tem estroma expandido; não descartado.
9. **GSE224356 n=3 suficiente?** INSUFICIENTE como prova estatística forte; é **consistência direcional**, não validação de grande poder.
10. **Três coortes independentes?** SIM — TCGA/GTEx (EUA), GSE224356 (BioProject PRJNA930711), GSE33630/GSE60542 (Ucrânia/Bélgica): pacientes, plataformas e centros distintos.
11. **Viés de plataforma favorecendo ITGA2?** Nenhuma evidência.
12. **ITGA2 em tecidos normais (risco de targeting)?** **SIM, RISCO REAL** — ITGA2 é expresso em plaquetas, endotélio, fibroblastos e múltiplos tecidos; **não** é tireoide-específico. Penalização de especificidade é insuficiente.
13. **Proteína acessível na superfície?** SIM (transmembrana), mas no contexto tumoral relevante não demonstrado.
14. **Evidência de internalização?** NENHUMA.
15. **Binding → internalização útil p/ delivery?** NÃO DEMONSTRADO.
16. **Evidência funcional em câncer tireoidiano?** NENHUMA.
17. **Evidência PTC-específica?** NENHUMA.
18. **Diferença PTC vs ATC?** NÃO ANALISADA.

## 3. Re-análise single-cell (pseudobulk)

ITGA2 % positivas por paciente: Thy01 3,2% · Thy04 5,8% · Thy05 1,2% · Thy06 3,4% · Thy10 3,4% ·
Thy15 2,1% · N3-GEX 2,0% → **consistentemente baixa (não é dropout)**.
FN1: 1,6–22,3% (mais variável), max=1982 → suporta **stromal/ECM**.

## 4. GSE232237

Decisão: **não resolvido por falta de anotação processada / custo metodológico incompatível**.
Só `RAW.tar` disponível (sem cell-type labels) → exigiria clustering/anotação completa. Documentado
como incerteza, não como resultado negativo.

## 5. Correções identificadas (erros reais / exageros)

1. **Score nanomédico super-premiou "proteína de membrana"** (accessibility=1,0 por categoria funcional)
   sem evidência de especificidade/localização tumoral. → corrigido na interpretação: acessibilidade
   conceitual ≠ acessibilidade tumoral demonstrada.
2. **"Replicado em 3 coortes"** deve ser qualificado: 2 coortes microarray (n moderado) + 1 RNA-seq
   pareada (n=3, consistência direcional). Não é "replicação de grande poder".
3. **"CNV não participa"** → reformulado: "não identificada alteração de CNV suficiente para explicar".
4. **"Expressão transcricional"** → reformulado: "não identificada mutação direta que explique".
5. **Metilação:** NÃO analisada com dados primários → permanece PENDENTE (não RESOLVIDA).

## 6. Teste de rejeição (10 argumentos hostis) e resposta

| # | Argumento hostil | Veredito |
|---|---|---|
| 1 | Batch TCGA×GTEx confunde tudo | RESOLVIDO (triangulação 3 coortes independentes) |
| 2 | Sinal muscular = composição | RESOLVIDO (single-cell: genes musculares ausentes na tireoide) |
| 3 | GSE224356 n=3 fraco | LIMITAÇÃO (consistência direcional, não prova forte) |
| 4 | Localização tumoral não demonstrada | **NÃO RESOLVIDO** (principal lacuna) |
| 5 | ITGA2 não específico (endotélio/plaquetas) | **NÃO RESOLVIDO** (risco de targeting) |
| 6 | Sem evidência funcional | NÃO RESOLVIDO |
| 7 | RPPA tumor-only | LIMITAÇÃO |
| 8 | Nanomedicina exagerada | CORRIGIDO (linguagem auditada; acessibilidade qualificada) |
| 9 | FN1/ITGA2 = estroma | PARCIAL (FN1 claramente estromal; ITGA2 não resolvido) |
| 10 | Efeito modesto em microarray | PARCIALMENTE RESOLVIDO (direção consistente, magnitude variável entre plataformas) |

## 7. Decisão final

### ITGA2 — categoria sustentada: **C. candidato de targeting** (candidato translacional prioritário)
NÃO é "D. target tumoral" nem "E. target nanomédico validado".

### FN1 — categoria sustentada: **C. componente de ECM/microambiente**

### Frase final conservadora (para resumo/manuscrito)

> "Em três coortes independentes e por três métodos, ITGA2 mostrou superexpressão tumoral
> consistente, com detecção proteica por RPPA e baixa expressão na tireoide normal avaliada por
> single-cell, sustentando sua priorização como candidato translacional para investigação de
> targeted delivery; a localização celular no tumor, a especificidade e a evidência funcional
> permanecem não demonstradas."

## 8. Próximo experimento falsificável

Perfil single-cell de PTC (GSE232237) ou IHC multiplexada de ITGA2 (tumor epitelial vs endotélio vs
estroma) — para determinar se a superexpressão é célula tumoral ou microambiente. Um resultado
"predominantemente endotelial/estromal" **falsificaria** a hipótese de targeting tumoral.

## 9. Veredito da falsificação

A hipótese ITGA2 **sobreviveu à tentativa de falsificação**, mas foi **rebaixada na linguagem**:
de "hipótese molecular mais defensável" para **"candidato translacional prioritário (não alvo
validado)"**, com três incertezas centrais registradas (localização tumoral, especificidade,
evidência funcional) e com a correção do viés de acessibilidade no score nanomédico.

---

## ATUALIZAÇÃO PÓS-ANÁLISE SINGLE-CELL TUMORAL (GSE232237)

**Nova evidência (marker-based, 7 PTC + 5 ATC):** ITGA2 é **predominantemente expresso em células
epiteliais/tumorais** do PTC (52,7% das células; média 3,25 vs 0,59 endotélio vs 0,15 fibroblasto),
e FN1 também é forte em células tumorais (168,9; 89,3%).

### Mudança de classificação

| Elo | Antes | Depois |
|---|---|---|
| Localização tumoral de ITGA2 | NOT_TESTABLE | **SUPPORTED (tumor-epitelial)** |
| FN1 origem | stromal/ECM | **tumor + estroma (EMT)** |
| Classificação falsificação ITGA2 | sobreviveu com incerteza | **sobreviveu com suporte adicional** |

### O que NÃO mudou (incertezas remanescentes)

1. **Especificidade tumoral** de ITGA2 permanece **NÃO ESTABELECIDA** — ITGA2 é expresso em
   plaquetas, endotélio e fibroblastos de tecidos normais (risco de on-target/off-tumor).
2. **Evidência funcional** em PTC (knockdown/CRISPR) permanece ausente.
3. **Internalização/binding** para delivery permanece não demonstrada.
4. Classificação por marcadores (sem anotação do autor no RAW.tar).

### Veredito atualizado

A expressão tumoral-epitelial de ITGA2 em single-cell **fortalece** a hipótese de targeting (a
proteína está nas células tumorais, não apenas no estroma), mas a **especificidade** (ITGA2 em
tecidos normais) e a **funcionalidade** permanecem as principais barreiras. A hipótese sobrevive,
porém permanece **candidato translacional, não alvo validado**.
