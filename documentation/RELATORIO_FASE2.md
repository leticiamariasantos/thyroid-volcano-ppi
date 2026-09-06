# RELATÓRIO CIENTÍFICO — FASE 2 (TRANSCRIPTOMA GLOBAL)
## Carcinoma de tireoide (THCA vs normal GTEx) — caracterização global, vias, redes e priorização

**Data:** 2026-09-06 · **Ambiente:** R 4.6.1 · **Tipo de estudo:** exploratório, gerador de hipóteses
**Escopo:** transcriptoma global (não mais restrito à via hsa04919)

> Estrutura de evidência (aplicada em todo o relatório):
> **N1** observação direta · **N2** associação estatística · **N3** associação de rede · **N4** inferência funcional · **N5** hipótese nanomédica · **N6** hipótese translacional.

---

## 1. RESUMO EXECUTIVO (separado por natureza da afirmação)

**RESULTADOS OBSERVADOS (N1–N2).** Sobre uma matriz global de **58.581 genes × 783 amostras** (unidade estabelecida: **log₂(TPM+0,001)**; 504 THCA / 279 tireoide normal GTEx), após filtro de baixa expressão foram testados **20.376 genes**, dos quais **8.161 são diferencialmente expressos** (1.530↑, 6.631↓) a FDR<0,05 e |log₂FC|>1. O GSEA global (372 vias KEGG) mostrou 50 vias com padj<0,25 e 22 com padj<0,05, **dominadas por vias imunes** (apresentação de antígeno, doença tireoidiana autoimune, rejeição de aloenxerto) e por proteassoma, fosforilação oxidativa, replicação de DNA e p53. Do painel a priori de 10 vias, **apenas p53 foi significativa** (NES=1,87; padj=0,001).

**INFERÊNCIAS (N3–N4).** A rede PPI (177 genes DEG∩vias → 148 nós/727 arestas, 8 comunidades) concentra-se em **matriz extracelular/adesão** (FN1, LAMA2, HSPG2, ITGA2, ITGA2B, ELN), **citoesqueleto muscular** (MYH7, FLNC, TPM2, TCAP, ACTA1) e **apresentação de antígeno** (HLA-DPA1, CTSS).

**HIPÓTESES (N5–N6).** *Ver seção 7.* Nenhuma nanopartícula/terapia foi selecionada a priori.

---

## 2. DADOS E PROCESSAMENTO (RASTREABILIDADE)

| Item | Valor | Status |
|---|---|---|
| Dataset | `TcgaTargetGtex_rsem_gene_tpm.gz` (TOIL recompute, UCSC Xena) | verificado |
| Acesso | 2026-09-06; Last-Modified 2021-04-09; 1.323.254.426 bytes (íntegro) | verificado |
| Genes | 60.498 → **58.581** após colapso de símbolos (regra: maior média) | verificado |
| Amostras | 783 (504 THCA + 279 GTEx), IDs idênticos aos da Fase 1 | verificado |
| Unidade | **log₂(TPM+0,001)** (min −9,97; máx 17,54; soma TPM≈1e6) | verificado |
| Filtro baixa expressão | TPM>1 em ≥10% das amostras → 20.376 genes testados | verificado |
| DEG | limma (`eBayes`, `trend=TRUE`), contraste THCA−Normal, BH | verificado |
| GSEA | `fgsea` pre-ranked (estatística t moderada), 372 vias KEGG | verificado |

**Divergência documentada com a Fase 1:** a Fase 1 usou `log₂(norm_count+1)`; a matriz global acessível é `log₂(TPM+0,001)`. Os arquivos de expected/norm count retornaram **HTTP 403**. Valores absolutos não são diretamente comparáveis entre fases (direção/ranking sim). **INFORMAÇÃO A CONFIRMAR:** normalização exata da Fase 1.

---

## 3. RESULTADOS (NÍVEL 1–2)

### 3.1 DEG genoma-wide
- Testados: 20.376 · DEGs: **8.161** (1.530↑ / 6.631↓).
- Estável em sensibilidade (FDR 0,01/0,05/0,10 × |logFC| 0,5/1/1,5 → 5.555–12.800).
- Ranking (t moderada) é essencialmente idêntico a signed −log10(p) (Spearman = 1,000).

### 3.2 GSEA global (exploratório)
Top vias significativas (padj<0,05) — **NÃO = seleção, apenas observação**:
apresentação/processamento de antígeno (hsa04612), proteassoma (hsa03050), rejeição de aloenxerto, diabetes tipo I, doença enxerto-vs-hospedeiro, envelope cornificado, **doença tireoidiana autoimune (hsa05320)**, interação por moléculas de adesão celular (hsa04514), fosforilação oxidativa, artrite reumatoide, p53 (hsa04115), replicação de DNA (hsa03030), etc.

### 3.3 Painel a priori (10 vias — imutável)
| Via | NES | padj | Conclusão |
|---|---|---|---|
| hsa05216 Thyroid cancer | 1,01 | 1 | não significativa |
| hsa04919 Thyroid hormone signaling | −0,94 | 1 | não significativa |
| hsa04010 MAPK | −1,04 | 0,94 | não significativa |
| hsa04151 PI3K-Akt | −1,07 | 0,77 | não significativa |
| hsa04150 mTOR | −1,00 | 1 | não significativa |
| **hsa04115 p53** | **1,87** | **0,001** | **significativa (enriquecimento associado)** |
| hsa04210 Apoptosis | 1,21 | 0,34 | não significativa (tendência) |
| hsa04110 Cell cycle | NA | NA | p-valor não calculado (desbalanceado) |
| hsa04310 Wnt | −0,97 | 1 | não significativa |
| hsa04064 NF-kB | −1,00 | 1 | não significativa |

> **Interpretação permitida (N2):** enriquecimento associado à condição. **Não permitida:** "via ativada/inibida" (RNA ≠ atividade proteica). A ausência de sinal para MAPK/PI3K-Akt/mTOR é **compatível** com o fato de essas vias serem ativadas por mutação/fosforilação, não por abundância de mRNA — resultado negativo relevante, não omitido.

---

## 4. REDUNDÂNCIA E SELEÇÃO DE VIAS (N3)

22 vias significativas (padj<0,05) → matriz Jaccard → 6 módulos → **7 vias não-redundantes** (representante de menor padj por módulo + p53):
hsa04612 (apresentação de antígeno), hsa03050 (proteassoma), hsa04820 (citoesqueleto muscular, **down**), hsa03030 (replicação de DNA), hsa00510 (N-glicana), hsa01232 (metabolismo de nucleotídeos), hsa04115 (p53).

---

## 5. PPI E MÓDULOS (N3)

- Entrada: **177 genes** (DEG ∩ vias selecionadas) → STRING v12.0, score≥700 → **148 nós / 727 arestas**, 8 comunidades (modularidade 0,67).
- Módulos coerentes: **ECM/adesão** (FN1↑, ITGA2↑, LAMA2↓, HSPG2↓, ELN↓, ITGA2B↓); **citoesqueleto/contrátil** (MYH7↓, FLNC↓, TPM2↓, TCAP↓, ACTA1↓); **apresentação de antígeno** (HLA-DPA1↑, CTSS↑); **ciclo celular** (CCND1↑).
- **Hub ≠ alvo/biomarcador.** FLNC (betweenness 0,27), ITGA2B (0,22), FN1 (0,16) são hubs topológicos, nada além disso.

---

## 6. PRIORIZAÇÃO MULTICRITÉRIO (N2+N3)

Modelo (pré-definido): 60% evidência transcriptômica + 40% topologia. Tabela completa em `09_target_prioritization/candidate_ranking.tsv`.

Top candidatos: **MYL1↓, MYL2↓, MYH7↓, FLNC↓, MYL7↓, ACTA1↓, FN1↑, CKM↓, LAMA2↓, ITGA2B↓, ACTN2↓, DES↓, ITGA2↑, CTSS↑, … CCND1↑**.

**Estabilidade:** composite vs pesos iguais Jaccard(top30)=0,88; composite vs DE-only=0,71; DE-only vs rede-only=0,28. O ranking é **dominado pela evidência DE** (esperado: sinal forte) e moderadamente estável.

### ⚠ CAVEAT CRÍTICO (não omitir)
A dominância de genes **musculares/contráteis downregulated** (MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM) provavelmente reflete **diferença de composição tecidual**: o tecido tireoidiano normal (GTEx) pode conter músculo esquelético adjacente (strap muscle), ausente no tumor. **Não é evidência de regulação tumoral intrínseca.** Este é um confundidor de composição, NÃO um achado de biologia do câncer. (Limitação estrutural do desenho TCGA-vs-GTEx, já reconhecida.)

---

## 7. AVALIAÇÃO DE PLAUSIBILIDADE NANOMÉDICA (OBRIGATORIAMENTE POSTERIOR)

**Pergunta-guia:** "Quais características moleculares observadas poderiam, de forma racional e explicitamente hipotética, justificar investigação futura por nanomedicina?"

**Resposta honesta:** existe uma **conexão hipotética limitada**, baseada em características de localização/expressão — não em validação funcional.

| Alvo | Evidência molecular (dados) | Característica relevante | Aplicação possível | Nível |
|---|---|---|---|---|
| FN1 (fibronectina) | ↑ logFC +4,32; FDR 6e-96; hub (grau 30) | proteína de matriz extracelular, expressa no microambiente tumoral | reconhecimento/alvo de ECM (hipotético) | **N5** |
| ITGA2 (integrina α2) | ↑ logFC +2,47; FDR 2e-89; hub | receptor de superfície celular / adesão | direcionamento molecular (hipotético) | **N5** |
| HLA-DPA1 / CTSS | ↑; módulo de apresentação de antígeno | superfície/imunidade | assinatura de imunovigilância/contexto imune | **N5** |
| CCND1 (ciclina D1) | ↑ logFC +2,26; FDR 6e-145; via p53 | intracelular (não-superfície) | alvo intracelular — **não** adequado a targeting de superfície | **N5/N6** |

**O que os dados NÃO demonstram:** eficácia, segurança, validade clínica, internalização, especificidade tumoral, nem que qualquer alvo seja "tratável por nanopartícula". A comparação é cross-coorte (batch confundido) e a assinatura é dominada por confundidor de composição.

**Conclusão da camada nanomédica:** **NÃO há evidência suficiente para propor uma estratégia nanomédica específica.** As características de ECM/superfície de FN1 e ITGA2 (superexpressos no tumor) constituem **hipótese para investigação futura** — e somente isso. Nenhuma nanopartícula, nanocarreador, sistema magnético, hipertermia ou biossensor foi selecionado.

---

## 8. LIMITAÇÕES (explícitas)

1. **Batch confundido com condição** (TCGA=tumor, GTEx=normal) — limitação estrutural, não corrigível sem destruir o sinal.
2. **Composição celular/tecidual** — assinatura muscular provavelmente artefato de composição; sem deconvolução.
3. **TPM vs count** — DEG sobre log₂(TPM+0,001); normalização exata da Fase 1 **A CONFIRMAR**.
4. **RNA ≠ proteína/atividade** — nenhuma inferência de atividade proteica.
5. **Sem validação externa independente** (não foi localizado dataset externo metodologicamente comparável nesta execução).
6. **Rede PPI in silico** (STRING) — associações preditas/inferidas.
7. **8.161 DEGs** — magnitude inflada pelo confundidor de coorte; interpretar com cautela.

---

## 9. TRABALHOS FUTUROS

1. Re-obter expected/norm count para DE com `voom`/DESeq2 (mais rigoroso que TPM).
2. Deconvolução de composição celular (ex.: CIBERSORTx) para separar sinal tumoral de estroma/músculo.
3. Validação externa (GEO THCA independente; proteômica).
4. Análise de localização subcelular sistemática (UniProt/Gene Ontology) para os candidatos prioritários.
5. Estudos funcionais experimentais — antes de qualquer hipótese nanomédica concreta.

---

## 10. CONCLUSÃO (respondendo às 9 perguntas)

1. **Alterações identificadas:** 8.161 DEGs (1.530↑/6.631↓) em 20.376 genes testados.
2. **Vias com evidência de enriquecimento:** imunes (apresentação de antígeno, tireoidite autoimune), proteassoma, fosforilação oxidativa, replicação de DNA, **p53**.
3. **Mecanismos convergentes:** inflamação/imunidade, proliferação (proteassoma/ciclo/p53), remodelamento de ECM.
4. **Redes/módulos relevantes:** ECM/adesão e apresentação de antígeno.
5. **Moléculas priorizadas:** FN1, ITGA2 (ECM/superfície, ↑); CCND1, CTSS; (genes musculares ↓ = provável artefato de composição).
6. **Robustez:** DEG e ranking GSEA robustos; priorização moderadamente estável (dominada por DE); rede preserva núcleo a score 900.
7. **Estratificação molecular:** dados agregados tumor-vs-normal — **não** permitem personalização individual; apenas "priorização molecular" e "nanomedicina de precisão potencial".
8. **Fundamento para nanomedicina:** **insuficiente para proposta específica**; apenas hipótese futura (FN1/ITGA2 de superfície/ECM).
9. **Perguntas em aberto:** contribuição real do tumor vs composição; validação independente; atividade proteica.

**Nunca afirmado neste relatório:** que nanopartícula "funciona", que alvo é terapêutico, causalidade a partir de transcriptômica, ou personalização individual a partir de dados agregados.

---

## 11. APÊNDICE — ANÁLISE BASEADA EM CONTAGENS (SENSIBILIDADE, N2)

### 11.1 Dados de contagem (FASE 4/7)

O arquivo de contagens original do Xena (`TcgaTargetGtex_rsem_gene_count.gz` e
`..._norm_count.gz`) retornou **HTTP 403 (AccessDenied)**. Como fonte alternativa legítima
para as **mesmas amostras**, usou-se o **recount3** (Lieber Institute; STAR, anotação GENCODE
G026). Proveniência registrada em `data/global/MANIFEST_counts.tsv`:

- TCGA THCA (`project="THCA"`, `file_source="tcga"`) → 572 amostras; `raw_counts`.
- GTEx THYROID (`project="THYROID"`, `file_source="gtex"`) → 706 amostras; `raw_counts`.
- **Correspondência:** 504/504 TCGA + 278/279 GTEx → **782/783 amostras** (única ausente:
  `GTEX-SUCS-0226-SM-5CHQG`).
- Matriz final: **56.937 genes × 782 amostras**, contagens inteiras, não-negativas, sem NA.
- **Caveat:** as contagens recount3 são *gene sums* (soma de cobertura por base), library sizes
  ~5×10⁸–1,4×10¹⁰ — proporcionais a read counts (válidas para voom/DESeq2), mas **não** são
  "raw read counts" RSEM e usam quantificação/anotação distintas do TOIL (GENCODE v23).

### 11.2 DE por contagem (FASE 7)

Filtro pré-definido (comparável ao TPM): CPM > 1 em ≥10% das amostras → 22.118 genes testados.

| Método | Genes testados | DEGs | up | down |
|---|---|---|---|---|
| limma (log₂ TPM, seção 3) | 20.376 | 8.161 | 1.530 | 6.631 |
| limma-voom (counts) | 22.118 | 7.185 | 2.710 | 4.475 |
| DESeq2 (counts) | 22.118 | 6.976 | 2.888 | 4.088 |

### 11.3 Concordância TPM × voom × DESeq2 (FASE 8)

Genes em comum = 18.459.

| Comparação | logFC Pearson | logFC Spearman | direção | Jaccard DEG | ranking Spearman |
|---|---|---|---|---|---|
| TPM × voom | 0,895 | 0,922 | 0,813 | 0,602 | 0,891 |
| TPM × DESeq2 | 0,871 | 0,906 | 0,808 | 0,557 | 0,907 |
| voom × DESeq2 | 0,954 | 0,978 | 0,964 | 0,857 | 0,983 |

**Leitura:** voom e DESeq2 (ambos baseados em contagem) concordam fortemente entre si
(esperado) e **bem** com o TPM (Spearman ≈ 0,90–0,92). Nenhum método foi selecionado pelo
número de DEGs; a concordância de direção/ranking é o critério de robustez.

Classificação de estabilidade (genes comuns): **ROBUSTO** 4.293 · **PARCIALMENTE ROBUSTO** 1.218
· **MÉTODO-DEPENDENTE** 2.471 · **NÃO ROBUSTO** 10.477.

### 11.4 Composição tecidual (FASE 9/10)

| Categoria | marcadores testados | DEGs | up | down |
|---|---|---|---|---|
| muscle_skeletal | 31 | 24 | 0 | **24** |
| thyroid_epithelial | 13 | 7 | 1 | 6 |
| fibroblast_ecm | 16 | 9 | 1 | 8 |
| endothelial | 9 | 9 | 0 | 9 |
| immune | 18 | 5 | 3 | 2 |

- **Score muscular:** THCA = −2,19 vs Normal = **+1,52** (log₂ TPM; Mann-Whitney p = 8,2×10⁻¹¹⁴).
- **PC1 correlaciona −0,845** com a assinatura muscular (e −0,657 epitélio tireoidiano,
  −0,681 endotelial).
- Genes musculares (MYH7, MYL2, ACTA1, TNNT3, CKM) são **ROBUSTOS** nos 3 métodos, todos **down**
  no tumor — porém o sinal está **concentrado nos GTEx normais**, compatível com músculo
  esquelético adjacente (strap muscle). **Conclusão: artefato de composição, não biologia
  tumoral intrínseca.** (MYL1 não entrou no conjunto comum por ausência/filtro na anotação G026.)

### 11.5 Genes-chave (FASE 16)

| Gene | logFC TPM | logFC voom | logFC DESeq2 | Robustez |
|---|---|---|---|---|
| FN1 | +4,32 | +4,96 | +6,10 | ROBUSTO |
| ITGA2 | +2,47 | +2,86 | +3,36 | ROBUSTO |
| MYH7 | −7,01 | −7,04 | −5,91 | ROBUSTO |
| MYL2 | −8,69 | −9,43 | −6,29 | ROBUSTO |
| ACTA1 | −7,17 | −5,41 | −5,53 | ROBUSTO |
| TNNT3 | −4,77 | −4,72 | −3,74 | ROBUSTO |
| CKM | −6,92 | −6,65 | −5,81 | ROBUSTO |

**FN1 e ITGA2 permanecem** superexpressos no tumor de forma robusta — mas a interpretação
continua sendo **hipótese exploratória (N5)**, não alvo terapêutico validado. Os genes
musculares, embora robustos, **não** são candidatos (artefato composicional).

### 11.6 Limitações adicionais desta sensibilidade

1. Contagem alternativa (recount3/STAR G026) ≠ TOIL RSEM (GENCODE v23) do TPM.
2. 1 amostra GTEx ausente (782 vs 783).
3. Deconvolução formal (CIBERSORTx/assinatura de referência) **não executada** — usado
   marker-based scoring como estimativa/modelo.
4. Validação externa (GEO GSE33630/GSE60542) **pendente** — deve ocorrer após congelar o
   discovery (FASE 18/19).

---

## 12. APÊNDICE B — VALIDAÇÃO EXTERNA (GEO, FASE 18/19)

### 12.1 Datasets independentes (congelados antes da validação)

| Dataset | Papel | Plataforma | PTC | Normal |
|---|---|---|---|---|
| GSE33630 | primária | GPL570 (HG-U133 Plus 2.0) | 49 | 45 (patient-matched) |
| GSE60542 | secundária/robustez | GPL570 | 33 | 30 |

**Ponto-chave:** nos GEO, o "normal" é **tecido pareado/adjacente** (não GTEx), o que **controla a
composição** e permite testar a hipótese composicional do sinal muscular.

### 12.2 Genes congelados — resultado

| Gene | discovery logFC | GSE33630 (FDR) | GSE60542 (FDR) | Classificação |
|---|---|---|---|---|
| FN1 | +4,32 | +2,91 (2e-25) | +3,06 (1e-14) | INDEPENDENTLY REPLICATED |
| ITGA2 | +2,47 | +2,00 (2e-16) | +2,50 (2e-11) | INDEPENDENTLY REPLICATED |
| CTSS | +2,42 | +1,38 (5e-7) | +1,12 (3e-4) | INDEPENDENTLY REPLICATED |
| HLA-DPA1 | +1,22 | +1,06 (1e-4) | +0,81 (2e-3) | INDEPENDENTLY REPLICATED |
| CCND1 | +2,26 | +1,42 (1e-21) | +1,07 (7e-11) | INDEPENDENTLY REPLICATED |
| MYH7 | −7,01 | −0,01 (0,86) | +0,12 (0,08) | NOT REPLICATED |
| MYL1 | −8,30 | +0,05 (0,79) | +0,11 (0,46) | NOT REPLICATED |
| MYL2 | −8,69 | −0,11 (0,51) | +0,05 (0,79) | NOT REPLICATED |
| ACTA1 | −7,17 | −0,12 (0,64) | +0,11 (0,49) | NOT REPLICATED |
| TNNT3 | −4,77 | −0,07 (0,27) | −0,07 (0,38) | NOT REPLICATED |
| CKM | −6,92 | +0,07 (0,70) | +0,15 (0,13) | NOT REPLICATED |

Concordância global (4.584 genes discovery-DEG testados): direção 0,67 (GSE33630) / 0,72 (GSE60542);
Spearman logFC 0,42 / 0,47.

### 12.3 Conclusão da validação

1. **FN1, ITGA2, CTSS, HLA-DPA1, CCND1** replicam direção e significância em **duas coortes
   independentes** → sinal tumor-intrínseco genuíno (N2 corroborado externamente).
2. **Genes musculares NÃO replicam** nos GEO (normal pareado, sem strap muscle) → **confirma
   definitivamente** que o sinal muscular do discovery (TCGA-vs-GTEx) era **artefato de composição
   tecidual**, não biologia tumoral intrínseca.
3. A hipótese composicional da FASE 9/10 é **corroborada por validação externa independente**.

**FN1/ITGA2:** permanecem como hipóteses exploratórias (N5/N6) — superexpressos no tumor e
replicados, mas sem validação funcional. **Genes musculares: descartados** (composicionais).
