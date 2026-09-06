# ADVANCED VALIDATION REPORT — Validação e triangulação independente

**Data:** 2026-09-06 · **R:** 4.6.1
**Discovery (congelado):** TCGA-THCA (504) vs GTEx normal (279), log2(TPM+0.001); complemento recount3 (voom/DESeq2).

---

## Respostas objetivas (ETAPA 33)

### 1. A composição tecidual explica quanto do sinal?
Substancialmente o componente **muscular**. O score muscular difere −3,71 entre TCGA e GTEx
(Cohen's d = −3,75; FDR 4e-113). Ao incluir o score muscular como covariável no modelo, os genes
musculares **colapsam** (MYL2 −8,69→−1,92; ACTA1 −7,17→−1,03; TNNT3 −4,77→**+0,66**, inverte).
A correlação global de ranking entre modelo com/sem composição é alta (Spearman 0,95), mas o
componente muscular específico é quase inteiramente composicional.

### 2. A deconvolução formal confirma a composição?
A deconvolução formal (CIBERSORTx/EPIC/MuSiC/BisqueRNA) **não foi aplicável** (falta referência
single-cell com músculo esquelético; CIBERSORTx exige Docker). O marker-based scoring
(MCPcounter-like) **confirma** o sinal muscular no normal GTEx, e a **validação externa** (GEO com
normal pareado, sem strap muscle) confirma que o sinal desaparece (diff −3,71 → −0,1).

### 3. FN1 permanece robusto? — SIM
Robusto em TPM (+4,32), voom (+4,96), DESeq2 (+6,10), após composição (+5,53), e **replicado**
em GSE33630 (+2,91) e GSE60542 (+3,06).

### 4. ITGA2 permanece robusto? — SIM
Robusto em TPM (+2,47), voom (+2,86), DESeq2 (+3,36), após composição (+3,46), e **replicado**
em GSE33630 (+2,00) e GSE60542 (+2,50).

### 5. Expressos por células tumorais ou microambiente?
**Não determinado diretamente** (single-cell não executado). FN1 é classicamente de estroma/ECM
(fibroblastos), o que é **compatível** com a interpretação de que o aumento bulk de FN1 reflete o
**microambiente tumoral** (remodelamento de ECM), não necessariamente células tumorais. Esta
questão permanece **hipótese** — requer single-cell (ETAPA 12, não executada por limitação técnica).

### 6. Existe replicação independente em RNA-seq?
**Parcial.** A análise recount3 (voom/DESeq2) usa as mesmas amostras (re-quantificação, não
totalmente independente). A validação externa usou **microarray** (GSE33630/GSE60542). Uma coorte
**bulk RNA-seq independente** foi **pesquisada** no GEO (9 candidatos avaliados) e **nenhuma** foi
adequada (single-cell, amostras pequenas, ou foco em genes específicos). Registrada como limitação.

### 7. Hallmark confirma p53?
**Não diretamente.** KEGG p53 (hsa04115) foi significativa (NES 1,87), mas `HALLMARK_P53_PATHWAY`
**não** foi significativa (NES 0,29; padj NA). São gene sets diferentes; a discrepância é registrada
honestamente (o sinal p53 é capturado pelo gene set KEGG, não pelo Hallmark).

### 8. Reactome confirma os principais processos?
**Sim, em parte.** Reactome mostra forte sinal de **tradução/ribossomo** (EUKARYOTIC_TRANSLATION,
NES 4,0), **proteassoma** (PROTEASOME_ASSEMBLY) e **processamento de antígeno** — convergente com o
KEGG (proteassoma, apresentação de antígeno). Adiciona o eixo de tradução/ribossomo como processo
emergente.

### 9. Existe evidência proteica? — PARCIAL / NÃO AVALIADA
TCGA RPPA cobre ~200 proteínas; **FN1** está no painel RPPA (fibronectina), **ITGA2 não**. A
extração dos dados RPPA de THCA não foi executada nesta fase (requer TCGAbiolinks/API; registrado
como pendente). Sem evidência proteica concluída.

### 10. Existe evidência de IHC? — NÃO AVALIADA
Human Protein Atlas possui dados de FN1/ITGA2, mas a extração quantitativa não foi executada.
Registrado como **not testable nesta execução**.

### 11. Existe relação com BRAF/RAS/TP53? — NÃO ANALISADA DIRETAMENTE
A paisagem mutacional de THCA (BRAF V600E ~60%, RAS ~10–15%, TP53 raro em PTC) é conhecida da
literatura, mas **não foi extraída de dados primários** nesta fase. A associação de FN1/ITGA2/CCND1
a subtipos mutacionais permanece **não testada**.

### 12. Existe evidência epigenética? — NÃO AVALIADA
Metilação de FN1/ITGA2 não foi analisada. Registrado como pendente.

### 13. Existe heterogeneidade entre coortes?
**Sim, moderada.** Concordância global de direção 0,67 (GSE33630) e 0,72 (GSE60542); Spearman de
logFC 0,42–0,47. Os candidatos tumor-intrínsecos são consistentes; o sinal muscular é específico do
desenho TCGA/GTEx.

### 14. Quais resultados são robustos?
FN1, ITGA2, CTSS, HLA-DPA1, CCND1 (replicados); vias imunes/proteassoma/replicação/nucleotídeos;
PPI núcleo ECM/adesão (FN1/ITGA2 hubs estáveis).

### 15. Quais são dependentes de composição?
Genes musculares (MYH7/MYL1/MYL2/ACTA1/TNNT3/CKM) e via "Cytoskeleton in muscle cells".

### 16. Quais são dependentes do método?
Nenhum dos candidatos não-musculares; o sinal p53 é capturado por KEGG mas não por Hallmark
(dependência da definição do gene set, não do método DE).

### 17. Quais resultados são contraditos?
A via muscular (direção invertida nos GEO). Nenhum candidato não-muscular foi contradito.

### 18. Quais permanecem apenas como hipóteses?
FN1 e ITGA2 como alvos de reconhecimento/targeting (N5/N6); a origem celular (tumor vs estroma);
qualquer aplicação nanomédica.

### 19. Evidência suficiente para afirmação translacional?
**Não.** Há associação estatística robusta e replicação independente, mas **sem** validação funcional,
proteica concluída ou single-cell. Máximo permitido: hipótese molecular prioritária.

### 20. Estratégia nanomédica suportada?
**NÃO.** Nenhuma nanopartícula/nanocarreador/hipertermia/biossensor é suportado pelos dados atuais.

---

## Conclusão geral
A triangulação independente **corrobora** os candidatos tumor-intrínsecos (FN1, ITGA2, CTSS,
HLA-DPA1, CCND1) e **confirma** que o sinal muscular era composicional. FN1/ITGA2 permanecem como
hipóteses moleculares prioritárias (N5/N6), sem elevação para biomarcador/alvo terapêutico.

*Resultados: `results/gsea/`, `results/composition/composition_covariate_*`, `results/integration/EVIDENCE_MATRIX.tsv`.*
