# PHASE 2 (reboot) — Limitações

## 1. Confundimento estrutural `source ≡ condition`

TCGA ≡ tumor e GTEx ≡ normal estão perfeitamente confundidos. Nenhuma correção de batch
elimina esse problema sem remover o sinal biológico. **Não** se alega correção completa;
a abordagem foi triangulação (limma/voom/DESeq2) + validação externa (GSE33630/GSE60542
pareadas) + single-cell + sensibilidade de composição.

## 2. Composição celular (músculo estriado)

A assinatura muscular (MYH7, MYL1, MYL2, ACTA1, TNNT3, CKM) é **composicional** — presente
no tecido normal GTEx (músculo esquelético adjacente/strap muscle). Genes musculares não
foram classificados como "artefato" *a priori*, mas a sensibilidade mostrou que as vias
robustas (proliferação/imune/metabolismo) independem deles.

## 3. Escala log2(TPM+0.001) vs contagens

O dataset principal é TPM (limma). A análise de contagens (voom/DESeq2) é sensibilidade.
As contagens (recount3) têm 1 amostra GTEx ausente e anotação gene-level distinta (GENCODE),
o que reduz a comparabilidade gene-a-gene com a matriz de símbolos.

## 4. Comparação TCGA vs GTEx — desequilíbrio de magnitude

O grande número de DEGs (12.200 no limma; 47% dos testados) reflete diferenças teciduais
amplas, não necessariamente oncogênese. Down (9.715) ≫ Up (2.485) é compatível com
composição (músculo/estroma no normal).

## 5. GSEA — desequilíbrio de estatísticas

Para vias com genes quase todos Up (ex.: Interferon Signaling), o fgsea retornou p-valores
NA em alguns métodos (estatísticas desbalanceadas). Reportado como "indeterminado", não
forçado.

## 6. PPI esparso

A rede STRING dos DEGs de maior |logFC| é esparsa (34 arestas) e dominada por composição
muscular. Centralidade é topológica, não funcional, e não identifica alvo terapêutico.

## 7. RPPA

ITGA2 e FN1 não possuem anticorpo no painel RPPA do PanCan Atlas. CCND1 presente, mas
correlação RNA-proteína fraca (0,051) — coerente com regulação pós-transcricional da ciclina D1.
RPPA é tumor-only: **não** prova localização celular.

## 8. Mutação/CNV

Exoma não capta mutações de promotor (TERT). Ausência de mutação em ITGA2/CCND1 **não**
significa ausência de mecanismo regulatório (a expressão pode ser transcricional/epigenética).

## 9. Single-cell

Classificação de GSE232237 é **marker-based** (sem anotação autoral); não tratada como
verdade absoluta. FN1 expresso também em fibroblastos (estroma), o que impede atribuir
exclusividade tumoral.

## 10. Generalização

Resultados restritos ao carcinoma papilífero de tireoide (THCA) no contexto TCGA/GTEx;
extrapolação a outros subtipos requer validação independente. Estudo exploratório, não
confirmatório.

## 11. Nanotecnologia

Nenhuma nanopartícula/nanocarreador foi pré-selecionado. A nanotecnologia é apresentada
apenas como consequência translacional **hipotética** da identificação de um candidato,
nunca como premissa da seleção de genes/vias.
