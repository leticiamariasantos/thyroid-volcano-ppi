# TRANSPARÊNCIA — Datasets pesquisados, incluídos, excluídos (ETAPA 29)

**Data:** 2026-09-06

## 1. Datasets externos PESQUISADOS para validação RNA-seq independente

Busca programática via NCBI E-utilities (`esearch db=gds`, termos: "thyroid cancer RNA-seq normal
thyroid", "papillary thyroid carcinoma normal RNA-seq", "PTC normal thyroid bulk RNA-seq").

| Accession | Título (resumo) | Motivo exclusão |
|---|---|---|
| GSE331276 | Pan-cancer KCNQ3 (oncogênico) | foco em gene único; não PTC-vs-normal |
| GSE297369 | FAM111B knockdown em HCC | célula hepática; não tireoide |
| GSE263800 | SK4 splicing em PTC | 6 amostras; foco em splicing |
| GSE274674 | single-cell CD4+ T em PTC | single-cell; não bulk |
| GSE239976 | neutrófilos high-density em câncer tireoidiano | subpopulação imune; não bulk |
| GSE249361 | single-cell antígeno em metástase hepática de mama | single-cell; mama |
| GSE201545 | decabromodifeniléter em PTC | exposição química; 6 amostras |
| GSE159330 | small RNA-seq (miRNA) em PTC | miRNA, não mRNA |
| GSE61844 | DNase I hypersensitive sites | epigenômica, não expressão |

**Conclusão:** nenhuma coorte bulk RNA-seq PTC-vs-normal adequada foi encontrada no GEO. O bulk
RNA-seq de THCA é dominado pelo próprio TCGA (já usado como discovery). A validação independente
foi feita em **microarray** (GSE33630, GSE60542). Registrado como limitação.

## 2. Datasets INCLUÍDOS

| Dataset | Papel | Tipo |
|---|---|---|
| TcgaTargetGtex (Xena) | discovery | RNA-seq TPM (bulk) |
| recount3 (THCA tcga + THYROID gtex) | sensibilidade de contagem | RNA-seq counts |
| GSE33630 | validação primária | microarray GPL570 |
| GSE60542 | validação secundária | microarray GPL570 |

## 3. Análises NÃO executáveis (limitações técnicas)

| Análise | Motivo |
|---|---|
| Deconvolução formal (MuSiC/BisqueRNA/CIBERSORTx) | sem referência single-cell com músculo; CIBERSORTx exige Docker |
| Single-cell (Seurat/Scanpy) | memória limitada; Seurat não instalado |
| Proteína (TCGA RPPA) | extração via TCGAbiolinks/API não concluída |
| Metilação (minfi) | dados não baixados; sem pergunta específica fechada |
| Mutação (maftools) | MAF não baixado; paisagem citada de literatura com ressalva |

## 4. Resultados negativos (preservados)

- Genes musculares NÃO replicam nos GEO (efeito ≈ 0).
- Via "Cytoskeleton in muscle cells" INVERTE direção nos GEO.
- HALLMARK_P53_PATHWAY não significativa (discrepância com KEGG p53 registrada).
- Nenhuma coorte bulk RNA-seq independente adequada no GEO.

## 5. Métodos testados vs descartados

- DE: limma (TPM), limma-voom, DESeq2 — todos executados e comparados.
- Composição: marker-based scoring (usado); deconvolução formal (descartada — inviável).
- GSEA: fgsea (KEGG, Hallmark, Reactome) — executado.
- Meta-análise formal: NÃO executada (coortes heterogêneas microarray × RNA-seq; sem efeitos
  padronizados comparáveis justificáveis).
