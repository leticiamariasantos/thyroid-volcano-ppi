# LITERATURE TO PROJECT MAP — Mapa entre artigos e partes do projeto

Cada etapa do pipeline é justificada por referência primária verificada (PMID/DOI). UNVERIFIED =
confirmar antes de uso.

| Etapa do projeto | Por que o método é apropriado | Referência |
|---|---|---|
| DATA (TCGA-THCA + GTEx) | paisagem molecular do PTC | TCGA THCA, Cell 2014, PMID 25417114 |
| DATA (GTEx) | atlas de expressão em tecidos normais | GTEx Consortium (UNVERIFIED — confirmar) |
| QC / normalização | transformação log, filtragem | Ritchie limma 2015, PMID 25605792 |
| DE (limma TPM) | modelagem linear moderada | Ritchie limma 2015, PMID 25605792 |
| DE (voom) | pesos de precisão p/ counts | Law voom 2014, PMID 24485249 |
| DE (DESeq2) | modelagem de dispersão p/ counts | Love DESeq2 2014, PMID 25516281 |
| Confounding TCGA×GTEx | source≡condition, batch | Leek 2010, PMID 20838408 |
| GSEA (KEGG/Hallmark/Reactome) | enriquecimento pré-ranqueado | Subramanian GSEA 2005, PMID 16199517 |
| Composição celular (bulk) | bulk reflete composição | CIBERSORT 2015 (PMID 25822800); MuSiC 2019 (PMID 30670690) |
| Tumor purity | expressão depende da pureza | tumor purity 2021, PMID 33954576 |
| Single-cell (GSE182416/GSE232237) | localização celular | Seurat (UNVERIFIED); ITGA2 single-cell PTC 2025 (PMID 40985182) |
| Pseudobulk | evitar pseudorreplicação | (pseudobulk — complementar/UNVERIFIED) |
| Proteína (RPPA) | evidência proteica (tumor-only) | RPPA TCGA 2026, PMID 42589374 |
| Mutação/CNV | ausência de mutação/CNV explicativa | GISTIC (UNVERIFIED) |
| ITGA2 biologia | receptor α2β1 de colágeno, membrana | Desgrosellier 2010 (PMID 20029421); Seguin 2015 (PMID 25572304) |
| ITGA2 funcional | knockdown/invasão em tireoide | PMID 40917053; PMID 36420922 |
| FN1 (ECM) | fibronectina/microambiente | (CAF/FN1 — complementar) |
| CCND1 (proliferação) | ciclina D1 | (ciclina D1 tireoide — complementar) |
| Targeting de superfície | critérios: acessibilidade, internalização, especificidade | Desgrosellier 2010 (PMID 20029421) |
| Nanomedicina | consequência translacional (não premissa) | (active targeting — complementar) |
| Falsificação/robustez | graus de liberdade do analista | Breznau 2022, PMID 36306328 |
