# LITERATURE READING GUIDE — Mapa de leitura progressiva

**Regra:** referências com PMID/DOI foram verificadas via NCBI E-utilities. Itens UNVERIFIED devem
ser confirmados antes de uso.

---

## BLOCO 1 — Carcinoma de tireoide e biologia molecular
**Por quê:** entender a doença e a paisagem molecular (MAPK/BRAF/RAS).
- Integrated genomic characterization of papillary thyroid carcinoma (TCGA). Cell 2014. PMID 25417114. [NÍVEL 1]
- Conceitos: PTC, BRAF V600E, RAS, fusões RET/NTRK.

## BLOCO 2 — Transcriptômica e expressão diferencial
**Por quê:** justificar limma/voom/DESeq2 e a interpretação de DEGs.
- Ritchie et al. limma. NAR 2015. PMID 25605792. [NÍVEL 1]
- Love et al. DESeq2. Genome Biol 2014. PMID 25516281. [NÍVEL 1]
- Law et al. voom. Genome Biol 2014. PMID 24485249. [NÍVEL 1]
- Conceitos: normalização, modelagem de dispersão, FDR, múltiplos testes.

## BLOCO 3 — TCGA/GTEx e confounding
**Por quê:** fundamentar a limitação estrutural do desenho discovery.
- Leek et al. Batch effects. Nat Rev Genet 2010. PMID 20838408. [NÍVEL 1]
- GTEx Consortium atlas (Science 2015) — UNVERIFIED (confirmar).
- Conceitos: source-condition confounding, triangulação, não "correção mágica".

## BLOCO 4 — GSEA e análise de vias
**Por quê:** justificar GSEA e a interpretação cautelosa de vias.
- Subramanian et al. GSEA. PNAS 2005. PMID 16199517. [NÍVEL 1]
- fgsea (Korotkevich 2021, bioRxiv) — UNVERIFIED.
- Conceitos: NES, FDR, redundância, "via ativa por fosforilação sem aumento de mRNA".

## BLOCO 5 — Bulk RNA-seq e composição celular
**Por quê:** fundamentar a interpretação composicional (MYH7 etc.).
- Newman et al. CIBERSORT. Nat Methods 2015. PMID 25822800. [NÍVEL 1]
- Wang et al. MuSiC. Nat Commun 2019. PMID 30670690. [NÍVEL 1]
- tumor purity prediction. Brief Bioinform 2021. PMID 33954576. [NÍVEL 2]
- Conceitos: bulk ≠ célula tumoral; deconvolução; contaminação.

## BLOCO 6 — Single-cell RNA-seq e pseudobulk
**Por quê:** justificar a localização celular e evitar pseudorreplicação.
- Seurat (Hao 2021 Cell) — UNVERIFIED.
- Conceitos: cell-type annotation por marcadores (limitação), pseudobulk por paciente.

## BLOCO 7 — ITGA2 e integrinas (NÚCLEO)
**Por quê:** o candidato central.
- Desgrosellier & Cheresh. Integrins in cancer. Nat Rev Cancer 2010. PMID 20029421. [NÍVEL 1]
- Seguin et al. Integrins and cancer. Trends Cell Biol 2015. PMID 25572304. [NÍVEL 1]
- FAK signaling. J Cell Sci 2024. PMID 39034922. [NÍVEL 2]
- ITGA2 single-cell PTC. Adv Sci 2025. PMID 40985182. [NÍVEL 1]
- ITGA2 knockdown → piroptose tireoide. Front Biosci 2025. PMID 40917053. [NÍVEL 2]
- ITGA2 hipometilação → invasão tireoide. 2023. PMID 36420922. [NÍVEL 2]
- Integrin α2 cascata metastática. J Transl Med 2026. PMID 41803969. [evidência a ponderar]
- Conceitos: α2β1 receptor de colágeno; membrana; FAK/SRC; papel dual.

## BLOCO 8 — FN1 e microambiente
**Por quê:** interpretar FN1 como ECM/microambiente (tumor + estroma).
- (revisões de FN1/CAF — complementar; buscar específico conforme lacuna)
- Conceitos: fibronectina, CAF, remodelamento estromal, EMT.

## BLOCO 9 — CCND1 e proliferação
**Por quê:** interpretar CCND1 como marcador de proliferação (intracelular).
- (ciclo celular/ciclina D1 em tireoide — complementar)
- Conceitos: ciclina D1, ciclo celular, sinalização oncogênica.

## BLOCO 10 — RNA-proteína e RPPA
**Por quê:** fundamentar que proteína aumenta plausibilidade, não resolve célula de origem.
- Proteomic Analysis of TCGA (RPPA). Int J Mol Sci 2026. PMID 42589374. [NÍVEL 3]
- Conceitos: correlação RNA-proteína limitada; RPPA tumor-only.

## BLOCO 11 — Mutação/CNV
**Por quê:** fundamentar "não identificada mutação/CNV que explique a expressão".
- GISTIC2.0 (Mermel 2011) — UNVERIFIED.
- Conceitos: ausência de mutação ≠ origem regulatória comprovada.

## BLOCO 12 — Heterogeneidade e biomarcadores
**Por quê:** média de coorte ≠ consistência individual.
- (heterogeneidade intratumoral — complementar)
- Conceitos: heterogeneidade inter/intratumoral; alvo heterogêneo.

## BLOCO 13 — Targeting de superfície
**Por quê:** distinguir "proteína de membrana" de "alvo de superfície" de "alvo validado".
- Integrins in cancer (Nat Rev Cancer 2010) — reusar.
- Conceitos: acessibilidade, internalização, especificidade, expressão normal (off-tumor).

## BLOCO 14 — Targeted delivery e nanomedicina (SÓ DEPOIS)
**Por quê:** fundamentar a possível consequência translacional, sem premissa nanomédica.
- (revisões de active targeting/nanopartículas — complementar)
- Conceitos: active targeting, receptor-mediated delivery, ligante.

## BLOCO 15 — Reproducibilidade, robustez e falsificação
**Por quê:** justificar a auditoria de falsificação e a análise de sensibilidade.
- Breznau et al. Observing many researchers. PNAS 2022. PMID 36306328. [NÍVEL 1]
- Conceitos: graus de liberdade do analista, robustez, falsificação.

---

## 10 artigos para ler primeiro (resumo executivo)

1. TCGA THCA (Cell 2014) — a paisagem molecular do PTC.
2. Ritchie limma (2015) — o método DE principal.
3. Love DESeq2 (2014) — DE baseada em counts.
4. Law voom (2014) — voom.
5. Leek batch effects (2010) — o confounding.
6. Subramanian GSEA (2005) — enriquecimento.
7. Desgrosellier Integrins in cancer (2010) — integrina/targeting.
8. ITGA2 single-cell PTC (Adv Sci 2025, PMID 40985182) — a evidência central.
9. ITGA2 knockdown tireoide (2025, PMID 40917053) — funcional.
10. Breznau Observing many researchers (2022) — robustez/falsificação.
