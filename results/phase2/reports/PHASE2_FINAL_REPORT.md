# PHASE 2 (reboot) — Relatório Final (Resumo Executivo)

**Data da execução:** 2026-09-06
**Escopo:** reconstrução independente do pipeline (dados → DE → painel 30 vias → GSEA global
→ composição → PPI → validação externa → single-cell → proteína/mutação) com painel
pré-especificado de 30 vias.

## Resumo executivo (16 itens)

1. **Nº final de DEGs:** limma (principal) 12.200 DEGs de 26.011 testados (2.485 Up,
   9.715 Down, |log2FC|>1 & FDR<0,05); voom 19.485; DESeq2 17.737.

2. **Principais alterações:** superexpressão tumoral de programas de proliferação/tradução,
   imunidade, proteostase e metabolismo mitocondrial; assinatura muscular **Down** no tumor
   (composição do normal GTEx, não regulação tumoral).

3. **Resultado das 30 vias (painel):** 6 vias significativas no limma, todas Up no tumor —
   Proteasome, Antigen processing, OXPHOS, DNA replication (ADICIONAIS) e p53, Cell cycle
   (ORIGINAIS). As demais (MAPK, PI3K-Akt, mTOR, Wnt, NF-kB, Thyroid cancer, Thyroid hormone,
   Apoptosis, adesão/ECM) não atingiram FDR<0,05 no método principal.

4. **GSEA global:** emergem Ribosome/tradução, Proteasome, resposta imune/inflamatória
   (graft-vs-host, antigen processing), Cell cycle/DNA replication, OXPHOS e p53 (Up);
   Myogenesis/cardiomiopatia dilatada (Down, músculo do normal).

5. **Vias robustas após controle de composição:** as mesmas 6 vias (Proteasome, Antigen
   processing, OXPHOS, DNA replication, p53, Cell cycle) permanecem robustas após remoção
   dos 48 marcadores musculares — nenhuma via significativa desapareceu.

6. **Principais módulos de PPI:** rede STRING reconstruída excluindo marcadores musculares
   (70 nós / 65 arestas, 18 comunidades); hub principal FN1 (grau 9), com CCND1 e genes de
   ciclo celular/EMT; centralidade é propriedade topológica, não alvo terapêutico. PRKCA não
   é mais hub (era hub apenas na Fase 1 direcionada à via hsa04919).

7. **ITGA2:** Up no tumor (logFC +2,47, FDR 4,1e-89), consistente nos 3 métodos, replicado
   externamente (GSE33630/GSE60542/GSE224356), predominantemente epitelial/tumoral no
   single-cell. NÃO é core enrichment das vias significativas (adesão/ECM é NS no limma);
   sem mutação; fora do RPPA.

8. **FN1:** Up (logFC +4,32, FDR 1,8e-95), consistente e replicado externamente; epitelial
   E fibroblástico (estroma) — sem exclusividade tumoral. NÃO é core enrichment.

9. **CCND1:** Up (logFC +2,26, FDR 2,7e-144), consistente e replicado externamente; epitelial;
   RPPA presente mas correlação RNA-proteína fraca (0,051). Entre os três candidatos, é o
   único que figura no core enrichment — restrito às vias ORIGINAIS (Cell cycle, p53) e a
   Cellular senescence/Thyroid cancer; não figura no core enrichment das vias adicionais.

10. **Validação externa:** ITGA2/FN1/CCND1 **Up** em GSE33630 (+2,00/+2,88/+1,42),
    GSE60542 (+2,50/+3,06/+1,07) e GSE224356 (T1N1/T2N2/T3N3).

11. **Evidência single-cell:** ITGA2 e CCND1 predominantemente epiteliais/tumorais; FN1
    epitelial e fibroblástico (classificação marker-based).

12. **Evidência proteica:** CCND1 no painel RPPA (correlação RNA-proteína 0,051); ITGA2 e
    FN1 fora do painel RPPA (não testáveis).

13. **Principais limitações:** `source≡condition` (TCGA/GTEx confundidos); composição
    muscular do normal; PPI esparso; RPPA tumor-only; TERT promotor não captado; single-cell
    marker-based; estudo exploratório.

14. **Principal hipótese translacional:** ITGA2 priorizado como **candidato translacional
    para investigação de direcionamento molecular** (não alvo validado) — superexpressão
    tumoral replicada, localização epitelial/tumoral, mas sem evidência de especificidade
    absoluta, acessibilidade, internalização ou eficácia.

15. **Classificação dos resultados:**
    - **Confirmatórios (hipótese a priori):** p53 e Cell cycle.
    - **Robustos (pré-especificados + tripla validação):** Proteasome, Antigen processing,
      OXPHOS, DNA replication.
    - **Convergência molecular do painel (vias ORIGINAIS):** eixo ciclo celular/p53 (TP53,
      CCND1, CDK1, CDKN1A, CCNB1/2, CCND2, CDK4, CDKN2A, CHEK1, MDM2, SFN) + replicação
      (MCM2–MCM5, PCNA) + PSME1 — 18 genes no leading edge de ≥2 vias robustas. As vias
      adicionais convergem para proteassoma/apresentação de antígeno/OXPHOS/replicação de DNA.
    - **Exploratórios:** demais vias adicionais (Interferon, Focal adhesion, Hippo,
      Ferroptose, Autophagy, etc. — método-dependentes).

16. **Questões em aberto:** ITGA2/FN1/CCND1 não têm alteração genômica que explique a
    expressão (mecanismo regulatório desconhecido); ausência de proteômica para ITGA2/FN1;
    acessibilidade e internalização de ITGA2 não testadas; generalização a outros subtipos
    histológicos não avaliada.

## Conclusão

O painel pré-especificado de 30 vias converge para o **eixo ciclo
celular/p53** (TP53, CDK1, CDKN1A, CDKN2A, MDM2, CHEK1 e a maquinaria MCM/PCNA). **CCND1 é
o único dos três candidatos que figura no core enrichment — restrito às vias ORIGINAIS (Cell
cycle, p53)** e a Cellular senescence/Thyroid cancer; não figura no core enrichment das 20
vias adicionais. **ITGA2 e FN1 não figuram no core enrichment** — emergem como DEGs
individuais robustos e replicados dentro de programas de adesão/ECM pré-especificados. Os
programas tumorais dominantes são proliferação, imunidade, proteostase e metabolismo
mitocondrial. O estudo permanece **exploratório e gerador de hipóteses**; ITGA2 é reportado
como candidato translacional, não como alvo validado.
