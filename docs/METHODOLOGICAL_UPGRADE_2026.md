# Upgrade metodológico 2026 — protocolo a priori

**Versão-alvo:** 5.0.0  
**Data de congelamento:** 2026-09-08  
**Seed global:** 42  
**Natureza:** estudo exploratório, gerador de hipóteses.

Este documento foi criado antes da execução das análises adicionadas pelos scripts
`20_`–`30_`. Ele registra decisões, limiares, análises primárias, sensibilidades e
critérios de interrupção. O painel de 30 vias permanece exatamente como definido em
`results/phase2/pathways/PANEL_SELECTION_RATIONALE.md`.

## 1. Inventário inicial

| Item | Estado em 4.1.0 | Consequência para o upgrade |
|---|---|---|
| Git | `main` limpo em `d861258` | Base histórica preservada; mudanças novas em commits lógicos |
| Expressão TOIL | 504 TCGA + 279 GTEx | Todas as 504 TCGA são código 01; não há normal adjacente |
| Contagens recount3 | 504 TCGA + 278 GTEx | Apenas amostras alinhadas à matriz TOIL foram retidas |
| Confundimento | `source` é idêntico a `condition` | Batch TCGA×GTEx não é identificável no desenho atual |
| DE | limma TPM, voom e DESeq2 | Acrescentar edgeR-QL e voom com pesos/qualidade |
| GSEA | fgsea; Hallmark, KEGG, Reactome e painel | Acrescentar GO-BP, clusterProfiler, GSVA e robustez |
| Composição | remoção/score muscular | Substituir por deconvolução e ajuste multiassinatura |
| PPI | top 500, STRING combinado, Walktrap | Separar universos/evidências e adicionar Louvain/Leiden |
| Validação | GEO, single-cell, RPPA e mutação | Acrescentar meta-análise, sobrevivência e superfície |
| Orquestração | sequencial, sempre reexecuta | Acrescentar checksum, cache e retomada |

## 2. Hierarquia de comparações

1. **Primária:** tumor primário TCGA-THCA versus tecido sólido normal adjacente do
   TCGA-THCA, processados de modo uniforme no recount3. Quando houver tumor e normal
   do mesmo participante, o pareamento será incorporado ao desenho.
2. **Sensibilidade externa bruta:** tumor TCGA versus tireoide normal GTEx no TOIL/recount3.
3. **Sensibilidade corrigida por source:** somente no conjunto ampliado que contém
   normais TCGA, tornando `source` e `condition` parcialmente separáveis.
4. **Sensibilidade ajustada por composição:** desenho anterior acrescido das proporções
   celulares estimadas ou de componentes composicionais pré-especificados.

O contraste TCGA×GTEx sem normais TCGA nunca será apresentado como livre de batch.
Se a matriz do desenho tiver posto incompleto, a etapa deverá parar com diagnóstico,
em vez de produzir um resultado rotulado como corrigido.

## 3. Pré-processamento e batch

- Contagens: genes com CPM ≥ 1 no mínimo em 25% do menor grupo; TMM para edgeR/voom.
- Escala contínua: log2-CPM com offset e normalização TMM.
- Correção principal: ComBat-seq nas contagens quando o desenho `~ condition + source`
  for identificável; `condition` é preservada no grupo de interesse.
- Sensibilidades: SVA no modelo de DE; RUVg somente se controles negativos empíricos
  forem definidos sem usar os candidatos; `removeBatchEffect` apenas para visualização
  e análises baseadas em escala contínua, nunca como substituto do modelo de contagens.
- Risco residual: ComBat/SVA/RUV podem remover sinal biológico real; concordância com o
  contraste TCGA-matched será reportada e nenhuma versão será declarada “verdade”.

## 4. Deconvolução e ajuste de composição

- Métodos-alvo: xCell e EPIC por `immunedeconv`; resultados são cacheados com versão,
  checksum da entrada e parâmetros. CIBERSORTx/BayesPrism não são obrigatórios porque
  exigem credenciais, execução externa ou referência single-cell validada; sua ausência
  deve permanecer explícita.
- Compartimentos mínimos: epitelial/tumoral, fibroblasto, endotélio, imune e músculo.
- Consenso: mediana de escores padronizados dentro de cada compartimento; resultados
  de métodos com escalas incompatíveis não serão tratados como frações absolutas.
- Ajuste: covariáveis composicionais não colineares no desenho de DE. A matriz
  “composition-adjusted” conterá resíduos com o efeito de condição preservado, para
  visualização/GSEA; testes de DE serão feitos no modelo, não sobre resíduos isolados.
- Contribuição à variância: diferença de R² parcial entre modelos com condição e com
  condição + composição para os top DEGs, limitada ao intervalo 0–100%.

## 5. DE, GSEA e robustez

- Métodos: limma, DESeq2, voom, edgeR quasi-likelihood e voomWithQualityWeights.
- Ranking principal de GSEA: estatística t do voom/voom-QW; DESeq2 `stat` é sensibilidade.
- Coleções: Hallmark, KEGG, Reactome, GO Biological Process e painel congelado de 30 vias.
- Engines: fgsea (primária), clusterProfiler (replicação de implementação) e GSVA
  (scores por amostra).
- Leading edge: genes recorrentes e overlap com marcadores/composição serão reportados.
- Robustez: 100 bootstraps estratificados (configurável para CI) e repetição sem outliers
  previamente definidos por QC; uma via é robusta quando mantém direção em ≥80% e
  FDR < 0,05 em ≥70% das reamostragens.

## 6. Validação e critérios de ITGA2

Meta-análise de efeitos aleatórios (REML) para ITGA2, FN1, CCND1 e os 18 genes de
convergência. Efeitos serão log2FC com erro-padrão; estudos sem erro-padrão não serão
combinados quantitativamente e aparecerão como evidência direcional.

ITGA2 só poderá receber a expressão **“candidato prioritário para investigação de
surface targeting”** se satisfizer simultaneamente:

1. direção positiva em ≥80% dos datasets elegíveis e nenhuma inversão significativa;
2. meta-log2FC ≥ 1,0, FDR da meta-análise < 0,05 e heterogeneidade I² ≤ 75%;
3. anotação de membrana plasmática em fonte versionada e evidência proteica;
4. log2FC ajustado por composição ≥ 1,0, FDR < 0,05 e retenção ≥70% do efeito bruto;
5. associação de sobrevivência, se reportada, não contraditória após ajuste clínico;
6. nenhum critério de segurança/especificidade será inferido de centralidade PPI.

Se qualquer critério obrigatório não puder ser avaliado ou falhar, a classificação será
**“DEG robusto de interesse”**. Em nenhuma situação serão usados “alvo terapêutico
validado”, “biomarcador clínico validado” ou “nanomedicina pronta”.

## 7. PPI

Serão construídas três redes: top 500 genome-wide; genes do painel de 30 vias; e genes
de leading edge das vias robustas. Arestas exigirão evidência experimental STRING e
correlação de expressão no conjunto analisado (|rho| ≥ 0,30; FDR < 0,05). Grau,
betweenness, closeness e eigenvector são descritores topológicos. Walktrap, Louvain e
Leiden serão comparados; módulos receberão enriquecimento funcional com universo explícito.

## 8. Falhas, cache e rastreabilidade

- Respostas de STRING, cBioPortal, MSigDB/GEO e demais APIs serão armazenadas em
  `data/cache/<fonte>/<versão>/` com metadados e checksum.
- O orquestrador só pula uma etapa quando checksum de script + inputs coincide e todos
  os outputs declarados existem e são não vazios.
- Falha de API usa cache válido; sem cache, a etapa falha com mensagem acionável.
- Resultados 4.1.0 permanecem disponíveis; novos resultados usam subdiretórios
  `raw`, `batch_corrected`, `composition_adjusted` e `tcga_matched`.
- `documentation/` e `scripts/legacy/` são históricos e não serão reescritos.

## 9. Limitações residuais previstas

Normais adjacentes não são tecido saudável perfeito, têm tamanho amostral menor e podem
conter field effects. Deconvolução depende de assinaturas e não prova origem celular.
Correção de batch sob confundimento parcial pode sobrecorrigir. Meta-análise combina
plataformas e coortes heterogêneas. Sobrevivência em THCA tem poucos eventos. Dados de
superfície não demonstram internalização, especificidade tumoral ou janela terapêutica.
