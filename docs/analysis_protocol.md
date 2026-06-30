# Analysis Protocol

## THCA Transcriptomic Analysis — Volcano Plot & PPI Network

> **Study type:** Exploratory, hypothesis-generating
> **Scope:** Differential expression + PPI network of KEGG hsa04919 pathway genes
> **Limitation:** TCGA vs GTEx comparison without batch effect correction

---

## 1. Data Source

**UCSC Xena Browser**  
TCGA THCA + GTEx Thyroid gene expression RNA-seq (TOIL recompute pipeline)

- URL: https://xenabrowser.net/?bookmark=c486b845ee2e750c3a9d2fc5145c8426
- Gene expression values: log₂(RSEM expected_count + 1), upper quartile normalized
- Samples: TCGA Thyroid Carcinoma (THCA, n = 504) + GTEx normal thyroid (n = 279)
- Total: 783 samples
- Download date: June 2026

**Reference**: Goldman MJ et al. (2020) *Nature Biotechnology* 38:675–678.

### 1.1 Known Limitations

⚠️ **TCGA-GTEx batch effects:** Tumor and normal samples come from different cohorts with distinct sequencing protocols, demographic profiles, and processing pipelines. No batch effect correction (e.g., ComBat, RUVseq) is applied. Observed differential expression may partially reflect technical rather than biological variation. This limitation is explicitly stated in the manuscript and README.

⚠️ **Cellular composition:** Expression differences may reflect differences in cell-type composition between tumor and normal tissue rather than transcriptional regulation per se. Deconvolution analysis is not performed.

---

## 2. Data Preprocessing

### 2.1 Column standardization
- Remove Xena underscore prefix (`_`) from metadata columns
- Drop redundant `samples` column if present

### 2.2 Condition assignment
- `TCGA_GTEX_main_category == "GTEX Thyroid"` → **Normal**
- `TCGA_GTEX_main_category == "TCGA Thyroid Carcinoma"` → **THCA**
- Samples with undefined condition excluded

### 2.3 Expression scale validation
- Assert: max(expression) ≤ 30 (consistent with log₂ scale)
- Check: non-integer fraction ≥ 1% (excludes raw count data)

### 2.4 Low-expression filter
- Gene retained if mean expression > 0.5 in ≥ 10% of samples
- Applied before differential expression to reduce multiple testing burden

---

## 3. Differential Expression Analysis

### 3.1 Method
- **Package**: `limma` v3.56+
- **Design**: ~0 + condition (intercept-free model matrix)
- **Contrast**: THCA − Normal
- **Moderation**: Empirical Bayes (`eBayes`)

### 3.2 Multiple testing correction
- Benjamini-Hochberg (BH) procedure
- Adjusted p-value threshold: **FDR < 0.05**

### 3.3 DEG classification
- **Upregulated in THCA**: log₂FC > 1.0 AND adj.P.Val < 0.05
- **Downregulated in THCA**: log₂FC < −1.0 AND adj.P.Val < 0.05
- **Not significant**: all other genes

---

## 4. Volcano Plot

### 4.1 Visualization
- All tested genes plotted; DEGs color-coded
- Threshold lines: |log₂FC| = 1.0, FDR = 0.05
- Gene labels: objective selection criteria (see `docs/figure_specs.md`)
- KEGG hsa04919 pathway genes: open ring highlight

### 4.2 Output
- 600 dpi PNG, PDF, SVG at 180 × 150 mm

---

## 5. PPI Network Analysis

### 5.1 STRING database
- **Version**: 12.0
- **Species**: Homo sapiens (NCBI taxon 9606)
- **Interaction score**: combined_score ≥ 700 (high confidence)

### 5.2 Gene-to-protein mapping
- All DEGs (Up + Down) submitted to STRING REST API
- Gene symbols mapped to STRING protein IDs (`9606.ENSP...`)
- Unmapped genes excluded (reported in mapping table)

### 5.3 Network construction
- **Edges**: STRING interactions with combined_score ≥ 700
- **Graph**: undirected, weighted (igraph)
- **Component**: largest connected component retained
- **Isolated nodes**: removed from visualization

### 5.4 Centrality analysis
- **Degree**: number of direct interaction partners
- **Betweenness**: normalized betweenness centrality
- **Closeness**: normalized closeness centrality
- **Hub score**: Kleinberg hub centrality

### 5.5 Hub identification
- Betweenness ≥ 75th percentile of non-zero values
- Degree ≥ 2
- Hub status used for node highlighting and labeling

### 5.6 Visualization (ggraph)
- Layout: `nicely` (optimal force-directed)
- Nodes: fill = regulation color, size ∝ degree, border = hub status
- Edges: thin grey, width ∝ combined_score
- Labels: hubs + KEGG pathway genes

### 5.7 Output
- 600 dpi PNG, PDF, SVG at 180 × 165 mm

---

## 6. KEGG Pathway Annotation

**Pathway**: hsa04919 — Thyroid hormone signaling pathway  
**Source**: KEGG REST API via `KEGGREST` package

### Nota sobre a contagem de genes

O website do KEGG (https://www.kegg.jp/entry/hsa04919) exibe **78 genes** no diagrama da via — são os genes manualmente curados com posição explícita na ilustração. A função `keggGet("hsa04919")` do pacote KEGGREST retorna **121 genes**, porque a API inclui todos os genes humanos que compartilham KEGG Orthology (KO) com os genes do diagrama. Isso significa que parálogos e membros de famílias gênicas são incluídos mesmo que não apareçam individualmente no desenho da via.

**Exemplo:** O diagrama utiliza "PRKC" (proteína quinase C). A API retorna todos os membros humanos dessa família: PRKCA, PRKCB, PRKCG, PRKCD, PRKCE, PRKCH, PRKCQ, PRKCI, PRKCZ. O mesmo ocorre com outras famílias como PLC (fosfolipase C), ATP1A (ATPase), NOTCH, MED (mediador), etc.

Dos 121 genes retornados pela API, **119 intersectam** com a matriz de expressão TCGA-GTEx após o pré-processamento (2 são excluídos pelo filtro de baixa expressão: expressão média > 0,5 em ≥ 10% das amostras).

Esta abordagem é **conservadora e recomendada** em estudos transcriptômicos baseados em vias KEGG, pois:
- Evita a perda de genes funcionalmente relacionados que participam da via mas não foram desenhados no diagrama por limitações de espaço
- Garante cobertura completa da via no contexto de uma análise exploratória
- A documentação explícita da discrepância permite que o leitor compreenda a origem e a magnitude da diferença

### Uso dos genes KEGG neste estudo

KEGG pathway genes are used for:
1. Volcano plot highlight (open rings on DEGs in pathway)
2. PPI network label priority
3. Context-specific gene tables

KEGG genes are **not** used to filter DEGs; all genes passing the statistical thresholds are included in the analysis.

---

## 7. Reproducibility

| Measure | Implementation |
|---------|---------------|
| Random seed | `set.seed(42)` in `R/00_setup.R` |
| Parameter centralization | All thresholds in `THRESHOLD` list |
| Path portability | `here::here()` throughout |
| Session capture | `sessionInfo()` logged to `logs/` |
| Pipeline modularity | Independent R scripts in `R/` |
| Dependency documentation | Package list in `R/00_setup.R` and README |

---

## 8. Software Versions

**Development environment:** R v4.6.0 (2026-04-26, "Jitterbug Madness") executed via Visual Studio Code v1.125 with R extension on Windows 11 x64.

See `logs/session_info.txt` for complete version information after each run.

Key dependencies:
- R ≥ 4.1 (tested on v4.6.0)
- limma ≥ 3.56
- ggplot2 ≥ 3.5
- ggraph ≥ 2.1
- igraph ≥ 1.5
- STRING v12.0 (REST API)
- KEGG (REST API via KEGGREST)

## 9. Audit Trail

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-21 | v3.0.0 | Initial pipeline with volcano + PPI |
| 2026-06-24 | v3.1.0 | Scientific audit: hypothesis reformulation, limitations documentation, PRKCA focus |

Full audit trail: `results/tables/S3_pipeline_audit_trail.tsv`
