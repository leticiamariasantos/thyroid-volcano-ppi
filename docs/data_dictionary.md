# Data Dictionary — thyroid-volcano-ppi

## Input Data

### `data/raw/XENA_THCA.tsv`

Expression matrix from UCSC Xena Browser combining TCGA THCA and GTEx Thyroid samples.

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier (e.g., `TCGA-BJ-A0YZ-01A`) |
| `_study` | Source study (`TCGA` or `GTEX`) |
| `TCGA_GTEX_main_category` | Tissue category (`TCGA Thyroid Carcinoma` or `GTEX Thyroid`) |
| `_sample_type` | Sample type (`Primary Tumor`, `Normal Tissue`, etc.) |
| `_primary_site` | Anatomical site (`Thyroid`) |
| `GENE_SYMBOL` ... | Expression values in **log₂(norm_count + 1)** |

---

## Output Tables

### `T01_sample_composition.tsv`
Sample counts by condition, study, and tissue type.

| Column | Description |
|--------|-------------|
| `condition` | `Normal` (GTEx) or `THCA` (TCGA) |
| `study` | Source study |
| `sample_type` | Tissue classification |
| `n` | Number of samples |

### `T02_deg_summary.tsv`
Differential expression analysis parameters and result counts.

| Column | Description |
|--------|-------------|
| `parameter` | Parameter name |
| `value` | Parameter value or count |

### `T03_deg_full_results.tsv`
Complete limma differential expression results for all tested genes.

| Column | Description |
|--------|-------------|
| `gene_symbol` | HGNC gene symbol |
| `logFC` | log₂ fold change (THCA vs Normal) |
| `AveExpr` | Average expression across all samples |
| `t` | Moderated t-statistic |
| `P.Value` | Raw p-value |
| `adj.P.Val` | Benjamini-Hochberg adjusted p-value |
| `B` | B-statistic (log-odds of differential expression) |
| `regulation` | `Up`, `Down`, or `NS` |
| `in_kegg` | Gene is in KEGG hsa04919 pathway (TRUE/FALSE) |

### `T04_top20_degs.tsv`
Top 20 differentially expressed genes by absolute log₂ fold change.

### `T05_kegg_missing_genes.tsv`
KEGG hsa04919 pathway genes not detected in the expression matrix.

### `T06_hub_proteins.tsv`
Hub proteins identified in the PPI network (top 25% betweenness + degree ≥ 2).

| Column | Description |
|--------|-------------|
| `gene_symbol` | HGNC gene symbol |
| `regulation` | Up/Down in THCA |
| `degree` | Number of direct PPI partners |
| `betweenness` | Normalized betweenness centrality |
| `closeness` | Normalized closeness centrality |
| `hub_score` | Kleinberg hub score |

### `T07_kegg_degs_ppi.tsv`
KEGG pathway DEGs with PPI network centrality metrics (integrated table).

---

## Network Metadata

### `N01_string_mapping.tsv`
Gene symbol → STRING protein ID mapping.

### `N02_string_interactions.tsv`
Filtered PPI edges (combined_score ≥ 700).

### `N03_centrality_metrics.tsv`
Complete centrality metrics for all network nodes.

### `N04_network_summary.tsv`
Global network statistics.
