# Figure Specifications

## Cell Press Publication Standards

All figures follow the visual guidelines of Cell Press journals (Cell, Cell Reports, iScience, Med).

### Color Palette

| Element | Hex Code | Description |
|---------|----------|-------------|
| Upregulated | `#4477AA` | Medium blue — genes upregulated in THCA |
| Downregulated | `#AA4488` | Magenta — genes downregulated in THCA |
| Not Significant | `#C8C8C8` | Light grey — non-DE genes |
| PPI Edges | `#E0E0E0` | Very light grey — network edges |
| Highlight Border | `#222222` | Near-black — hub protein borders |

### Typography
- **Family**: Sans-serif (Arial/Helvetica)
- **Base size**: 14pt
- **Labels**: 5.2pt italic for gene symbols

### Dimensions
- **Volcano Plot**: 180 × 150 mm (double column)
- **PPI Network**: 180 × 180 mm (double column)
- **Resolution**: 600 dpi (PNG)

### Export Formats
- **PNG**: 600 dpi raster for publication
- **PDF**: Vector for editorial review
- **SVG**: Vector for figure editing

---

## Figure 1: Volcano Plot

**File**: `Fig1_Volcano_THCA_vs_Normal.{png,pdf,svg}`

**Content**:
- X-axis: log₂ fold change (THCA / Normal)
- Y-axis: -log₁₀ adjusted p-value
- Points: genes colored by regulation status (blue = upregulated, magenta = downregulated, grey = NS)
- Threshold lines: |log₂FC| = 1.0 (vertical dashed), FDR = 0.05 (horizontal dashed)
- Labels: top 15 DEGs by |log₂FC| + top 5 KEGG pathway DEGs by significance
- Highlight: KEGG hsa04919 pathway genes with open ring
- Legend: inset (upper-left), regulation colors only

**Label selection criteria**:
1. Top 15 DEGs ranked by absolute log₂ fold change
2. Union with top 5 KEGG pathway DEGs ranked by adjusted p-value
3. `ggrepel` with seed = 42 for reproducible label placement

---

## Figure 2: PPI Network

**File**: `Fig2_PPI_Network_THCA_DEGs.{png,pdf,svg}`

**Content**:
- Layout: Fruchterman-Reingold (force-directed)
- Nodes: filled circles, sized by degree centrality
- Node colors: regulation-based (vermilion = Up, blue = Down)
- Hub nodes: thick dark border highlight (betweenness ≥ 70th percentile + degree ≥ 2, or degree ≥ 3)
- Edges: thin light grey, width proportional to STRING combined_score
- Labels: hub proteins + KEGG pathway genes only
- Legend: inset annotation (upper-left)

**Hub definition**:
- Betweenness centrality ≥ 75th percentile of non-zero betweenness values
- Degree ≥ 2
- Both criteria must be met

**Label strategy**:
- All hub proteins
- All KEGG hsa04919 pathway genes present in the network
- Union of both sets → `ggraph::geom_node_text` with repel
