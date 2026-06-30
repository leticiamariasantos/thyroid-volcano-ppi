# Data Directory

## data/raw/

Place the `XENA_THCA.tsv` file here before running the pipeline.

### Download

```bash
Rscript scripts/download_data.R
```

Or manually from Xena Browser:
1. Go to: https://xenabrowser.net/?bookmark=c486b845ee2e750c3a9d2fc5145c8426
2. Click "Download" → "Download current visualization data"
3. Save as: `data/raw/XENA_THCA.tsv`

### Format

- **Source**: UCSC Xena Browser (TCGA THCA + GTEx Thyroid)
- **Values**: log₂(norm_count + 1)
- **Rows**: samples
- **Columns**: sample metadata (1-5) + gene symbols (6+)

## data/processed/

Intermediate files generated during pipeline execution. Not versioned.

## data/string_cache/

STRING database cache files. Auto-populated on first run. Not versioned.
