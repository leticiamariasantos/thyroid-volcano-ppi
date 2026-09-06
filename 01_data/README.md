# Dados (01_data) — manifestos

Os dados brutos/processados estão em `data/global/` e `data/raw/`. Esta pasta mantém
os manifestos de proveniência.

| Arquivo | Conteúdo |
|---|---|
| `MANIFEST_TPM.tsv` | Matriz TPM (TOIL Xena) — 58.581 genes × 783 amostras, log2(TPM+0.001) |
| `MANIFEST_counts.tsv` | Matriz de contagens (recount3) — 56.937 genes × 782 amostras, raw counts |

**Arquivos de dados reais:** `data/global/TCGA_GTEx_thyroid_tpm.tsv` (352 MB),
`data/global/TCGA_GTEx_thyroid_counts.tsv` (184 MB),
`data/raw/XENA_THCA.tsv` (Fase 1, 121 genes).
