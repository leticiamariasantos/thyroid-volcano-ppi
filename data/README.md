# Data Directory

Dados de entrada e proveniência do projeto (Fase 2 — reboot).

## data/global/

Matrizes globais de expressão (TCGA-THCA + GTEx tireoide) e o arquivo bruto:

| Arquivo | Conteúdo |
|---|---|
| `TcgaTargetGtex_rsem_gene_tpm.gz` | Matriz TOIL global (60.498 genes × 19.131 amostras, todos os tecidos) |
| `TCGA_GTEx_thyroid_tpm.tsv` | Subconjunto tireoide (58.581 genes × 783 amostras), log2(TPM+0.001) |
| `TCGA_GTEx_thyroid_counts.tsv` | Contagens recount3 (56.937 genes × 782 amostras) |
| `MANIFEST.tsv` / `MANIFEST_counts.tsv` | Manifestos de proveniência |
| `gencode.v23.annotation.gene.probemap.tsv` | Mapa de sondas GENCODE v23 |

## data/external/

Dados externos usados na validação, single-cell e gene sets:

- `GSE33630/`, `GSE60542/` — matrizes GEO (microarray GPL570)
- `GSE224356/` — listas de DEGs pré-computadas (xlsx)
- `GSE232237/` — single-cell de PTC (contagens e resultados marker-based)
- `GSE182416/` — single-cell de tireoide normal
- `reactome/ReactomePathways.gmt` — gene sets Reactome

## data/raw/

`XENA_THCA.tsv` — arquivo legado da Fase 1 (via direcionada hsa04919). Não é usado
pela Fase 2 (reboot).

## data/manifests/

Manifestos de proveniência da Fase 2 (`MANIFEST_TPM.tsv`, `MANIFEST_counts.tsv`).

## Download

O arquivo bruto global pode ser obtido em:
https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz
