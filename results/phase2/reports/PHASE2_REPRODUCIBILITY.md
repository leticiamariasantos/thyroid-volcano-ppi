# PHASE 2 (reboot) — Reprodutibilidade

## Ambiente

- **R:** 4.6.1 (2026-06-24, ucrt).
- **Pacotes:** versões registradas em `renv.lock`; sessão em `logs/session_info.txt` (a
  regenerar a cada execução).
- **Seed:** 42 (todos os scripts com aleatoriedade).
- **Caminhos:** `here::here()` (raiz = repositório); nenhum caminho absoluto.
- **Data da nova execução:** 2026-09-06 (registrada em `00_config.R` e nos relatórios).

## Execução

Scripts em `scripts/phase2/` (numerados 00–16), executados na ordem:

```
00_config.R 01_audit_data.R 02_qc.R 03_de.R 04_panel.R 05_gsea_global.R
06_gsea_panel.R 07_redundancy.R 08_composition.R 09_ppi.R 10_validation.R
11_singlecell.R 12_protein_mutation.R 13_robustness.R 14_convergence.R
15_figures.R 16_validate.R
```

## Dependências externas

- **KEGG:** `KEGGREST::keggGet` (gene sets do painel).
- **Reactome:** GMT oficial `ReactomePathways.gmt` (copiado para `data/external/reactome/`).
- **MSigDB:** `msigdbr` (KEGG_LEGACY, Reactome, Hallmark) — baixado on-the-fly.
- **STRING:** API REST (PPI, escore≥700).
- **cBioPortal:** API (mutações e RPPA, THCA PanCan Atlas 2018).
- **GEO:** matrizes GSE33630/GSE60542 preservadas em `data/external/`.

## Caching / reprodutibilidade determinística

- Matrizes intermediárias em `results/phase2/data/*.rds`.
- Sem dependência de arquivos temporários externos.
- `fgseaMultilevel` com `seed=42` (via `set.seed`).

## Validação técnica

`Rscript scripts/phase2/16_validate.R` — verifica existência, dimensões, consistência e
integridade dos outputs. Objetivo: **0 failures**.

## Rastreabilidade

- Resultados antigos removidos da camada analítica (diretórios 02–14 e `results/` antigos).
- Documentação histórica preservada em `documentation/` e `docs/` (rastreabilidade).
- `results/phase2/` contém exclusivamente outputs da nova fase.
