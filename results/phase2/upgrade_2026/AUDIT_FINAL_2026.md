# Auditoria Final — thyroid-volcano-ppi v5.0.0

**Data de conclusão:** 2026-09-12 (America/Sao_Paulo)
**Protocolo congelado:** 2026-09-08 (docs/METHODOLOGICAL_UPGRADE_2026.md)
**Resultado:** `Phase 2 + upgrade 2026 complete (0 failures).` — validação técnica `0 failure(s)`.

---

## 1. Status das etapas 20–31

| Etapa | Script | Status | Observação |
|---|---|---|---|
| 20 | 20_upgrade_config.R | ✅ concluída | congela configuração 5.0.0 |
| 21 | 21_tcga_adjacent_normal.R | ✅ concluída | 505 tumores + 59 normais adjacentes (recount3) |
| 22 | 22_batch_correction.R | ✅ concluída | ComBat-seq + SVA |
| 23 | 23_deconvolution.R | ✅ concluída | xCell + EPIC + consenso (51 não-convergências EPIC divulgadas) |
| 24 | 24_multiverse_de.R | ✅ concluída | 5 variantes × 5 métodos; gate de equivalência formal `passed=TRUE` |
| 25 | 25_multiverse_gsea.R | ✅ concluída | 68.784 resultados fgsea + clusterProfiler + GSVA |
| 26 | 26_composition_variance.R | ✅ concluída | R² parcial de 100 top DEGs |
| 27 | 27_validation_meta_survival_surface.R | ✅ concluída | meta REML + OS/PFS/DSS + superfície + gate ITGA2 |
| 28 | 28_ppi_multinetwork.R | ✅ concluída | 3 redes PPI (STRING experimental + coexpressão) |
| 29 | 29_gsea_robustness.R | ✅ concluída | 100 bootstraps; 14 vias robustas |
| 30 | 30_render_upgrade_report.R | ✅ concluída | HTML self-contained com figuras embutidas |
| 31 | 31_validate_upgrade.R | ✅ concluída | 42/42 checks `TRUE` |

### Resultados-chave

- **DE:** `DEG_summary_all_variants.tsv` com 25 combinações variante×método. Ex.: `tcga_matched/voom_quality_weights` = 4.015 DEGs (1.961 Up / 2.054 Down); `composition_adjusted/voom_quality_weights` = 2.592 DEGs.
- **Equivalência acelerada (etapa 24):** `DE_acceleration_equivalence.tsv` → `passed=TRUE`; diferenças máximas ~1e-13 (logFC, estatística, p, FDR, pesos).
- **GSEA robustez:** 14 vias robustas (entre elas Proteassoma `hsa03050`, Ciclo celular `hsa04110`, p53 `hsa04115`, Replicação de DNA `hsa03030`, Adesão focal `hsa04510`, ECM-receptor `hsa04512`, TNF `hsa04668`, Câncer de tireoide `hsa05216`).
- **Meta-análise ITGA2:** meta_logFC = 2,42 (IC 95% 1,94–2,91); FDR = 5,6e-22; I² = 66,7.
- **Gate ITGA2 (pré-especificado):** 5/6 critérios `TRUE`; **`surface_and_protein` = FALSE** → classificação final **"DEG robusto de interesse"** (nunca "alvo validado").
- **PPI:** 3 universos concluídos — genome_wide_top500 (34 nós), panel30 (1.333 nós), robust_leading_edge (1.109 nós); centralidade reportada como descritor topológico.

---

## 2. Checagens de integridade (etapa 31)

`VALIDATION_SUMMARY.tsv` registra **42 checks, todos `TRUE`**, incluindo:

- Matrizes críticas não-vazias e legíveis; checksums TCGA (`c001dfd1…`) e raw filtrada (`8614d59b…`) conferem.
- Design full-rank; CPM filter pré-especificado satisfeito.
- 5 métodos DE e 5 variantes presentes; 5 coleções GSEA presentes.
- Correção de batch não-catastrófica (rho=0,899; top500_direction=0,998; median_abs_delta=0,097); sinal de source reduzido (raw 0,096 → corrigido 0,064).
- EPIC: 51 amostras não-convergentes divulgadas.
- Gate ITGA2 completo e vocabulário restrito (sem "alvo terapêutico validado"/"nanomedicina pronta").

Manifesto de checksums completo: `results/phase2/upgrade_2026/validation/OUTPUT_FILE_CHECKSUMS.tsv` (125 artefatos).

---

## 3. Desvios do protocolo a priori (todos corretivos, sem mudança de critérios)

1. **Etapa 24 — checagem de identidade.** Causa-raiz: `limma::voomWithQualityWeights` armazena pesos como coluna de `data.frame`, descartando os nomes das amostras; o gate comparava `names()` (vazio) vs legado (564). Correção: restaurar nomes via `setNames(..., colnames(dge))` + alinhar por nome, e adicionar `gate_sha256` à `validation_key`.
2. **Etapa 27 — mutação cBioPortal.** A API retorna `entrezGeneId` (sem `hugoGeneSymbol`); usar `gene_name` (conhecido da query).
3. **Etapa 27 — SURFY.** Tabela S3 tem cabeçalho em duas linhas; ler com `skip=1` e detectar a coluna "UniProt gene".
4. **Etapa 28 — `igraph::simplify`.** Mascarado pelo genérico S4 `simplify` de `clusterProfiler`; qualificar com `igraph::simplify`.
5. **Etapa 28 — STRING >2000 proteínas.** API limita a 2000; `panel30` (2.477) e `leading_edge` (3.590) excedem. Dividir em lotes ≤1900 e mesclar (preserva a definição a priori dos universos).
6. **Etapa 31 — escopo `data.table` aninhado.** `cmp[gene_symbol %in% head(raw[order(-abs(statistic)), ...])]` falha por escopo; pré-computar `top_genes` antes do subset.
7. **Relatório (etapa 30) — figuras.** `include_graphics` com `output_dir` ≠ diretório do Rmd quebrava o caminho; embutir via `knitr::image_uri` + `<img>` (HTML self-contained).

**Nenhuma mudança** no painel de 30 vias nem nos critérios pré-especificados de ITGA2.

---

## 4. Arquivos gerados (resumo)

- `results/phase2/upgrade_2026/matrices/` — 9 matrizes (counts, logcpm, SV, manifestos, metadados).
- `results/phase2/upgrade_2026/differential_expression/` — 31 arquivos (5 variantes × 5 métodos + `DEG_summary_all_variants.tsv`).
- `results/phase2/upgrade_2026/gsea/` — fgsea, clusterProfiler, GSVA, bootstrap, robustez, leading-edge overlap.
- `results/phase2/upgrade_2026/deconvolution/` — composição, convergência EPIC, R² parcial.
- `results/phase2/upgrade_2026/validation/` — meta-análise, sobrevivência, superfície, gate ITGA2, `VALIDATION_SUMMARY.tsv`, `OUTPUT_FILE_CHECKSUMS.tsv`.
- `results/phase2/upgrade_2026/ppi/` — 3 redes (arestas, centralidade, módulos, enriquecimento, grafos).
- `results/phase2/reports/METHODOLOGICAL_UPGRADE_2026.html` — relatório self-contained.
- `results/phase2/figures/Fig23_Composition_Variance.png` e `Fig24_MetaAnalysis_Forest.png`.
- `data/cache/multiverse_de/v5.0.0/` — caches de equivalência e de métodos por variante.
- `results/phase2/upgrade_2026/pipeline/` — state files, `pipeline_benchmark.tsv`, `DE_acceleration_equivalence.tsv`, proveniências.

---

## 5. Reprodução

### A partir da etapa 24 (após as correções)

```powershell
powershell -NoProfile -Command "& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' scripts/run_phase2.R --from=24"
```

### Retomada incremental (qualquer ponto)

```powershell
powershell -NoProfile -Command "& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' scripts/run_phase2.R --from=N"
```

### Validação completa

```powershell
powershell -NoProfile -Command "& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' scripts/phase2/31_validate_upgrade.R"
```

**⚠️ Windows/MSYS2:** executar `Rscript.exe` diretamente pelo bash MSYS2 causa *segfault* em `data.table::fwrite/saveRDS`; sempre lance via **PowerShell** (ver `resume_phase2.ps1` e `HANDOFF_RESUME_2026-09-11.md`).

---

## 6. Rastreabilidade

- Versões de pacotes: `renv.lock` (R 4.6.1; limma 3.68.0; edgeR 4.10.0; DESeq2 1.51.7; fgsea 1.37.4; igraph 2.3.0; data.table 1.18.4).
- Seed: 42 (`SEED` em `00_config.R`).
- `logs/phase2_resume_*.log` preservam o stdout/stderr de cada retomada.
