# Handoff de retomada — thyroid-volcano-ppi (v5.0.0)

**Data do checkpoint:** 2026-09-11 ~00:03 (America/Sao_Paulo)
**Ponto de parada:** Etapa 25 (`25_multiverse_gsea.R`) interrompida após início (nenhum output escrito).
**Próximo passo:** retomar a partir da etapa 25.

---

## 1. Resumo do estado

| Etapa | Status | Observação |
|---|---|---|
| 20–23 | ✅ concluídas (anteriores) | state files em `pipeline/` |
| 24 (`24_multiverse_de.R`) | ✅ **concluída** | checkpoint salvo 23:00:35; 31 arquivos em `differential_expression/` |
| 25 (`25_multiverse_gsea.R`) | ⏸ interrompida no início | nenhum output em `gsea/` escrito; reinicia limpo (cache msigdb já salvo) |
| 26–31 | ⏳ pendentes | dependem de 25 |

### Auditoria da etapa 24
- **Gate de equivalência formal:** `DE_acceleration_equivalence.tsv` → `passed=TRUE`,
  todas as 9 checagens `TRUE`, diferenças máximas ~1e-13.
- **`DEG_summary_all_variants.tsv`:** 25 combinações variante×método (5 variantes × 5 métodos).
- **Proveniência:** `tcga_matched` = `official_reference_cache` (legacy validado);
  `raw`/`batch_corrected`/`composition_adjusted`/`tcga_paired` = `equivalent_linear_algebra`.
- **Tempo da etapa 24:** ~8,2 h (dominado pelo DESeq2 de `composition_adjusted`, design
  ~19 covariáveis; e `tcga_paired`, design ~61 covariáveis). Todos os caches por método
  foram persistidos → a retomada **não** recalcula a etapa 24.

---

## 2. Causa-raiz corrigida (checagem de identidade da etapa 24)

O `limma::voomWithQualityWeights` armazena os pesos como coluna de `data.frame`
(`$targets$sample.weights`), o que **descarta os nomes das amostras**. O gate comparava
`names()` do vetor recomputado (vazio) contra o cache legado (564 nomes) → `samples=FALSE`
→ `passed=FALSE`.

**Correção mínima em `scripts/phase2/24_multiverse_de.R`** (função `validate_accelerated_reference`):
1. Restaurar nomes com `setNames(vqw$targets$sample.weights, colnames(dge))` e alinhar por
   nome antes de comparar — espelhando o que `run_variant()` já faz.
2. Adicionar `gate_sha256 = file_sha256("scripts/phase2/24_multiverse_de.R")` à
   `validation_key`, para que mudanças no próprio gate invalidem o cache de validação.

Os caches obsoletos de `data/cache/multiverse_de/v5.0.0/reference_validation/` (13 colunas,
`passed=FALSE`) foram removidos e regenerados corretamente.

---

## 3. Como retomar

### Comando recomendado (PowerShell — obrigatório no Windows, ver §4)

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File resume_phase2.ps1 -From 25
```

Ou, diretamente:

```powershell
powershell -NoProfile -Command "& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' scripts/run_phase2.R --from=25"
```

Monitorar progresso:

```powershell
Get-Content logs\phase2_resume_*_stdout.log -Tail 20
# ou: tail -c 2000 logs/phase2_resume_*_stdout.log  (via bash MSYS2)
```

### Retomada incremental (se alguma etapa 25–31 falhar)
Após corrigir, retome do ponto de falha: `--from=<N>` (ex.: `--from=27`).

> `--force` só deve ser usado se a assinatura antiga estiver irrecuperavelmente
> inconsistente. As etapas 20–24 NÃO precisam de `--force`; elas serão puladas pelo cache.

---

## 4. Aviso operacional crítico (Windows + git-bash/MSYS2)

Executar `Rscript.exe` **diretamente** a partir do bash MSYS2 causa **segfault** em
`data.table::fwrite`/`saveRDS` (data.table 1.18.4 compilado para R 4.6.0 rodando em R 4.6.1,
agravado pelo ambiente MSYS2). Sempre lance via **PowerShell** (ou cmd/RStudio).

Verificado: `powershell -Command "& 'Rscript.exe' ..."` funciona; invocação direta do bash falha.

---

## 5. Mudanças não commitadas (working tree)

- `scripts/phase2/24_multiverse_de.R` — correção do gate + `gate_sha256`.
- `results/phase2/upgrade_2026/pipeline/DE_acceleration_equivalence.tsv` — regenerado (passed=TRUE).
- `results/phase2/upgrade_2026/pipeline/pipeline_benchmark.tsv` — entrada `completed` da etapa 24.
- `logs/session_info.txt` — atualizado pelo `00_config.R`.
- Diversos outputs novos em `results/phase2/upgrade_2026/differential_expression/` e `pipeline/DE_*`.
- `resume_phase2.ps1` — helper de retomada.

Nada foi commitado. Revise e commite quando desejar.

---

## 6. Etapas restantes (25 → 31)

| Etapa | Descrição | Output principal |
|---|---|---|
| 25 | fgsea + clusterProfiler + GSVA (5 variantes × 5 coleções) | `gsea/GSEA_fgsea_all_variants.tsv` |
| 26 | R² parcial por compartimento | `deconvolution/top_deg_celltype_variance.tsv` |
| 27 | Meta-análise REML + sobrevivência + superfície + gate ITGA2 | `validation/ITGA2_prespecified_gate.tsv` |
| 28 | Três redes PPI (STRING + coexpressão) | `ppi/PPI_upgrade_summary.tsv` |
| 29 | 100 bootstraps GSEA + remoção de outliers | `gsea/GSEA_robustness_summary.tsv` |
| 30 | Renderiza relatório HTML | `reports/METHODOLOGICAL_UPGRADE_2026.html` |
| 31 | Validação técnica integral + checksums | `validation/VALIDATION_SUMMARY.tsv` |

Etapas 27 e 28 dependem de internet (cBioPortal/HPA/SURFY/UniProt e STRING), com caches
locais versionados; sem rede e sem cache válido, falham explicitamente.
