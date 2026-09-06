# FINAL AUDIT REPORT — Publication-grade audit (ETAPA 35)

**Data:** 2026-09-06

| Domínio | Status | Observação |
|---|---|---|
| Estatística (DE, FDR, thresholds) | PASS | limma/voom/DESeq2, BH, thresholds pré-definidos |
| Código (scripts R) | PASS | scripts 01–26 + orquestradores; idempotentes |
| Dados (matrizes) | PASS | TPM (58581×783), counts (56937×782) verificadas |
| Provenance (manifestos) | PASS | MANIFEST, MANIFEST_counts, provenance.tsv |
| Versões (R/pacotes) | WARNING | R 4.6.1 difere de renv.lock (R 4.6.0); renv.lock desatualizado |
| Hashes (SHA256) | PASS | discovery_freeze/SHA256SUMS.txt + provenance |
| Figuras | PASS | 14_figures + results/figures |
| Tabelas | PASS | results/* e numbered dirs |
| Documentação | PASS | FREEZE/LEAKAGE/REPORTS/TRACEABILITY/TRANSPARENCY |
| Containers (Docker) | WARNING | Dockerfile presente (Fase 1, R 4.4.0), não atualizado p/ Fase 2 |
| CI (GitHub Actions) | PASS | `.github/workflows/ci.yml` criado (structure check + smoke test + lint) |
| Seeds | PASS | set.seed(42) em GSEA/DESeq2/figuras |
| Batch/confound | WARNING | TCGA×GTEx estruturalmente confundido; documentado, não modelado formalmente |
| Composition | PASS | marker-based + covariável + stress test GEO |
| Data leakage | PASS | VALIDATION_LEAKAGE_AUDIT.md (nenhum vazamento) |
| Metadata/sample IDs | PASS | metadata.tsv GSE33630/GSE60542 |
| Gene annotation | WARNING | mapeamento probe 1:many (regra pré-registrada); MYL1 ausente em G026 |
| External validation | PASS | GSE33630 + GSE60542 (microarray) |
| Negative results | PASS | musculares não replicam; via muscular inverte |
| Single-cell | WARNING | não executado (limitação técnica) |
| Proteína/IHC/metilação/mutação | WARNING | não extraídos de dados primários (pendente) |

## Veredito

**NÃO "publication-ready"** — há **WARNINGs** remanescentes (renv.lock, Docker, single-cell,
proteína/mutação/metilação, batch modeling). O núcleo científico (discovery + composição +
validação externa) é sólido e rastreável, e o CI foi implementado, mas a reprodutibilidade total
(Docker/renv) e a triangulação multi-ômica (proteína/single-cell/mutação) permanecem pendentes.

## Ações corretivas pendentes (prioridade)

1. Atualizar `renv.lock` para R 4.6.1 + pacotes usados na Fase 2.
3. Atualizar Dockerfile para o ambiente da Fase 2.
4. Extrair TCGA RPPA (FN1) e paisagem mutacional (BRAF/RAS/TP53) de dados primários.
5. Referência single-cell (tireoide + músculo) para deconvolução formal e origem de FN1/ITGA2.
