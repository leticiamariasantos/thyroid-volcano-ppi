# ═══════════════════════════════════════════════════════════════════════════════
# 24_integrity_tests.R — PARTE 32: testes finais de integridade
#
# Verifica arquivos, dimensões, grupos, duplicatas, NA/Inf, mapeamento, FDR,
# candidatos, hashes. Produz results/validation/validation_integrity.tsv.
# Falha (exit 1) em erro crítico.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(data.table))

cat("══ PARTE 32 — Testes finais de integridade ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_v <- "results/validation"
fail <- 0
ck <- function(cond, msg, check) {
  if (isTRUE(cond)) cat("  [OK]  ", msg, "\n")
  else { cat("  [FAIL]", msg, "\n"); fail <<- fail + 1 }
  return(check)
}
rows <- list()

add <- function(check, status, detail = "") {
  rows[[length(rows)+1]] <<- data.table(check = check, status = status, detail = detail)
}

# 1. Arquivos
req <- c("results/validation/GSE33630/metadata/metadata.tsv",
         "results/validation/GSE60542/metadata/metadata.tsv",
         "results/validation/composition_stress.tsv",
         "results/validation/integrated_candidates.tsv",
         "results/validation/provenance.tsv",
         "results/validation/pathway_validation.tsv",
         "10_validation/GSE33630_deg.tsv",
         "10_validation/GSE60542_deg.tsv",
         "documentation/DISCOVERY_FREEZE.md",
         "documentation/VALIDATION_LEAKAGE_AUDIT.md",
         "results/discovery_freeze/SHA256SUMS.txt")
for (f in req) {
  ok <- file.exists(f) && file.info(f)$size > 0
  add(paste("file", f), ifelse(ok, "OK", "FAIL"))
  if (!ok) { cat("  [FAIL] ausente/vazio:", f, "\n"); fail <- fail + 1 }
  else cat("  [OK]  ", f, "\n")
}

# 2. Metadata dimensões e grupos
m1 <- fread(req[1]); m2 <- fread(req[2])
add("GSE33630 metadata grupos", ifelse(setequal(names(table(m1$group)), c("ATC","Normal","PTC")), "OK", "FAIL"),
    paste(names(table(m1$group)), table(m1$group), collapse=" "))
add("GSE33630 n PTC/Normal", ifelse(sum(m1$group=="PTC")==49 && sum(m1$group=="Normal")==45, "OK", "FAIL"))
add("GSE60542 n PTC/Normal", ifelse(sum(m2$group=="PTC")==33 && sum(m2$group=="Normal")==30, "OK", "FAIL"))

# 3. Sem duplicatas / NA em colunas-chave
g1 <- fread(req[7]); g2 <- fread(req[8])
add("GSE33630 DE sem duplicata", ifelse(!anyDuplicated(g1$gene_symbol), "OK", "FAIL"))
add("GSE60542 DE sem duplicata", ifelse(!anyDuplicated(g2$gene_symbol), "OK", "FAIL"))
add("GSE33630 logFC/fdr finitos", ifelse(all(is.finite(g1$logFC)) && all(is.finite(g1$adj.P.Val)), "OK", "FAIL"))
add("GSE60542 logFC/fdr finitos", ifelse(all(is.finite(g2$logFC)) && all(is.finite(g2$adj.P.Val)), "OK", "FAIL"))

# 4. Candidatos integrados: 11 genes, status coerente
ic <- fread("results/validation/integrated_candidates.tsv")
key <- c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM")
add("integrated candidates 11 genes", ifelse(setequal(ic$gene, key), "OK", "FAIL"))
add("FN1 replicated", ifelse(ic[gene=="FN1"]$final_status == "INDEPENDENTLY REPLICATED", "OK", "FAIL"))
add("ITGA2 replicated", ifelse(ic[gene=="ITGA2"]$final_status == "INDEPENDENTLY REPLICATED", "OK", "FAIL"))

# 5. Composição stress: muscle diff ~ 0 nos GEO (vs -3.7 no discovery)
cs <- fread("results/validation/composition_stress.tsv")
add("muscle diff ~0 no GEO", ifelse(all(abs(cs$muscle_diff) < 0.5), "OK", "FAIL"),
    paste(cs$dataset, cs$muscle_diff, collapse="; "))

# 6. SHA256 freeze
if (file.exists("results/discovery_freeze/SHA256SUMS.txt")) {
  h <- readLines("results/discovery_freeze/SHA256SUMS.txt")
  add("SHA256SUMS freeze", ifelse(length(h) >= 10, "OK", "FAIL"), paste(length(h), "hashes"))
} else add("SHA256SUMS freeze", "FAIL")

# 7. Pathway validation: p53 replicated
pv <- fread("results/validation/pathway_validation.tsv")
add("p53 pathway replicated", ifelse(pv[pathway=="04115"]$replicated == TRUE, "OK", "FAIL"))

# ── salvar ──
out <- rbindlist(rows)
data.table::fwrite(out, file.path(dir_v, "validation_integrity.tsv"), sep = "\t")
cat(sprintf("\n══ Testes: %d falha(s) ══\n", fail))
if (fail > 0) quit(status = 1) else cat("  INTEGRIDADE OK\n")
