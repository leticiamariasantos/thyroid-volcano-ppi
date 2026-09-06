#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════════════
# 07_prioritization.R — FASE 2: Priorização multicritério de candidatos
# thyroid-volcano-ppi
#
# Modelo (definido ANTES de ver o ranking final, pesos justificados):
#   composite = 0.30*z(|logFC|) + 0.30*z(-log10 FDR) + 0.25*z(betweenness) + 0.15*z(degree)
#   → 60% evidência transcriptômica (Nível 2) + 40% topologia de rede (Nível 3).
# Sensibilidade: DE-only, rede-only, pesos iguais.
# HUB NÃO É ALVO: centralidade é apenas um dos componentes.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

PROJECT_ROOT <- here::here()
dir_pri <- file.path(PROJECT_ROOT, "09_target_prioritization")
dir.create(dir_pri, recursive = TRUE, showWarnings = FALSE)

cat("══ FASE 2 — Priorização multicritério ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

deg <- data.table::fread(file.path(PROJECT_ROOT, "04_differential_expression", "DEG_full_results.tsv"))
cen <- data.table::fread(file.path(PROJECT_ROOT, "08_ppi", "PPI_centrality.tsv"))
sel <- data.table::fread(file.path(PROJECT_ROOT, "07_pathway_redundancy", "pathways_selected.tsv"),
                         colClasses = c(id = "character"))
kegg <- readRDS(file.path(PROJECT_ROOT, "07_pathway_redundancy", "kegg_pathways.rds"))

# vias a que cada gene pertence (das selecionadas)
gene2path <- data.table(gene = unlist(kegg[sel$id], use.names = FALSE),
                        pathway = rep(sel$id, lengths(kegg[sel$id])))
gene2path <- gene2path[, .(n_pathways = .N, pathways = paste(sort(unique(pathway)), collapse = ";")), by = gene]

# ── 1. Juntar evidências ───────────────────────────────────────────────────────
d <- merge(cen[, .(gene_symbol, degree, betweenness, community)],
           deg[, .(gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)], by = "gene_symbol")
d <- merge(d, gene2path, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
d[is.na(n_pathways), n_pathways := 0]

d[, logFC_mag := abs(logFC)]
d[, sig := -log10(pmax(adj.P.Val, .Machine$double.xmin))]

z <- function(x) as.numeric(scale(x))
d[, z_lfc := z(logFC_mag)]
d[, z_sig := z(sig)]
d[, z_bet := z(betweenness)]
d[, z_deg := z(degree)]
d[, z_bet := ifelse(is.na(z_bet), 0, z_bet)]
d[, z_deg := ifelse(is.na(z_deg), 0, z_deg)]

# composite (60% DE, 40% rede)
d[, composite := 0.30*z_lfc + 0.30*z_sig + 0.25*z_bet + 0.15*z_deg]
# variantes de sensibilidade
d[, score_DE_only := 0.5*z_lfc + 0.5*z_sig]
d[, score_net_only := 0.6*z_bet + 0.4*z_deg]
d[, score_equal := 0.25*z_lfc + 0.25*z_sig + 0.25*z_bet + 0.25*z_deg]

d <- d[order(-composite)]

# ── 2. Tabela final de candidatos (campos obrigatórios; NA = não disponível) ──
cand <- d[, .(
  gene = gene_symbol,
  logFC = round(logFC, 3),
  FDR = signif(adj.P.Val, 3),
  pathways = pathways,
  n_pathways = n_pathways,
  degree = degree,
  betweenness = betweenness,
  community = community,
  # localização celular: NÃO inferida automaticamente (sem fonte validada)
  localization = "A CONFIRMAR",
  external_evidence = NA_character_,
  composite = round(composite, 3)
)]

data.table::fwrite(cand, file.path(dir_pri, "candidate_ranking.tsv"), sep = "\t")
data.table::fwrite(d, file.path(dir_pri, "prioritization_scores.tsv"), sep = "\t")

cat("  Top 25 candidatos (composite):\n")
print(head(cand, 25), row.names = FALSE)

# ── 3. Estabilidade do ranking entre os modelos ────────────────────────────────
r1 <- d$gene[order(-d$composite)]
r2 <- d$gene[order(-d$score_DE_only)]
r3 <- d$gene[order(-d$score_net_only)]
r4 <- d$gene[order(-d$score_equal)]
top30 <- function(r) head(r, 30)
overlap <- function(a, b) length(intersect(a, b)) / length(union(a, b))
stab <- data.frame(
  comparison = c("composite vs DE-only","composite vs net-only","composite vs equal",
                 "DE-only vs net-only"),
  jaccard_top30 = c(overlap(top30(r1), top30(r2)), overlap(top30(r1), top30(r3)),
                    overlap(top30(r1), top30(r4)), overlap(top30(r2), top30(r3)))
)
data.table::fwrite(stab, file.path(dir_pri, "ranking_stability.tsv"), sep = "\t")
cat("\n  Estabilidade (Jaccard top-30):\n")
print(stab, row.names = FALSE)

cat("══ FASE 2 — Priorização CONCLUÍDA ══\n")
