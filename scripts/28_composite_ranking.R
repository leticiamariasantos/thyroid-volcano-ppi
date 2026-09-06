# ═══════════════════════════════════════════════════════════════════════════════
# 28_composite_ranking.R — Ranking composto multidimensional (Rank A molecular,
# Rank B translacional, Rank C nanomédico). Transparente, auditável, com
# sensibilidade de pesos, penalizações e tiers. Sem circularidade.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({library(data.table)})

cat("══ Ranking composto multidimensional ══\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dir_r <- "results/prioritization"
dir.create(dir_r, recursive = TRUE, showWarnings = FALSE)

# ── 1. Universo de candidatos (congelado, não-circular) ───────────────────────
# = genes com contexto (pathway/PPI) + top DE + controles negativos (musculares)
tpm  <- fread("04_differential_expression/DEG_full_results.tsv")
voo  <- fread("results/counts_deg/voom_full_results.tsv")
dds  <- fread("results/counts_deg/deseq2_full_results.tsv")
cen  <- fread("08_ppi/PPI_centrality.tsv")
g1   <- fread("10_validation/GSE33630_deg.tsv")
g2   <- fread("10_validation/GSE60542_deg.tsv")
pri  <- fread("09_target_prioritization/prioritization_scores.tsv")

muscle <- c("MYH1","MYH2","MYH7","MYH6","MYL1","MYL2","MYL3","MYL7","ACTA1","ACTC1",
            "TNNT1","TNNT3","TNNC1","TNNI1","TNNI2","TPM1","TPM2","TPM3","CKM","CKMT2",
            "DES","MB","ENO3","MYBPC1","MYBPC2","TTN","NEB","MYOM1","MYOM2","LDB3",
            "TCAP","MYOT","FLNC","ATP2A1","CASQ1","CASQ2")
epi <- c("TG","TPO","TSHR","PAX8","NKX2-1","FOXE1","SLC5A5","DIO1","DIO2","CALCA","TFF3","KRT7","KRT19")
immune_m <- c("PTPRC","CD3D","CD3E","CD8A","CD4","CD19","MS4A1","CD68","CD14","NCAM1",
              "FCGR3A","ITGAX","NKG7","GNLY","GZMB","PRF1","CD163","CD79A")
fibro <- c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","FN1","DCN","LUM","VIM","ACTA2","FAP","PDGFRB","POSTN","MMP2","LOX","BGN")
endo <- c("PECAM1","VWF","CDH5","ENG","KDR","FLT1","EMCN","CLDN5","TEK")

universe <- unique(c(cen$gene_symbol, tpm$gene_symbol[order(-abs(tpm$logFC))][1:60], muscle,
                     c("FN1","ITGA2","CTSS","HLA-DPA1","CCND1","SOX10","TFF1","DIO3")))
universe <- universe[!is.na(universe)]

# RPPA (proteína): genes com dados proteicos TCGA
rppa_present <- c("FN1","ITGA2","CCND1","BRAF","TP53","NRAS")

# ── 2. Componentes por gene ───────────────────────────────────────────────────
z <- function(x){ s <- as.numeric(scale(x)); s[is.na(s)] <- 0; s }

build <- function(g) {
  t <- tpm[gene_symbol==g]; v <- voo[gene_symbol==g]; d <- dds[gene_symbol==g]
  a <- g1[gene_symbol==g]; b <- g2[gene_symbol==g]; c <- cen[gene_symbol==g]; p <- pri[gene_symbol==g]

  # A. robustez transcriptômica
  lfc_t <- if(nrow(t)) t$logFC[1] else NA_real_
  lfc_v <- if(nrow(v)) v$logFC[1] else NA_real_
  lfc_d <- if(nrow(d)) d$log2FoldChange[1] else NA_real_
  fdr_t <- if(nrow(t)) t$adj.P.Val[1] else NA_real_
  fdr_v <- if(nrow(v)) v$adj.P.Val[1] else NA_real_
  fdr_d <- if(nrow(d)) d$padj[1] else NA_real_
  sigs <- c(fdr_t<0.05 & abs(lfc_t)>1, fdr_v<0.05 & abs(lfc_v)>1, fdr_d<0.05 & abs(lfc_d)>1)
  dirs <- sign(c(lfc_t,lfc_v,lfc_d))
  dir_cons <- length(unique(dirs[!is.na(dirs)]))==1
  method_conc <- sum(sigs, na.rm=TRUE)/3
  de_mag <- if(!is.na(lfc_t)) abs(lfc_t) else 0
  de_sig <- if(!is.na(fdr_t)) -log10(pmax(fdr_t,1e-300)) else 0

  # B. replicação externa
  rep_score <- 0; rep_n <- 0
  for (geo in list(a,b)) {
    if (nrow(geo)) {
      rep_n <- rep_n + 1
      dir_ok <- sign(geo$logFC[1])==sign(lfc_t) & !is.na(lfc_t)
      if (dir_ok & geo$adj.P.Val[1] < 0.05) rep_score <- rep_score + 1
      else if (dir_ok) rep_score <- rep_score + 0.5
    }
  }
  rep_score <- if(rep_n>0) rep_score/rep_n else NA_real_

  # C. proteína
  prot <- if(g %in% rppa_present) 1 else NA_real_  # NA = não medido (não penalizar)

  # D. composição (flag)
  is_muscle <- g %in% muscle
  is_epi <- g %in% epi
  is_fibro <- g %in% fibro
  is_immune <- g %in% immune_m
  is_endo <- g %in% endo

  # E. pathway context
  n_path <- if(nrow(p)) p$n_pathways[1] else 0

  # F. PPI
  deg <- if(nrow(c)) c$degree[1] else 0
  bet <- if(nrow(c)) c$betweenness[1] else 0

  data.table(gene=g, logFC_TPM=lfc_t, logFC_voom=lfc_v, logFC_DESeq2=lfc_d,
             FDR_TPM=fdr_t, method_concordance=method_conc, dir_consistent=dir_cons,
             de_mag=de_mag, de_sig=de_sig, replication=rep_score, protein=prot,
             is_muscle=is_muscle, is_epi=is_epi, is_fibro=is_fibro, is_immune=is_immune,
             is_endo=is_endo, n_pathways=n_path, degree=deg, betweenness=bet)
}
U <- rbindlist(lapply(universe, build))

# ── 3. Normalização 0–1 + z ───────────────────────────────────────────────────
U[, de_mag_n := pmin(abs(de_mag)/max(abs(de_mag),na.rm=TRUE),1)]
U[, de_sig_n := pmin(de_sig/max(de_sig,na.rm=TRUE),1)]
U[, ppi_deg_n := pmin(degree/max(degree),1)]
U[, ppi_bet_n := pmin(betweenness/max(betweenness),1)]
U[, pathway_n := pmin(n_pathways/2,1)]
U[is.na(replication), replication := NA_real_]

# scores de dimensão
U[, TRANSCRIPTOMIC_ROBUSTNESS := 0.35*de_mag_n + 0.35*de_sig_n + 0.20*method_concordance + 0.10*dir_consistent]
U[, EXTERNAL_REPLICATION := replication]
U[, PROTEIN_SUPPORT := protein]  # NA = not measured
U[, PATHWAY_CONTEXT := pathway_n]
U[, PPI_SUPPORT := 0.5*ppi_deg_n + 0.5*ppi_bet_n]
U[, COMPOSITION_FLAG := fifelse(is_muscle, -1, fifelse(is_epi | is_endo, -0.3, fifelse(is_fibro, -0.2, 0)))]

# ── 4. MOLECULAR_EVIDENCE_SCORE (sem proteína quando NA) ──────────────────────
U[, prot_impute := fifelse(is.na(PROTEIN_SUPPORT), 0.5, PROTEIN_SUPPORT)]  # neutral quando ausente
U[, MOLECULAR := 0.30*TRANSCRIPTOMIC_ROBUSTNESS + 0.20*ifelse(is.na(EXTERNAL_REPLICATION),0.5,EXTERNAL_REPLICATION)
              + 0.15*prot_impute + 0.15*PATHWAY_CONTEXT + 0.10*PPI_SUPPORT + 0.10*de_mag_n]

# ── 5. TRANSLATIONAL (especificidade + localização) ───────────────────────────
# especificidade: penaliza muscular/epitelial/endotelial/fibroblasto (não-tumor específico)
U[, specificity := fifelse(is_muscle, 0, fifelse(is_immune, 0.6, fifelse(is_fibro, 0.5, fifelse(is_endo,0.5, fifelse(is_epi,0.7, 0.8)))))]
U[, TRANSLATIONAL := 0.25*TRANSCRIPTOMIC_ROBUSTNESS + 0.20*ifelse(is.na(EXTERNAL_REPLICATION),0.5,EXTERNAL_REPLICATION)
                    + 0.15*prot_impute + 0.15*specificity + 0.10*PATHWAY_CONTEXT + 0.10*PPI_SUPPORT + 0.05*de_mag_n]

# ── 6. NANOMEDICINE (acessibilidade: extracelular/superfície > intracelular) ──
# acessibilidade estimada por categoria funcional (ECM/superfície = alta; intracelular = baixa)
surface_ecm <- c("FN1","ITGA2","ITGA2B","ITGA5","ITGA7","ITGA10","LAMA2","HSPG2","COL6A1","COL6A2","COL6A3","ELN","DCN","BGN","FBLN1","SDC1","SDC4","THBS1","THBS2","THBS3")
access <- fifelse(U$gene %in% surface_ecm, 1.0,
           fifelse(U$gene %in% c("CTSS","HLA-DPA1","HLA-DQB2"), 0.8,  # secretado/superfície imune
           fifelse(U$gene %in% c("CCND1","SOX10","TP53"), 0.3,  # intracelular/nuclear
           0.5)))  # desconhecido = neutro
U[, TARGET_ACCESSIBILITY := access]
U[, NANOMEDICINE := 0.20*TRANSCRIPTOMIC_ROBUSTNESS + 0.20*ifelse(is.na(EXTERNAL_REPLICATION),0.5,EXTERNAL_REPLICATION)
                  + 0.15*prot_impute + 0.15*TARGET_ACCESSIBILITY + 0.15*specificity
                  + 0.10*PPI_SUPPORT + 0.05*de_mag_n]

# ── 7. Penalizações explícitas ────────────────────────────────────────────────
U[, PENALTY := 0]
U[is_muscle==TRUE, PENALTY := PENALTY + 0.35]        # composição muscular
U[!is.na(replication) & replication < 0.5, PENALTY := PENALTY + 0.15]  # replicação fraca
U[specificity < 0.5, PENALTY := PENALTY + 0.10]      # baixa especificidade
U[, FINAL_MOLECULAR := pmax(0, MOLECULAR - PENALTY)]
U[, FINAL_TRANSLATIONAL := pmax(0, TRANSLATIONAL - PENALTY)]
U[, FINAL_NANOMEDICINE := pmax(0, NANOMEDICINE - PENALTY)]

# ── 8. Sensibilidade de pesos (1000 combinações) ──────────────────────────────
`%||%` <- function(a,b) if(is.null(a)) b else a
set.seed(42)
combos <- 1000
top10_freq <- list()
for (i in 1:combos) {
  w <- runif(6); w <- w/sum(w)
  sc <- w[1]*U$TRANSCRIPTOMIC_ROBUSTNESS + w[2]*ifelse(is.na(U$EXTERNAL_REPLICATION),0.5,U$EXTERNAL_REPLICATION) +
        w[3]*U$prot_impute + w[4]*U$PATHWAY_CONTEXT + w[5]*U$PPI_SUPPORT + w[6]*U$de_mag_n
  r <- U$gene[order(-sc)][1:10]
  for (g in r) top10_freq[[g]] <- (top10_freq[[g]] %||% 0) + 1
}
stab <- data.table(gene = names(top10_freq), top10_freq = unlist(top10_freq))
U <- merge(U, stab, by="gene", all.x=TRUE)
U[is.na(top10_freq), top10_freq := 0]
U[, RANKING_STABILITY := top10_freq/combos]

# ── 9. Tiers ──────────────────────────────────────────────────────────────────
U[, TIER := fifelse(FINAL_MOLECULAR >= quantile(FINAL_MOLECULAR, 0.9), "TIER 1",
             fifelse(FINAL_MOLECULAR >= quantile(FINAL_MOLECULAR, 0.75), "TIER 2",
             fifelse(FINAL_MOLECULAR >= quantile(FINAL_MOLECULAR, 0.5), "TIER 3",
             fifelse(is_muscle==TRUE, "TIER 5 (composicional)", "TIER 4"))))]
U[, NANO_TIER := fifelse(FINAL_NANOMEDICINE >= 0.5 & !is_muscle, "NANO-TIER B",
                  fifelse(is_muscle, "NANO-TIER E",
                  fifelse(FINAL_NANOMEDICINE >= 0.35, "NANO-TIER C", "NANO-TIER D")))]

# ── 10. Ordenar e salvar ──────────────────────────────────────────────────────
U <- U[order(-FINAL_MOLECULAR)]
cols <- c("gene","logFC_TPM","logFC_voom","logFC_DESeq2","FDR_TPM","method_concordance",
          "replication","protein","is_muscle","n_pathways","degree","TARGET_ACCESSIBILITY",
          "TRANSCRIPTOMIC_ROBUSTNESS","EXTERNAL_REPLICATION","PROTEIN_SUPPORT","PATHWAY_CONTEXT",
          "PPI_SUPPORT","MOLECULAR","TRANSLATIONAL","NANOMEDICINE","PENALTY",
          "FINAL_MOLECULAR","FINAL_TRANSLATIONAL","FINAL_NANOMEDICINE","RANKING_STABILITY",
          "TIER","NANO_TIER")
data.table::fwrite(U[, ..cols], file.path(dir_r, "integrated_gene_ranking.tsv"), sep="\t")
data.table::fwrite(U[order(-FINAL_MOLECULAR)], file.path(dir_r, "molecular_evidence_ranking.tsv"), sep="\t")
data.table::fwrite(U[order(-FINAL_TRANSLATIONAL)], file.path(dir_r, "translational_priority_ranking.tsv"), sep="\t")
data.table::fwrite(U[order(-FINAL_NANOMEDICINE)], file.path(dir_r, "nanomedicine_hypothesis_ranking.tsv"), sep="\t")

cat("── Top 15 — evidência molecular ──\n")
print(U[order(-FINAL_MOLECULAR)][1:15, .(gene, FINAL_MOLECULAR, TIER, replication, is_muscle, degree)])
cat("\n── Top 15 — prioridade translacional ──\n")
print(U[order(-FINAL_TRANSLATIONAL)][1:15, .(gene, FINAL_TRANSLATIONAL, specificity, TARGET_ACCESSIBILITY, is_muscle)])
cat("\n── Top 15 — hipótese nanomédica ──\n")
print(U[order(-FINAL_NANOMEDICINE)][1:15, .(gene, FINAL_NANOMEDICINE, NANO_TIER, TARGET_ACCESSIBILITY, replication)])
cat("\n── Controles negativos (musculares) devem ter baixa prioridade ──\n")
print(U[gene %in% c("MYH7","MYL1","MYL2","ACTA1","TNNT3","CKM"), .(gene, FINAL_MOLECULAR, FINAL_TRANSLATIONAL, FINAL_NANOMEDICINE, TIER, NANO_TIER)])

cat("\n══ Ranking composto CONCLUÍDO ══\n")
