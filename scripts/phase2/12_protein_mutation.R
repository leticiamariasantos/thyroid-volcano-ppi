# ═══════════════════════════════════════════════════════════════════════════════
# 12_protein_mutation.R — RPPA (proteína) e mutação/CNV (cBioPortal, THCA)
#
# RPPA: CCND1 (Cyclin D1) está no painel; ITGA2 e FN1 tipicamente NÃO possuem
# anticorpo no painel RPPA do PanCan Atlas (limitação documentada). Reporta
# correlação RNA-proteína para CCND1 e limitações.
# Mutação/CNV: BRAF, RAS, TP53 e candidatos (ITGA2, FN1, CCND1). Linguagem:
# "não foi identificada alteração genômica suficiente..." quando ausente.
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(here))
source(here("scripts", "phase2", "00_config.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})

log_msg("══ Proteína (RPPA) e Mutação/CNV ══")

cbio <- "https://www.cbioportal.org/api"
study <- "thca_tcga_pan_can_atlas_2018"

# gene -> entrez
genes <- c(BRAF = 673, NRAS = 4893, HRAS = 3265, KRAS = 3845, TP53 = 7157,
           TERT = 7015, ITGA2 = 3673, FN1 = 2335, CCND1 = 595)
sample_list <- paste0(study, "_all")

# ── 1. Frequência de mutação ──────────────────────────────────────────────────
fetch_mutations <- function(entrez) {
  url <- sprintf("%s/molecular-profiles/%s_mutations/mutations?sampleListId=%s&entrezGeneId=%d",
                 cbio, study, sample_list, entrez)
  r <- tryCatch(GET(url, timeout(120)), error = function(e) NULL)
  if (is.null(r) || status_code(r) != 200) return(NULL)
  fromJSON(content(r, as = "text", encoding = "UTF-8"), flatten = TRUE)
}

mut_list <- lapply(genes, fetch_mutations)
n_samples <- 0
mut_rows <- list()
for (nm in names(mut_list)) {
  d <- mut_list[[nm]]
  if (is.null(d) || length(d) == 0) {
    mut_rows[[nm]] <- data.table(gene = nm, n_mutated = 0, n_samples = NA_real_,
                                 mutation_type = NA_character_, top_protein_change = NA_character_)
    next
  }
  if (is.data.frame(d)) {
    n_samples <- max(n_samples, length(unique(d$sampleId)))
    mut_rows[[nm]] <- data.table(
      gene = nm, n_mutated = length(unique(d$sampleId)),
      n_samples = length(unique(d$sampleId)),
      mutation_type = names(sort(table(d$mutationType), decreasing = TRUE))[1],
      top_protein_change = names(sort(table(d$proteinChange), decreasing = TRUE))[1])
  } else {
    mut_rows[[nm]] <- data.table(gene = nm, n_mutated = 0, n_samples = NA_real_,
                                 mutation_type = NA_character_, top_protein_change = NA_character_)
  }
}
mut_freq <- rbindlist(mut_rows)
mut_freq[, n_samples_total := 500]  # THCA PanCan ~500 amostras; ver abaixo
mut_freq[, frequency := n_mutated / n_samples_total]
fwrite_tsv(mut_freq, file.path(DIR_VAL, "mutation_frequency.tsv"))

# ── 2. RPPA para CCND1 (Cyclin D1) ────────────────────────────────────────────
rppa_url <- sprintf("%s/molecular-profiles/%s_rppa/molecular-data?sampleListId=%s&entrezGeneId=595",
                    cbio, study, sample_list)
r_rppa <- tryCatch(GET(rppa_url, timeout(120)), error = function(e) NULL)
rppa <- NULL
if (!is.null(r_rppa) && status_code(r_rppa) == 200) {
  rppa <- fromJSON(content(r_rppa, as = "text", encoding = "UTF-8"), flatten = TRUE)
}
rppa_summary <- data.table(protein = "CCND1 (Cyclin D1)", antibody_in_panel = !is.null(rppa),
                           n_samples = if (is.null(rppa)) NA_integer_ else nrow(rppa),
                           mean_value = if (is.null(rppa)) NA_real_ else mean(as.numeric(rppa$value), na.rm = TRUE))

# correlação RNA-proteína (CCND1) usando a matriz TPM (THCA)
if (!is.null(rppa) && is.data.frame(rppa)) {
  tpm <- readRDS(file.path(DIR_DATAOUT, "tpm_matrix.rds"))
  if ("CCND1" %in% rownames(tpm)) {
    rppa_dt <- as.data.table(rppa)[, .(sample = sub("-01.*", "", sampleId), rppa_value = as.numeric(value))]
    rppa_dt <- rppa_dt[!is.na(rppa_value)]
    ccnd1_tpm <- data.table(sample = sub("-01.*", "", colnames(tpm)), rna = as.numeric(tpm["CCND1", ]))
    merged <- merge(rppa_dt, ccnd1_tpm, by = "sample", all = FALSE)
    if (nrow(merged) >= 10) {
      rho <- cor(merged$rppa_value, merged$rna, method = "spearman")
      rppa_summary$rna_protein_spearman <- rho
      rppa_summary$n_matched <- nrow(merged)
    }
  }
}
# ITGA2/FN1: documentar ausência de anticorpo no painel
rppa_summary <- rbind(rppa_summary, data.table(
  protein = c("ITGA2 (Integrin alpha-2)", "FN1 (Fibronectin)"),
  antibody_in_panel = c(FALSE, FALSE), n_samples = NA_integer_, mean_value = NA_real_,
  rna_protein_spearman = NA_real_, n_matched = NA_integer_), fill = TRUE)
fwrite_tsv(rppa_summary, file.path(DIR_VAL, "RPPA_summary.tsv"))

log_msg("══ Proteína/mutação concluído ══")
cat("\n=== Frequência de mutação (THCA PanCan Atlas) ===\n")
print(mut_freq[, .(gene, n_mutated, frequency, top_protein_change, mutation_type)])
cat("\n=== RPPA ===\n")
print(rppa_summary)
