# Random-effects meta-analysis, TCGA survival, surface evidence and gated prioritization.
suppressPackageStartupMessages({
  library(here); library(data.table); library(httr); library(jsonlite); library(metafor)
  library(survival); library(ggplot2)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

conv <- fread(file.path(DIR_VAL, "panel_convergence_genes.tsv"))
candidates <- unique(c("ITGA2", "FN1", "CCND1", head(conv[order(-n_robust_pathways_LE, -n_pathways_LE), gene], 18L)))
external <- fread(file.path(DIR_VAL, "validation_candidates.tsv"))[gene %in% candidates]
external[, `:=`(effect = logFC_external,
  se = abs(logFC_external / qnorm(pmax(pval, 1e-300) / 2, lower.tail = FALSE)), source_type = "GEO")]
external[!is.finite(se) | se <= 0, se := NA_real_]
effects <- external[is.finite(effect) & is.finite(se), .(gene, dataset, effect, se, source_type)]

tcga <- fread(file.path(DIR_DE_MULTI, "tcga_matched", "voom_qw_full_results.tsv"))[gene_symbol %in% candidates]
tcga_eff <- tcga[, .(gene = gene_symbol, dataset = "TCGA_THCA_matched", effect = logFC,
                     se = abs(logFC / statistic), source_type = "TCGA_matched")]
effects <- rbind(effects, tcga_eff[is.finite(se) & se > 0])
fwrite_tsv(effects, file.path(DIR_META, "meta_analysis_effects.tsv"))

meta_rows <- list(); forest_rows <- list()
for (g in candidates) {
  x <- effects[gene == g]
  if (nrow(x) < 2L) {
    meta_rows[[g]] <- data.table(gene = g, k = nrow(x), meta_logFC = NA_real_, meta_se = NA_real_,
                                 ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, I2 = NA_real_)
    next
  }
  fit <- metafor::rma(yi = x$effect, sei = x$se, method = "REML")
  meta_rows[[g]] <- data.table(gene = g, k = fit$k, meta_logFC = as.numeric(fit$b), meta_se = fit$se,
    ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval, I2 = fit$I2)
  forest_rows[[g]] <- rbind(x[, .(gene, label = dataset, estimate = effect,
                                   ci_low = effect - 1.96 * se, ci_high = effect + 1.96 * se, type = "study")],
    data.table(gene = g, label = "REML random effects", estimate = as.numeric(fit$b),
               ci_low = fit$ci.lb, ci_high = fit$ci.ub, type = "pooled"))
}
meta <- rbindlist(meta_rows); meta[, FDR := p.adjust(p, "BH")]
fwrite_tsv(meta, file.path(DIR_META, "meta_analysis_REML.tsv"))
forest <- rbindlist(forest_rows, fill = TRUE)
fwrite_tsv(forest, file.path(DIR_META, "forest_plot_data.tsv"))
pforest <- ggplot(forest[gene %in% c("ITGA2", "FN1", "CCND1", "TP53")],
  aes(estimate, reorder(label, estimate), xmin = ci_low, xmax = ci_high, shape = type)) +
  geom_vline(xintercept = 0, colour = "grey60") + geom_pointrange() + facet_wrap(~gene, scales = "free_y") +
  labs(x = "log2FC (95% CI)", y = NULL, title = "Meta-análise de efeitos aleatórios (REML)") + theme_minimal()
ggsave(file.path(DIR_FIG, "Fig24_MetaAnalysis_Forest.png"), pforest, width = 10, height = 7, dpi = 300)

# cBioPortal clinical and mutation data, cached and reshaped to one row per patient.
cbio <- "https://www.cbioportal.org/api"; study <- "thca_tcga_pan_can_atlas_2018"
fetch_json <- function(url) {
  r <- GET(url, add_headers(Accept = "application/json"), timeout(180)); stop_for_status(r)
  fromJSON(content(r, as = "text", encoding = "UTF-8"), flatten = TRUE)
}
clinical_cache <- file.path(DIR_CACHE, "cBioPortal", study, "patient_clinical.rds")
clinical_long <- tryCatch(safe_api_cache(clinical_cache, function() fetch_json(sprintf(
  "%s/studies/%s/clinical-data?clinicalDataType=PATIENT&projection=DETAILED&pageSize=100000", cbio, study)),
  validate = function(x) is.data.frame(x) && nrow(x) > 100L), error = function(e) NULL)

surv_results <- data.table(gene = character(), endpoint = character(), stratum_type = character(),
  stratum = character(), n = integer(), events = integer(), HR = numeric(),
  ci_low = numeric(), ci_high = numeric(), p = numeric(), model_type = character(),
  adjustment = character())
mutation_class <- data.table(participant = character(), mutation_class = character())
if (!is.null(clinical_long)) {
  cl <- as.data.table(clinical_long)
  clinical <- dcast(cl, patientId ~ clinicalAttributeId, value.var = "value", fun.aggregate = function(x) x[1])
  mutation_ids <- c(BRAF = 673, NRAS = 4893, HRAS = 3265, KRAS = 3845)
  mut <- rbindlist(Map(function(gene_name, entrez) {
    u <- sprintf("%s/molecular-profiles/%s_mutations/mutations?sampleListId=%s_all&entrezGeneId=%d",
                 cbio, study, study, entrez)
    mutation_cache <- file.path(DIR_CACHE, "cBioPortal", study,
      paste0("mutation_", gene_name, "_", entrez, ".rds"))
    z <- tryCatch(safe_api_cache(mutation_cache, function() fetch_json(u),
      validate = function(x) is.data.frame(x)), error = function(e) NULL)
    if (is.null(z) || !nrow(z)) return(data.table())
    # cBioPortal mutation responses for this study carry entrezGeneId but not
    # hugoGeneSymbol; the gene is already known from the query (gene_name).
    as.data.table(z)[, .(participant = substr(sampleId, 1L, 12L), gene = gene_name)]
  }, names(mutation_ids), mutation_ids), fill = TRUE)
  if (nrow(mut)) mutation_class <- mut[, .(mutation_class = if ("BRAF" %in% gene) "BRAF" else "RAS"), by = participant]

  expr <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
  md <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))[condition == "Tumor" & source == "TCGA"]
  ids <- intersect(md$sample, colnames(expr)); genes <- intersect(candidates, rownames(expr))
  ex <- as.data.table(t(expr[genes, ids, drop = FALSE]), keep.rownames = "sample")
  ex[, patientId := substr(sample, 1L, 12L)]
  ex <- ex[, lapply(.SD, mean), by = patientId, .SDcols = genes]
  dat <- merge(clinical, ex, by = "patientId", all = FALSE)
  dat <- merge(dat, mutation_class, by.x = "patientId", by.y = "participant", all.x = TRUE)
  dat[is.na(mutation_class), mutation_class := "WT_or_other"]
  hist_col <- intersect(c("HISTOLOGICAL_DIAGNOSIS", "HISTOLOGICAL_SUBTYPE", "SUBTYPE"), names(dat))
  if (length(hist_col)) hist_col <- hist_col[1]
  endpoints <- list(OS = c("OS_MONTHS", "OS_STATUS"), PFS = c("PFS_MONTHS", "PFS_STATUS"),
                    DSS = c("DSS_MONTHS", "DSS_STATUS"))
  fit_one <- function(d, gene, endpoint, stratum_type, stratum) {
    tm <- suppressWarnings(as.numeric(d[[endpoints[[endpoint]][1]]]))
    ev <- grepl("^(1:|DECEASED|PROGRESSED|YES|TRUE)", toupper(d[[endpoints[[endpoint]][2]]]))
    dd <- data.frame(time = tm, event = ev, expression = as.numeric(d[[gene]]))
    age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "DIAGNOSIS_AGE"), names(d))
    sex_col <- intersect(c("SEX", "GENDER"), names(d))
    stage_col <- intersect(c("AJCC_PATHOLOGIC_TUMOR_STAGE", "PATHOLOGIC_STAGE", "TUMOR_STAGE"), names(d))
    if (length(age_col)) dd$age <- suppressWarnings(as.numeric(d[[age_col[1]]]))
    if (length(sex_col)) dd$sex <- factor(d[[sex_col[1]]])
    if (length(stage_col)) {
      stage <- toupper(as.character(d[[stage_col[1]]]))
      dd$stage_advanced <- ifelse(grepl("(^| )STAGE (III|IV)|^(III|IV)", stage), 1,
        ifelse(grepl("(^| )STAGE (I|II)|^(I|II)", stage), 0, NA))
    }
    base_ok <- is.finite(dd$time) & dd$time >= 0 & is.finite(dd$expression)
    if (sum(base_ok) < 20L || sum(dd$event[base_ok]) < 5L) return(data.table())
    available <- intersect(c("age", "stage_advanced", "sex"), names(dd))
    available <- available[vapply(available, function(v) {
      z <- dd[[v]][base_ok]
      mean(!is.na(z)) >= 0.80 && length(unique(z[!is.na(z)])) > 1L
    }, logical(1))]
    max_covariates <- max(0L, floor(sum(dd$event[base_ok]) / 10L) - 1L)
    covariates <- head(available, max_covariates)
    repeat {
      vars <- c("time", "event", "expression", covariates)
      model_data <- dd[complete.cases(dd[, vars, drop = FALSE]), vars, drop = FALSE]
      if (!length(covariates) || (sum(model_data$event) >= 10L * (1L + length(covariates)) && nrow(model_data) >= 20L)) break
      covariates <- head(covariates, -1L)
    }
    if (nrow(model_data) < 20L || sum(model_data$event) < 5L) return(data.table())
    form <- as.formula(paste("Surv(time, event) ~ expression",
      if (length(covariates)) paste("+", paste(covariates, collapse = " + ")) else ""))
    fit <- tryCatch(coxph(form, data = model_data), error = function(e) NULL)
    if (is.null(fit) || !"expression" %in% rownames(summary(fit)$coef)) return(data.table())
    s <- summary(fit)
    data.table(gene = gene, endpoint = endpoint, stratum_type = stratum_type, stratum = stratum,
      n = nrow(model_data), events = sum(model_data$event), HR = s$coef["expression", "exp(coef)"],
      ci_low = s$conf.int["expression", "lower .95"], ci_high = s$conf.int["expression", "upper .95"],
      p = s$coef["expression", "Pr(>|z|)"],
      model_type = if (length(covariates)) "clinically_adjusted" else "expression_only_low_events",
      adjustment = if (length(covariates)) paste(covariates, collapse = "+") else "none")
  }
  surv_results <- rbindlist(lapply(genes, function(g) rbindlist(lapply(names(endpoints), function(ep) {
    parts <- list(fit_one(dat, g, ep, "all", "all"))
    mut_parts <- split(dat, dat$mutation_class)
    parts <- c(parts, lapply(names(mut_parts), function(s)
      fit_one(mut_parts[[s]], g, ep, "mutation", s)))
    if (length(hist_col)) {
      hist_parts <- split(dat, dat[[hist_col]])
      parts <- c(parts, lapply(names(hist_parts), function(s)
        fit_one(hist_parts[[s]], g, ep, "histology", s)))
    }
    rbindlist(parts, fill = TRUE)
  }))), fill = TRUE)
}
fwrite_tsv(surv_results, file.path(DIR_META, "TCGA_survival_stratified.tsv"))

# Human Protein Atlas downloadable table; inability to retrieve is an explicit failed criterion.
hpa_cache <- file.path(DIR_CACHE, "HPA", "2026-06-03", "proteinatlas.rds")
hpa <- tryCatch(safe_api_cache(hpa_cache, function() {
  z <- tempfile(fileext = ".zip"); download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", z, mode = "wb", quiet = TRUE)
  member <- unzip(z, list = TRUE)$Name[1]; as.data.table(read.delim(unz(z, member), check.names = FALSE))
}, validate = function(x) is.data.frame(x) && nrow(x) > 10000L), error = function(e) NULL)

surface <- data.table(gene = candidates, HPA_available = FALSE, plasma_membrane = FALSE,
                      protein_evidence = FALSE, antibody_localization_evidence = FALSE,
                      surfaceome_evidence = FALSE,
                      uniprot_membrane_annotation = FALSE)
if (!is.null(hpa)) {
  gene_col <- intersect(c("Gene", "Gene name"), names(hpa))
  if (length(gene_col)) gene_col <- gene_col[1]
  loc_cols <- grep("Subcellular|Main location|Additional location", names(hpa), value = TRUE)
  ev_cols <- grep("Reliability|Protein evidence|Evidence", names(hpa), value = TRUE)
  for (g in candidates) {
    z <- if (length(gene_col)) hpa[get(gene_col) == g] else hpa[0]
    if (nrow(z)) {
      loc <- paste(unlist(z[, ..loc_cols]), collapse = ";")
      ev <- paste(unlist(z[, ..ev_cols]), collapse = ";")
      surface[gene == g, `:=`(HPA_available = TRUE,
        plasma_membrane = grepl("plasma membrane|cell junction", loc, ignore.case = TRUE),
        protein_evidence = nzchar(ev) && !grepl("not detected", ev, ignore.case = TRUE),
        antibody_localization_evidence = grepl("plasma membrane|cell junction", loc, ignore.case = TRUE) &&
          nzchar(ev) && !grepl("not detected", ev, ignore.case = TRUE))]
    }
  }
}
# In silico Human Surfaceome (SURFY; Bausch-Fluck et al., 2018).
surfy_cache <- file.path(DIR_CACHE, "Surfaceome_SURFY", "2018", "table_S3_surfaceome.rds")
surfy <- tryCatch(safe_api_cache(surfy_cache, function() {
  z <- tempfile(fileext = ".xlsx")
  download.file("https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx", z, mode = "wb", quiet = TRUE)
  # SURFY Table S3 has a two-line header (title row + column names); skip the
  # title so "UniProt gene" (the gene-symbol column) becomes a real column name.
  as.data.table(readxl::read_excel(z, sheet = 1, skip = 1))
}, validate = function(x) is.data.frame(x) && nrow(x) > 2000L), error = function(e) NULL)
if (!is.null(surfy)) {
  gene_col <- grep("gene.*name|hgnc|symbol|uniprot.*gene", names(surfy), ignore.case = TRUE, value = TRUE)
  if (length(gene_col)) gene_col <- gene_col[1]
  if (length(gene_col)) {
    surfy_genes <- unique(unlist(strsplit(as.character(surfy[[gene_col]]), "[;, ]+")))
    surface[, surfaceome_evidence := gene %in% surfy_genes]
  }
}

# Independent reviewed membrane/topology annotation from UniProt.
for (g in candidates) {
  u <- paste0("https://rest.uniprot.org/uniprotkb/search?query=(gene_exact:", g,
              ")+AND+(organism_id:9606)+AND+(reviewed:true)&format=tsv&fields=gene_names,cc_subcellular_location,ft_topo_dom")
  uniprot_cache <- file.path(DIR_CACHE, "UniProt", "reviewed_human", paste0(g, ".rds"))
  txt <- tryCatch(safe_api_cache(uniprot_cache, function() {
    response <- GET(u, timeout(60)); stop_for_status(response)
    content(response, as = "text", encoding = "UTF-8")
  }, validate = function(x) is.character(x) && length(x) == 1L && nzchar(x)), error = function(e) "")
  surface[gene == g, uniprot_membrane_annotation := grepl("Cell membrane|topological domain", txt, ignore.case = TRUE)]
}
fwrite_tsv(surface, file.path(DIR_META, "surface_accessibility_evidence.tsv"))

# Quantitative gate; missing evidence is failure, not neutral evidence.
it_meta <- meta[gene == "ITGA2"]
raw_it <- fread(file.path(DIR_DE_MULTI, "raw", "voom_qw_full_results.tsv"))[gene_symbol == "ITGA2"]
adj_it <- fread(file.path(DIR_DE_MULTI, "composition_adjusted", "voom_qw_full_results.tsv"))[gene_symbol == "ITGA2"]
dir_it <- effects[gene == "ITGA2"]
ok <- function(x) length(x) == 1L && isTRUE(x)
surv_it <- surv_results[gene == "ITGA2" & stratum_type == "all" & model_type == "clinically_adjusted"]
criteria <- data.table(criterion = c("direction_consistency", "meta_effect", "heterogeneity",
  "surface_and_protein", "composition_independence", "survival_not_contradictory"), passed = c(
  nrow(dir_it) >= 2L && mean(dir_it$effect > 0) >= 0.8 && !any(dir_it$effect < 0 & abs(dir_it$effect / dir_it$se) > 1.96),
  ok(nrow(it_meta) == 1L && it_meta$meta_logFC >= ITGA2_MIN_LOGFC && it_meta$FDR < FDR_THRESH),
  ok(nrow(it_meta) == 1L && it_meta$I2 <= META_I2_MAX),
  ok(surface[gene == "ITGA2", plasma_membrane & protein_evidence & antibody_localization_evidence &
    surfaceome_evidence & uniprot_membrane_annotation]),
  ok(nrow(raw_it) == 1L && nrow(adj_it) == 1L && adj_it$logFC >= ITGA2_MIN_LOGFC &&
    adj_it$adj.P.Val < FDR_THRESH && adj_it$logFC / raw_it$logFC >= ITGA2_MIN_RETENTION),
  nrow(surv_it) > 0L && !nrow(surv_it[p < FDR_THRESH & HR < 1])
))
classification <- if (all(criteria$passed)) "candidato prioritário para investigação de surface targeting" else "DEG robusto de interesse"
criteria[, classification := classification]
fwrite_tsv(criteria, file.path(DIR_META, "ITGA2_prespecified_gate.tsv"))
log_msg("Validation complete; ITGA2 classification:", classification)
