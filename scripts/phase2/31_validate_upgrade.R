suppressPackageStartupMessages({ library(here); library(data.table) })
source(here("scripts", "phase2", "20_upgrade_config.R"))

checks <- list()
add <- function(name, ok, detail = "") checks[[length(checks)+1L]] <<- data.table(check=name, passed=isTRUE(ok), detail=as.character(detail))
required <- c(
  file.path(DIR_VARIANTS,"tcga_thca_primary_normal_counts.rds"), file.path(DIR_VARIANTS,"raw_counts.rds"),
  file.path(DIR_VARIANTS,"batch_corrected_counts.rds"), file.path(DIR_VARIANTS,"composition_adjusted_logcpm.rds"),
  file.path(DIR_DECONV,"epic_convergence_status.tsv"),
  file.path(DIR_DE_MULTI,"DEG_summary_all_variants.tsv"), file.path(DIR_GSEA_MULTI,"GSEA_fgsea_all_variants.tsv"),
  file.path(DIR_DECONV,"top_deg_celltype_variance.tsv"), file.path(DIR_META,"meta_analysis_REML.tsv"),
  file.path(DIR_META,"TCGA_survival_stratified.tsv"), file.path(DIR_META,"surface_accessibility_evidence.tsv"),
  file.path(DIR_META,"ITGA2_prespecified_gate.tsv"), file.path(DIR_PPI_MULTI,"PPI_upgrade_summary.tsv"),
  file.path(DIR_REP,"METHODOLOGICAL_UPGRADE_2026.html"))
for (f in required) add(paste("nonempty", basename(f)), file.exists(f) && file.info(f)$size > 0, f)

binary_files <- list.files(DIR_UPGRADE, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
binary_ok <- vapply(binary_files, function(f) !inherits(try(readRDS(f), silent = TRUE), "try-error"), logical(1))
add("all upgrade RDS readable", length(binary_files) > 0L && all(binary_ok),
    paste(basename(binary_files[!binary_ok]), collapse = ","))
tabular_files <- list.files(DIR_UPGRADE, pattern = "[.]tsv$", recursive = TRUE, full.names = TRUE)
tabular_ok <- vapply(tabular_files, function(f) !inherits(try(fread(f), silent = TRUE), "try-error"), logical(1))
add("all upgrade TSV files readable", length(tabular_files) > 0L && all(tabular_ok),
    paste(basename(tabular_files[!tabular_ok]), collapse = ","))

md <- fread(file.path(DIR_VARIANTS,"analysis_metadata.tsv")); design <- model.matrix(~ condition + source, md)
add("batch design full rank", qr(design)$rank == ncol(design), paste(qr(design)$rank,ncol(design),sep="/"))
tm <- fread(file.path(DIR_VARIANTS,"tcga_thca_primary_normal_metadata.tsv"))
add("adjacent normals >=20", sum(tm$condition=="Normal") >= 20, sum(tm$condition=="Normal"))
add("paired participants >=20", uniqueN(tm[paired==TRUE,participant]) >= 20, uniqueN(tm[paired==TRUE,participant]))
tcga_counts <- readRDS(file.path(DIR_VARIANTS,"tcga_thca_primary_normal_counts.rds"))
raw_counts <- readRDS(file.path(DIR_VARIANTS,"raw_counts.rds"))
batch_counts <- readRDS(file.path(DIR_VARIANTS,"batch_corrected_counts.rds"))
add("TCGA counts metadata aligned", identical(colnames(tcga_counts), tm$sample), paste(dim(tcga_counts), collapse="x"))
add("augmented counts metadata aligned", identical(colnames(raw_counts), md$sample), paste(dim(raw_counts), collapse="x"))
add("batch counts dimensions preserved", identical(dim(raw_counts), dim(batch_counts)) &&
      identical(dimnames(raw_counts), dimnames(batch_counts)), paste(dim(batch_counts), collapse="x"))
add("sample identifiers unique and nonmissing", !anyDuplicated(md$sample) && all(nzchar(md$sample)), nrow(md))
add("pre-specified CPM filter satisfied", all(keep_by_cpm(raw_counts, md$condition,
    min_cpm = MIN_EXPR_CPM, fraction = EXPR_FRAC)),
    sprintf("CPM>=%s; fraction=%s; smaller_group=%s", MIN_EXPR_CPM, EXPR_FRAC, min(table(md$condition))))

tcga_manifest <- fread(file.path(DIR_VARIANTS,"tcga_thca_primary_normal_manifest.tsv"))
tcga_expected <- tcga_manifest[field == "matrix_fingerprint_sha256", value][1]
batch_manifest <- fread(file.path(DIR_VARIANTS,"batch_correction_manifest.tsv"))
raw_expected <- batch_manifest[field == "raw_matrix_fingerprint_sha256", value][1]
add("TCGA matrix checksum", identical(matrix_checksum(tcga_counts), tcga_expected), tcga_expected)
add("raw filtered matrix checksum", identical(matrix_checksum(raw_counts), raw_expected), raw_expected)

ds <- fread(file.path(DIR_DE_MULTI,"DEG_summary_all_variants.tsv"))
add("five DE methods", all(c("limma_continuous","voom","voom_quality_weights","edgeR_QL","DESeq2") %in% ds$method), paste(unique(ds$method),collapse=","))
add("five variants", all(c("raw","batch_corrected","composition_adjusted","tcga_matched","tcga_paired") %in% ds$variant), paste(unique(ds$variant),collapse=","))

deconv <- fread(file.path(DIR_DECONV,"deconvolution_status.tsv"))
add("formal xCell and EPIC completed", all(deconv[method %in% c("xCell","EPIC"), completed]),
    paste(deconv$method, deconv$completed, collapse=","))
epic_status <- fread(file.path(DIR_DECONV,"epic_convergence_status.tsv"))
add("EPIC convergence status aligned", nrow(epic_status) == nrow(md) &&
    identical(epic_status$sample, md$sample) && !anyNA(epic_status$epic_converged),
    sprintf("samples=%s; nonconverged=%s", nrow(epic_status), sum(!epic_status$epic_converged)))
add("EPIC nonconvergence disclosed", deconv[method == "EPIC", nonconverged_samples] ==
    sum(!epic_status$epic_converged), deconv[method == "EPIC", nonconverged_samples])

gsea <- fread(file.path(DIR_GSEA_MULTI,"GSEA_fgsea_all_variants.tsv"))
add("GSEA five collections", all(c("Hallmark","KEGG","Reactome","GO_BP","Panel30") %in% gsea$collection),
    paste(unique(gsea$collection), collapse=","))
add("GSEA five variants", all(c("raw","batch_corrected","composition_adjusted","tcga_matched","tcga_paired") %in% gsea$variant),
    paste(unique(gsea$variant), collapse=","))
cp <- fread(file.path(DIR_GSEA_MULTI,"GSEA_clusterProfiler_primary.tsv"))
add("clusterProfiler replication nonempty", nrow(cp) > 0L, nrow(cp))

raw <- fread(file.path(DIR_DE_MULTI,"raw","voom_qw_full_results.tsv"))
bc <- fread(file.path(DIR_DE_MULTI,"batch_corrected","voom_qw_full_results.tsv"))
cmp <- merge(raw[,.(gene_symbol,raw=logFC)],bc[,.(gene_symbol,corrected=logFC)],by="gene_symbol")
rho <- cor(cmp$raw,cmp$corrected,method="spearman",use="complete.obs")
top <- cmp[gene_symbol %in% head(raw[order(-abs(statistic)),gene_symbol],500)]
direction <- mean(sign(top$raw)==sign(top$corrected))
median_delta <- median(abs(cmp$raw-cmp$corrected),na.rm=TRUE)
drastic <- rho < 0.50 || direction < 0.70 || median_delta > 2.5
add("batch correction not catastrophic", !drastic, sprintf("rho=%.3f; top500_direction=%.3f; median_abs_delta=%.3f",rho,direction,median_delta))
fwrite_tsv(data.table(metric=c("spearman_logFC","top500_direction","median_abs_delta","drastic_change"),value=c(rho,direction,median_delta,drastic)), file.path(DIR_META,"batch_regression_diagnostics.tsv"))

bd <- fread(file.path(DIR_VARIANTS,"batch_diagnostics.tsv"))
raw_source <- bd[matrix == "raw", mean(source_R2, na.rm=TRUE)]
corrected_source <- bd[matrix == "batch_corrected", mean(source_R2, na.rm=TRUE)]
add("source signal reduced after correction", is.finite(raw_source) && is.finite(corrected_source) && corrected_source < raw_source,
    sprintf("raw_mean_R2=%.3f; corrected_mean_R2=%.3f", raw_source, corrected_source))

surv <- fread(file.path(DIR_META,"TCGA_survival_stratified.tsv"))
add("survival endpoints analyzed", all(c("OS","PFS","DSS") %in% surv$endpoint), paste(unique(surv$endpoint),collapse=","))
add("survival mutation stratification", any(surv$stratum_type == "mutation"), sum(surv$stratum_type == "mutation"))
add("survival histology stratification", any(surv$stratum_type == "histology"), sum(surv$stratum_type == "histology"))
add("survival clinical adjustment", any(surv$model_type == "clinically_adjusted"), paste(unique(surv$model_type),collapse=","))
surface <- fread(file.path(DIR_META,"surface_accessibility_evidence.tsv"))
add("surface evidence sources represented", all(c("HPA_available","antibody_localization_evidence",
    "surfaceome_evidence","uniprot_membrane_annotation") %in% names(surface)) &&
    any(surface$HPA_available) && any(surface$surfaceome_evidence) && any(surface$uniprot_membrane_annotation), nrow(surface))

ppi <- fread(file.path(DIR_PPI_MULTI,"PPI_upgrade_summary.tsv"))
add("three PPI universes audited", setequal(ppi$network, c("genome_wide_top500","panel30","robust_leading_edge")),
    paste(ppi$network, ppi$status, collapse=","))
centrality_files <- list.files(DIR_PPI_MULTI, pattern="_centrality[.]tsv$", full.names=TRUE)
centrality_ok <- length(centrality_files) > 0L && all(vapply(centrality_files, function(f)
  all(fread(f)$interpretation == "topological descriptor; not therapeutic evidence"), logical(1)))
add("PPI centrality language constrained", centrality_ok, length(centrality_files))

gate <- fread(file.path(DIR_META,"ITGA2_prespecified_gate.tsv"))
add("ITGA2 gate complete", nrow(gate)==6L && all(!is.na(gate$passed)), unique(gate$classification))
add("ITGA2 vocabulary", all(unique(gate$classification) %in% c("DEG robusto de interesse","candidato prioritário para investigação de surface targeting")))
communication_files <- c(list.files(DIR_UPGRADE, pattern = "[.](tsv|txt|md|html)$", recursive = TRUE,
  full.names = TRUE), file.path(DIR_REP,"METHODOLOGICAL_UPGRADE_2026.html"))
communication_text <- paste(unlist(lapply(communication_files[file.exists(communication_files)],
  readLines, warn = FALSE)), collapse = "\n")
add("no prohibited therapeutic-readiness claims",
    !grepl("alvo terapêutico validado|nanomedicina pronta", communication_text, ignore.case = TRUE))

res <- rbindlist(checks); fwrite_tsv(res,file.path(DIR_META,"VALIDATION_SUMMARY.tsv"))
fail <- sum(!res$passed); cat("Upgrade validation:",fail,"failure(s)\n")
if (fail) { print(res[passed==FALSE]); quit(status=1) }

checksum_path <- file.path(DIR_META, "OUTPUT_FILE_CHECKSUMS.tsv")
artifact_files <- c(list.files(DIR_UPGRADE, recursive = TRUE, full.names = TRUE),
  file.path(DIR_REP, "METHODOLOGICAL_UPGRADE_2026.html"),
  file.path(DIR_FIG, c("Fig23_Composition_Variance.png", "Fig24_MetaAnalysis_Forest.png")))
artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]
artifact_files <- setdiff(normalizePath(artifact_files, winslash = "/", mustWork = TRUE),
  normalizePath(checksum_path, winslash = "/", mustWork = FALSE))
root_prefix <- paste0(normalizePath(here(), winslash = "/"), "/")
relative_paths <- ifelse(startsWith(artifact_files, root_prefix),
  substring(artifact_files, nchar(root_prefix) + 1L), artifact_files)
checksums <- data.table(
  path = relative_paths,
  bytes = file.info(artifact_files)$size,
  sha256 = vapply(artifact_files, digest::digest, character(1), file = TRUE, algo = "sha256")
)
fwrite_tsv(checksums, checksum_path)
