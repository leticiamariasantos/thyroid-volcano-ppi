# Resumable Phase 2 + 2026 upgrade orchestrator.
# Usage: Rscript scripts/run_phase2.R [--force] [--from=20] [--adopt-existing]
suppressPackageStartupMessages({ library(digest); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
adopt <- "--adopt-existing" %in% args
from_arg <- grep("^--from=", args, value = TRUE)
from_num <- if (length(from_arg)) as.integer(sub("^--from=", "", from_arg[1])) else 0L
RS <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
state_dir <- file.path("results", "phase2", "upgrade_2026", "pipeline")
dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
benchmark_file <- file.path(state_dir, "pipeline_benchmark.tsv")
benchmark_history <- if (file.exists(benchmark_file)) fread(benchmark_file) else data.table()
benchmark_rows <- list()
invocation_id <- paste(format(Sys.time(), "%Y%m%dT%H%M%S"), Sys.getpid(), sep = "-")

step <- function(script, inputs, outputs) list(script = script, inputs = inputs, outputs = outputs)
P <- file.path
steps <- list(
  step(P("scripts","phase2","00_config.R"), c("renv.lock"), "logs/session_info.txt"),
  step(P("scripts","phase2","01_audit_data.R"), c("data/global/TCGA_GTEx_thyroid_tpm.tsv","data/global/TCGA_GTEx_thyroid_counts.tsv"), c("results/phase2/data/tpm_matrix.rds","results/phase2/data/counts_matrix.rds")),
  step(P("scripts","phase2","02_qc.R"), "results/phase2/data/tpm_matrix.rds", "results/phase2/preprocessing/QC_summary.txt"),
  step(P("scripts","phase2","03_de.R"), c("results/phase2/data/tpm_matrix.rds","results/phase2/data/counts_matrix.rds"), "results/phase2/differential_expression/DEG_summary.tsv"),
  step(P("scripts","phase2","04_panel.R"), character(), "results/phase2/pathways/PANEL_30_PATHWAYS.tsv"),
  step(P("scripts","phase2","05_gsea_global.R"), "results/phase2/data/limma_rank_t.rds", "results/phase2/gsea/GSEA_GLOBAL_ALL.tsv"),
  step(P("scripts","phase2","06_gsea_panel.R"), c("results/phase2/data/panel30_genesets.rds","results/phase2/data/limma_rank_t.rds"), "results/phase2/pathways/PANEL_30_GSEA_RESULTS.tsv"),
  step(P("scripts","phase2","07_redundancy.R"), "results/phase2/data/panel30_genesets.rds", "results/phase2/pathways/PATHWAY_REDUNDANCY.tsv"),
  step(P("scripts","phase2","08_composition.R"), c("results/phase2/data/tpm_matrix.rds","results/phase2/data/panel30_genesets.rds"), "results/phase2/sensitivity/gsea_composition_sensitivity.tsv"),
  step(P("scripts","phase2","09_ppi.R"), "results/phase2/differential_expression/limma_full_results.tsv", "results/phase2/ppi/PPI_summary.tsv"),
  step(P("scripts","phase2","10_validation.R"), character(), "results/phase2/validation/validation_candidates.tsv"),
  step(P("scripts","phase2","11_singlecell.R"), character(), "results/phase2/validation/singlecell_predominant_compartment.tsv"),
  step(P("scripts","phase2","12_protein_mutation.R"), "results/phase2/data/tpm_matrix.rds", "results/phase2/validation/mutation_frequency.tsv"),
  step(P("scripts","phase2","13_robustness.R"), "results/phase2/pathways/PANEL_30_GSEA_RESULTS.tsv", "results/phase2/pathways/PATHWAY_ROBUSTNESS_MATRIX.tsv"),
  step(P("scripts","phase2","14_convergence.R"), c("results/phase2/pathways/PANEL_30_GSEA_RESULTS.tsv","results/phase2/validation/validation_candidates.tsv"), "results/phase2/validation/convergence_candidates.tsv"),
  step(P("scripts","phase2","15_figures.R"), "results/phase2/validation/convergence_candidates.tsv", "results/phase2/figures/Fig10_Conceptual_Model.png"),
  step(P("scripts","phase2","16_validate.R"), character(), character()),
  step(P("scripts","phase2","20_upgrade_config.R"), "docs/METHODOLOGICAL_UPGRADE_2026.md", character()),
  step(P("scripts","phase2","21_tcga_adjacent_normal.R"), character(), c("results/phase2/upgrade_2026/matrices/tcga_thca_primary_normal_counts.rds","results/phase2/upgrade_2026/matrices/tcga_thca_primary_normal_metadata.tsv")),
  step(P("scripts","phase2","22_batch_correction.R"), c("results/phase2/upgrade_2026/matrices/tcga_thca_primary_normal_counts.rds","results/phase2/data/counts_matrix.rds"), c("results/phase2/upgrade_2026/matrices/batch_corrected_counts.rds","results/phase2/upgrade_2026/matrices/batch_corrected_logcpm.rds")),
  step(P("scripts","phase2","23_deconvolution.R"), "results/phase2/upgrade_2026/matrices/batch_corrected_logcpm.rds", c("results/phase2/upgrade_2026/deconvolution/composition_consensus.tsv","results/phase2/upgrade_2026/deconvolution/epic_convergence_status.tsv","results/phase2/upgrade_2026/matrices/composition_adjusted_logcpm.rds")),
  step(P("scripts","phase2","24_multiverse_de.R"), c("results/phase2/upgrade_2026/matrices/batch_corrected_counts.rds","results/phase2/upgrade_2026/matrices/composition_adjusted_logcpm.rds"), "results/phase2/upgrade_2026/differential_expression/DEG_summary_all_variants.tsv"),
  step(P("scripts","phase2","25_multiverse_gsea.R"), "results/phase2/upgrade_2026/differential_expression/DEG_summary_all_variants.tsv", "results/phase2/upgrade_2026/gsea/GSEA_fgsea_all_variants.tsv"),
  step(P("scripts","phase2","26_composition_variance.R"), c("results/phase2/upgrade_2026/deconvolution/composition_consensus.tsv","results/phase2/upgrade_2026/differential_expression/raw/voom_qw_full_results.tsv"), "results/phase2/upgrade_2026/deconvolution/top_deg_celltype_variance.tsv"),
  step(P("scripts","phase2","27_validation_meta_survival_surface.R"), c("results/phase2/upgrade_2026/differential_expression/tcga_matched/voom_qw_full_results.tsv","results/phase2/validation/validation_candidates.tsv"), c("results/phase2/upgrade_2026/validation/meta_analysis_REML.tsv","results/phase2/upgrade_2026/validation/ITGA2_prespecified_gate.tsv")),
  step(P("scripts","phase2","28_ppi_multinetwork.R"), c("results/phase2/upgrade_2026/gsea/GSEA_fgsea_all_variants.tsv","results/phase2/upgrade_2026/differential_expression/batch_corrected/voom_qw_full_results.tsv"), "results/phase2/upgrade_2026/ppi/PPI_upgrade_summary.tsv"),
  step(P("scripts","phase2","29_gsea_robustness.R"), c("results/phase2/upgrade_2026/matrices/composition_adjusted_logcpm.rds","results/phase2/data/panel30_genesets.rds"), "results/phase2/upgrade_2026/gsea/GSEA_robustness_summary.tsv"),
  step(P("scripts","phase2","30_render_upgrade_report.R"), c("results/phase2/upgrade_2026/differential_expression/DEG_summary_all_variants.tsv","results/phase2/upgrade_2026/validation/ITGA2_prespecified_gate.tsv"), "results/phase2/reports/METHODOLOGICAL_UPGRADE_2026.html"),
  step(P("scripts","phase2","31_validate_upgrade.R"), character(), "results/phase2/upgrade_2026/validation/VALIDATION_SUMMARY.tsv")
)

file_sig <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return("none")
  digest(list(path = paths, md5 = unname(tools::md5sum(paths))), algo = "sha256")
}
valid_outputs <- function(paths) !length(paths) || all(file.exists(paths) & file.info(paths)$size > 0)
record_benchmark <- function(script, status, started, signature, outputs) {
  benchmark_rows[[length(benchmark_rows) + 1L]] <<- data.table(
    invocation_id = invocation_id, step = basename(script), status = status,
    started_at = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    finished_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    input_signature_sha256 = signature,
    output_signature_sha256 = file_sig(outputs),
    outputs = paste(outputs, collapse = "|"))
  current <- benchmark_rows |> rbindlist(fill = TRUE)
  fwrite(rbindlist(list(benchmark_history, current), fill = TRUE), benchmark_file, sep = "\t", quote = FALSE)
}

for (s in steps) {
  num <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", basename(s$script))))
  if (!is.na(num) && num < from_num) next
  if (!file.exists(s$script)) stop("Missing step: ", s$script)
  missing_inputs <- s$inputs[!file.exists(s$inputs)]
  if (length(missing_inputs)) stop("Missing inputs for ", s$script, ": ", paste(missing_inputs, collapse = ", "))
  signature <- file_sig(c(s$script, s$inputs, "renv.lock"))
  state_file <- file.path(state_dir, paste0(basename(s$script), ".rds"))
  old <- if (file.exists(state_file)) readRDS(state_file) else NULL
  skip <- !force && valid_outputs(s$outputs) && !is.null(old) && identical(old$signature, signature)
  if (adopt && is.null(old) && valid_outputs(s$outputs)) skip <- TRUE
  if (skip) {
    started <- Sys.time()
    cat("[cache]", s$script, "\n")
    if (is.null(old)) saveRDS(list(signature = signature, completed = Sys.time(), adopted = TRUE), state_file)
    record_benchmark(s$script, if (is.null(old)) "adopted" else "cached", started, signature, s$outputs)
    next
  }
  cat("[run]", s$script, "\n")
  started <- Sys.time()
  status <- system2(RS, s$script)
  if (status != 0L || !valid_outputs(s$outputs)) {
    record_benchmark(s$script, "failed", started, signature, s$outputs)
    stop("Step failed or outputs invalid: ", s$script)
  }
  saveRDS(list(signature = signature, completed = Sys.time(), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))), state_file)
  record_benchmark(s$script, "completed", started, signature, s$outputs)
}
cat("Phase 2 + upgrade 2026 complete (0 failures).\n")
