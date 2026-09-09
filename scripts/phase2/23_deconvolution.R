# xCell + EPIC formal deconvolution; auditable multi-signature fallback.
suppressPackageStartupMessages({
  library(here); library(data.table); library(GSVA); library(limma)
  # EPIC 1.1.8 resolves its bundled TRef object through "package:EPIC".
  # Attaching the locked package is therefore an explicit runtime contract.
  library(EPIC); library(xCell)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

expr_log <- readRDS(file.path(DIR_VARIANTS, "batch_corrected_logcpm.rds"))
meta <- fread(file.path(DIR_VARIANTS, "analysis_metadata.tsv"))
stopifnot(identical(colnames(expr_log), meta$sample))
expr_linear <- pmax(2^expr_log - 2, 0)

run_xcell <- function() {
  # xCell 1.1.0 ships its reference objects as package data, but the exported
  # function refers to `xCell.data` without loading it. Load the locked package
  # data explicitly and pass every reference object to avoid hidden global-state
  # dependence (and the runtime error "object 'xCell.data' not found").
  xcell_env <- new.env(parent = emptyenv())
  utils::data("xCell.data", package = "xCell", envir = xcell_env)
  if (!exists("xCell.data", envir = xcell_env, inherits = FALSE)) {
    stop("Locked xCell installation does not contain the xCell.data reference set")
  }
  xcell_ref <- get("xCell.data", envir = xcell_env, inherits = FALSE)
  required <- c("signatures", "genes", "spill")
  if (!all(required %in% names(xcell_ref))) {
    stop("xCell.data is incomplete; missing: ", paste(setdiff(required, names(xcell_ref)), collapse = ", "))
  }
  x <- xCell::xCellAnalysis(
    expr_linear,
    signatures = xcell_ref$signatures,
    genes = xcell_ref$genes,
    spill = xcell_ref$spill,
    parallel.sz = 2
  )
  x <- as.data.table(x, keep.rownames = "cell_type")
  melt(x, id.vars = "cell_type", variable.name = "sample", value.name = "estimate")[,
    `:=`(method = "xcell", status = "formal_relative_enrichment")]
}
run_epic <- function() {
  convergence_warnings <- character()
  fit <- withCallingHandlers(
    EPIC::EPIC(bulk = expr_linear, reference = "TRef"),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("optimization didn't fully converge", msg, fixed = TRUE)) {
        convergence_warnings <<- c(convergence_warnings, msg)
        invokeRestart("muffleWarning")
      }
    }
  )
  fractions <- fit$cellFractions
  if (!setequal(rownames(fractions), meta$sample)) {
    stop("EPIC output sample identifiers do not match analysis metadata")
  }
  fractions <- fractions[meta$sample, , drop = FALSE]
  warning_text <- paste(convergence_warnings, collapse = "\n")
  nonconverged <- meta$sample[vapply(meta$sample, function(id) {
    grepl(id, warning_text, fixed = TRUE)
  }, logical(1))]
  if (length(convergence_warnings) && !length(nonconverged)) {
    stop("EPIC reported nonconvergence but affected sample IDs could not be parsed")
  }
  convergence <- data.table(
    sample = meta$sample,
    epic_converged = !(meta$sample %in% nonconverged),
    status = fifelse(meta$sample %in% nonconverged,
                     "formal_fraction_nonconverged", "formal_fraction_converged")
  )
  x <- as.data.table(t(fractions), keep.rownames = "cell_type")
  estimates <- melt(x, id.vars = "cell_type", variable.name = "sample", value.name = "estimate")
  estimates <- merge(estimates, convergence[, .(sample, status)], by = "sample", all.x = TRUE,
                     sort = FALSE)
  estimates[, method := "epic"]
  list(estimates = estimates, convergence = convergence)
}
if (!requireNamespace("xCell", quietly = TRUE) || !requireNamespace("EPIC", quietly = TRUE))
  stop("Formal deconvolution requires locked packages xCell and EPIC")
xcell_result <- run_xcell()
epic_result <- run_epic()
fwrite_tsv(epic_result$convergence, file.path(DIR_DECONV, "epic_convergence_status.tsv"))
formal <- rbindlist(list(xcell_result, epic_result$estimates), fill = TRUE)

signatures <- list(
  thyroid_epithelial = c("TG", "TPO", "SLC5A5", "TSHR", "PAX8", "FOXE1", "KRT8", "KRT18"),
  fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2"),
  immune = c("PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "LYZ", "FCER1G"),
  muscle = c("ACTA1", "MYH7", "MYL1", "MYL2", "TNNT3", "CKM", "DES")
)
signatures <- lapply(signatures, intersect, y = rownames(expr_log))
if (utils::packageVersion("GSVA") >= "2.0.0") {
  param <- GSVA::ssgseaParam(expr_log, signatures, normalize = TRUE)
  fallback_scores <- GSVA::gsva(param, verbose = FALSE)
} else {
  fallback_scores <- GSVA::gsva(expr_log, signatures, method = "ssgsea", verbose = FALSE)
}
fallback <- as.data.table(t(fallback_scores), keep.rownames = "sample")
fallback_long <- melt(fallback, id.vars = "sample", variable.name = "cell_type", value.name = "estimate")
fallback_long[, `:=`(method = "multi_signature_ssgsea", status = "fallback_reference_score")]
all_long <- rbindlist(list(formal, fallback_long), fill = TRUE)
fwrite_tsv(all_long, file.path(DIR_DECONV, "cell_composition_long.tsv"))

compartment <- function(x) {
  fcase(
    grepl("fibro|strom|CAF", x, ignore.case = TRUE), "fibroblast",
    grepl("endothel", x, ignore.case = TRUE), "endothelial",
    grepl("muscle|myocyte", x, ignore.case = TRUE), "muscle",
    grepl("thyroid|epithelial", x, ignore.case = TRUE), "thyroid_epithelial",
    grepl("otherCells|cancer|keratino|sebocyte|hepatocyte|neuron", x, ignore.case = TRUE), "other_parenchymal",
    grepl("(^|[._ -])(B|T|NK)[._ -]?cells?|immune|lymph|monocyte|macrophage|dendritic|DC$|neutrophil|eosinophil|basophil|mast|plasma|megakaryocyte", x, ignore.case = TRUE), "immune",
    default = NA_character_
  )
}
all_long[, compartment := compartment(cell_type)]
mapping <- unique(all_long[, .(method, cell_type, compartment)])
fwrite_tsv(mapping, file.path(DIR_DECONV, "celltype_compartment_mapping.tsv"))
all_long <- all_long[!is.na(compartment)]
all_long[, z := as.numeric(scale(estimate)), by = .(method, cell_type)]
cons <- all_long[is.finite(z), .(score = median(z)), by = .(sample, compartment)]
cons_wide <- dcast(cons, sample ~ compartment, value.var = "score")
comp_mat <- as.data.frame(cons_wide[, -1]); rownames(comp_mat) <- cons_wide$sample
saveRDS(comp_mat, file.path(DIR_DECONV, "composition_consensus.rds"))
fwrite_tsv(cons_wide, file.path(DIR_DECONV, "composition_consensus.tsv"))

adjusted <- residualize_preserving_condition(expr_log, meta, comp_mat)
saveRDS(adjusted, file.path(DIR_VARIANTS, "composition_adjusted_logcpm.rds"))
method_status <- data.table(method = c("xCell", "EPIC", "multi_signature_ssgsea"),
  completed = c(any(formal$method == "xcell"), any(formal$method == "epic"), TRUE),
  samples = nrow(meta),
  nonconverged_samples = c(NA_integer_, sum(!epic_result$convergence$epic_converged), NA_integer_),
  interpretation = c("relative enrichment", "estimated fractions; convergence audited per sample",
                     "reference score; not a fraction"))
fwrite_tsv(method_status, file.path(DIR_DECONV, "deconvolution_status.tsv"))
log_msg("Deconvolution complete. Formal rows:", nrow(formal),
        "; EPIC nonconverged samples:", sum(!epic_result$convergence$epic_converged))
