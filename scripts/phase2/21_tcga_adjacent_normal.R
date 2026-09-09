# Acquire all uniformly processed recount3 TCGA-THCA samples and retain codes 01/11.
suppressPackageStartupMessages({
  library(here); library(data.table); library(recount3); library(SummarizedExperiment)
})
source(here("scripts", "phase2", "20_upgrade_config.R"))

log_msg("══ TCGA-THCA tumor/normal adjacente (recount3) ══")
cache <- file.path(DIR_CACHE, "recount3", "v1.22.0", "TCGA_THCA_all_gene_counts.rds")

rse <- safe_api_cache(cache, function() {
  projects <- recount3::available_projects()
  hit <- projects[projects[["project"]] == "THCA" & projects[["file_source"]] == "tcga", ]
  if (nrow(hit) != 1L) stop("Expected one recount3 TCGA/THCA project; found ", nrow(hit))
  recount3::create_rse(hit)
}, validate = function(x) inherits(x, "SummarizedExperiment") && ncol(x) >= 550L)

barcodes <- as.character(colData(rse)[["tcga.tcga_barcode"]])
if (is.null(barcodes)) stop("recount3 metadata lacks tcga.tcga_barcode")
codes <- sample_code(barcodes)
keep <- codes %in% c("01", "11")
counts <- assay(rse, "raw_counts")[, keep, drop = FALSE]
symbols <- rowData(rse)[["gene_name"]]
if (is.null(symbols)) symbols <- rowData(rse)[["gene_id"]]
counts <- collapse_counts(counts, symbols)
colnames(counts) <- tcga_key(barcodes[keep])

meta <- data.table(
  sample = colnames(counts),
  participant = substr(colnames(counts), 1L, 12L),
  sample_code = codes[keep],
  condition = factor(ifelse(codes[keep] == "01", "Tumor", "Normal"), levels = c("Normal", "Tumor")),
  source = "TCGA",
  cohort = "TCGA_THCA_recount3"
)
if (meta[condition == "Normal", .N] < 20L) stop("Too few adjacent normals: ", meta[condition == "Normal", .N])
if (meta[condition == "Tumor", .N] < 450L) stop("Too few primary tumors: ", meta[condition == "Tumor", .N])
if (anyDuplicated(meta$sample)) stop("Duplicated TCGA sample keys after barcode normalization")

paired_ids <- meta[, .(has_tumor = any(condition == "Tumor"), has_normal = any(condition == "Normal")), by = participant][has_tumor & has_normal, participant]
meta[, paired := participant %in% paired_ids]

saveRDS(counts, file.path(DIR_VARIANTS, "tcga_thca_primary_normal_counts.rds"))
fwrite_tsv(meta, file.path(DIR_VARIANTS, "tcga_thca_primary_normal_metadata.tsv"))
write_manifest(file.path(DIR_VARIANTS, "tcga_thca_primary_normal_manifest.tsv"), list(
  source = "recount3 TCGA THCA", assay = "raw_counts", annotation = "GENCODE G026",
  tumor_code = "01", normal_code = "11", n_tumor = meta[condition == "Tumor", .N],
  n_normal = meta[condition == "Normal", .N], n_paired_participants = length(paired_ids),
  matrix_fingerprint_sha256 = matrix_checksum(counts), acquired = as.character(Sys.Date())
))
log_msg("TCGA-matched:", meta[condition == "Tumor", .N], "tumores;",
        meta[condition == "Normal", .N], "normais;", length(paired_ids), "participantes pareados")
