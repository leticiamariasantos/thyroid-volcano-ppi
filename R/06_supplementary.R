# ═══════════════════════════════════════════════════════════════════════════════
# R/06_supplementary.R — Supplementary Tables S1–S4
# thyroid-volcano-ppi
#
# Automatically generates supplementary tables for scientific submission.
# Run after the main pipeline (00_setup.R + 01_functions.R already loaded):
#   source("R/06_supplementary.R")
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n── M5: Supplementary Tables ──\n")

out_dir <- here::here("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_sup <- function(x, name) {
  f <- file.path(out_dir, name)
  readr::write_tsv(x, f)
  cat(sprintf("  → %s (%d×%d)\n", basename(f), nrow(x), ncol(x)))
}

# ── S1: Computational Environment ─────────────────────────────────────────────
pkg_info <- tibble::tribble(
  ~package,          ~main_function,                           ~source,        ~required,
  "limma",           "Differential expression analysis (DEG)", "Bioconductor", "Yes",
  "ggplot2",         "Visualization (Volcano Plot)",           "CRAN",         "Yes",
  "ggrepel",         "Repulsive text labels",                  "CRAN",         "Yes",
  "igraph",          "Network construction and analysis",      "CRAN",         "Yes",
  "ggraph",          "Network visualization",                  "CRAN",         "Yes",
  "dplyr",           "Data manipulation",                      "CRAN",         "Yes",
  "tidyr",           "Data tidying",                           "CRAN",         "Yes",
  "readr",           "TSV read/write",                         "CRAN",         "Yes",
  "tibble",          "Data structures",                        "CRAN",         "Yes",
  "stringr",         "String manipulation",                    "CRAN",         "Yes",
  "purrr",           "Functional programming",                 "CRAN",         "Yes",
  "here",            "Path portability",                       "CRAN",         "Yes",
  "KEGGREST",        "KEGG API (signaling pathway)",           "Bioconductor", "Yes",
  "org.Hs.eg.db",    "Human gene annotation",                  "Bioconductor", "Yes",
  "AnnotationDbi",   "Annotation interface",                   "Bioconductor", "Yes",
  "httr",            "HTTP requests (STRING API)",             "CRAN",         "Yes",
  "jsonlite",        "JSON parsing (STRING API)",              "CRAN",         "Yes"
)

inst <- as.data.frame(utils::installed.packages()[, c("Package", "Version")],
                      stringsAsFactors = FALSE)
names(inst) <- c("package", "version")

s1 <- pkg_info |>
  dplyr::left_join(inst, by = "package") |>
  dplyr::select(package, version, main_function, source, required)

write_sup(s1, "S1_computational_environment.tsv")
cat("  S1: Computational environment —", nrow(s1), "packages\n")

# Also write R version and platform info
si <- sessionInfo()
s1_env <- tibble::tibble(
  property = c("R version", "Platform", "OS", "System"),
  value = c(
    paste0(si$R.version$major, ".", si$R.version$minor),
    si$platform,
    si$running,
    Sys.info()[["sysname"]]
  )
)
write_sup(s1_env, "S1b_system_info.tsv")
cat("  S1b: System info written\n")

# ── S2: External Databases ────────────────────────────────────────────────────
org_ver <- inst$version[inst$package == "org.Hs.eg.db"]
if (length(org_ver) == 0) org_ver <- "N/A"
ann_ver <- inst$version[inst$package == "AnnotationDbi"]
if (length(ann_ver) == 0) ann_ver <- "N/A"

s2 <- tibble::tribble(
  ~database,               ~identifier,                          ~version,      ~access_date,  ~description,
  "KEGG PATHWAY",          "hsa04919",                           "112.0",       as.character(Sys.Date()), "Thyroid hormone signaling pathway",
  "STRING",                "Homo sapiens (9606)",                "12.0",        as.character(Sys.Date()), "PPI network — high confidence (score ≥700)",
  "TCGA + GTEx (Xena)",    "c486b845ee2e750c3a9d2fc5145c8426",  "TOIL RSEM",   as.character(Sys.Date()), "Gene expression log2(norm+1), 783 samples",
  "org.Hs.eg.db",          "Homo sapiens",                       org_ver,       as.character(Sys.Date()), "Entrez/Symbol gene annotations",
  "AnnotationDbi",         "—",                                  ann_ver,       as.character(Sys.Date()), "Bioconductor annotation interface"
)

write_sup(s2, "S2_external_databases.tsv")
cat("  S2: External databases —", nrow(s2), "databases\n")

# ── S3: Pipeline Audit Trail ──────────────────────────────────────────────────
s3 <- tibble::tribble(
  ~step,              ~script,            ~input,                        ~output,                                    ~dependencies,          ~validation,
  "Setup",            "00_setup.R",       "—",                          "Parameters, directories, packages",         "here, BiocManager",    "sessionInfo()",
  "Import & QC",      "02_import.R",      "XENA_THCA.tsv",              "Expression matrix + metadata",              "readr, dplyr",         "log2 scale, NAs, duplicates",
  "DEG (limma)",      "03_deg.R",         "Expression matrix",          "DEG results (119 genes)",                   "limma",                "eBayes, BH, thresholds",
  "Volcano Plot",     "04_volcano.R",     "DEG results",                "Fig1_Volcano_THCA_vs_Normal.png",           "ggplot2, ggrepel",     "seed=42, parameters",
  "PPI Network",      "05_ppi.R",         "DEG list",                   "Fig2_PPI_Network_THCA_DEGs.png",            "igraph, ggraph, httr", "STRING≥700, walktrap",
  "Supplementary",    "06_supplementary.R","sessionInfo() + results",    "S1–S4 (TSV)",                               "dplyr, readr",         "Manual cross-check"
)

write_sup(s3, "S3_pipeline_audit_trail.tsv")
cat("  S3: Audit trail —", nrow(s3), "steps\n")

# ── S4: AI-Assisted Tasks ─────────────────────────────────────────────────────
# Declaração conforme Portaria CNPq nº 2.664/2026
s4 <- tibble::tribble(
  ~task,                       ~tool,              ~human_participation,              ~validation_method,              ~status,
  "R code auditing",           "DeepSeek-v4-pro",  "Line-by-line review",             "Manual pipeline execution",     "Compliant",
  "Namespace corrections",     "DeepSeek-v4-pro",  "Verification and approval",       "Zero-warning execution",        "Compliant",
  "Statistical function review","DeepSeek-v4-pro", "Methodological validation",       "limma documentation cross-check","Compliant",
  "Script generation/debugging","Codex",           "Code review and testing",         "Pipeline execution validation", "Compliant",
  "Technical documentation",   "Codex",            "Author review and editing",       "Manual proofreading",           "Compliant",
  "Textual review",            "ChatGPT 5.5",      "Writing, review, and approval",   "Author review",                 "Compliant",
  "Documentation structure",   "ChatGPT 5.5",      "Content organization review",     "Cross-referencing validation",  "Compliant",
  "Reproducibility suggestions","ChatGPT 5.5",     "Implementation review",           "renv + Docker validation",      "Compliant",
  "Exploratory analysis",      "Grok (xAI)",       "Review and approval",             "Statistical re-validation",     "Compliant",
  "Volcano Plot prototyping",  "Grok (xAI)",       "Design review and final approval","Visual inspection + parameters","Compliant",
  "PPI Network prototyping",   "Grok (xAI)",       "Design review and final approval","Visual inspection + metrics",   "Compliant",
  "Methodological support",    "Grok (xAI)",       "Review and approval",             "Protocol cross-check",          "Compliant",
  "Supplementary tables",      "Grok (xAI)",       "Content review",                  "Cross-referencing",             "Compliant"
)

write_sup(s4, "S4_ai_assisted_tasks.tsv")
cat("  S4: AI-assisted tasks —", nrow(s4), "tasks\n")

cat("  Supplementary tables complete.\n")
