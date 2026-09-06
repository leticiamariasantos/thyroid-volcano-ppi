suppressPackageStartupMessages({ library(recount3); library(data.table) })

# Explora a estrutura dos RSE recount3 para TCGA THCA e GTEx THYROID
ap <- available_projects()

# TCGA THCA
pi_tcga <- ap[ap$project == "THCA" & ap$file_source == "tcga", ]
cat("=== TCGA THCA ===\n"); print(pi_tcga)
rse_t <- create_rse(pi_tcga)
cat("TCGA dims:", dim(rse_t), "\n")
cat("assays:", assayNames(rse_t), "\n")
cat("colData cols:", paste(head(colnames(colData(rse_t)), 40), collapse=", "), "\n")
cat("rowData cols:", paste(head(colnames(rowData(rse_t)), 15), collapse=", "), "\n")
cat("sample barcode head:\n"); print(head(colData(rse_t)$tcga.tcga_barcode))

cat("\n=== GTEx THYROID ===\n")
pi_gt <- ap[ap$project == "THYROID" & ap$file_source == "gtex", ]
print(pi_gt)
rse_g <- create_rse(pi_gt)
cat("GTEx dims:", dim(rse_g), "\n")
cat("colData cols:", paste(head(colnames(colData(rse_g)), 40), collapse=", "), "\n")
cat("sample sampid head:\n"); print(head(colData(rse_g)$gtex.sampid))
