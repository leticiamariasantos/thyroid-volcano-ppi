# Script para inicializar renv
if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv", repos = "https://cloud.r-project.org")

renv::init(bare = TRUE, force = TRUE, restart = FALSE)

# Instalar pacotes CRAN
renv::install(c(
  "here", "dplyr", "stringr", "tibble", "tidyr", "purrr", "readr",
  "ggplot2", "ggrepel", "igraph", "ggraph", "httr", "jsonlite"
), prompt = FALSE)

# Instalar Bioconductor e pacotes
if (!requireNamespace("BiocManager", quietly = TRUE))
  renv::install("BiocManager", prompt = FALSE)
BiocManager::install(
  c("limma", "KEGGREST", "org.Hs.eg.db", "AnnotationDbi"),
  update = FALSE, ask = FALSE
)

# Salvar lockfile
renv::snapshot(prompt = FALSE)
cat("renv configurado com sucesso.\n")
