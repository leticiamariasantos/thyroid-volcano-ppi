# ═══════════════════════════════════════════════════════════════════════════════
# Dockerfile — Ambiente computacional reprodutível (Fase 2)
# thyroid-volcano-ppi
#
# Atualizado para R 4.6.x + Bioconductor com os pacotes da Fase 2
# (limma, edgeR, DESeq2, fgsea, msigdbr, KEGGREST, org.Hs.eg.db, data.table,
#  igraph, ggrepel, ggplot2, recount3, hgu133plus2.db, GEOquery, SummarizedExperiment).
#
# Construir:
#   docker build -t thyroid-volcano-ppi .
# Executar pipeline completo da FASE 2 (reboot, painel de 30 vias):
#   docker run --rm -v "$(pwd):/work" -w /work thyroid-volcano-ppi Rscript scripts/run_phase2.R
# ═══════════════════════════════════════════════════════════════════════════════

FROM rocker/r-ver:4.6.1

LABEL maintainer="thyroid-volcano-ppi" \
      description="Ambiente reprodutível thyroid-volcano-ppi (Fase 2)" \
      version="2.0.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libfontconfig1-dev && apt-get clean && rm -rf /var/lib/apt/lists/*

# CRAN
RUN install2.r --error --skipinstalled \
    here data.table ggplot2 ggrepel igraph httr jsonlite dplyr tidyr readr \
    stringr tibble purrr pheatmap corrplot

# Bioconductor (Fase 2)
RUN R -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(c("limma","edgeR","DESeq2","fgsea","msigdbr","KEGGREST","org.Hs.eg.db","AnnotationDbi","hgu133plus2.db"), update=FALSE, ask=FALSE)'

WORKDIR /work
COPY . .

CMD ["Rscript", "scripts/run_phase2.R"]
