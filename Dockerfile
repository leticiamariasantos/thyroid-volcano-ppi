# ═══════════════════════════════════════════════════════════════════════════════
# Dockerfile — Ambiente computacional reprodutível
# thyroid-volcano-ppi
#
# Construir:
#   docker build -t thyroid-volcano-ppi .
#
# Executar:
#   docker run --rm -v "$(pwd)/results:/home/rstudio/results" thyroid-volcano-ppi
#
# A imagem inclui R 4.6.0 com todos os pacotes CRAN e Bioconductor necessários.
# ═══════════════════════════════════════════════════════════════════════════════

FROM rocker/r-ver:4.4.0

# Metadados
LABEL maintainer="Ryan de Paulo Santos <ryan.santos@example.com>" \
      description="Ambiente reprodutível para thyroid-volcano-ppi" \
      version="1.0.0"

# Instalar dependências de sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Instalar pacotes CRAN
RUN install2.r --error --skipinstalled \
    here \
    dplyr \
    stringr \
    tibble \
    tidyr \
    purrr \
    readr \
    ggplot2 \
    ggrepel \
    igraph \
    ggraph \
    httr \
    jsonlite

# Instalar Bioconductor e pacotes
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("limma", "KEGGREST", "org.Hs.eg.db", "AnnotationDbi"), update = FALSE, ask = FALSE)'

# Criar diretórios do projeto
RUN mkdir -p /home/rstudio/thyroid-volcano-ppi/data/raw \
    /home/rstudio/thyroid-volcano-ppi/results/figures \
    /home/rstudio/thyroid-volcano-ppi/results/tables \
    /home/rstudio/thyroid-volcano-ppi/results/network \
    /home/rstudio/thyroid-volcano-ppi/logs

WORKDIR /home/rstudio/thyroid-volcano-ppi

# Copiar código-fonte
COPY R/         R/
COPY scripts/   scripts/
COPY run_pipeline.R .

# Instruções
CMD ["Rscript", "run_pipeline.R"]
