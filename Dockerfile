# ═══════════════════════════════════════════════════════════════════════════════
# Dockerfile — Ambiente computacional reprodutível (Fase 2 — reboot)
# thyroid-volcano-ppi
#
# Restaura o ambiente exato via renv::restore() (renv.lock com 165 pacotes),
# garantindo as mesmas versões usadas na análise.
#
# Construir:
#   docker build -t thyroid-volcano-ppi .
# Executar pipeline completo da FASE 2 (painel de 30 vias):
#   docker run --rm -v "$(pwd):/work" -w /work thyroid-volcano-ppi Rscript scripts/run_phase2.R
# ═══════════════════════════════════════════════════════════════════════════════

FROM rocker/r-ver:4.6.1

LABEL maintainer="thyroid-volcano-ppi" \
      description="Ambiente reprodutível thyroid-volcano-ppi (Fase 2, reboot)" \
      version="4.1.0"

# Dependências de sistema para compilação de pacotes (graphics, xml, ssl, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libfontconfig1-dev && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Restaurar o ambiente exato a partir do lockfile (versões fixadas)
COPY renv.lock renv.lock
COPY renv/ renv/
RUN R -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv"); renv::restore()'

# Copiar o restante do repositório
COPY . .

CMD ["Rscript", "scripts/run_phase2.R"]
