suppressPackageStartupMessages({ library(here); library(rmarkdown) })
source(here("scripts", "phase2", "20_upgrade_config.R"))
out <- rmarkdown::render(here("reports", "methodological_upgrade_2026.Rmd"),
  output_file = "METHODOLOGICAL_UPGRADE_2026.html", output_dir = DIR_REP,
  envir = new.env(parent = globalenv()), quiet = TRUE)
if (!file.exists(out) || file.info(out)$size == 0) stop("Report rendering failed")
log_msg("Upgrade report rendered:", out)
