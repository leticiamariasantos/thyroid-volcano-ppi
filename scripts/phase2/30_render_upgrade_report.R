suppressPackageStartupMessages({ library(here); library(rmarkdown) })
source(here("scripts", "phase2", "20_upgrade_config.R"))
if (!rmarkdown::pandoc_available()) {
  for (p in c(file.path(Sys.getenv("LOCALAPPDATA"), "Pandoc"), "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")) {
    if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); if (rmarkdown::pandoc_available()) break }
  }
}
out <- rmarkdown::render(here("reports", "methodological_upgrade_2026.Rmd"),
  output_file = "METHODOLOGICAL_UPGRADE_2026.html", output_dir = DIR_REP,
  envir = new.env(parent = globalenv()), quiet = TRUE)
if (!file.exists(out) || file.info(out)$size == 0) stop("Report rendering failed")
log_msg("Upgrade report rendered:", out)
