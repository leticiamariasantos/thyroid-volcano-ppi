# ═══════════════════════════════════════════════════════════════════════════════
# tests/testthat.R — Test runner (testthat)
# thyroid-volcano-ppi
#
# Usage:
#   testthat::test_dir("tests/testthat")
# ═══════════════════════════════════════════════════════════════════════════════

library(testthat)
library(here)

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_functions.R"))

test_dir(here::here("tests", "testthat"))
