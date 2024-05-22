library(here)
library(dplyr)
library(qs)
library(future)
library(furrr)
library(purrr)
library(lubridate)


cl <- makeClusterPSOCK(
  availableWorkers(),
  rscript_args = c("-e", shQuote('setwd("/home/kelp2501/EezyPeezyISIOxy")')),
  rscript_libs = .libPaths(),
  # revtunnel = FALSE,
  homogeneous = TRUE,
  setup_strategy = "parallel"
)

plan(cluster, workers = cl)


source(here("R/tidy_data.R"))

tidy_annual_min_o2()