# ---- libraries --------------
library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(lubridate)
library(here)
library(purrr)
library(lubridate)
library(pracma)
library(ggpubr)
library(future)
library(furrr)
library(ncdf4)
library(tictoc)

source(here("R/thermal_script.R"))
source(here("R/thermocline_helper.R"))

# parallel processing - comment out next line for sequential run
plan(multisession, workers = 5)


# ---- script -----------------

lake_folder <- "ObsDOTest"
lakes <- list.files(here(lake_folder), full.names = T)

lakes %>% 
  future_walk(~{
    lake <- .x
    
    # writes thermal_info into lake folder
    tictoc::tic(msg = basename(lake))
    lake %>% thermal_info()
    tictoc::toc()
})

