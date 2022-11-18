library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(lubridate)
library(here)
library(purrr)
library(lubridate)
library(pracma)
library(tictoc)
library(ncdf4)


source(here("R/thermal_script.R"))
source(here("R/thermocline_helper.R"))

lake_folder <- "InputTest"
lakes <- list.files(here(lake_folder), full.names = T)

lakes %>% 
  walk(~{
    lake <- .x
    
    # writes thermal_info.csv into lake folder
    tic(msg = basename(lake))
    lake %>% thermal_info()
    toc()
})

