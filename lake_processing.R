library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(lubridate)
library(here)
library(purrr)
library(ncdf4)

source(here("R/thermal_script.R"))

lake_folder <- "InputTest"
lakes <- list.files(here(lake_folder), full.names = T)


lakes %>% 
  walk(~{
    lake <- .x
    
    # writes thermal_info.csv into lake folder
    lake %>% thermal_info()
    
  })
