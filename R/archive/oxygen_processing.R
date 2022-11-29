rm(list = ls())

library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(deSolve)
library(pracma)
library(ncdf4)
library(here)
library(lubridate)

source('R/oxygen_helper.R')
source("R/thermal_script.R")
source('R/thermocline_helper.R')

working_folder = 'ObsDOTest'
lake_id = 'ob13083'
trophy = 'eutro'
iterations = 3
method = 'patankar-rk4'

hypsography_data <- get_hypsography(lake_id = lake_id,
                               working_folder = working_folder)

# thermal_data <- get_thermal_data(lake_id = lake_id,
#                                  working_folder = working_folder,
#                                  hypsography_data = hypsography_data)

thermal <- thermal_info(here(working_folder, lake_id))

thermal_to_model <- tibble(data = thermal %>%
         filter(!is.na(stratified)) %>%
         filter(duration > 2) %>%
         group_split(strat_id))
# count = stratification_batches

oxygen_output <- run_oxygen_model(thermal_data = thermal_to_model,
                                  method = method,
                                  trophy = trophy,
                                  iterations = iterations)

save_model_output(lake_id = lake_id,
                  working_folder = working_folder,
                  oxygen_output = oxygen_output)

create_plot(lake_id = lake_id,
            working_folder = working_folder,
            oxygen_output = oxygen_output)
