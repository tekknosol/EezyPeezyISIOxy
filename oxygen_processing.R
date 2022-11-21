rm(list = ls())

library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(deSolve)
library(pracma)
library(lubridate)

source('R/oxygen_helper.R')

working_folder = 'ObsDOTest'
lake_id = 'ob490'


hypsography_data <- get_hypsography(lake_id = lake_id,
                               working_folder = working_folder)

thermal_data <- get_thermal_data(lake_id = lake_id,
                                 working_folder = working_folder,
                                 hypsography_data = hypsography_data)

oxygen_output <- run_oxygen_model(thermal_data = thermal_data,
                                  method = 'rk4')

save_model_output(lake_id = lake_id,
                  working_folder = working_folder,
                  oxygen_output = oxygen_output)

create_plot(lake_id = lake_id,
            working_folder = working_folder,
            oxygen_output = oxygen_output)
