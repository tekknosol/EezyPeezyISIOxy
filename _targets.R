# Load packages required to define the pipeline:
library(targets)
library(tidyr)
library(tibble)
library(here)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c(
    "tibble", 
    "ncdf4", 
    "readr", 
    "dplyr", 
    "rLakeAnalyzer", 
    "pracma", 
    "LakeMetabolizer", 
    "deSolve",
    "tidyr",
    "lubridate",
    "ggplot2",
    "stringr",
    "here"
  ), # packages that your targets need to run
  # format = "rds" # default storage format
  format = "qs" # default storage format
)

tar_config_set(store = "~/scratch/isioxy/") # Folder for target's internal data storage
tar_option_set(storage = "worker", retrieval = "worker") # let workers store and retrieve data directely

# Slurm configs
tar_option_set(
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(memory = "2G", time = "01:00:00"))
  )
)

# tar_make_clustermq() configuration:
# options(clustermq.scheduler = "multicore") # parallel processing on local machine
options(clustermq.scheduler = "slurm") # Slurm on HPC
options(clustermq.template = "clustermq.tmpl") # Slurm sbatch template

# source required functions
source("R/thermal_script.R")
source("R/thermocline_helper.R")
source("R/oxygen_helper.R")

# settings for computations:
lake_folder <- "ObsDOTest" # Folder containing isimip results
numit <- 1000 # number of iterations for oxygen model
stratification_batches <- 5 # Number of batches of stratification events per lake

# Total number of targets for computation: lakes * batches

# target list:
# lakes <- expand_grid(
#  lake_id = list.files(here(lake_folder), full.names = F)[1],
#  trophy = c("oligo", "eutro")
# )

lakes <- tibble(
  lake_id = list.files(here(lake_folder), full.names = F)
)

glob_trophy <- tar_target(trophy, c("oligo", "eutro"))
glob_params <- tar_target(oxy_params, get_prior(trophy, n = numit), pattern = trophy)



targets <- tar_map(
  values = lakes,

  #process thermal module per lake
  tar_target(thermal, thermal_info(here(lake_folder, lake_id))),

  # batch computation based on stratification events. Number of batches defined by 'stratification_batches'
  tar_group_count(
    thermal_to_model,
    tibble(data = thermal %>%
      filter(!is.na(stratified)) %>%
      filter(duration > 2) %>%
      group_split(strat_id)),
    count = stratification_batches
  ),

  # run oxygen model
  tar_target(oxygen,
             run_oxygen_model(
               thermal_to_model,
                method = "rk4",
                trophy = trophy,
                iterations = numit,
               params = oxy_params
              ),
      pattern = cross(thermal_to_model, trophy)
    ),

  # store results of oxygen model in results folder
  tar_target(write_oxygen, save_model_output(oxygen, lake_id), format = "file"),
  # load observations (Abby's lakes)
  tar_target(observations, read_observations(lake_id)),
  # Create QC plots for oxygen
  tar_target(plot_qc_oxy, save_qc_plot_oxygen(oxygen, lake_id, observations), format = "file"),
  # Create plots of temperature and thermocline depth
  tar_target(plot_thermal, create_plots_thermal(thermal, lake_id, lake_folder), format = "file")
)

list(glob_trophy, glob_params, targets)