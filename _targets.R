# Load packages required to define the pipeline:
library(targets)
library(tibble)
library(here)
library(tarchetypes) # Load other packages as needed. # nolint

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

tar_config_set(store = "~/scratch/isioxy/")
tar_option_set(storage = "worker", retrieval = "worker")

tar_option_set(
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(memory = "2G", time = "01:00:00"))
  )
)

# tar_make_clustermq() configuration:
# options(clustermq.scheduler = "multicore")
options(clustermq.scheduler = "slurm")
options(clustermq.template = "clustermq.tmpl")

source("R/thermal_script.R")
source("R/thermocline_helper.R")
source("R/oxygen_helper.R")

# target list:
lake_folder <- "ObsDOTest"

lakes <- tibble(
 lake_id = list.files(here(lake_folder), full.names = F)[1]
)

numit <- 1000

# tarobservations <- tar_target(observations, "data/observed.rds", format = "file")

targets <- tar_map(
  values = lakes,
  tar_target(thermal, thermal_info(here("ObsDOTest", lake_id))),
  # tar_target(oxygen, run_oxygen_model(thermal, method = 'rk4', trophy = 'oligo',
                                      # iterations = numit)),
  # tar_target(plot_oxygen, create_plot(oxygen, lake_id, "ObsDOTest"), format = "file"),
  tar_group_count(
    thermal_to_model,
    tibble(data = thermal %>%
      filter(!is.na(stratified)) %>%
      filter(duration > 2) %>% 
      group_split(strat_id)),
    count = 10
    # pattern = map(index_batch)
    #   group_by(strat_id) %>%
    #   tar_group(),
    # iteration = "group"
  ),
  tar_target(oxygen,
             run_oxygen_model(
               thermal_to_model,
                method = "rk4",
                trophy = "oligo",
                iterations = numit
              ),
      pattern = map(thermal_to_model)
    ),
  tar_target(write_oxygen, save_model_output(oxygen, lake_id), format = "file"),
  tar_target(observations, read_observations(lake_id)),
  tar_target(plot_qc_oxy, save_qc_plot_oxygen(oxygen, lake_id, observations), format = "file"),
  tar_target(plot_thermal, create_plots_thermal(thermal, lake_id, "ObsDOTest"), format = "file")
)

list(targets)
