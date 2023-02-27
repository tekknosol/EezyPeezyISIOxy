# Load packages required to define the pipeline:
library(targets)
library(tidyr)
library(tibble)
library(here)
library(tarchetypes)
library(compiler)
library(readr)
library(dplyr)
library(stringr)

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
    "here", 
    "purrr",
    "hydroGOF",
    "matrixStats",
    "oxypatankr"
  ), # packages that your targets need to run
  # format = "rds" # default storage format
  format = "qs" # default storage format
)

tar_config_set(store = "~/scratch/isioxy/") # Folder for target's internal data storage
tar_option_set(storage = "worker", retrieval = "worker") # let workers store and retrieve data directely

# Slurm configs
tar_option_set(
  resources = tar_resources(
    clustermq = tar_resources_clustermq(template = list(memory = "4G", time = "04:00:00"))
  )
)

# tar_make_clustermq() configuration:
options(clustermq.scheduler = "multicore") # parallel processing on local machine
# options(clustermq.scheduler = "slurm") # Slurm on HPC
options(clustermq.template = "clustermq.tmpl") # Slurm sbatch template

# source required functions from R subfolder
tar_source()

# settings for computations:
lake_folder <- "ObsDOTest" # Folder containing isimip results
numit <- 1000 # number of iterations for oxygen model
stratification_batches <- 1 # Number of batches of stratification events per lake

# Total number of targets for computation: lakes * batches

# target list:
# lakes <- expand_grid(
#  lake_id = list.files(here(lake_folder), full.names = F)[1],
#  trophy = c("oligo", "eutro")
# )

rafalakes <- read_delim(file = here("data/IDs_1000_lakes.csv"), delim = " ")
id_lookup <- read_csv(here("data/coord_area_depth.csv"))

id_lookup <- id_lookup %>% 
  mutate(isimip_id = str_split(id, "_", simplify = T)[,2])

rafalakes <- rafalakes %>% 
  left_join(
    id_lookup %>% 
      select(hydrolakes, isimip_id), 
    by = c("ID_vector" = "hydrolakes")
  )

lakes <- tibble(
  # lake_id = list.files(here(lake_folder), full.names = F)[1:1]
  lake_id = rafalakes$isimip_id[1:1]
)

glob_trophy <- tar_target(trophy, c("oligo", "eutro"), deployment = "main")
# glob_methods <- tar_target(methods, c("rk4", "rk4_zero", "patankar-rk2"), deployment = "main")
# glob_methods <- tar_target(methods, c("patankar-rk2_c", "lsoda_event", "lsoda"), deployment = "main")
glob_methods <- tar_target(methods, c("patankar-rk2_c"), deployment = "main")
# glob_methods <- tar_target(methods, c("rk4", "rk4_zero"), deployment = "main")
glob_params <- tar_target(oxy_params, get_prior(trophy, n = numit), pattern = trophy, deployment = "main")



targets <- tar_map(
  unlist = FALSE, # Return a nested list from tar_map()
  values = lakes,

  #process thermal module per lake
  tar_target(thermal, thermal_info(lake_id)),

  # batch computation based on stratification events. Number of batches defined by 'stratification_batches'
  # tar_group_count(
  #   thermal_to_model,
  #   tibble(data = thermal %>%
  #     filter(!is.na(stratified)) %>%
  #     filter(duration > 2) %>%
  #     group_split(strat_id)),
  #   count = stratification_batches
  # ),

  # run oxygen model
  tar_target(oxygen,
             run_oxygen_model(
               thermal,
               lake_id,
                method = methods, # rk4, patankar-rk2
                trophy = trophy,
                iterations = numit,
               params = oxy_params
              ),
      pattern = cross(trophy, methods)
    ),

  # store results of oxygen model in results folder
  tar_target(write_oxygen, save_model_output(oxygen, lake_id), format = "file"),
  # load observations (Abby's lakes)
  tar_target(observations, read_observations(lake_id, thermal)),
  # Model quality
  # tar_target(oxy_quality, oxy_qa(oxygen, observations), error = "null"),
  # tar_target(oxy_scatter, oxy_qa_full(oxygen, observations), error = "null"),
  # Create QC plots for oxygen
  tar_target(plot_qc_oxy, save_qc_plot_oxygen(oxygen, lake_id, observations), format = "file", error = "null"),
  # Create plots of temperature and thermocline depth
  tar_target(plot_thermal, create_plots_thermal(thermal, lake_id, lake_folder), format = "file")
)

combined <- list(
  tar_combine(
    runtimes,
    targets[[2]],
    command = combine_lakes(!!!.x)
  )

  # tar_combine(
  #   combined_oxy_qa,
  #   targets[[5]],
  #   command = bind_rows(!!!.x)
  # ),
  # 
  # tar_target(plot_oxy_qa, plot_oxygen_qa(combined_oxy_qa), format = "file"),
  # 
  # tar_combine(
  #   full_qa,
  #   targets[[6]],
  #   command = plot_full_qa(!!!.x),
  #   format = "file"
  # )
)

glob_runtime <- tar_target(plot_runtime, plot_runtimes(runtimes), format = "file")

list(glob_trophy, glob_methods, glob_params, targets, combined, glob_runtime)
# list(glob_trophy, glob_methods, glob_params, targets, glob_runtime)