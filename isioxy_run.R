args = commandArgs(trailingOnly=TRUE)

library(qs)
library(tidyr)
library(tibble)
library(here)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(ncdf4)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(pracma)
library(future)
library(furrr)
library(oxypatankr)

# plan(multisession, workers = 5)

message("Start setup")

cl <- makeClusterPSOCK(
  availableWorkers(),
  rscript_args = c("-e", shQuote('setwd("/home/kelp2501/EezyPeezyISIOxy")')),
  rscript_libs = .libPaths(),
  # revtunnel = FALSE,
  homogeneous = TRUE,
  setup_strategy = "parallel"
)

plan(cluster, workers = cl)


walk(list.files(here("R/"), full.names = T), source)

path_store <-  "~/scratch/isi" # Folder for target's internal data storage

lake_folder <- "data/isimip/20CRv3-ERA5" # Folder containing isimip results
numit <- 1000 # number of iterations for oxygen model
# stratification_batches <- 1 # Number of batches of stratification events per lake
lake_batches <- 9

fullrun <- TRUE
run_thermal <- FALSE
run_oxygen <- TRUE

isimip_lakes <- list.files(here(lake_folder), full.names = F)
isimip_lakes <- isimip_lakes %>% str_split_i("_", 5)

lakes <- tibble(
  # lake_id = list.files(here(lake_folder), full.names = F)[1:1]
  # lake_id = rafalakes$isimip_id
  # lake_id = od_ids$isimip_id  #lakes with observations
  # lake_id = 18005
  lake_id = isimip_lakes
) %>% 
  arrange(lake_id)

lakes <- lakes %>% 
  mutate(grp = (row_number()-1) %/% (n()/lake_batches))

lakes <- lakes %>% 
  filter(grp == as.numeric(args[1]))

trophy <- c("oligo", "eutro")

oxy_params <- trophy %>% map_dfr(get_prior, n = numit)

message("Start model")

lakes %>% 
  split(.$lake_id) %>% 
  future_walk(~{
    lake <- .x$lake_id
    if(run_thermal & (fullrun | !file.exists(paste0(path_store, "/thermal/thermal_", lake, ".qs")))){
      lake %>% thermal_walk(gcm = args[3], type = args[2]) 
    } else{message("skip thermal lake ",lake)}
    
    trophy %>%
      walk(~{
        trophy <- .x
        if(run_oxygen & (fullrun | !file.exists(paste0(path_store,"/oxygen/oxygen_",trophy, "_", lake,".qs")))){
          oxygen_walk(
            lake,
            gcm = args[3],
            type = args[2],
            trophy = trophy,
            method = "patankar-rk2_c",
            iterations = numit,
            params = oxy_params
          )
        } else{message("skip oxygen ", trophy, " lake ",lake)}
      })
      
  }, future.seed = TRUE)
