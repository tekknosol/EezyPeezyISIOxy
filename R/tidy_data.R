
tidy_annual_min_o2 <- function(path = "~/scratch/isi/oxygen", type = "obsclim", gcm = "20CRv3-ERA5"){
  path <- here(path, type, gcm)
  files <- list.files(path)
  
  annual_min_o2 <- files %>% 
    future_map_dfr(~{
      file <- .x
      
      tryCatch({
        oxy <- qread(here(path,file)) %>% 
          select(lake_id, datetime, oxygen_mean, trophic_state) %>% 
          group_by(trophic_state, year = year(datetime)) %>% 
          summarise(
            oxygen_min = min(oxygen_mean, na.rm = T), 
            lake_id = first(lake_id),
            .groups = "drop"
          )
        
      }, error = function(e){
        NULL
      })
    })
  
  path_out <- here("results", paste0(paste(type, tolower(gcm), "annual_min_o2", sep = "_"), ".qs"))
  qs::qsave(annual_min_o2, path_out)
}

tidy_observations <- function(){
  observations <- lakes %>% 
    split(.$lake_id) %>% 
    map_dfr(~{
      lake <- .x
      read_observations(lake$lake_id)
    })
  
  qs::qsave(observations, "data/observations.qs")
}