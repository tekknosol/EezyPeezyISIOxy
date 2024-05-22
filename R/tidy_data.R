
tidy_annual_min_o2 <- function(path = "~/scratch/isi/oxygen"){
  files <- list.files(here(path))
  
  annual_min_o2 <- files %>% 
    map_dfr(~{
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
  
  qs::qsave(annual_min_o2, "results/annual_min_o2.qs")
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