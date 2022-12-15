combine_lakes <- function(...){
  # current <- tibble(...)[[1]]
  dots <- rlang::list2(...)
  dots %>% 
    map_dfr(~{
      .x %>% group_by(lake_id, trophic_state, method, strat_id) %>%
        summarise(runtime = first(runtime)) %>% 
        summarise(runtime = sum(runtime))
      
    })
  
}

oxy_qa <- function(oxygen, observations){
  df <- oxygen %>% 
    mutate(datetime = as_date(datetime)) %>% 
    # select(lake_id, datetime, oxygen_mean, trophic_state) %>% 
    left_join(
      observations
    ) 
  
  tryCatch({
    df %>% 
      group_by(lake_id, method, trophic_state) %>% 
      summarise(RMSE = hydroGOF::rmse(oxygen_mean, DO_mgL), R2 = hydroGOF::rPearson(oxygen_mean, DO_mgL)^2) %>% 
      pivot_longer(any_of(c("RMSE", "R2")), values_to = "gof", names_to = "names")
    # group_modify(~{
    #   gof <- gof(.x$oxygen_mean, .x$DO_mgL)
    #   tibble(names = rownames(gof), gof = gof[,1])
    # })
  }, error = function(e){return(NULL)})
}

oxy_qa_full <- function(oxygen, observations){
  df <- oxygen %>% 
    mutate(datetime = as_date(datetime)) %>% 
    select(lake_id, datetime, trophic_state, method, oxygen_mean) 
  
  df %>% 
    pivot_wider(names_from = method, values_from = oxygen_mean) %>% 
    left_join(observations)
}
