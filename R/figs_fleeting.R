
fig_ts <- function(path = "results/plots"){
  annual_min_o2 <- qread(here("results/annual_min_o2.qs"))
  
  annual_min_o2 <- annual_min_o2 %>% 
    mutate("0.5" = oxygen_min < 0.5) %>% 
    mutate("1" = oxygen_min < 1) %>% 
    mutate("2" = oxygen_min < 2) %>% 
    pivot_longer(any_of(c("0.5", "1", "2")), names_to = "oxy_lvl", values_to = "anoxic")
  
  # PLot Time series (Fig 4)
  plot_df <- annual_min_o2 %>% 
    filter(anoxic == TRUE) %>% 
    group_by(year, trophic_state, oxy_lvl) %>% 
    count()
  
  fig <- ggplot(plot_df, aes(year, n))+
    geom_line()+
    facet_wrap(trophic_state~oxy_lvl, scales = "free_y")+
    ylab("Number of lakes")+
    theme_bw()
  
  path <- here(path, "fig_ts.jpg")
  ggsave(path, fig, width = 9.5, height = 6)
}


fig_obs_ridge <- function(path = ""){
  df_obs <- observations %>% 
    group_by(year = year(datetime), lake_id) %>%
    summarise(across(everything(), min)) %>%
    ungroup() %>% 
    mutate(trophic_state = "xobs") %>% 
    rename(oxygen_min = DO_mgL)
  
  
  v <- annual_min_o2 %>% bind_rows(df_obs)
  
  sort <- v %>% 
    filter(trophic_state == "xobs") %>% 
    group_by(lake_id, trophic_state) %>% 
    summarise(sort = mean(oxygen_min, na.rm = T)) %>% 
    ungroup() %>% 
    select(-trophic_state)
  
  v <- v %>% left_join(sort)
  
  w <- v %>% 
    mutate(oxygen_min = ifelse(trophic_state == "xobs", oxygen_min, oxygen_min)) %>% 
    mutate(lake_id = factor(lake_id)) %>% 
    mutate(lake_id = forcats::fct_reorder(lake_id, sort, .na_rm = FALSE))
  
  fig <- ggplot(w, aes(x = oxygen_min, fill = trophic_state, y = factor(lake_id)))+ 
    geom_density_ridges(scale = 5, alpha = .8)+
    theme_ridges(center_axis_labels = TRUE)+
    scale_fill_discrete("", labels = c("Eutrophic model", "Oligotrophic model", "Observed"))+
    labs(x = "Annual mimimum DO", y = "Lake ID")+
    theme(legend.position = "bottom")
  
  ggsave(path, fig)
}