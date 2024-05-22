
fig_ts <- function(path = "results/plots", type = "obsclim", gcm = "20CRv3-ERA5"){
  path_data <- here("results", paste0(paste(type, tolower(gcm), "annual_min_o2", sep = "_"), ".qs"))
  annual_min_o2 <- qread(path_data)
  
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
  
  path <- here(path, paste0(paste("fig",type, tolower(gcm), "ts", sep = "_"), ".jpg"))
  ggsave(path, fig, width = 9.5, height = 6)
}

fig_map <- function(path = "results/plots", type = "obsclim", gcm = "20CRv3-ERA5"){
  path_data <- here("results", paste0(paste(type, tolower(gcm), "annual_min_o2", sep = "_"), ".qs"))
  annual_min_o2 <- qread(path_data)
  
  
  annual_min_o2 <- annual_min_o2 %>% 
    mutate("0.5" = oxygen_min < 0.5) %>% 
    mutate("1" = oxygen_min < 1) %>% 
    mutate("2" = oxygen_min < 2) %>% 
    pivot_longer(any_of(c("0.5", "1", "2")), names_to = "oxy_lvl", values_to = "anoxic")
  
  id_lookup <- read_csv(here("data/coord_area_depth.csv"))
  
  id_lookup <- id_lookup %>% 
    mutate(isimip_id = str_split(id, "_", simplify = T)[,2])
  
  hlakes <- read_sf("data/HydroLAKES/HydroLAKES_points_v10.shp")
  hlakes <- hlakes %>% 
    left_join(
      id_lookup %>% 
        select(hydrolakes, isimip_id),
      by = c("Hylak_id" = "hydrolakes")
    )
  
  plot_map <- annual_min_o2 %>% 
    group_by(trophic_state, oxy_lvl, lake_id) %>% 
    summarise(anoxic = max(anoxic)) %>% 
    mutate(anoxic = as.logical(anoxic)) %>% 
    left_join(
      hlakes,
      by = c("lake_id" = "isimip_id")
    ) %>% 
    st_as_sf()
  
  map_bg <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf") %>% 
    filter(continent != "Antarctica")
  
  # map_bg <- map_bg %>%  group_by(continent) %>% summarise()
  
  fig <- ggplot(plot_map %>% filter(Depth_avg > 6.5) %>% filter(!is.na(anoxic)))+
    geom_sf(data = map_bg)+
    geom_sf(aes(color = anoxic))+
    coord_sf(crs = "ESRI:53030")+
    facet_wrap(oxy_lvl~trophic_state, ncol = 2)+
    scale_color_discrete(direction = -1)+
    theme_bw()+
    theme(legend.position = "bottom")
  
  path <- here(path, paste0(paste("fig",type, tolower(gcm), "map", sep = "_"), ".jpg"))
  ggsave(path, fig, width = 10, height = 8.5)
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