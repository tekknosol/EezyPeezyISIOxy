thermal_info <- function(lake){
  # Read Temp & Depth
  
  Temp_raw <- read_tsv(here(lake, 'output_temp.txt'), skip = 8, show_col_types = F)
  
  Depth_raw <- read_tsv(here(lake, 'output_z.txt'), skip = 8, show_col_types = F)
  Depth <- as.numeric(Depth_raw[1, -1]) * (-1)
  
  Temp <- Temp_raw
  colnames(Temp) <- c('Datetime', paste0('wtr_', Depth))
  Temp <- data.frame('datetime' = Temp$Datetime, rev(Temp[, 2:ncol(Temp)]))
  
  # Read Hyspograph
  hypso <- read_delim(here(lake, "hypsograph.dat"), skip = 1, delim = " ", show_col_types = F)
  colnames(hypso) <- c("depths", "areas")
  hypso$depths <- (hypso$depths - hypso$depths[length(hypso$depths)]) * -1
  
  # Calculate Density difference and Stratifcation
  Density.Diff <- water.density(Temp[, ncol(Temp)]) - water.density(Temp[, 2])
  Stratification <- ifelse(Density.Diff >= 0.1, 1, NA)
  
  Thermocline.depth <- ts.thermo.depth(wtr = Temp)
  
  thermocline_smooth <- Thermocline.depth %>% 
    tibble() %>% 
    na.omit() %>%
    mutate(rn = row_number()) %>% 
    group_by(year = year(datetime)) %>% 
    # filter(year == 2015) %>% 
    mutate(smooth = predict(loess(thermo.depth ~ rn))) %>% 
    ungroup() %>% 
    select(datetime, thermocline_depth_smooth = smooth, thermocline_depth = thermo.depth)
  
  
  thermal_info <- tibble(
    datetime = Temp$datetime,
    bottom_temperature = Temp[,length(Temp)],
    stratified = Stratification
  ) %>% 
    left_join(thermocline_smooth)
  
  hypo_temp <- ts.layer.temperature(wtr = Temp, top = thermal_info$thermocline_depth, bottom = max(Depth), bathy = hypso)
  
  thermal_info <- thermal_info %>% 
    left_join(hypo_temp %>% rename(hypo_temp = layer.temp))
  
  write_csv(thermal_info, file = here(lake, "/thermal_info.csv"))
}