thermal_info <- function(lake){
  # Read Temp & Depth
  
  Temp <- lake %>% read_temp()

  hypso <- lake %>% read_hypso()
  
  # Calculate Density difference and Stratifcation
  Density.Diff <- water.density(Temp[, ncol(Temp)]) - water.density(Temp[, 2])
  Stratification <- ifelse(Density.Diff >= 0.1, 1, NA)
  
  # Calculate temperature difference and Stratifcation
  # Temp.Diff <- abs(Temp[, ncol(Temp)] - Temp[, 2])
  # Stratification <- ifelse(Temp.Diff < 1, NA, 1)
  
  Thermocline.depth <- ts_thermo_depth(wtr = Temp, seasonal = F, mixed.cutoff = 0)
  
  # Calculate Schmidt stability
  # schmidt <- ts.schmidt.stability(wtr = Temp, bathy = hypso)
  
  # Calculate length of stratified periods
  strat_duration <- rle(Stratification)
  strat_duration <- rep(strat_duration$values*strat_duration$lengths,strat_duration$lengths)
  
  # calculate length of non-stratified periods
  strat_duration_inv <- rle(is.na(Stratification))
  strat_duration_inv <- rep(strat_duration_inv$values*strat_duration_inv$lengths,strat_duration_inv$lengths)
  
  # combine length vectors
  strat_duration <- ifelse(is.na(strat_duration), strat_duration_inv, strat_duration)
  
  thermal_info <- tibble(
    datetime = Temp$datetime,
    bottom_temperature = Temp[,length(Temp)],
    stratified = Stratification
  ) %>% 
    left_join(Thermocline.depth %>% rename(thermocline_depth = thermo.depth), by = "datetime") %>% 
    bind_cols(duration = strat_duration)
    # left_join(schmidt, by = "datetime") %>% 
    # mutate(schmidt.stability = ifelse(schmidt.stability > 30, 1, NA))
  
  # Define stratification based on schmidt
  # thermal_info <- thermal_info %>% 
  #   mutate(thermocline_depth = ifelse(schmidt.stability > 30 , thermocline_depth, NA))
  
  # thermal_info <- thermal_info %>% 
  #   mutate(thermocline_depth = ifelse(is.na(stratified) ,NA, thermocline_depth))
  
  thermal_info <- thermal_info %>% 
    mutate(rn = row_number()) %>%
    group_by(year = year(datetime)) %>% 
    mutate(thermocline_depth_smooth = predict(loess(thermocline_depth ~ rn, na.action = na.exclude))) %>% 
    ungroup() %>% 
    select(-rn, -year)
  
  hypo_temp <- ts.layer.temperature(wtr = Temp, top = thermal_info$thermocline_depth, bottom = max(hypso$depths), bathy = hypso)
  
  thermal_info <- thermal_info %>% 
    left_join(hypo_temp %>% rename(hypo_temp = layer.temp), by = "datetime")
  
  # Calculate Area and Volume from hypsograph
  Area_interp <- approx(hypso$depths, hypso$areas, seq(0, max(hypso$depths), 0.1))$y
  Area_interp <- data.frame('Depth' = seq(0, max(hypso$depths), 0.1), 'Area' =
                              Area_interp)
  #Area_interp$Volume <- rev(cumsum(Area_interp$Depth * Area_interp$Area))
  Area_interp$Volume <- rev(cumtrapz(Area_interp$Depth , Area_interp$Area))
  
  Volume <- approx(Area_interp$Depth, Area_interp$Volume, thermal_info$thermocline_depth_smooth,
                   rule =2)$y # VOLUME FROM AREA DEPENDING ON THERMOCLINE DEPTH
  Area <- approx(Area_interp$Depth, Area_interp$Area, thermal_info$thermocline_depth_smooth,
                 rule = 2)$y # ACTIVE AREA FOR SEDIMENT FLUX
  
  thermal_info <- thermal_info %>% 
    bind_cols(volume = Volume, area = Area)
  
  # write_csv(thermal_info, file = here(lake, "/thermal_info.csv.gz"))
  write_rds(thermal_info, file = here(lake, "/thermal_info.rds"), compress = "gz")
}

read_temp <- function(lake){
  Temp_raw <- read_tsv(here(lake, 'output_temp.txt'), skip = 8, show_col_types = F)
  
  Depth_raw <- read_tsv(here(lake, 'output_z.txt'), skip = 8, show_col_types = F, n_max = 5)
  Depth <- as.numeric(Depth_raw[1, -1]) * (-1)
  
  Temp <- Temp_raw
  colnames(Temp) <- c('Datetime', paste0('wtr_', Depth))
  Temp <- data.frame('datetime' = Temp$Datetime, rev(Temp[, 2:ncol(Temp)]))
  
  Temp
}

read_hypso <- function(lake){
  # Read Hyspograph
  hypso <- read_delim(here(lake, "hypsograph.dat"), skip = 1, delim = " ", show_col_types = F, col_names = F)
  colnames(hypso) <- c("depths", "areas")
  hypso$depths <- (hypso$depths - hypso$depths[length(hypso$depths)]) * -1
  hypso <- hypso %>% arrange(depths)
  
  hypso
}