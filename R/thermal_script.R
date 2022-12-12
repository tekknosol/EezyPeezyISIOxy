thermal_info <- function(lake){
  # Read Temp & Depth
  
  Temp <- lake %>% read_temp_nc()

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
  
  
  thermal_info <- tibble(
    datetime = Temp$datetime,
    bottom_temperature = Temp[,length(Temp)],
    stratified = Stratification
  ) %>% 
    left_join(Thermocline.depth %>% rename(thermocline_depth = thermo.depth), by = "datetime")
  
    # left_join(schmidt, by = "datetime") %>% 
    # mutate(schmidt.stability = ifelse(schmidt.stability > 30, 1, NA))
  
  # Define stratification based on schmidt
  # thermal_info <- thermal_info %>% 
  #   mutate(thermocline_depth = ifelse(schmidt.stability > 30 , thermocline_depth, NA))
  
  # thermal_info <- thermal_info %>% 
  #   mutate(thermocline_depth = ifelse(is.na(stratified) ,NA, thermocline_depth))
  
  thermal_info <- thermal_info %>% 
    mutate(s = ifelse(!is.na(stratified), stratified, 2)) %>% 
    mutate(strat_id = cumsum(c(TRUE, diff(s) != 0))) %>% 
    group_by(strat_id) %>% 
    mutate(thermocline_depth_smooth = plssmooth(thermocline_depth, lambda = 150)) %>% 
    mutate(duration = n()) %>% 
    ungroup() %>% 
    select(-s)
    
    # mutate(rn = row_number()) %>%
    # group_by(year = year(datetime)) %>% 
    # mutate(thermocline_depth_smooth = predict(loess(thermocline_depth ~ rn, na.action = na.exclude))) %>% 
    # ungroup() %>% 
    # select(-rn, -year)
  
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
    bind_cols(hypo_volume = Volume, hypo_area = Area)
  
  # write_csv(thermal_info, file = here(lake, "/thermal_info.csv.gz"))
  # write_rds(thermal_info, file = here(lake, "/thermal_info.rds"), compress = "gz")
  return(thermal_info)
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

read_temp_nc <- function(lake){
  ncin <- nc_open(here(lake, "output.nc"))
  Temp_raw_nc <- ncvar_get(ncin, "temp")
  Temp_raw_nc <- data.frame(t(Temp_raw_nc))
  tunits <- ncatt_get(ncin, "time", "units")
  seq_date <- seq(from = as.Date(unlist(strsplit(tunits$value, "seconds since "))[2]), 
                  length.out = dim(Temp_raw_nc)[1], by = "days")
  Temp_raw_nc <- cbind(date=seq_date, Temp_raw_nc)
  
  Depth_raw_nc <- ncvar_get(ncin, "z")
  Depth <- as.numeric(-Depth_raw_nc[,1])
  
  Temp <- Temp_raw_nc
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

plssmooth <- function(x, lambda=1000){
  m <- length(x)
  dd <- diag(m)
  D <- diff(dd)
  a <- dd+lambda*t(D)%*%D
  rs <- solve(a,x)
  return(rs)
}

reduce_thermal <- function(thermal){
  thermal <- thermal %>% 
    filter(!is.na(stratified)) %>% 
    filter(duration > 2)
  return(thermal)
}

create_plots_thermal <- function(thermal_data, lake_id, lake_folder){
  # year <- min(c(2010, year(max(thermal_data$datetime))-1))
  nyears <- 9
  years <- rev(unique(year(thermal_data$datetime)))[1:min(nyears, length(unique(year(thermal_data$datetime))))]
  
  plot_df <- thermal_data %>% 
    # filter(year(datetime) %in% year) %>%
    filter(year(datetime) %in% years) %>% 
    pivot_longer(any_of(c("thermocline_depth", "thermocline_depth_smooth")))
  
  temp <- read_temp_nc(here(lake_folder, lake_id))
  
  plot_temp <- temp %>%
    filter(year(datetime) %in% years) %>%
    pivot_longer(2:last_col()) %>%
    mutate(name = as.numeric(str_sub(name, 5, nchar(name))))
  
  plot <- ggplot(plot_temp, aes(datetime, name))+
    geom_raster(aes(fill = value))+
    # scale_fill_viridis_c("Temp", option = "A")+
    scale_fill_distiller("Temp", palette = "Spectral")+
    scale_color_manual("Thermocline", labels = c("Depth raw", "Depth smooth"), values = c("grey50", "black"))+
    geom_line(data = plot_df, aes(datetime, value, color = name))+
    labs(x = "Date", y = "Depth")+
    scale_y_reverse(expand = c(0,0))+
    scale_x_date(expand=c(0,0), date_labels = "%b %Y")+
    facet_wrap(~year(datetime), scales = "free")
  
  filename1 <- paste0('results/plots/thermo/', lake_id,'_thermocline_recent.jpg')
  ggsave(filename =  filename1, width = 15,
         height = 10, units = 'in', bg = "white")
  
  
  c(filename1)
}

hypo_oxy <- function(thermo, oxy, depths, bthA, bthD){
  
  dz <- 0.1
  
  thermo <- ifelse(is.na(thermo), 0, thermo)
  
  if (thermo >= max(depths)){
    return(NA)
  }
  
  if (length(oxy) < 2){
    return(NA)
  }
  
  layerD <- seq(thermo, max(bthD), dz)
  # layerD <- seq(thermo, max(depths), dz)
  layerO <- stats::approx(depths, oxy, layerD, rule = 2)$y
  layerA <- stats::approx(bthD, bthA, layerD)$y
  weightedO <- layerA * layerO * dz
  aveOxygen <- sum(weightedO)/(sum(layerA))/dz
  return(aveOxygen)
}

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

plot_runtimes <- function(runtimes){
  
  
  ggplot(runtimes, aes(method, runtime/60))+
    geom_boxplot(aes(fill = trophic_state), outlier.shape = NA)+
    # stat_summary(aes(color = trophic_state), geom = "pointrange", fun = mean, fun.min = min, fun.max = max, position = position_dodge(width = .2))+
    stat_summary(aes(color = trophic_state), fun.y = max, geom = "point", size = 2, position = position_dodge(width = .75)) +
    stat_summary(aes(color = trophic_state), fun.y = min, geom = "point", size = 2, position = position_dodge(width = .75)) +
    labs(y = "Runtime (minutes)", x = "Method")+
    theme_bw()
  
  filename1 <- paste0('results/plots/qc/runtime.jpg')
  ggsave(filename = filename1,
         width = 15, height = 10, units = 'in', bg = "white")
  
  filename1
}

plot_oxygen_qa <- function(oxygen_qa){
  
  criteria <- c("RMSE", "R2")
  
  plot_df <- oxygen_qa %>% 
    filter(names %in% criteria)
  
  ggplot(plot_df, aes(trophic_state, gof, fill = method))+
    geom_boxplot()+
    facet_grid(names~., scales = "free")+
    # labs(y = "Runtime (minutes)", x = "Method")+
    theme_bw()
  
  filename1 <- paste0('results/plots/qc/gof.jpg')
  ggsave(filename = filename1,
         width = 15, height = 10, units = 'in', bg = "white")
  
  filename1
}

plot_full_qa <- function(...){
  dots <- rlang::list2(...)
  plot_df <- bind_rows(dots)
  
  plot_df2 <- plot_df %>% 
    select(any_of(c("rk4", "rk4_zero", "patankar-rk2", "trophic_state")))
  
  plot1 <- GGally::ggpairs(plot_df2, columns = 1:(ncol(plot_df2)-1), diag = "blank", aes(color = trophic_state))+
    theme_bw()
  
  # plot_df <- plot_df %>% 
  #   select
  # 
  # a <- ggplot(plot_df, aes(rk4, `patankar-rk2`, color = trophic_state))+
  #   geom_point()+
  #   geom_abline()+
  #   theme_bw()
  # 
  # b <- ggplot(plot_df, aes(rk4_zero, `patankar-rk2`, color = trophic_state))+
  #   geom_point()+
  #   geom_abline()+
  #   theme_bw()
  # 
  # c <- ggplot(plot_df, aes(rk4, rk4_zero, color = trophic_state))+
  #   geom_point()+
  #   geom_abline()+
  #   theme_bw()
  # 
  # plot1 <- ggpubr::ggarrange(
  #   a,b,c, 
  #   labels = "auto", 
  #   ncol = 3, 
  #   common.legend = T,
  #   legend = "bottom"
  # )
  
  plot_df <- plot_df %>%
    pivot_longer(any_of(c("rk4", "rk4_zero", "patankar-rk2")))

  a <- ggplot(plot_df, aes(DO_mgL, value))+
    geom_point(color = "grey")+
    stat_summary(geom = "pointrange", aes(x = plyr::round_any(DO_mgL, 1)))+
    facet_wrap_equal(name~trophic_state, ncol = 2, nrow = 3)+
    geom_abline()+
    theme_bw()+
    labs(x = "Observed DO (mg L⁻¹)", y = "Modelled DO (mg L⁻¹)")
  
  
  
  filename1 <- paste0('results/plots/qc/gof_scatter.jpg')
  ggsave(plot1, filename = filename1,
         width = 8, height = 8, units = 'in', bg = "white")
  
  filename2 <- paste0('results/plots/qc/gof_scatter2.jpg')
  ggsave(a, filename = filename2,
         width = 8, height = 10, units = 'in', bg = "white")
  
  # c(filename1)
  c(filename1, filename2)
}


# ----- facets ----------
library(ggplot2)
FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)

facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}