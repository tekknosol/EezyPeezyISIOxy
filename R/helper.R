
plssmooth <- function(x, lambda=1000){
  if (all(is.na(x))) return(NA)
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

read_observations <- function(lake_id, thermal){
  
  # bathy <- read_hypso(here("ObsDOTest", lake_id))
  hypso_pattern <- paste0("h_", lake_id, ".dat")
  bathy <- read_hypso(here("data/isimip/hypsography", hypso_pattern))
  
  # lake_id <- as.numeric(str_sub(lake_id, 3, nchar(lake_id)))
  
  # obs_lookup <- read_csv(here("data/ID_ODtest.csv"))
  
  hylakid <- id_lookup %>% filter(isimip_id == lake_id) %>% pull(hydrolakes)
  
  obs <- read_rds("data/observed.rds") 
  obs <- obs %>% 
    filter(hylak_id == hylakid) %>%
    filter(Depth_m >= 0) %>% 
    filter(DO_mgL >= 0) %>% 
    filter(DO_mgL < 30)
  
  obs <- obs %>% 
    select(datetime = Date, Depth_m, DO_mgL) %>% 
    left_join(
      thermal %>% 
        select(datetime, thermocline_depth_smooth)
    ) 
  
  obs %>% 
    # filter(datetime == ymd("2010-01-03")) %>% 
    na.exclude() %>% 
    group_by(datetime) %>% 
    group_modify(~{
      
      .x %>% 
        na.exclude() %>% 
        summarise(DO_mgL = hypo_oxy(first(thermocline_depth_smooth), DO_mgL, Depth_m, bathy$areas, bathy$depths))
      
    }) %>% 
    ungroup() %>% 
    mutate(lake_id = lake_id)
}

get_hypsography <- function(lake_id, working_folder){
  
  Area_raw <- read_delim(paste0(working_folder, '/', lake_id, '/hypsograph.dat'), skip = 1, delim = ' ',
                         col_names = F)
  
  Area_list <- data.frame('Depth' = (Area_raw[, 1] - max(Area_raw[, 1])) * (-1), 'Area' = Area_raw[, 2])
  Area_list <- apply(Area_list, 2, rev)
  
  Hypsography <- data.frame(Area_list)
  
  colnames(Hypsography) <- c('depth', 'area')
  
  message('Finished reading hypsography.')
  
  return(Hypsography)
}

get_thermal_data <- function(lake_id, working_folder, hypsography_data){
  
  # get thermal information
  thermo_information <- readRDS(paste0(working_folder, '/', lake_id, '/thermal_info.rds'))
  
  # Strat.Period <- (na.contiguous(thermo_information$stratified))
  #
  # thermo_information <- thermo_information[attributes(Strat.Period)$tsp[1]:
  #                                    attributes(Strat.Period)$tsp[2],]
  
  Area_interp <- approx(hypsography_data$depth, hypsography_data$area, seq(0, max(hypsography_data$depth), 0.1))$y
  Area_interp <- data.frame('Depth' = seq(0, max(hypsography_data$depth), 0.1), 'Area' =
                              Area_interp)
  
  Area_interp$Volume <- rev(cumtrapz(Area_interp$Depth , Area_interp$Area))
  
  
  # VOLUME FROM AREA DEPENDING ON THERMOCLINE DEPTH
  Volume <- approx(Area_interp$Depth, Area_interp$Volume, thermo_information$thermocline_depth_smooth,
                   rule =2)$y
  # ACTIVE AREA FOR SEDIMENT FLUX
  Area <- approx(Area_interp$Depth, Area_interp$Area, thermo_information$thermocline_depth_smooth,
                 rule = 2)$y
  
  thermo_information$hypo_area <- Area
  thermo_information$hypo_volume <- Volume
  
  message('Finished reading thermal information.')
  return(thermo_information)
}

summarize_oxygen_matrix <- function(oxygen){
  df <- tibble(
    oxygen_mean = rowMeans(oxygen),
    oxygen_median = rowMedians(oxygen),
    oxygen_sd = rowSds(oxygen),
    oxygen_upperPercentile = rowQuantiles(oxygen, probs = c(0.975)),
    oxygen_lowerPercentile = rowQuantiles(oxygen, probs = c(0.025))
  )
  return(df)
}

create_oxygen_output <- function(datetime, strat_id, oxygen){
  df <- tibble(
    datetime = datetime,
    strat_id = thermal_data$strat_id
  ) %>% 
    bind_cols(
      summarize_oxygen_matrix(oxygen)
    )
  
  return(df)
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