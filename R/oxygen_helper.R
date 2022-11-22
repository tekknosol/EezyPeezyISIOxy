
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

run_oxygen_model <- function(thermal_data, method = 'rk4', trophy = 'oligo',
                             iterations = 100){

  id_na <- which(is.na(thermal_data$stratified))
  thermal_data_reduced <- thermal_data[-id_na, ]

  oxygen_df <- c()

  for (i in unique(thermal_data_reduced$strat_id)){

    data_to_model <- thermal_data_reduced[which(thermal_data_reduced$strat_id == i), ]

    if (mean(data_to_model$duration <= 2)){
      next
    }

    model_output <- consume_oxygen(thermal_subset = data_to_model, method = method, strat_id = i,
                                   trophy = trophy, iterations = iterations)

    oxygen_df <- rbind(oxygen_df, model_output)

    message(paste0('Finished running ',match(i,unique(thermal_data_reduced$strat_id)),' out of ', length(unique(thermal_data_reduced$strat_id)),'.'))
  }

  return(oxygen_df)
}

get_prior <- function(trophy){
  if (trophy == 'oligo'){

    Flux <- rnorm(1, mean = -0.32, sd = 0.096)  # (g / m2 / d)
    Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(1, mean = 1.07, sd = 0.03)

  } else if (trophy == 'eutro'){

    Flux <- rnorm(1, mean = -0.32, sd = 0.096)  # (g / m2 / d)
    Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
    Theta <- rnorm(1, mean = 1.07, sd = 0.03)

  }
  return(c(Flux = Flux, Khalf = Khalf, Theta = Theta))
}

consume_oxygen <- function(thermal_subset, method, trophy,
                           iterations){
  strat_id <- thermal_subset$strat_id
  
  Time_linear <- seq(1, nrow(thermal_subset), 1)
  Area_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_area, method = "linear", rule = 2)
  Volume_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_volume, method = "linear", rule = 2)
  Temp_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_temp, method = "linear", rule = 2)

  yini <- c(cO2 = o2.at.sat.base(temp = thermal_subset$hypo_temp[1], altitude = 500)) # g/m3

  Output = c(NULL)
  for (k in 1:iterations){
    parameters <- get_prior(trophy = trophy)

    Output_ode <- ode(times = Time_linear, y = yini, func = o2_model,
                      parms = parameters, method = 'rk4',
                      Area_linear = Area_linear, Volume_linear = Volume_linear,
                      Temp_linear = Temp_linear)

    Output <- cbind(Output, Output_ode[, 2])
  }

  Output_df = data.frame('Time' = thermal_subset$datetime, Output)

  m.Output_df = pivot_longer(Output_df, 2:last_col())
  m.Output_df <- m.Output_df %>%
    mutate(Time = lubridate::as_datetime(Time))

  df <- m.Output_df %>%
    group_by(Time) %>%
    summarise(oxygen_mean = mean(value),
              oxygen_median = median(value),
              oxygen_sd = sd(value),
              oxygen_upperPercentile = quantile(value, probs = c(0.975)),
              oxygen_lowerPercentile = quantile(value, probs = c(0.025)),
              tropic_state = 'oligotrophic') %>%
    rename(datetime = Time)

  df$strat_id <- strat_id

  return(df)
}

# Model code
o2_model <- function(Time, State, Pars, Area_linear, Temp_linear, Volume_linear) {
  with(as.list(c(State, Pars)), {
    # cO2 <- y

    SedimentFlux    <- Area_linear(Time) * Flux # m2 * g/m2/d
    MichaelisMenten   <- ((cO2) / (Khalf + cO2)) # g/m3 / g/m3
    ArrheniusCorrection <- Theta^(Temp_linear(Time) - 20) # -

    dcO2        <-  SedimentFlux * MichaelisMenten * ArrheniusCorrection / Volume_linear(Time)
    # m2 g/m2/d g/m3 / g/m3 m3 / m3 / m3 = g/m3/d

    return(list(c(dcO2)))
  })
}

save_model_output <- function(oxygen_output, lake_id){

  oxygen_output$lake_id <- lake_id

  # write_csv(file = paste0(working_folder, '/', lake_id, '/oxygen_info_meta.csv'), x = oxygen_output)
  filename <- paste0("results/oxygen_", lake_id,".rds")
  saveRDS(object = oxygen_output, file = filename)

  return(filename)
  # message('Output saved.')
}

create_plot <- function(oxygen_output, lake_id, working_folder){

  oxygen_output$year <- year(oxygen_output$datetime)
  ggplot(oxygen_output) +
    geom_point(aes(datetime, oxygen_mean, col = 'mean'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_median, col = 'median'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_mean - oxygen_sd, col = 'sd_lower'), linetype = 'dashed', lwd = 1.5) +
    geom_point(aes(datetime, oxygen_mean + oxygen_sd, col = 'sd_upper'), linetype = 'dashed', lwd = 1.5) +
    geom_point(aes(datetime, oxygen_upperPercentile, col = '97.5'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_lowerPercentile, col = '2.5'), lwd = 1.5) +
    ylab('DO conc. (g/m3)') + xlab("") +
    facet_wrap(~ year, scales = 'free') +
    theme_bw()

  filename1 <- paste0('results/plots/oxygen/', lake_id, '_oxygen_plot.jpg')
  ggsave(filename = filename1,
         width = 10, height = 10, units = 'in', bg = "white")

  ggplot(subset(oxygen_output, year >= 2010)) +
    geom_point(aes(datetime, oxygen_mean, col = 'mean'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_median, col = 'median'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_mean - oxygen_sd, col = 'sd_lower'), linetype = 'dashed', lwd = 1.5) +
    geom_point(aes(datetime, oxygen_mean + oxygen_sd, col = 'sd_upper'), linetype = 'dashed', lwd = 1.5) +
    geom_point(aes(datetime, oxygen_upperPercentile, col = '97.5'), lwd = 1.5) +
    geom_point(aes(datetime, oxygen_lowerPercentile, col = '2.5'), lwd = 1.5) +
    ylab('DO conc. (g/m3)') + xlab("") +
    theme_bw()

  filename2 <- paste0('results/plots/oxygen/', lake_id,'_oxygen_plot_recent.jpg')
  ggsave(filename =  filename2, width = 7,
         height = 3, units = 'in', bg = "white")
  
  
  c(filename1, filename2)
  
  # lake_id

  # message('Plots saved.')
}
