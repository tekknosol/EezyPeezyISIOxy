


run_oxygen_model <- function(thermal_data, lake_id, method = 'rk4', trophy = 'oligo',
                             iterations = 100, params = NULL){

  
  
  # This filtering is done directly in the target workflow  
  # id_na <- which(is.na(thermal_data$stratified))
  # thermal_data_reduced <- thermal_data[-id_na, ]
  
  # thermal_data is a tibble where each row contains a list of stratification event.
  # This is necessary because of the target workflow
  # bind_row re-structures it into a tibble
  thermal_data_reduced <- thermal_data$data %>% bind_rows()

  oxygen_df <- c()

  for (i in unique(thermal_data_reduced$strat_id)){

    data_to_model <- thermal_data_reduced[which(thermal_data_reduced$strat_id == i), ]

    if (mean(data_to_model$duration <= 2)){
      next
    }

    model_output <- consume_oxygen(thermal_subset = data_to_model, method = method, strat_id = i,
                                   trophy = trophy, iterations = iterations, params = params)

    oxygen_df <- rbind(oxygen_df, model_output)

    # message(paste0('Finished running ',match(i,unique(thermal_data_reduced$strat_id)),' out of ', length(unique(thermal_data_reduced$strat_id)),'.'))
  }
  
  oxygen_df$lake_id <- lake_id

  return(oxygen_df)
}



consume_oxygen <- function(thermal_subset, method, trophy,
                           iterations, strat_id, params = NULL){
  # strat_id <- thermal_subset$strat_id
  
  Time_linear <- seq(1, nrow(thermal_subset), 1)
  Area_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_area, method = "linear", rule = 2)
  Volume_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_volume, method = "linear", rule = 2)
  Temp_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_temp, method = "linear", rule = 2)

  yini <- c(cO2 = o2.at.sat.base(temp = thermal_subset$hypo_temp[1], altitude = 500)) # g/m3

  Output = c(NULL)
  if (!is.null(params)) {params <- params %>% filter(trophic_state == trophy)}
  start <- Sys.time()
  
  if(method == "patankar-rk2_c"){
    Output_ode <- o2_model_patankarrk2_cpp(
      times = Time_linear, 
      y = as.numeric(yini['cO2']),
      Flux = as.numeric(params[['Flux']]),
      Khalf = as.numeric(params[['Khalf']]),
      Theta = as.numeric(params[['Theta']]),
      Area_linear = Area_linear(Time_linear), 
      Volume_linear = Volume_linear(Time_linear),
      Temp_linear = Temp_linear(Time_linear),
      iterations = iterations
    )
    end <- Sys.time()
    runtime <- round(as.numeric(end - start),2)
    
    Output_ode <- tibble(
        datetime = lubridate::as_datetime(thermal_subset$datetime)
      ) %>% 
      bind_cols(Output_ode) %>% 
      bind_cols(
        tibble(
          trophic_state = trophy, 
          method = method,
          strat_id = strat_id,
          runtime = runtime
        )
      ) %>% 
      ungroup()
    
    
    
    return(Output_ode)
    
  }else{
    
  
  
  for (k in 1:iterations){
    # parameters <- get_prior(trophy = trophy)
    parameters <- params[k,]
    
    if (method == 'rk4' | method == "lsoda"){
      Output_ode <- ode(times = Time_linear, y = yini, func = o2_model_rk4,
                        parms = parameters, method = method,
                        Area_linear = Area_linear, Volume_linear = Volume_linear,
                        Temp_linear = Temp_linear)
      
      Output <- cbind(Output, Output_ode[, 2])
    } else if (method == 'rk4_zero'){
      Output_ode <- ode(times = Time_linear, y = yini, func = o2_model_rk4_zero,
                        parms = parameters, method = 'rk4',
                        Area_linear = Area_linear, Volume_linear = Volume_linear,
                        Temp_linear = Temp_linear)
      
      Output <- cbind(Output, Output_ode[, 2])
    } else if (method == 'patankar-rk2'){
      
      Output_ode <- o2_model_patankarrk2(times = Time_linear, y = yini,
                        parms = parameters,
                        Area_linear = Area_linear, Volume_linear = Volume_linear,
                        Temp_linear = Temp_linear)
      
      Output <- cbind(Output, Output_ode)
      
    } else if (method == 'patankar-rk2_c'){
      
      Output_ode <- o2_model_patankarrk2_cpp(
          times = Time_linear, 
          y = as.numeric(yini['cO2']),
          Flux = as.numeric(params['Flux']),
          Khalf = as.numeric(params['Khalf']),
          Theta = as.numeric(params['Theta']),
          Area_linear = Area_linear(Time_linear), 
          Volume_linear = Volume_linear(Time_linear),
          Temp_linear = Temp_linear(Time_linear),
          iterations = iterations
        )
      
      Output <- cbind(Output, Output_ode)
      
      
    }else {
      
      warning('No numerical scheme selected! Please choose either rk4 or patankar-rk2.')

      }

  }
  }
  end <- Sys.time()
  runtime <- round(as.numeric(end - start),2)
    
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
              trophic_state = trophy, method = method) %>%
    rename(datetime = Time)

  df$strat_id <- strat_id
  df$runtime <- runtime

  return(df)
}


# Model code
dcdt <- function(SedimentFlux, MichaelisMenten, ArrheniusCorrection, Volume){
  dc <- SedimentFlux * MichaelisMenten * ArrheniusCorrection / Volume
  
  return(dc)
}

o2_model_rk4 <- function(Time, State, Pars, Area_linear, Temp_linear, Volume_linear) {
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

o2_model_rk4_zero <- function(Time, State, Pars, Area_linear, Temp_linear, Volume_linear) {
  with(as.list(c(State, Pars)), {
    # cO2 <- y
    
    SedimentFlux    <- Area_linear(Time) * Flux # m2 * g/m2/d
    MichaelisMenten   <- ((cO2) / (Khalf + cO2)) # g/m3 / g/m3
    ArrheniusCorrection <- Theta^(Temp_linear(Time) - 20) # -
    
    MichaelisMenten <- max(c(MichaelisMenten, 0))
    
    dcO2        <-  SedimentFlux * MichaelisMenten * ArrheniusCorrection / Volume_linear(Time)
    
    # m2 g/m2/d g/m3 / g/m3 m3 / m3 / m3 = g/m3/d
    
    return(list(c(dcO2)))
  })
}

o2_model_patankarrk2 <- function(times, y, parms, Area_linear, Temp_linear, Volume_linear,
                                 dt = 1) {
    # cO2 <- y
  
    # parms <- get_prior(trophy = 'eutro')
    
    Flux = parms['Flux']
    Khalf = parms['Khalf']
    Theta = parms['Theta']
    
    dcO2 <- rep(NA, length(times))
    dcO2[1] <- y
    
    len_y0 = length(c(y))
    eye = diag(len_y0)
    eye = ifelse(eye == 1, TRUE, FALSE)
    eyetilde = ifelse(eye == 1, FALSE, TRUE)
    avec = eye * 0
    r = avec[,1] * 0
    
    for (i in times[2: length(times)]){
      
      SedimentFlux    <- Area_linear(i) * Flux # m2 * g/m2/d
      MichaelisMenten   <- ((dcO2[i - 1]) / (Khalf + dcO2[i - 1])) # g/m3 / g/m3
      ArrheniusCorrection <- Theta^(Temp_linear(i) - 20) # -
      
      p0 = 1e-10
      d0 = abs(SedimentFlux * MichaelisMenten * ArrheniusCorrection / Volume_linear(i))
      
      ydat <- dcO2[i - 1]
      
      avec[eye] <- dt * d0 / ydat + 1
      
      c_rep = ydat
      
      avec[eyetilde] = -dt * p0 / c_rep[eyetilde]
      
      r =  ydat + dt * p0[eye]
      
      cproxy = solve(as.numeric(avec), r)
      
      
      SedimentFlux    <- Area_linear(i) * Flux # m2 * g/m2/d
      MichaelisMenten   <- ((cproxy) / (Khalf + cproxy)) # g/m3 / g/m3
      ArrheniusCorrection <- Theta^(Temp_linear(i) - 20) # -
      
      p1 = 1e-10
      d1 =  abs(SedimentFlux * MichaelisMenten * ArrheniusCorrection / Volume_linear(i))
      
      p = 0.5 * (p0 + p1)
      d = 0.5 * (d0 + d1)
      
      avec[eye] = dt * d / cproxy + 1
      
      c_rep = cproxy
      
      avec[eyetilde] = -dt * p / c_rep[eyetilde]
      
      r = ydat + dt * p[eye]
      
      dcO2[match(i, times)] = solve(as.numeric(avec), r)
    }
    
    # m2 g/m2/d g/m3 / g/m3 m3 / m3 / m3 = g/m3/d
    
    return(dcO2)
}

save_model_output <- function(oxygen_output, lake_id){

  oxygen_output$lake_id <- lake_id

  # write_csv(file = paste0(working_folder, '/', lake_id, '/oxygen_info_meta.csv'), x = oxygen_output)
  filename <- paste0("results/oxygen/oxygen_", lake_id,".rds")
  saveRDS(object = oxygen_output, file = filename)

  return(filename)
  # message('Output saved.')
}



