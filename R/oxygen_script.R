
oxygen_walk <- function(
    lake, 
    method = 'patankar-rk2_c', 
    trophy = 'oligo',
    iterations = 100, 
    params = NULL
){
  thermal_path <- "~/scratch/isi/thermal/"
  
  oxygen_path <- "~/scratch/isi/oxygen/"
  if(!dir.exists(oxygen_path)) dir.create(oxygen_path, recursive = TRUE)
  
  tryCatch({    
      td_subset <- qs::qread(paste0(thermal_path,"thermal_", lake, ".qs"))
      
      run_oxygen_model(td_subset, lake, method, trophy, iterations, params) %>% 
        qs::qsave(paste0(oxygen_path,"oxygen_",trophy,"_",lake,".qs"))
    
  }, error = function(err) {
    
  }, finally = {
    
  }) # END tryCatch
}

run_oxygen_model <- function(thermal_data, lake_id, method = "patankar-rk2_c", trophy = 'oligo',
                             iterations = 100, params = NULL){

  
  start <- Sys.time()
  # This filtering is done directly in the target workflow  
  # id_na <- which(is.na(thermal_data$stratified))
  # thermal_data_reduced <- thermal_data[-id_na, ]
  
  # thermal_data is a tibble where each row contains a list of stratification event.
  # This is necessary because of the target workflow
  # bind_row re-structures it into a tibble
  # thermal_data_reduced <- thermal_data$data %>% bind_rows()
  
  thermal_data_reduced <- thermal_data %>%
    filter(!is.na(stratified)) %>%
    filter(duration > 2)
  
  
  output = oxygen_model_switch(thermal_data_reduced, method = method, trophy = trophy, iterations = iterations, params = params, lakeid = lake_id)
  end <- Sys.time()
  runtime <- round(as.numeric(end) - as.numeric(start),2)
  
  output <- output %>% 
    bind_cols(
      tibble(
        runtime = runtime,
        lake_id = lake_id
      )
    ) %>% 
    ungroup()
  return(output)
  
  
}


oxygen_model_switch <- function(thermal_data_reduced, method, trophy,
                           iterations, params = NULL, lakeid = NULL){
  
  if(method == "lsoda_event"){
    output <-  consume_oxygen_event(thermal_data_reduced, trophy = trophy, method = method, iterations = iterations, params = params)
  }else if (method == "lsoda_event_f"){
    output <-  cmp_consume_oxygen_event_f(thermal_data_reduced, trophy = trophy, method = method, iterations = iterations, params = params)
  } else if (method == "patankar-rk2_c"){
    output <-  consume_oxygen_patankar(thermal_data_reduced, method = method, trophy = trophy, iterations = iterations, params = params, lakeid = lakeid)
  } else {
    output <-  consume_oxygen_per_stratevent(thermal_data_reduced, method = method, trophy = trophy, iterations = iterations, params = params)
  }
  
  return(output)
}

consume_oxygen_per_stratevent <- function(thermal_subset, method, trophy,
                                          iterations, params = NULL){
  oxygen_df <- c(NULL)
  for (i in unique(thermal_subset$strat_id)){
    
    data_to_model <- thermal_subset[which(thermal_subset$strat_id == i), ]
    
    if (mean(data_to_model$duration <= 2)){
      next
    }
    
    
    model_output <- consume_oxygen(thermal_subset = data_to_model, method = method, strat_id = i,
                                           trophy = trophy, iterations = iterations, params = params)
    
    
    oxygen_df <- rbind(oxygen_df, model_output)
    
    # message(paste0('Finished running ',match(i,unique(thermal_data_reduced$strat_id)),' out of ', length(unique(thermal_data_reduced$strat_id)),'.'))
  }
  
  return(oxygen_df)
  
}

consume_oxygen_patankar <- function(thermal_subset, method, trophy,
                                    iterations, params = NULL, lakeid = NULL){
  
  thermal_subset <- thermal_subset %>% 
    mutate(onset = (strat_id - lag(strat_id, default = 9999))) %>% 
    mutate(onset = ifelse(onset == 0,0,1)) %>% 
    mutate(O2_sat = LakeMetabolizer::o2.at.sat.base(temp = hypo_temp, altitude = 500))
  
  
  Time_linear <- seq(1, nrow(thermal_subset), 1)
  Area_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_area, method = "linear", rule = 2)
  Volume_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_volume, method = "linear", rule = 2)
  Temp_linear <- approxfun(x = Time_linear, y = thermal_subset$hypo_temp, method = "linear", rule = 2)
  
  Output = c(NULL)
  if (!is.null(params)) {params <- params %>% filter(trophic_state == trophy)}

  Output_ode <- o2_model_patankarrk2_cpp(
    times = Time_linear,
    y = thermal_subset$O2_sat,
    onset = thermal_subset$onset,
    Flux = as.numeric(params[['Flux']]),
    Khalf = as.numeric(params[['Khalf']]),
    Theta = as.numeric(params[['Theta']]),
    Area_linear = Area_linear(Time_linear),
    Volume_linear = Volume_linear(Time_linear),
    Temp_linear = Temp_linear(Time_linear),
    iterations = iterations
  )
  
  
  # qs::qsave(
  #   Output_ode %>% 
  #     as_tibble() %>% 
  #     mutate(datetime = lubridate::as_datetime(thermal_subset$datetime))
  #     , 
  #   file = here("results/validation/", paste0(lakeid,"_",trophy, ".qs")))
  
  # Output_ode <- o2_summarize_matrix(Output_ode)
  # Output_ode <- summarize_oxygen_matrix(Output_ode)
  
  

  Output_ode <- tibble(
    datetime = lubridate::as_datetime(thermal_subset$datetime)
  ) %>%
    bind_cols(Output_ode) %>%
    bind_cols(
      tibble(
        trophic_state = trophy,
        method = method,
        strat_id = thermal_subset$strat_id
      )
    ) %>%
    ungroup()



  return(Output_ode)
  
  
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
    
    
    Output_ode <- tibble(
        datetime = lubridate::as_datetime(thermal_subset$datetime)
      ) %>% 
      bind_cols(Output_ode) %>% 
      bind_cols(
        tibble(
          trophic_state = trophy, 
          method = method,
          strat_id = strat_id
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
      
    } else {
      
      warning('No numerical scheme selected! Please choose either rk4 or patankar-rk2.')

      }

  }
  }
  
  df <- tibble(
    datetime = thermal_subset$datetime,
    strat_id = thermal_data$strat_id,
    trophic_state = trophy, 
    method = method
  ) %>% 
    bind_cols(
      summarize_oxygen_matrix(oxygen)
    )
  
  return(df)
    
}

consume_oxygen_event_f <- function(thermal_data, trophy, method, iterations , params){
  
  times <- seq(1, nrow(thermal_data), 1)
  n_times <- length(times)
  
  events <- thermal_data %>% 
    mutate(time = times) %>% 
    group_by(strat_id) %>% 
    summarise(across(any_of(c("hypo_temp", "time")), first)) %>% 
    mutate(
      var = "cO2", value = o2.at.sat.base(temp = hypo_temp, altitude = 500), 
      method = "replace"
    ) %>% 
    select(-strat_id, -hypo_temp) %>% 
    relocate(var) %>%
    data.frame()
  
  # events <- tibble(
  #   var = events$var,
  #   time = events$time,
  #   value = events$value,
  #   method = events$method
  # ) %>%
  #   data.frame()
  
  # Area_linear <- approxfun(x = times, y = thermal_data$hypo_area, method = "linear", rule = 2)
  # Volume_linear <- approxfun(x = times, y = thermal_data$hypo_volume, method = "linear", rule = 2)
  # Temp_linear <- approxfun(x = times, y = thermal_data$hypo_temp, method = "linear", rule = 2)
  
  # compiled version
  # times1 <- times[c(1,n_times)]
  # area <- matrix(c(times,Area_linear(times)), ncol=2)
  # volume <- matrix(c(times,Volume_linear(times)), ncol=2)
  # temp <- matrix(c(times,Temp_linear(times)), ncol=2)
  
  # times1 <- times[c(1,n_times)]
  area <- matrix(c(times,thermal_data$hypo_area), ncol=2)
  volume <- matrix(c(times,thermal_data$hypo_volume), ncol=2)
  temp <- matrix(c(times,thermal_data$hypo_temp), ncol=2)
  
  DLLname <- "o2model2"
  recompile <- TRUE
  
  if (!is.null(params)) {params <- params %>% filter(trophic_state == trophy)}
  
  # yini <- c(cO2 = o2.at.sat.base(temp = thermal_data$hypo_temp[1], altitude = 500)) # g/m3  
  yini <- c(cO2 = events$value[1])
            
  Output = matrix(nrow = n_times, ncol = iterations)
  
  # if(recompile) system(paste0("R CMD SHLIB ","inst/ext/", DLLname,".f"))
  dyn.load(paste0("inst/ext/", DLLname, .Platform$dynlib.ext))
  for (k in 1:iterations){
    # parameters <- get_prior(trophy = trophy)
    # parameters <- params[k,]
    Pars <- c(Flux = params[k,]$Flux, Khalf = params[k,]$Khalf, Theta = params[k,]$Theta)
    
    Output_ode <- ode(
      y = yini, 
      times = times, 
      func = "simderivs",
      parms = Pars, 
      dllname = DLLname,
      initforc = "simforc", 
      forcings = list(area, volume, temp),
      initfunc = "siminit", 
      nout = 0, 
      nspec = 1,
      outnames = c("dcO2"),
      events = list(data=events)
    )
    
    Output[,k] <- Output_ode[, 2]
  
  }
  dyn.unload(paste0("inst/ext/",DLLname, .Platform$dynlib.ext))
  
  for(row in events$time){
    Output[row,] <- NA
  }
  
  df <- tibble(
    datetime = thermal_data$datetime,
    strat_id = thermal_data$strat_id,
    trophic_state = trophy, 
    method = method
  ) %>% 
    bind_cols(
      summarize_oxygen_matrix(Output)
    )
  
  return(df)
}

# cmp_consume_oxygen_event_f <- cmpfun(consume_oxygen_event_f)

consume_oxygen_event <- function(thermal_data, trophy, method, iterations , params){

  
  events <- thermal_data %>% 
    mutate(time = seq(1, nrow(thermal_data), 1)) %>% 
    group_by(strat_id) %>% 
    summarise(across(any_of(c("hypo_temp", "time")), first)) %>% 
    mutate(
      var = "cO2", value = o2.at.sat.base(temp = hypo_temp, altitude = 500), 
      method = "replace"
    ) %>% 
    select(-strat_id, -hypo_temp)
  
  events <- tibble(
    var = events$var,
    time = events$time, 
    value = events$value,
    method = events$method
  ) %>% 
    data.frame()
  
  
  Time_linear <- seq(1, nrow(thermal_data), 1)
  # Time_linear <- as.numeric(thermal_data$datetime)
  Area_linear <- approxfun(x = Time_linear, y = thermal_data$hypo_area, method = "linear", rule = 2)
  Volume_linear <- approxfun(x = Time_linear, y = thermal_data$hypo_volume, method = "linear", rule = 2)
  Temp_linear <- approxfun(x = Time_linear, y = thermal_data$hypo_temp, method = "linear", rule = 2)
  
  yini <- c(cO2 = o2.at.sat.base(temp = thermal_data$hypo_temp[1], altitude = 500)) # g/m3  
  
  Output = c(NULL)
  
  if (!is.null(params)) {params <- params %>% filter(trophic_state == trophy)}
  
  for (k in 1:iterations){
    # parameters <- get_prior(trophy = trophy)
    parameters <- params[k,]
    
    
    # out <- ode(y = yini, times = Time_linear, func = model, parms = parameters, 
    #          events = list(data = events))
    
    Output_ode <- ode(times = Time_linear, y = yini, func = o2_model_rk4,
                      parms = parameters, events = list(data = events),
                      Area_linear = Area_linear, Volume_linear = Volume_linear,
                      Temp_linear = Temp_linear)
    
    Output <- cbind(Output, Output_ode[, 2])
    
  }
  
  for(row in events$time){
    Output[row,] <- NA
  }
  
  df <- tibble(
    datetime = thermal_data$datetime,
    strat_id = thermal_data$strat_id,
    trophic_state = trophy, 
    method = method
  ) %>% 
    bind_cols(
      summarize_oxygen_matrix(Output)
    )
  
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



