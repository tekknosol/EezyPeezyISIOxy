library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)

Temp_raw <- read_delim('InputTest/temperate/output_temp.txt', skip = 8, delim = '\t')
str(Temp_raw)
head(Temp_raw)

Depth_raw <- read_delim('InputTest/temperate/output_z.txt', skip = 8, delim = '\t')
Depth <- as.numeric(Depth_raw[1, -1]) * (-1)

Area_raw <- read_delim('InputTest/temperate/hypsograph.dat', skip = 1, delim = ' ',
                       col_names = F)
# dO2/dt = Flux * (O2 / (Khalf + O2)) * Theta^(Temp - 20) * Area
# g/day = g/m2/day  * m2

Output <- c(NULL)

for (k in 1:10){
  # Normal distributions for parameter assumptions
  Flux <- rnorm(1, mean = -0.32, sd = 0.096)  # (g / m2 / d)
  Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
  Theta <- rnorm(1, mean = 1.07, sd = 0.03) # (-)
  
  dt = 1
  
  
  Area_list = data.frame('Depth' = (Area_raw[, 1] - max(Area_raw[, 1])) * (-1), 'Area' = Area_raw[, 2])
  Area_list = apply(Area_list, 2, rev)
  Area_list <- as.data.frame(Area_list)
  
  Temp <- Temp_raw
  colnames(Temp) <- c('Datetime', paste0('wtr_', Depth))
  Temp <- data.frame('datetime' = Temp$Datetime, rev(Temp[, 2:ncol(Temp)]))
  head(Temp)
  
  Temp <- Temp %>%
    dplyr::filter(datetime >= as.POSIXct('2000-01-01') &
                    datetime <= as.POSIXct('2000-12-31'))
  
  Stratification <- data.frame('Datetime' = Temp$datetime, 'Density.Diff' = water.density(Temp[, ncol(Temp)]) - water.density(Temp[, 2]))
  Stratification$Stratif.Check <- ifelse(Stratification$Density.Diff >= 0.1, 1, NA)
  ggplot(Stratification) + geom_line(aes(Datetime, Stratif.Check))
  
  Strat.Period <- (na.contiguous(Stratification$Stratif.Check))
  
  Temp <- Temp[attributes(Strat.Period)$tsp[1]: 
               attributes(Strat.Period)$tsp[2],]
  
  Thermocline.depth <- ts.center.buoyancy(wtr = Temp)
  str(Thermocline.depth)
  ggplot(subset(Thermocline.depth, datetime > as.Date('2000-01-01') &
                  datetime < as.Date('2000-12-31')), 
         aes(datetime, cent.n2)) + 
    geom_path() +
    geom_smooth() + 
    scale_y_reverse()
  
  Thermocline.PK <- read.csv('InputTest/temperate/thermo_information.csv')
  
  Strat.Period <- (na.contiguous(Thermocline.PK$stratified))
  
  Thermocline.PK <- Thermocline.PK[attributes(Strat.Period)$tsp[1]: 
                           attributes(Strat.Period)$tsp[2],]
  
  Area_interp <- approx(Area_list$X1, Area_list$X2, seq(0, max(Area_list$X1), 0.1))$y
  Area_interp <- data.frame('Depth' = seq(0, max(Area_list$X1), 0.1), 'Area' =
                            Area_interp)
  Area_interp$Volume <- rev(cumsum(Area_interp$Depth * Area_interp$Area))
  
  Volume <- approx(Area_interp$Depth, Area_interp$Volume, Thermocline.PK$thermocline_depth_smooth,
                   rule =2)$y # VOLUME FROM AREA DEPENDING ON THERMOCLINE DEPTH
  Area <- approx(Area_interp$Depth, Area_interp$Area, Thermocline.PK$thermocline_depth_smooth,
                 rule = 2)$y # ACTIVE AREA FOR SEDIMENT FLUX
  
  Temp <- Thermocline.PK$hypo_temp
  
  Oxygen <- rep(NA, length(Volume))
  
  Oxygen[1] <- o2.at.sat.base(temp = Temp[1], altitude = 500) * Volume[1] # g
    
  for (i in 2:length(Oxygen)){
    
    Oxygen[i] <- (Oxygen[i - 1] + 
                    Area[i - 1] * (Flux * 
                    ((Oxygen[i - 1]/ Volume[i -1 ]) / (Khalf + Oxygen[i - 1] / Volume[i - 1] )) * 
                    Theta^(Temp[i-1] - 20)) * dt ) * Volume[i] / Volume[i-1]
    
  }
  
  
  
  Output <- cbind(Output, Oxygen / Volume)

}

Output_df = data.frame('Time' = Thermocline.PK$datetime, Output)

m.Output_df = pivot_longer(Output_df, 2:last_col())
m.Output_df <- m.Output_df %>% 
  mutate(Time = lubridate::as_datetime(Time))

ggplot(m.Output_df) +
  geom_line(aes(Time, value, group = name)) +
  theme_bw() + xlab('Date') + ylab('O2 [mg/L]') #+
  # scale_x_datetime(date_labels = "%M-")



library(deSolve)

Time_linear <- seq(1, length(Temp), 1)
Area_linear <- approxfun(x = Time_linear, y = Area, method = "linear", rule = 2)
Volume_linear <- approxfun(x = Time_linear, y = Volume, method = "linear", rule = 2)
Temp_linear <- approxfun(x = Time_linear, y = Temp, method = "linear", rule = 2)

# Model code
o2_model <- function(t, y, parms){
  dO2 <- (Area_linear(t) * (Flux *
        ((y/ Volume_linear(t)) / (Khalf + y/ Volume_linear(t))) *
        Theta^(Temp_linear(t) - 20)) ) * Volume_linear(t) / Volume_linear(t-1)

  return(list(dO2))
}

# Define parameters for model
parameters <- c(Flux, Khalf, Theta)

yini <- o2.at.sat.base(temp = Temp[1], altitude = 500) * Volume[1] # g

# Runge-Kutta 4th-order model solver
Output_ode <- ode(times = Time_linear, y = yini, func = o2_model, 
                  parms = parameters, method = 'rk4')

plot(Output_ode[, 2]/Volume)

