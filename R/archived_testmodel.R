library(tidyverse)
library(rLakeAnalyzer)
library(LakeMetabolizer)
library(deSolve)
library(pracma)

# Temp_raw <- read_delim('InputTest/temperate/output_temp.txt', skip = 8, delim = '\t')
# str(Temp_raw)
# head(Temp_raw)
# 
# Depth_raw <- read_delim('InputTest/temperate/output_z.txt', skip = 8, delim = '\t')
# Depth <- as.numeric(Depth_raw[1, -1]) * (-1)
# 
Area_raw <- read_delim('InputTest/temperate/hypsograph.dat', skip = 1, delim = ' ',
                       col_names = F)
# # dO2/dt = Flux * (O2 / (Khalf + O2)) * Theta^(Temp - 20) * Area
# # g/day = g/m2/day  * m2
# 
# 
Area_list = data.frame('Depth' = (Area_raw[, 1] - max(Area_raw[, 1])) * (-1), 'Area' = Area_raw[, 2])
Area_list = apply(Area_list, 2, rev)
Area_list <- as.data.frame(Area_list)
# 
# Temp <- Temp_raw
# colnames(Temp) <- c('Datetime', paste0('wtr_', Depth))
# Temp <- data.frame('datetime' = Temp$Datetime, rev(Temp[, 2:ncol(Temp)]))
# head(Temp)
# 
# Temp <- Temp %>%
#   dplyr::filter(datetime >= as.POSIXct('2000-01-01') &
#                   datetime <= as.POSIXct('2000-12-31'))
# 
# Stratification <- data.frame('Datetime' = Temp$datetime, 'Density.Diff' = water.density(Temp[, ncol(Temp)]) - water.density(Temp[, 2]))
# Stratification$Stratif.Check <- ifelse(Stratification$Density.Diff >= 0.1, 1, NA)
# ggplot(Stratification) + geom_line(aes(Datetime, Stratif.Check))
# 
# Strat.Period <- (na.contiguous(Stratification$Stratif.Check))
# 
# Temp <- Temp[attributes(Strat.Period)$tsp[1]: 
#                attributes(Strat.Period)$tsp[2],]
# 
# Thermocline.depth <- ts.center.buoyancy(wtr = Temp)
# str(Thermocline.depth)
# ggplot(subset(Thermocline.depth, datetime > as.Date('2000-01-01') &
#                 datetime < as.Date('2000-12-31')), 
#        aes(datetime, cent.n2)) + 
#   geom_path() +
#   geom_smooth() + 
#   scale_y_reverse()

lake_id = 'temperate'

Thermocline.PK <- read.csv('InputTest/temperate/thermo_information.csv')

Strat.Period <- (na.contiguous(Thermocline.PK$stratified))

Thermocline.PK <- Thermocline.PK[attributes(Strat.Period)$tsp[1]: 
                                   attributes(Strat.Period)$tsp[2],]

Area_interp <- approx(Area_list$X1, Area_list$X2, seq(0, max(Area_list$X1), 0.1))$y
Area_interp <- data.frame('Depth' = seq(0, max(Area_list$X1), 0.1), 'Area' =
                            Area_interp)
#Area_interp$Volume <- rev(cumsum(Area_interp$Depth * Area_interp$Area))
Area_interp$Volume <- rev(cumtrapz(Area_interp$Depth , Area_interp$Area))



Volume <- approx(Area_interp$Depth, Area_interp$Volume, Thermocline.PK$thermocline_depth_smooth,
                 rule =2)$y # VOLUME FROM AREA DEPENDING ON THERMOCLINE DEPTH
Area <- approx(Area_interp$Depth, Area_interp$Area, Thermocline.PK$thermocline_depth_smooth,
               rule = 2)$y # ACTIVE AREA FOR SEDIMENT FLUX

Temp <- Thermocline.PK$hypo_temp

###################################

Output <- c(NULL)

for (k in 1:10){
  # Normal distributions for parameter assumptions
  Flux <- rnorm(1, mean = -0.32, sd = 0.096)  # (g / m2 / d) 
  # 0.32 g/m2/d / 32 g/mol = 10 mmol O2/m2/d
  Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
  Theta <- rnorm(1, mean = 1.07, sd = 0.03) # (-)
  
  dt = 1
  
  Oxygen <- rep(NA, length(Volume))
  
  Oxygen[1] <- o2.at.sat.base(temp = Temp[1], altitude = 500) #* Volume[1] # g
    
  for (i in 2:length(Oxygen)){
    
    # EULER SCHEME
    Oxygen[i] <- Oxygen[i - 1] + 
                    ((Area[i - 1] * Flux * 
                    ((Oxygen[i - 1]) / (Khalf + Oxygen[i - 1])) * 
                    Theta^(Temp[i-1] - 20))) * dt  / Volume[i-1]
    
    F_Oxygen = (((Area[i - 1] * Flux *
                    ((Oxygen[i - 1]) / (Khalf + Oxygen[i - 1])) *
                    Theta^(Temp[i-1] - 20)))) / Volume[i-1] * (-1)
    
    # PATANKAR EULER SCHEME
    Oxygen[i]  <- (Oxygen[i - 1] +  dt * 1e-10) / (1 + dt   * F_Oxygen/Oxygen[i - 1])
  }
  
  
  
  Output <- cbind(Output, Oxygen)

}

Output_df = data.frame('Time' = Thermocline.PK$datetime, Output)

m.Output_df = pivot_longer(Output_df, 2:last_col())
m.Output_df <- m.Output_df %>% 
  mutate(Time = lubridate::as_datetime(Time))

ggplot(m.Output_df) +
  geom_line(aes(Time, value, group = name)) +
  theme_bw() + xlab('Date') + ylab('O2 [mg/L]') #+
  # scale_x_datetime(date_labels = "%M-")


###################################

Time_linear <- seq(1, length(Temp), 1)
Area_linear <- approxfun(x = Time_linear, y = Area, method = "linear", rule = 2)
Volumediff_linear <- approxfun(x = Time_linear, y = (Volume)/lag(Volume), method = "linear", rule = 2)
Volume_linear <- approxfun(x = Time_linear, y = Volume, method = "linear", rule = 2)
Temp_linear <- approxfun(x = Time_linear, y = Temp, method = "linear", rule = 2)

# Model code
o2_model <- function(Time, State, Pars) {
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

# Define parameters for model
parameters <- c(Flux = Flux, Khalf = Khalf, Theta = Theta)

yini <- c(cO2 = o2.at.sat.base(temp = Temp[1], altitude = 500)) # g/m3

Output = c(NULL)
for (k in 1:1000){
  # Normal distributions for parameter assumptions
  Flux <- rnorm(1, mean = -0.32, sd = 0.096)  # (g / m2 / d) 
  Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
  Theta <- rnorm(1, mean = 1.07, sd = 0.03) # (-)
  
  parameters <- c(Flux = Flux, Khalf = Khalf, Theta = Theta)
  
  Output_ode <- ode(times = Time_linear, y = yini, func = o2_model, 
                    parms = parameters, method = 'rk4')
  
  Output <- cbind(Output, Output_ode[, 2])
}

Output_df = data.frame('Time' = Thermocline.PK$datetime, Output)

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

ggplot(df) +
  geom_line(data = m.Output_df, aes(Time, value, group = name), col = 'grey', alpha = 0.5) +
  geom_line(aes(datetime, oxygen_mean, col = 'mean'), lwd = 1.5) +
  geom_line(aes(datetime, oxygen_median, col = 'median'), lwd = 1.5) +
  geom_line(aes(datetime, oxygen_mean - oxygen_sd, col = 'sd_lower'), linetype = 'dashed', lwd = 1.5) +
  geom_line(aes(datetime, oxygen_mean + oxygen_sd, col = 'sd_upper'), linetype = 'dashed', lwd = 1.5) +
  geom_line(aes(datetime, oxygen_upperPercentile, col = '97.5'), lwd = 1.5) +
  geom_line(aes(datetime, oxygen_lowerPercentile, col = '2.5'), lwd = 1.5) +
  ylab('DO conc. (g/m3)') + xlab("") +
  theme_minimal()

write_csv(file = paste0('InputTest/oxygen_info_meta_',lake_id), x = df)
