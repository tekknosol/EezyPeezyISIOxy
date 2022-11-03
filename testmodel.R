library(tidyverse)

# dO2/dt = Flux * (O2 / (Khalf + O2)) * Theta^(Temp - 20) * Area
# g/day = g/m2/day  * m2


# Normal distributions for parameter assumptions
Flux <- rnorm(1, mean = -0.32, sd = -0.096)  # (g / m2 / d)
Khalf <- rnorm(1, mean = 0.224, sd = 0.032)   # (g / m3)
Theta <- rnorm(1, mean = 1.07, sd = 0.03) # (-)

dt = 1

library(rLakeAnalyzer)
library(LakeMetabolizer)

Temp_raw <- read_delim('InputTest/temperate/output_temp.txt', skip = 8, delim = '\t')
str(Temp_raw)
head(Temp_raw)

Depth_raw <- read_delim('InputTest/temperate/output_z.txt', skip = 8, delim = '\t')
Depth <- as.numeric(Depth_raw[1, -1]) * (-1)

Area_raw <- read_delim('InputTest/temperate/hypsograph.dat', skip = 1, delim = ' ',
                       col_names = F)

Area = data.frame('Depth' = (Area_raw[, 1] - max(Area_raw[, 1])) * (-1), 'Area' = Area_raw[, 2])
Area = apply(Area, 2, rev)

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


Volume <- # VOLUME FROM AREA DEPENDING ON THERMOCLINE DEPTH
Area <- # ACTIVE AREA FOR SEDIMENT FLUX

Oxygen <- rep(NA, length(SOMETHING))

Oxygen[1] <- o2.at.sat.base(temp = Temp[1], altitude = 500) * Volume[1] # g
  
for (i in 2:SOMETHING){
  
  Oxygen[i] <- (Oxygen[i - 1] + Area[i - 1] * (Flux * 
                  (output[i - 1] / (Khalf + output[i - 1] )) * 
                  Theta^(Temp - 20)) * dt ) * Volume[n] / Volume[n-1]
  
}







# Model code
o2_model <- function(t, y, parms){
  dO2 <- Flux * (O2 / (Khalf + O2)) * Theta^(Temp - 20)
}

# Define parameters for model
parameters <- c(Flux, Khalf, Theta)

# Runge-Kutta 4th-order model solver
ode(times = times, y = yini, func = o2_model, parms = parameters, method = 'rk4')

