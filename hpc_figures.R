args = commandArgs(trailingOnly=TRUE)

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(qs)
library(ggplot2)
library(ggridges)
library(sf)
library(rnaturalearth)
library(stringr)
library(ggbeeswarm)

source(here("R/figs_fleeting.R"))

# fig_ts(type = args[1], gcm = args[2])

fig_boxplots()
fig_boxplots(type = "counterclim")