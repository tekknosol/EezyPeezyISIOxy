
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


save_qc_plot_oxygen <- function(oxygen_data, lake_id, observed){
  
  nyears <- 6
  
  years <- rev(sort(unique(year(observed$datetime[which(year(observed$datetime)<2020)]))))
  years <-  years[1:min(nyears, length(years))]
  
  plot_df <- oxygen_data %>% 
    filter(year(datetime) %in% years) %>% 
    mutate(datetime = as_date(datetime))
  
  plot_df2 <- observed %>% 
    filter(year(datetime) %in% years) %>% 
    filter(!is.na(DO_mgL)) 
  # group_by(datetime) %>% 
  # filter(Depth_m == max(Depth_m)) %>% 
  # ungroup() %>% 
  # select(datetime, DO_mgL)
  
  ggplot()+
    geom_ribbon(data = plot_df, aes(datetime, ymin = oxygen_lowerPercentile, ymax = oxygen_upperPercentile, fill = "Percentile", group = interaction(strat_id, trophic_state)))+
    geom_line(data = plot_df2, aes(datetime, DO_mgL, color = "XObserved"))+
    geom_point(data = plot_df2, aes(datetime, DO_mgL, color = "XObserved"))+
    geom_line(data = plot_df, aes(datetime, oxygen_mean, group = interaction(strat_id, trophic_state), color = trophic_state))+
    scale_fill_manual(name = NULL, values = c("grey"))+
    scale_color_brewer(name = NULL, palette = "Set1", labels = c("Model Eutrophic", "Model Oligotrophic", "Observed"))+
    facet_wrap(~method, scales = "free_y", ncol =  1)+
    scale_x_date(date_labels = "%m-%Y", date_breaks = "1 year")+
    labs(x = "Date", y = "Oxygen (mg L⁻¹)")+
    theme_bw()+
    theme(legend.position = "bottom")
  
  
  filename1 <- paste0('results/plots/qc/', 'oxygen_', lake_id,'.jpg')
  ggsave(filename = filename1,
         width = 15, height = 10, units = 'in', bg = "white")
  
  filename1
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



plot_runtimes <- function(runtimes){
  
  
  ggplot(runtimes, aes(method, runtime))+
    geom_boxplot(aes(fill = trophic_state), outlier.shape = NA)+
    geom_point(aes(color = trophic_state), position = position_jitterdodge(jitter.width = .2), alpha = .6)+
    # stat_summary(aes(color = trophic_state), geom = "pointrange", fun = mean, fun.min = min, fun.max = max, position = position_dodge(width = .2))+
    # stat_summary(aes(color = trophic_state), fun.y = max, geom = "point", size = 2, position = position_dodge(width = .75)) +
    # stat_summary(aes(color = trophic_state), fun.y = min, geom = "point", size = 2, position = position_dodge(width = .75)) +
    scale_color_manual(values = c("black", "black"))+
    labs(y = "Runtime per lake for 100 iterations (seconds)", x = "Method")+
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
    select(any_of(c("lsoda","lsoda_event", "lsoda_event_f", "rk4", "rk4_zero", "patankar-rk2", "patankar-rk2_c", "trophic_state")))
  
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
    na.exclude() %>% 
    pivot_longer(any_of(c("lsoda","lsoda_event", "lsoda_event_f", "rk4", "rk4_zero", "patankar-rk2", "patankar-rk2_c")))
  
  a <- ggplot(plot_df, aes(DO_mgL, value))+
    geom_point(color = "grey")+
    stat_summary(geom = "pointrange", aes(x = plyr::round_any(DO_mgL, 1)))+
    facet_wrap_equal(name~trophic_state, ncol = 2, nrow = 3, scales = "free")+
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
