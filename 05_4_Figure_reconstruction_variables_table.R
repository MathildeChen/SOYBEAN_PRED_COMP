
library(tidyverse)

# > name variables and abbrevations
vars_names <- data.frame(clim.var = c("max_2m_temperature", "min_2m_temperature", 
                                      "et0", "surface_net_solar_radiation", 
                                      "total_precipitation", "vapor_pressure_deficit")) %>% 
  mutate(clim.var_abb = recode(clim.var, 
                               "min_2m_temperature"         ="min_temp",
                               "max_2m_temperature"         ="max_temp",
                               "et0"                        ="et0",
                               "surface_net_solar_radiation"="rad",
                               "total_precipitation"        ="prec",
                               "vapor_pressure_deficit_1"   ="vpd"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature",
                               "max_2m_temperature"         ="Maximum temperature",
                               "et0"                        ="Evapotranspiration ref",
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))

path_to_dim_red <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/"

# --------------------------------------------------------------------------------

exp_var_M_world <- rbind(
  data.frame(clim.var=rep("et0", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(3.42, 3.55, 0.71, 0.03, 0.02, 0.00, 0.00), 
             cumulative.variance.percent=c(3.42, 6.98, 7.69, 7.72, 7.73, 7.73, 7.73)),
  data.frame(clim.var=rep("max_temp", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(21.21,  0.70,  1.09,  0.04,  0.02,  0.00,  0.00), 
             cumulative.variance.percent=c(21.21, 21.91, 23.01, 23.05, 23.07, 23.07, 23.07)),
  data.frame(clim.var=rep("min_temp",7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(20.89,  1.15,  0.26,  0.42,  0.05,  0.00,  0.00), 
             cumulative.variance.percent=c(20.89, 22.04, 22.30, 22.72, 22.77, 22.77, 22.77)),
  data.frame(clim.var=rep("prec", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(18.05,  8.00,  0.42,  0.07,  0.01,  0.00,  0.00), 
             cumulative.variance.percent=c(18.05, 26.05, 26.47, 26.54, 26.55, 26.55, 26.55)),
  data.frame(clim.var=rep("rad", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(28.58,  0.50,  0.28,  0.05,  0.03,  0.00,  0.00), 
             cumulative.variance.percent=c(28.58, 29.09, 29.37, 29.42, 29.45, 29.45, 29.45)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(1.37, 0.33, 0.05, 0.03, 0.02, 0.00, 0.00), 
             cumulative.variance.percent=c(1.37, 1.69, 1.74, 1.77, 1.79, 1.79, 1.79))) %>% 
  mutate(model = "FPLSR",
         data_type = "Monthly averages", 
         type_data = "M",
         country   = "WORLD")

exp_var_D_world <- rbind(
  data.frame(clim.var=rep("et0", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(0.21,  4.63,  1.80,  0.93,  0.65,  0.34,  0.59,  0.56,  0.56,  0.59,  0.27,  0.33,  0.33, 0.24,  0.12,  0.11,  0.09 , 0.06,  0.02,  0.01), 
             cumulative.variance.percent=c(0.21,  4.83,  6.63,  7.57,  8.22,  8.56,  9.15,  9.71, 10.28, 10.87, 11.14, 11.47, 11.80, 12.04, 12.16, 12.27, 12.36, 12.42, 12.44, 12.44)),
  data.frame(clim.var=rep("max_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(21.09,  0.31,  1.45,  0.39,  1.20,  0.45,  0.66,  0.31,  0.39,  0.60,  0.55,  0.36,  0.29, 0.18,  0.20,  0.16,  0.09,  0.07,  0.03,  0.00), 
             cumulative.variance.percent=c(21.09, 21.40, 22.85, 23.24, 24.44, 24.89, 25.55, 25.85, 26.25, 26.85, 27.40, 27.77, 28.05, 28.23, 28.43, 28.60, 28.69, 28.76, 28.78, 28.79)),
  data.frame(clim.var=rep("min_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(20.23,  1.03,  0.47,  1.85,  1.29,  0.35,  0.68,  0.47,  0.43,  0.50,  0.43,  0.32,  0.28, 0.17,  0.19,  0.15,  0.11,  0.05,  0.03,  0.01), 
             cumulative.variance.percent=c(20.23, 21.26, 21.73, 23.58, 24.87, 25.22, 25.90, 26.37, 26.80, 27.30, 27.73, 28.05, 28.33, 28.50, 28.68, 28.83, 28.94, 28.99, 29.02, 29.03)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(14.97,  7.48,  2.93,  0.98,  0.36,  0.31,  0.57,  0.27,  0.31,  0.27,  0.29,  0.20,  0.16, 0.21,  0.14,  0.07,  0.05,  0.04,  0.02,  0.00), 
             cumulative.variance.percent=c(14.97, 22.46, 25.38, 26.37, 26.73, 27.04, 27.61, 27.88, 28.18, 28.46, 28.75, 28.95, 29.11, 29.32, 29.46, 29.53, 29.58, 29.61, 29.63, 29.63)),
  data.frame(clim.var=rep("rad", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(29.11,  0.12,  0.22,  0.17,  0.59,  0.29,  0.43,  0.24,  0.32,  0.25,  0.19,  0.19,  0.15, 0.21,  0.09,  0.06,  0.05,  0.03,  0.01,  0.00), 
             cumulative.variance.percent=c(29.11, 29.23, 29.45, 29.62, 30.21, 30.50, 30.93, 31.17, 31.49, 31.74, 31.94, 32.12, 32.28, 32.49, 32.58, 32.64, 32.69, 32.71, 32.73, 32.73)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(1.29,  0.24,  0.14,  0.39,  0.70,  0.32,  0.23,  0.43,  0.24,  0.49,  0.27,  0.26,  0.10, 0.16,  0.13,  0.12,  0.07,  0.05,  0.02,  0.00), 
             cumulative.variance.percent=c(1.29,  1.53,  1.68,  2.06,  2.76,  3.08,  3.31,  3.74,  3.98,  4.47,  4.74,  5.00,  5.09, 5.25,  5.38,  5.50,  5.57,  5.62,  5.63,  5.64))) %>% 
  mutate(model = "FPLSR",
         data_type = "Daily averages", 
         type_data = "D",
         country   = "WORLD")

# --------------------------------------------------------------------------------

exp_var_M_usa <- rbind(
  data.frame(clim.var=rep("et0", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(20.47,  3.97,  1.06,  0.21,  0.02,  0.00,  0.00), 
             cumulative.variance.percent=c(20.47, 24.44, 25.50, 25.71, 25.73, 25.73, 25.73)),
  data.frame(clim.var=rep("max_temp", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(30.00,  1.20,  0.20,  0.06,  0.00,  0.00,  0.00), 
             cumulative.variance.percent=c(30.00, 31.20, 31.40, 31.46, 31.46, 31.46, 31.46)),
  data.frame(clim.var=rep("min_temp",7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(33.85,  0.88,  0.23,  0.06,  0.00,  0.00,  0.00), 
             cumulative.variance.percent=c(33.85, 34.72, 34.96, 35.01, 35.02, 35.02, 35.02)),
  data.frame(clim.var=rep("prec", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(50.69,  1.01,  0.05,  0.00,  0.00,  0.00,  0.00), 
             cumulative.variance.percent=c(50.69, 51.70, 51.75, 51.75, 51.75, 51.75, 51.75)),
  data.frame(clim.var=rep("rad", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(29.93,  2.90,  0.37,  0.18,  0.03,  0.00,  0.00), 
             cumulative.variance.percent=c(29.93, 32.82, 33.19, 33.37, 33.40, 33.40, 33.40)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(1.34, 1.47, 0.88, 0.64, 0.09, 0.00, 0.00), 
             cumulative.variance.percent=c(1.34, 2.80, 3.68, 4.31, 4.40, 4.40, 4.40))) %>% 
  mutate(model = "FPLSR",
         data_type = "Monthly averages", 
         type_data = "M",
         country   = "USA")

exp_var_D_usa <- rbind(
  data.frame(clim.var=rep("et0", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(0.79, 18.73,  6.23,  0.65,  0.74,  0.74,  0.30,  0.60,  0.40,  0.52,  0.60,  0.42,  0.49,  0.41,  0.21,  0.23,  0.05,  0.04,  0.02,  0.01), 
             cumulative.variance.percent=c(0.79, 19.52, 25.75, 26.39, 27.14, 27.88, 28.18, 28.78, 29.18, 29.70, 30.29, 30.71, 31.21, 31.62, 31.83, 32.05, 32.10, 32.14, 32.17, 32.17)),
  data.frame(clim.var=rep("max_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(29.23,  1.58,  0.41,  0.73,  0.74,  0.54,  0.79,  0.48,  0.62,  0.36,  0.51,  0.67,  0.38,  0.52,  0.25,  0.15,  0.08,  0.07,  0.04,  0.01), 
             cumulative.variance.percent=c(29.23, 30.81, 31.22, 31.95, 32.69, 33.23, 34.01, 34.49, 35.11, 35.47, 35.98, 36.65, 37.03, 37.54, 37.79, 37.94, 38.03, 38.09, 38.13, 38.15)),
  data.frame(clim.var=rep("min_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(33.47,  0.93,  0.36,  1.46,  0.85,  0.92,  0.48,  0.79,  0.33,  0.31,  0.53,  0.28,  0.39,  0.41,  0.25,  0.32,  0.16,  0.08,  0.05,  0.01), 
             cumulative.variance.percent=c(33.47, 34.40, 34.76, 36.22, 37.07, 37.99, 38.48, 39.27, 39.59, 39.90, 40.43, 40.71, 41.10, 41.50, 41.75, 42.07, 42.23, 42.31, 42.36, 42.37)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(49.05,  1.51,  0.56,  1.11,  0.49,  0.47,  0.29,  0.61,  0.48,  0.51,  0.57,  0.53,  0.17,  0.22,  0.17,  0.11,  0.11,  0.05,  0.03,  0.01), 
             cumulative.variance.percent=c(49.05, 50.56, 51.12, 52.23, 52.71, 53.18, 53.47, 54.08, 54.56, 55.07, 55.64, 56.17, 56.35, 56.56, 56.73, 56.83, 56.94, 56.99, 57.03, 57.03)),
  data.frame(clim.var=rep("rad", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(30.19,  1.32,  1.94,  0.31,  0.58,  0.70,  1.16,  0.70,  0.39,  0.82,  0.31,  0.46,  0.38,  0.28,  0.20,  0.26,  0.13,  0.05,  0.04,  0.01), 
             cumulative.variance.percent=c(30.19, 31.52, 33.46, 33.76, 34.34, 35.04, 36.20, 36.90, 37.29, 38.11, 38.42, 38.88, 39.25, 39.53, 39.73, 39.99, 40.11, 40.16, 40.20, 40.21)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(1.21,  0.68,  1.12,  1.37,  0.86,  0.40,  0.36,  0.58,  0.42,  0.77,  0.46,  0.32,  0.37,  0.29,  0.30,  0.23,  0.11,  0.08,  0.04,  0.01), 
             cumulative.variance.percent=c(1.21,  1.89,  3.02,  4.39,  5.25,  5.65,  6.02,  6.60,  7.02,  7.78,  8.24,  8.56,  8.94,  9.23,  9.53,  9.76,  9.86,  9.94,  9.97,  9.99 ))) %>% 
  mutate(model = "FPLSR",
         data_type = "Daily averages", 
         type_data = "D",
         country   = "USA")

# --------------------------------------------------------------------------------

exp_var_M_bra <- rbind(
  data.frame(clim.var=rep("et0", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(9.04, 24.32,  4.06,  0.11,  0.02,  0.00,  0.00), 
             cumulative.variance.percent=c(9.04, 33.36, 37.42, 37.53, 37.55, 37.55, 37.55)),
  data.frame(clim.var=rep("max_temp", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(40.54, 24.33,  0.52,  0.10,  0.02,  0.00,  0.00), 
             cumulative.variance.percent=c(40.54, 64.87, 65.39, 65.49, 65.51, 65.51, 65.51)),
  data.frame(clim.var=rep("min_temp",7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(46.97, 21.60,  1.05,  0.13,  0.01,  0.00,  0.00), 
             cumulative.variance.percent=c(46.97, 68.57, 69.63, 69.75, 69.76, 69.76, 69.76)),
  data.frame(clim.var=rep("prec", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(55.79,  5.56,  0.32,  0.02,  0.00,  0.00,  0.00), 
             cumulative.variance.percent=c(55.79, 61.35, 61.67, 61.69, 61.70, 61.70, 61.70)),
  data.frame(clim.var=rep("rad", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7),  
             variance.percent           =c(50.20,  9.70,  0.85,  0.21,  0.05,  0.00,  0.00), 
             cumulative.variance.percent=c(50.20, 59.90, 60.75, 60.96, 61.02, 61.02, 61.02)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 7), 
             PC=1:7, 
             eigenvalue=rep(NA, 7), 
             variance.percent           =c(7.43, 14.10,  4.37,  0.41,  0.03,  0.00,  0.00), 
             cumulative.variance.percent=c(7.43, 21.53, 25.90, 26.32, 26.35, 26.35, 26.35))) %>% 
  mutate(model = "FPLSR",
         data_type = "Monthly averages", 
         type_data = "M",
         country   = "BRA")

exp_var_D_bra <- rbind(
  data.frame(clim.var=rep("et0", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(5.50, 18.69,  5.25,  0.62,  0.75,  0.26,  0.31,  0.65,  0.53,  0.75,  0.41,  0.32,  0.27,  0.29,  0.18,  0.13,  0.06,  0.04,  0.02,  0.00), 
             cumulative.variance.percent=c(5.50, 24.19, 29.44, 30.06, 30.81, 31.07, 31.39, 32.04, 32.56, 33.32, 33.73, 34.05, 34.32, 34.61, 34.79, 34.91, 34.97, 35.02, 35.04, 35.04)),
  data.frame(clim.var=rep("max_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(36.06, 20.69,  0.28,  0.77,  0.29,  0.43,  0.30,  0.19,  0.14,  0.36,  0.37,  0.24,  0.20,  0.21,  0.23,  0.12,  0.04  ,0.04,  0.02,  0.01), 
             cumulative.variance.percent=c(36.06, 56.74, 57.03, 57.80, 58.09, 58.52, 58.82, 59.01, 59.15, 59.52, 59.89, 60.13, 60.32, 60.53, 60.76, 60.88, 60.92 ,60.96, 60.98, 60.99)),
  data.frame(clim.var=rep("min_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(42.77, 14.71,  0.51,  0.29,  0.33,  0.38,  0.26,  0.38,  0.30,  0.24,  0.25,  0.15,  0.13,  0.12,  0.11,  0.10,  0.06,  0.04,  0.01,  0.00), 
             cumulative.variance.percent=c(42.77, 57.48, 57.99, 58.28, 58.61, 58.99, 59.24, 59.62, 59.92, 60.16, 60.40, 60.55, 60.68, 60.80, 60.90, 61.00, 61.06, 61.10, 61.12, 61.12)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(53.57,  4.70,  2.24,  1.62,  1.21,  0.87,  0.73,  0.37,  0.95,  0.57,  0.62,  0.40,  0.52,  0.36,  0.18,  0.16,  0.10,  0.05,  0.02  ,0.01), 
             cumulative.variance.percent=c(53.57, 58.28, 60.52, 62.14, 63.35, 64.21, 64.94, 65.31, 66.26, 66.83, 67.46, 67.85, 68.37, 68.73, 68.92, 69.08, 69.18, 69.23, 69.25 ,69.25)),
  data.frame(clim.var=rep("rad", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20),  
             variance.percent           =c(29.05, 12.55,  0.60,  2.18,  0.52,  0.51,  0.96,  0.56,  0.83,  0.71,  0.91,  0.56,  0.37,  0.55,  0.33,  0.12,  0.10,  0.07,  0.03,  0.00), 
             cumulative.variance.percent=c(29.05, 41.60, 42.20, 44.38, 44.90, 45.41, 46.37, 46.92, 47.75, 48.46, 49.38, 49.94, 50.31, 50.86, 51.19, 51.31, 51.41, 51.48, 51.51, 51.52)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(5.02,  8.79,  3.64,  1.34,  0.56,  0.29,  0.33,  0.48,  0.33,  0.27,  0.28,  0.37,0.35  ,0.31  ,0.14  ,0.24,  0.14,  0.04,  0.02,  0.00), 
             cumulative.variance.percent=c(5.02, 13.81, 17.45, 18.79, 19.35, 19.64, 19.97, 20.45, 20.78, 21.05, 21.33, 21.70,22.05 ,22.36 ,22.50 ,22.75, 22.89, 22.92, 22.94, 22.95))) %>% 
  mutate(model = "FPLSR",
         data_type = "Daily averages", 
         type_data = "D",
         country   = "BRA")
# --------------------------------------------------------------------------------
exp_var_fplsr <- rbind(exp_var_M_world, exp_var_D_world, exp_var_M_usa, exp_var_D_usa, exp_var_M_bra, exp_var_D_bra) %>% 
  left_join(., vars_names, by = c("clim.var" = "clim.var_abb")) %>%
  mutate(clim.var_lab = recode(clim.var_lab, 
                               'Maximum temperature'=        'Maximum\ntemperature (°C)', 
                               'Minimum temperature'=        'Minimum\ntemperature (°C)', 
                               'Evapotranspiration ref'=       'Evapotransp.\nref (mm/day)', 
                               'Solar radiations'=    'Net solar\nradiations (MJ/m²)', 
                               'Precipitation'=       'Total\nprecipitations (mm)', 
                               'vapor_pressure_deficit'=          'Vapor pressure\ndeficit')) %>% 
  mutate(clim.var_lab = factor(clim.var_lab, levels = c('Maximum\ntemperature (°C)',
                                                        'Minimum\ntemperature (°C)', 
                                                        'Evapotransp.\nref (mm/day)', 
                                                        'Net solar\nradiations (MJ/m²)', 
                                                        'Total\nprecipitations (mm)', 
                                                        'Vapor pressure\ndeficit', 
                                                        'All variables\n(only for MFPCA)')))

# Check
exp_var_fplsr %>%
  group_by(country, type_data, clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent)) %>% 
  mutate(diff=sum_var-last_cum_var) %>% 
  arrange(desc(diff))


# --------------------------------------------------------------------------------

load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/tab_var_dim_red.rda")
tab_var_full <- rbind(tab_var %>% dplyr::select(-clim.var) %>% 
                        dplyr::rename("clim.var"="clim.var_abb") %>% dplyr::select(country, type_data, data_type, model, clim.var, clim.var_lab, PC, eigenvalue, variance.percent, cumulative.variance.percent),
      exp_var_fplsr %>% dplyr::select(country, type_data, data_type, model, clim.var, clim.var_lab, PC, eigenvalue, variance.percent, cumulative.variance.percent))

tab_var_full %>% 
  group_by(country, model, clim.var_lab, data_type) %>% count() %>% View(.)

save(tab_var_full, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/tab_var_dim_red_full.rda")
