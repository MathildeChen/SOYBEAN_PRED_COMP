# -------------------------------------------------------------------------
# 
#           DIMENSION REDUCTION IN CLIMATIC PREDICTORS 
#           GLOBAL - RECONSTRUCTION FIGURE
# 
# -------------------------------------------------------------------------

# TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS
# TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
# TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
# TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION
# TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION

# --> all techniques applied on daily and monthly data 

# ----------------------------------
# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
# > plots
library(cowplot) ; 
# > maps 
library(terra) ; library(raster) ; library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(rmapshaper) ; library(tidygeocoder) ; 
# > exploratory analyses
library(arsenal) ; library(mgcv) ; library(multcomp) ; library(corrplot)
# > PCA
library(FactoMineR) ; library(factoextra) 
# > machine learning
library(parallel) ; library(doParallel); library(foreach)
library(caret) ; library(ranger) ; library(fastshap) ; library(boot) ; library(coxed)
# > functional analysis 
library(fda) ; library(MFPCA) ; library(fda.usc)
# > others
library(hydroGOF)

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
                               "vapor_pressure_deficit"     ="vapor_pressure_deficit"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature (°C)",
                               "max_2m_temperature"         ="Maximum temperature (°C)",
                               "et0"                        ="Evapotransp. ref (mm/day)",
                               "surface_net_solar_radiation"="Net solar radiations (MJ/m²)",
                               "total_precipitation"        ="Total precipitations (mm)",
                               "vapor_pressure_deficit"     ="Vapor pressure deficit")) %>% 
  mutate(clim.var_lab2 = recode(clim.var, 
                                "min_2m_temperature"         ="minimum temperature\n(°C)",
                                "max_2m_temperature"         ="maximum temperature\n(°C)",
                                "et0"                        ="evapotranspiration of reference\n(mm/day)",
                                "surface_net_solar_radiation"="net solar radiations\n(MJ/m²)",
                                "total_precipitation"        ="total precipitations\n(mm)",
                                "vapor_pressure_deficit"     ="vapor pressure deficit\n"))

# ----------------------------------
# > Path to files
path_to_dim_red <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/"
patterns <- c("M_usa", "D_usa", "M_bra", "D_bra", "M_world", "D_world")

exp_var <- list()
# > Load and read
for(p in patterns)
{
  
  # > Files to read 
  tab_files <- data.frame(filename = list.files(path_to_dim_red, pattern = paste0(p, ".rda$")))  %>% 
    # > add variable 
    mutate(technique = case_when(
      str_detect(filename, "mfpca") == T ~ "MFPCA",
      str_detect(filename, "fpca")  == T ~ "FPCA",
      str_detect(filename, "pca")   == T ~ "PCA",
      str_detect(filename, "fplsr") == T ~ "FPLSR",
      str_detect(filename, "plsr")  == T ~ "PLSR"
    ))
  
  exp_var_p <- list()
  # > Load and retrieve the scores from FPLS
  for(i in 1:length(unique(tab_files$filename)))
  {
    
    # Load 
    dim_red_i <- loadRDa(file = paste0(path_to_dim_red, tab_files$filename[i]))
    technique_i <- tab_files$technique[i]
    
    
    # > PCA 
    if(technique_i == "PCA")
    {
      exp_var_i <- dim_red_i$list_pca_per_variable %>% 
        map_dfr(., ~{ .x$pca.eig %>% 
            mutate(PC = 1:nrow(.x$pca.eig))
        }, .id="clim.var")%>% 
        dplyr::select(clim.var, PC, eigenvalue, variance.percent, cumulative.variance.percent) %>% 
        mutate(model = technique_i)
      
      exp_var_p[[paste0("PCA")]] <- exp_var_i
      
    }
    
    # > FPCA
    if(technique_i == "FPCA")
    {
      exp_var_i <- dim_red_i$list_fpca_per_variable %>% 
        map_dfr(., ~{
          
          data.frame(PC = 1:length(.x$values),
                     eigenvalue = .x$values) %>% 
            # > total eigen value
            mutate(sum_eig=sum(eigenvalue)) %>% 
            group_by(PC) %>%
            # > compute explained variance for each component
            mutate(variance.percent=100*(eigenvalue/sum_eig)) %>% 
            ungroup() %>% 
            # > compute cumulated explained variance
            mutate(cumulative.variance.percent = cumsum(variance.percent))
        }, .id="clim.var") %>% 
        dplyr::select(clim.var, PC, eigenvalue, variance.percent, cumulative.variance.percent) %>% 
        mutate(model = technique_i)
      
      exp_var_p[[paste0("FPCA")]] <- exp_var_i
      
    }
    
    # > MFPCA
    if(technique_i == "MFPCA")
    { 
      # > MFPCA - 1 table for all variables
      exp_var_i <- data.frame(PC = 1:length(dim_red_i$list_mfpca$values),
                              eigenvalue = dim_red_i$list_mfpca$values) %>% 
        # > total eigen value
        mutate(sum_eig=sum(eigenvalue)) %>% 
        group_by(PC) %>%
        # > compute explained variance for each component
        mutate(variance.percent=100*(eigenvalue/sum_eig)) %>%  
        ungroup() %>% 
        # > compute cumulated explained variance
        mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
        dplyr::select(PC, eigenvalue, variance.percent, cumulative.variance.percent) %>% 
        mutate(clim.var = "All variables\n(only for MFPCA)") %>% 
        dplyr::select(clim.var, PC, eigenvalue, variance.percent, cumulative.variance.percent) %>% 
        mutate(model = technique_i)
      
      exp_var_p[[paste0("MFPCA")]] <- exp_var_i
    
    }
    
    # > PLSR
    if(technique_i == "PLSR")
    {
    
      exp_var_i <- dim_red_i$list_pls_per_variable %>% 
        map_dfr(., ~{
          
          data.frame(PC = 1:length(.x$Xvar),
                     eigenvalue = .x$Xvar) %>% 
            # > total eigen value
            mutate(sum_eig=.x$Xtotvar) %>% 
            group_by(PC) %>%
            # > compute explained variance for each component
            mutate(variance.percent=100*(eigenvalue/sum_eig)) %>% 
            ungroup() %>% 
            # > compute cumulated explained variance
            mutate(cumulative.variance.percent = cumsum(variance.percent))
        }, .id = "clim.var") %>% 
        dplyr::select(clim.var, PC, eigenvalue, variance.percent, cumulative.variance.percent) %>% 
        mutate(model = technique_i)
      
      exp_var_p[[paste0("PLSR")]] <- exp_var_i
    }
    
    # > FPLSR
    if(technique_i == "FPLSR")
    {
    
      exp_var_p[[paste0("FPLSR")]] <- NULL
      
    }
    
    tab_exp_var_p <- plyr::ldply(exp_var_p, data.frame, .id = "model_check") 
    testthat::expect_equivalent(as.character(tab_exp_var_p$model), as.character(tab_exp_var_p$model_check))
    
    exp_var[[paste0(p)]] <- tab_exp_var_p
    
  }
  
}


# ---------------------------------------
# EXTRACT SCORES FOR FPLSR for global dataset
# > on monthly data 
# > Path to files
path_to_fplsr_M <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/fplsr_M_world/"
files_M <- list.files(path_to_fplsr_M, pattern = ".rda$")

#[1] "fplsr_M_world_var_et0.rda"                   
#[2] "fplsr_M_world_var_max_temp.rda"              
#[3] "fplsr_M_world_var_min_temp.rda"              
#[4] "fplsr_M_world_var_prec.rda"                  
#[5] "fplsr_M_world_var_rad.rda"                   
#[6] "fplsr_M_world_var_vapor_pressure_deficit.rda"

fplsr_et0 <- loadRDa(file = paste0(path_to_fplsr_M, files_M[1]))
summary(fplsr_et0$list_fpls_per_variable$et0, draw=F)
#- R^2 by component (%)
#PLS1 PLS2 PLS3 PLS4 PLS5 
#3.42 3.54 0.70 0.03 0.02 
#- Cumulative R^2 (%)
#PLS1 PLS2 PLS3 PLS4 PLS5 
#3.42 6.96 7.66 7.69 7.71 

fplsr_max_temp <- loadRDa(file = paste0(path_to_fplsr_M, files_M[2]))
summary(fplsr_max_temp$list_fpls_per_variable[[1]], draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#21.21  0.70  1.09  0.04  0.02 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5 
#21.21 21.91 23.00 23.04 23.06

fplsr_min_temp <- loadRDa(file = paste0(path_to_fplsr_M, files_M[3]))
summary(fplsr_min_temp$list_fpls_per_variable[[1]], draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#20.89  1.15  0.26  0.41  0.05 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5 
#20.89 22.04 22.30 22.71 22.76 

fplsr_prec <- loadRDa(file = paste0(path_to_fplsr_M, files_M[4]))
summary(fplsr_prec$list_fpls_per_variable[[1]], draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4
#17.96  7.76  0.40  0.04
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4 
#17.96 25.73 26.13 26.17

fplsr_rad <- loadRDa(file = paste0(path_to_fplsr_M, files_M[5]))
summary(fplsr_rad$list_fpls_per_variable[[1]], draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#28.55  0.50  0.26  0.05  0.03
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5 
#28.55 29.05 29.32 29.37 29.40 

fplsr_vpd <- loadRDa(file = paste0(path_to_fplsr_M, files_M[6]))
summary(fplsr_vpd$list_fpls_per_variable[[1]], draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#1.37 0.33 0.05 0.03 0.02
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5 
#1.37 1.69 1.74 1.77 1.79

exp_var_M_world <- rbind(
  data.frame(clim.var=rep("et0", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(3.42, 3.54, 0.70, 0.03, 0.02), 
             cumulative.variance.percent=c(3.42, 6.96, 7.66, 7.69, 7.71),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("max_temp", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(21.21, 0.70, 1.09, 0.04, 0.02), 
             cumulative.variance.percent=c(21.21, 21.91, 23.00, 23.04, 23.06),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("min_temp", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(20.89, 1.15, 0.26, 0.41, 0.05), 
             cumulative.variance.percent=c(20.89, 22.04, 22.30, 22.71, 22.76),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("prec", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(17.96, 7.76, 0.40, 0.04), 
             cumulative.variance.percent=c(17.96, 25.73, 26.13, 26.17),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("rad", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(28.55, 0.50, 0.26, 0.05, 0.03), 
             cumulative.variance.percent=c(28.55, 29.05, 29.32, 29.37, 29.40),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(1.37,  0.33,  0.05,  0.03,  0.02), 
             cumulative.variance.percent=c(1.37,  1.69,  1.74,  1.77,  1.79),
             model = rep("FPLSR", 5))
)

# check if equal in all variable
exp_var_M_world %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
         last_cum_var = last(cumulative.variance.percent))

# > Add to other data 
exp_var$M_world <- rbind(exp_var$M_world %>% dplyr::select(-model_check), exp_var_M_world)
rm(fplsr_et0, fplsr_max_temp, fplsr_min_temp, fplsr_prec, fplsr_rad, fplsr_vpd, exp_var_M_world)

# > on daily data 
# > Path to files
path_to_fplsr_D <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/fplsr_D_world/"
files_D <- list.files(path_to_fplsr_D, pattern = ".rda$") ; files_D
#[1] "fplsr_D_world_var_et0.rda"                   
#[2] "fplsr_D_world_var_max_temp.rda"              
#[3] "fplsr_D_world_var_min_temp.rda"              
#[4] "fplsr_D_world_var_prec.rda"                  
#[5] "fplsr_D_world_var_rad.rda"                   
#[6] "fplsr_D_world_var_vapor_pressure_deficit.rda"

fplsr_et0 <- loadRDa(file = paste0(path_to_fplsr_D, files_D[1]))
summary(fplsr_et0$list_fpls_per_variable$et0, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#0.21  4.63  1.80  0.93  0.65  0.34  0.59  0.56  0.56  0.59  0.27  0.33  0.33  0.24  0.12  0.11  0.09  0.06  0.02  0.01
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 0.21  4.83  6.63  7.57  8.22  8.56  9.15  9.71 10.28 10.87 11.14 11.47 11.80 12.04 12.16 12.27 12.36 12.42 12.44 12.44 

fplsr_max_temp <- loadRDa(file = paste0(path_to_fplsr_D, files_D[2]))
summary(fplsr_max_temp$list_fpls_per_variable$max_temp, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#21.09  0.31  1.45  0.39  1.20  0.45  0.66  0.31  0.39  0.60  0.55  0.36  0.29  0.18  0.20  0.16  0.09  0.07  0.03  0.00
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#21.09 21.40 22.85 23.24 24.44 24.89 25.55 25.85 26.25 26.85 27.40 27.77 28.05 28.23 28.43 28.60 28.69 28.76 28.78 28.79

fplsr_min_temp <- loadRDa(file = paste0(path_to_fplsr_D, files_D[3]))
summary(fplsr_min_temp$list_fpls_per_variable$min_temp, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#20.23  1.03  0.47  1.85  1.29  0.35  0.68  0.47  0.43  0.50  0.43  0.32  0.28  0.17  0.19  0.15  0.11  0.05  0.03  0.01
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#20.23 21.26 21.73 23.58 24.87 25.22 25.90 26.37 26.80 27.30 27.73 28.05 28.33 28.50 28.68 28.83 28.94 28.99 29.02 29.03

fplsr_prec <- loadRDa(file = paste0(path_to_fplsr_D, files_D[4]))
summary(fplsr_prec$list_fpls_per_variable$prec, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#14.97  7.48  2.93  0.98  0.36  0.31  0.57  0.27  0.31  0.27  0.29  0.20  0.16  0.21  0.14  0.07  0.05  0.04  0.02  0.00
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#14.97 22.46 25.38 26.37 26.73 27.04 27.61 27.88 28.18 28.46 28.75 28.95 29.11 29.32 29.46 29.53 29.58 29.61 29.63 29.63

fplsr_rad <- loadRDa(file = paste0(path_to_fplsr_D, files_D[5]))
summary(fplsr_rad$list_fpls_per_variable$rad, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#29.11  0.12  0.22  0.17  0.59  0.29  0.43  0.24  0.32  0.25  0.19  0.19  0.15  0.21  0.09  0.06  0.05  0.03  0.01  0.00
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20  
#29.11 29.23 29.45 29.62 30.21 30.50 30.93 31.17 31.49 31.74 31.94 32.12 32.28 32.49 32.58 32.64 32.69 32.71 32.73 32.73 

fplsr_vpd <- loadRDa(file = paste0(path_to_fplsr_D, files_D[6]))
summary(fplsr_vpd$list_fpls_per_variable$vapor_pressure_deficit, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#1.29  0.24  0.14  0.39  0.70  0.32  0.23  0.43  0.24  0.49  0.27  0.26  0.10  0.16  0.13  0.12  0.07  0.05  0.02  0.00
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 LS5 
#1.29  1.53  1.68  2.06  2.76  3.08  3.31  3.74  3.98  4.47  4.74  5.00  5.09  5.25  5.38  5.50  5.57  5.62  5.63  5.64

exp_var_D_world <- rbind(
  data.frame(clim.var=rep("et0", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(0.21,  4.63,  1.80,  0.93,  0.65,  0.34,  0.59,  0.56,  0.56,  0.59,  0.27,  0.33,  0.33,  0.24,  0.12,  0.11,  0.09,  0.06,  0.02,  0.01), 
             cumulative.variance.percent=c(0.21,  4.83,  6.63,  7.57,  8.22,  8.56,  9.15,  9.71, 10.28, 10.87, 11.14, 11.47, 11.80, 12.04, 12.16, 12.27, 12.36, 12.42, 12.44, 12.44),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("max_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(21.09,  0.31,  1.45,  0.39,  1.20,  0.45,  0.66,  0.31,  0.39,  0.60, 0.55, 0.36, 0.29, 0.18, 0.20, 0.16, 0.09, 0.07, 0.03, 0.00), 
             cumulative.variance.percent=c(21.09, 21.40, 22.85, 23.24, 24.44, 24.89, 25.55, 25.85, 26.25, 26.85, 27.40, 27.77, 28.05, 28.23, 28.43, 28.60, 28.69, 28.76, 28.78, 28.79),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("min_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(20.23, 1.03, 0.47, 1.85, 1.29, 0.35, 0.68, 0.47, 0.43, 0.50, 0.43, 0.32, 0.28, 0.17, 0.19, 0.15, 0.11, 0.05, 0.03, 0.01), 
             cumulative.variance.percent=c(20.23, 21.26, 21.73, 23.58, 24.87, 25.22, 25.90, 26.37, 26.80, 27.30, 27.73, 28.05, 28.33, 28.50, 28.68, 28.83, 28.94, 28.99, 29.02, 29.03),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(14.97, 7.48, 2.93, 0.98, 0.36, 0.31, 0.57, 0.27, 0.31, 0.27, 0.29, 0.20, 0.16, 0.21, 0.14, 0.07, 0.05, 0.04, 0.02, 0.00), 
             cumulative.variance.percent=c(14.97, 22.46, 25.38, 26.37, 26.73, 27.04, 27.61, 27.88, 28.18, 28.46, 28.75, 28.95, 29.11, 29.32, 29.46, 29.53, 29.58, 29.61, 29.63, 29.63),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("rad", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(29.11, 0.12, 0.22, 0.17, 0.59, 0.29, 0.43, 0.24, 0.32, 0.25, 0.19, 0.19, 0.15, 0.21, 0.09, 0.06, 0.05, 0.03, 0.01, 0.00), 
             cumulative.variance.percent=c(29.11, 29.23, 29.45, 29.62, 30.21, 30.50, 30.93, 31.17, 31.49, 31.74, 31.94, 32.12, 32.28, 32.49, 32.58, 32.64, 32.69, 32.71, 32.73, 32.73),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(1.29, 0.24, 0.14, 0.39, 0.70, 0.32, 0.23, 0.43, 0.24, 0.49, 0.27, 0.26, 0.10, 0.16, 0.13, 0.12, 0.07, 0.05, 0.02, 0.00), 
             cumulative.variance.percent=c(1.29,  1.53,  1.68,  2.06,  2.76,  3.08,  3.31,  3.74,  3.98,  4.47,  4.74,  5.00,  5.09,  5.25,  5.38,  5.50,  5.57,  5.62,  5.63,  5.64),
             model = rep("FPLSR", 20))
  
)

# check if equal in all variable
exp_var_D_world %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent))

exp_var$D_world <- rbind(exp_var$D_world %>% dplyr::select(-model_check), exp_var_D_world)
rm(fplsr_et0, fplsr_max_temp, fplsr_min_temp, fplsr_prec, fplsr_rad, fplsr_vpd, exp_var_D_world)

# ----------------------------------
# EXTRACT SCORES FOR FPLSR for USA

# > on monthly data 
fplsr_M_usa <- loadRDa(file = paste0(path_to_dim_red, "fplsr_M_usa.rda"))
summary(fplsr_M_usa$list_fpls_per_variable$et0)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4 
#20.35  3.95  1.04  0.14 
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4 
#20.35 24.30 25.34 25.48 

summary(fplsr_M_usa$list_fpls_per_variable$max_temp, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4 
#30.00  1.19  0.20  0.06 
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4 
#30.00 31.18 31.38 31.44

summary(fplsr_M_usa$list_fpls_per_variable$min_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4 
#33.85  0.86  0.22  0.05 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4 
#33.85 34.70 34.92 34.97

summary(fplsr_M_usa$list_fpls_per_variable$prec, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  
#49.83  0.94  0.04
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  
#49.83 50.76 50.80

summary(fplsr_M_usa$list_fpls_per_variable$rad, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4 PLS5
#29.89  2.73  0.36  0.18  0.03
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4 PLS5
#29.89 32.63 32.99 33.17 33.20

summary(fplsr_M_usa$list_fpls_per_variable$vapor_pressure_deficit, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5
#1.34  1.46  0.87  0.62  0.08
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5
#1.34  2.79  3.67  4.28  4.36

exp_var_fplsr_M_usa <- rbind(
  data.frame(clim.var=rep("et0", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(20.35, 3.95, 1.04, 0.14), 
             cumulative.variance.percent=c(20.35, 24.30, 25.34, 25.48),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("max_temp", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(30.00,  1.19,  0.20,  0.06), 
             cumulative.variance.percent=c(30.00, 31.18, 31.38, 31.44),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("min_temp", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(33.85,  0.86,  0.22,  0.05), 
             cumulative.variance.percent=c(33.85, 34.70, 34.92, 34.97),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("prec", 3), 
             PC=1:3, 
             eigenvalue=rep(NA, 3), 
             variance.percent           =c(49.83,  0.94,  0.04), 
             cumulative.variance.percent=c(49.83, 50.76, 50.80),
             model = rep("FPLSR", 3)),
  data.frame(clim.var=rep("rad", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(29.89,  2.73,  0.36,  0.18,  0.03), 
             cumulative.variance.percent=c(29.89, 32.63, 32.99, 33.17, 33.20),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(1.34, 1.46, 0.87, 0.62, 0.08), 
             cumulative.variance.percent=c(1.34, 2.79, 3.67, 4.28, 4.36),
             model = rep("FPLSR", 5))
)

# check if equal in all variable
exp_var_fplsr_M_usa %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent))

# > Add to other data 
exp_var$M_usa <- rbind(exp_var$M_usa %>% dplyr::select(-model_check), exp_var_fplsr_M_usa)
rm(fplsr_M_usa, exp_var_fplsr_M_usa)

# > on daily data 
fplsr_D_usa <- loadRDa(file = paste0(path_to_dim_red, "fplsr_D_usa.rda"))
summary(fplsr_D_usa$list_fpls_per_variable$et0)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 0.79 18.73  6.23  0.65  0.74  0.74  0.30  0.60  0.40  0.52  0.60  0.42  0.49  0.41  0.21  0.23  0.05  0.04  0.02  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 0.79 19.52 25.75 26.39 27.14 27.88 28.18 28.78 29.18 29.70 30.29 30.71 31.21 31.62 31.83 32.05 32.10 32.14 32.17 32.17

summary(fplsr_D_usa$list_fpls_per_variable$max_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#29.23  1.58  0.41  0.73  0.74  0.54  0.79  0.48  0.62  0.36  0.51  0.67  0.38  0.52  0.25  0.15  0.08  0.07  0.04  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#29.23 30.81 31.22 31.95 32.69 33.23 34.01 34.49 35.11 35.47 35.98 36.65 37.03 37.54 37.79 37.94 38.03 38.09 38.13 38.15

summary(fplsr_D_usa$list_fpls_per_variable$min_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20  
#33.47  0.93  0.36  1.46  0.85  0.92  0.48  0.79  0.33  0.31  0.53  0.28  0.39  0.41  0.25  0.32  0.16  0.08  0.05  0.01  
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20  
#33.47 34.40 34.76 36.22 37.07 37.99 38.48 39.27 39.59 39.90 40.43 40.71 41.10 41.50 41.75 42.07 42.23 42.31 42.36 42.37

summary(fplsr_D_usa$list_fpls_per_variable$prec, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#49.05  1.51  0.56  1.11  0.49  0.47  0.29  0.61  0.48  0.51  0.57  0.53  0.17  0.22  0.17  0.11  0.11  0.05  0.03  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#49.05 50.56 51.12 52.23 52.71 53.18 53.47 54.08 54.56 55.07 55.64 56.17 56.35 56.56 56.73 56.83 56.94 56.99 57.03 57.03

summary(fplsr_D_usa$list_fpls_per_variable$rad, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#30.19  1.32  1.94  0.31  0.58  0.70  1.16  0.70  0.39  0.82  0.31  0.46  0.38  0.28  0.20  0.26  0.13  0.05  0.04  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#30.19 31.52 33.46 33.76 34.34 35.04 36.20 36.90 37.29 38.11 38.42 38.88 39.25 39.53 39.73 39.99 40.11 40.16 40.20 40.21

summary(fplsr_D_usa$list_fpls_per_variable$vapor_pressure_deficit, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 1.21  0.68  1.12  1.37  0.86  0.40  0.36  0.58  0.42  0.77  0.46  0.32  0.37  0.29  0.30  0.23  0.11  0.08  0.04  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 1.21  1.89  3.02  4.39  5.25  5.65  6.02  6.60  7.02  7.78  8.24  8.56  8.94  9.23  9.53  9.76  9.86  9.94  9.97  9.99

exp_var_fplsr_D_usa <- rbind(
  data.frame(clim.var=rep("et0", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(0.79, 18.73,  6.23,  0.65,  0.74,  0.74,  0.30,  0.60,  0.40,  0.52,  0.60,  0.42,  0.49,  0.41,  0.21,  0.23,  0.05,  0.04,  0.02,  0.01), 
             cumulative.variance.percent=c(0.79, 19.52, 25.75, 26.39, 27.14, 27.88, 28.18, 28.78, 29.18, 29.70, 30.29, 30.71, 31.21, 31.62, 31.83, 32.05, 32.10, 32.14, 32.17, 32.17),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("max_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(29.23,  1.58,  0.41,  0.73,  0.74,  0.54,  0.79,  0.48,  0.62,  0.36,  0.51,  0.67,  0.38,  0.52,  0.25,  0.15,  0.08,  0.07,  0.04,  0.01), 
             cumulative.variance.percent=c(29.23, 30.81, 31.22, 31.95, 32.69, 33.23, 34.01, 34.49, 35.11, 35.47, 35.98, 36.65, 37.03, 37.54, 37.79, 37.94, 38.03, 38.09, 38.13, 38.15),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("min_temp", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(33.47,  0.93,  0.36,  1.46,  0.85,  0.92,  0.48,  0.79,  0.33,  0.31,  0.53,  0.28,  0.39,  0.41,  0.25,  0.32,  0.16,  0.08,  0.05,  0.01), 
             cumulative.variance.percent=c(33.47, 34.40, 34.76, 36.22, 37.07, 37.99, 38.48, 39.27, 39.59, 39.90, 40.43, 40.71, 41.10, 41.50, 41.75, 42.07, 42.23, 42.31, 42.36, 42.37),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("rad", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(30.19,  1.32,  1.94,  0.31,  0.58,  0.70,  1.16,  0.70,  0.39,  0.82,  0.31,  0.46,  0.38,  0.28,  0.20,  0.26,  0.13,  0.05,  0.04,  0.01), 
             cumulative.variance.percent=c(30.19, 31.52, 33.46, 33.76, 34.34, 35.04, 36.20, 36.90, 37.29, 38.11, 38.42, 38.88, 39.25, 39.53, 39.73, 39.99, 40.11, 40.16, 40.20, 40.21),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(49.05,  1.51,  0.56,  1.11,  0.49,  0.47,  0.29,  0.61,  0.48,  0.51,  0.57,  0.53,  0.17,  0.22,  0.17,  0.11,  0.11,  0.05,  0.03,  0.01), 
             cumulative.variance.percent=c(49.05, 50.56, 51.12, 52.23, 52.71, 53.18, 53.47, 54.08, 54.56, 55.07, 55.64, 56.17, 56.35, 56.56, 56.73, 56.83, 56.94, 56.99, 57.03, 57.03),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(1.21,  0.68,  1.12,  1.37,  0.86,  0.40,  0.36,  0.58,  0.42,  0.77,  0.46,  0.32,  0.37,  0.29,  0.30,  0.23,  0.11,  0.08,  0.04,  0.01), 
             cumulative.variance.percent=c(1.21,  1.89,  3.02,  4.39,  5.25,  5.65,  6.02,  6.60,  7.02,  7.78,  8.24,  8.56,  8.94,  9.23,  9.53,  9.76,  9.86,  9.94,  9.97,  9.99),
             model = rep("FPLSR", 20))
)

# check if equal in all variable
exp_var_fplsr_D_usa %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent))

# > Add to other data 
exp_var$D_usa <- rbind(exp_var$D_usa %>% dplyr::select(-model_check), exp_var_fplsr_D_usa)
rm(fplsr_D_usa, exp_var_fplsr_D_usa)

# ----------------------------------
# EXTRACT SCORES FOR FPLSR for BRA

# > on monthly data 
fplsr_M_bra <- loadRDa(file = paste0(path_to_dim_red, "fplsr_M_bra.rda"))
summary(fplsr_M_bra$list_fpls_per_variable$et0)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4 
#9.04 24.04  4.05  0.08 
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4 
#9.04 33.08 37.13 37.21 

summary(fplsr_M_bra$list_fpls_per_variable$max_temp, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5
#40.54 24.32  0.52  0.10  0.02 
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#40.54 64.85 65.37 65.47 65.49

summary(fplsr_M_bra$list_fpls_per_variable$min_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4 
#46.96 21.57  1.05  0.12 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4 
#46.96 68.53 69.58 69.70

summary(fplsr_M_bra$list_fpls_per_variable$prec, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4 PLS5
#55.50  5.55  0.31  0.02
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4 PLS5
#55.50 61.05 61.36 61.38

summary(fplsr_M_bra$list_fpls_per_variable$rad, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#50.06  9.67  0.82  0.21  0.05 
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5 
#50.06 59.73 60.55 60.76 60.81

summary(fplsr_M_bra$list_fpls_per_variable$vapor_pressure_deficit, draw=F)
#- R^2 by component (%)
#PLS1  PLS2  PLS3  PLS4  PLS5
#7.43 13.90  4.35  0.29
#- Cumulative R^2 (%)
#PLS1  PLS2  PLS3  PLS4  PLS5
#7.43 21.33 25.68 25.98

exp_var_fplsr_M_bra <- rbind(
  data.frame(clim.var=rep("et0", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(9.04, 24.04,  4.05,  0.08), 
             cumulative.variance.percent=c(9.04, 33.08, 37.13, 37.21),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("max_temp", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(40.54, 24.32,  0.52,  0.10,  0.02), 
             cumulative.variance.percent=c(40.54, 64.85, 65.37, 65.47, 65.49),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("min_temp", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(46.96, 21.57,  1.05,  0.12), 
             cumulative.variance.percent=c(46.96, 68.53, 69.58, 69.70),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("prec", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(55.50,  5.55,  0.31,  0.02), 
             cumulative.variance.percent=c(55.50, 61.05, 61.36, 61.38),
             model = rep("FPLSR", 4)),
  data.frame(clim.var=rep("rad", 5), 
             PC=1:5, 
             eigenvalue=rep(NA, 5), 
             variance.percent           =c(50.06,  9.67,  0.82,  0.21,  0.05), 
             cumulative.variance.percent=c(50.06, 59.73, 60.55, 60.76, 60.81),
             model = rep("FPLSR", 5)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 4), 
             PC=1:4, 
             eigenvalue=rep(NA, 4), 
             variance.percent           =c(7.43, 13.90,  4.35,  0.29), 
             cumulative.variance.percent=c(7.43, 21.33, 25.68, 25.98),
             model = rep("FPLSR", 4))
)

# check if equal in all variable
exp_var_fplsr_M_bra %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent))

# > Add to other data 
exp_var$M_bra <- rbind(exp_var$M_bra %>% dplyr::select(-model_check), exp_var_fplsr_M_bra)
rm(fplsr_M_bra, exp_var_fplsr_M_bra)

# > on daily data 
fplsr_D_bra <- loadRDa(file = paste0(path_to_dim_red, "fplsr_D_bra.rda"))
summary(fplsr_D_bra$list_fpls_per_variable$et0)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15  
# 5.50 18.69  5.25  0.62  0.75  0.26  0.30  0.59  0.47  0.63  0.27  0.19  0.12  0.10  0.03  
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15  
# 5.50 24.19 29.44 30.06 30.81 31.06 31.36 31.96 32.43 33.06 33.33 33.52 33.64 33.74 33.77 

summary(fplsr_D_bra$list_fpls_per_variable$max_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 
#36.06 20.69  0.28  0.77  0.29  0.42  0.29  0.16  0.11  0.23  0.21  0.12  0.07  0.06  0.04  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 
#36.06 56.74 57.03 57.80 58.08 58.50 58.79 58.95 59.06 59.29 59.50 59.62 59.69 59.74 59.78 59.79

summary(fplsr_D_bra$list_fpls_per_variable$min_temp, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 
#42.77 14.71  0.51  0.29  0.32  0.37  0.25  0.35  0.26  0.20  0.17  0.08  0.06  0.04  0.02  0.00 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 
#42.77 57.48 57.99 58.28 58.60 58.97 59.21 59.56 59.83 60.03 60.19 60.28 60.33 60.38 60.40 60.40

summary(fplsr_D_bra$list_fpls_per_variable$prec, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#53.57  4.70  2.24  1.62  1.21  0.87  0.73  0.37  0.95  0.57  0.62  0.40  0.52  0.36  0.18  0.16  0.10  0.05  0.02  0.01 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#53.57 58.28 60.52 62.14 63.35 64.21 64.94 65.31 66.26 66.83 67.46 67.85 68.37 68.73 68.92 69.08 69.18 69.23 69.25 69.25

summary(fplsr_D_bra$list_fpls_per_variable$rad, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#29.05 12.55  0.60  2.18  0.52  0.51  0.96  0.56  0.83  0.71  0.91  0.56  0.37  0.55  0.33  0.12  0.10  0.07  0.03  0.00 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
#29.05 41.60 42.20 44.38 44.90 45.41 46.37 46.92 47.75 48.46 49.38 49.94 50.31 50.86 51.19 51.31 51.41 51.48 51.51 51.52

summary(fplsr_D_bra$list_fpls_per_variable$vapor_pressure_deficit, draw=F)
#- R^2 by component (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 5.02  8.79  3.64  1.34  0.56  0.29  0.33  0.48  0.33  0.27  0.28  0.37  0.35  0.31  0.14  0.24  0.14  0.04  0.02  0.00 
#- Cumulative R^2 (%)
# PLS1  PLS2  PLS3  PLS4  PLS5  PLS6  PLS7  PLS8  PLS9 PLS10 PLS11 PLS12 PLS13 PLS14 PLS15 PLS16 PLS17 PLS18 PLS19 PLS20 
# 5.02 13.81 17.45 18.79 19.35 19.64 19.97 20.45 20.78 21.05 21.33 21.70 22.05 22.36 22.50 22.75 22.89 22.92 22.94 22.95

exp_var_fplsr_D_bra <- rbind(
  data.frame(clim.var=rep("et0", 15), 
             PC=1:15, 
             eigenvalue=rep(NA, 15), 
             variance.percent           =c(5.50, 18.69,  5.25,  0.62,  0.75,  0.26,  0.30,  0.59,  0.47,  0.63,  0.27,  0.19,  0.12,  0.10,  0.03), 
             cumulative.variance.percent=c(5.50, 24.19, 29.44, 30.06, 30.81, 31.06, 31.36, 31.96, 32.43, 33.06, 33.33, 33.52, 33.64, 33.74, 33.77),
             model = rep("FPLSR", 15)),
  data.frame(clim.var=rep("max_temp", 16), 
             PC=1:16, 
             eigenvalue=rep(NA, 16), 
             variance.percent           =c(36.06, 20.69,  0.28,  0.77,  0.29,  0.42,  0.29,  0.16,  0.11,  0.23,  0.21,  0.12,  0.07,  0.06,  0.04,  0.01), 
             cumulative.variance.percent=c(36.06, 56.74, 57.03, 57.80, 58.08, 58.50, 58.79, 58.95, 59.06, 59.29, 59.50, 59.62, 59.69, 59.74, 59.78, 59.79),
             model = rep("FPLSR", 16)),
  data.frame(clim.var=rep("min_temp", 4), 
             PC=1:16, 
             eigenvalue=rep(NA, 16), 
             variance.percent           =c(42.77, 14.71,  0.51,  0.29,  0.32,  0.37,  0.25,  0.35,  0.26,  0.20,  0.17,  0.08,  0.06,  0.04,  0.02,  0.00), 
             cumulative.variance.percent=c(42.77, 57.48, 57.99, 58.28, 58.60, 58.97, 59.21, 59.56, 59.83, 60.03, 60.19, 60.28, 60.33, 60.38, 60.40, 60.40),
             model = rep("FPLSR", 16)),
  data.frame(clim.var=rep("prec", 20), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(53.57,  4.70,  2.24,  1.62,  1.21,  0.87,  0.73,  0.37,  0.95,  0.57,  0.62,  0.40,  0.52,  0.36,  0.18,  0.16,  0.10,  0.05,  0.02,  0.01), 
             cumulative.variance.percent=c(53.57, 58.28, 60.52, 62.14, 63.35, 64.21, 64.94, 65.31, 66.26, 66.83, 67.46, 67.85, 68.37, 68.73, 68.92, 69.08, 69.18, 69.23, 69.25, 69.25),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("rad", 5), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(29.05, 12.55,  0.60,  2.18,  0.52,  0.51,  0.96,  0.56,  0.83,  0.71,  0.91,  0.56,  0.37,  0.55,  0.33,  0.12,  0.10,  0.07,  0.03,  0.00), 
             cumulative.variance.percent=c(29.05, 41.60, 42.20, 44.38, 44.90, 45.41, 46.37, 46.92, 47.75, 48.46, 49.38, 49.94, 50.31, 50.86, 51.19, 51.31, 51.41, 51.48, 51.51, 51.52),
             model = rep("FPLSR", 20)),
  data.frame(clim.var=rep("vapor_pressure_deficit", 4), 
             PC=1:20, 
             eigenvalue=rep(NA, 20), 
             variance.percent           =c(5.02,  8.79,  3.64,  1.34,  0.56,  0.29,  0.33,  0.48,  0.33,  0.27,  0.28,  0.37,  0.35,  0.31,  0.14,  0.24,  0.14,  0.04,  0.02,  0.00), 
             cumulative.variance.percent=c(5.02, 13.81, 17.45, 18.79, 19.35, 19.64, 19.97, 20.45, 20.78, 21.05, 21.33, 21.70, 22.05, 22.36, 22.50, 22.75, 22.89, 22.92, 22.94, 22.95),
             model = rep("FPLSR", 20))
)

# check if equal in all variable
exp_var_fplsr_D_bra %>% 
  group_by(clim.var) %>% 
  summarise(sum_var = sum(variance.percent),
            last_cum_var = last(cumulative.variance.percent))

# > Add to other data 
exp_var$D_bra <- rbind(exp_var$D_bra %>% dplyr::select(-model_check), exp_var_fplsr_D_bra)
rm(fplsr_D_bra, exp_var_fplsr_D_bra)


tab_var <- plyr::ldply(exp_var, data.frame, .id = "type_data") %>% 
  separate(col = "type_data", into = c("type_data", "country"), remove = T) %>% 
  # > formatting
  mutate(country = recode(country, 
                          "usa"="USA",
                          "bra"="BRA",
                          "world"="WORLD")) %>% 
  mutate(data_type = recode(type_data,
                            "M"="Monthly averages",
                            "D"="Daily averages")) %>% 
  mutate(model = factor(model, levels=c("PCA", "FPCA", "MFPCA", "PLSR", "FPLSR"))) %>%
  dplyr::rename("clim.var_abb"="clim.var") %>% 
  left_join(., vars_names, by ="clim.var_abb") %>% 
  mutate(clim.var = if_else(clim.var_abb == "vapor_pressure_deficit", "vpd", clim.var),
         clim.var_lab = if_else(clim.var_abb == "vapor_pressure_deficit", "Vapor pressure deficit", clim.var_lab),
         clim.var_lab2 = if_else(clim.var_abb == "vapor_pressure_deficit", "Vapor pressure deficit\n", clim.var_lab2)) %>% 
  mutate(clim.var_lab = if_else(is.na(clim.var_lab), "All variables\n(only for MFPCA)", clim.var_lab)) %>% 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               'Maximum temperature (°C)'=        'Maximum\ntemperature (°C)', 
                               'Minimum temperature (°C)'=        'Minimum\ntemperature (°C)', 
                               'Evapotransp. ref (mm/day)'=       'Evapotransp.\nref (mm/day)', 
                               'Net solar radiations (MJ/m²)'=    'Net solar\nradiations (MJ/m²)', 
                               'Total precipitations (mm)'=       'Total\nprecipitations (mm)', 
                               'Vapor pressure deficit'=          'Vapor pressure\ndeficit', 
                               'All variables\n(only for MFPCA)'='All variables\n(only for MFPCA)'   )) %>% 
  mutate(clim.var_lab = factor(clim.var_lab, levels = c('Maximum\ntemperature (°C)',
                                                        'Minimum\ntemperature (°C)', 
                                                        'Evapotransp.\nref (mm/day)', 
                                                        'Net solar\nradiations (MJ/m²)', 
                                                        'Total\nprecipitations (mm)', 
                                                        'Vapor pressure\ndeficit', 
                                                        'All variables\n(only for MFPCA)')))

tab_var %>%
  filter(PC == 1) %>% 
  group_by(type_data, country, model) %>% 
  count() %>% 
  spread(key=model, value=n)

save(tab_var, file = paste0(data_path, "tab_var_dim_red.rda"))


stop()

# ---------------------------------------
# Link the performance of the model with the cumulative explained variance 
# SUPPLEMENTARY FIGURES
# > Load table with models performances
load(paste0(data_path, "tab_perf_models.rda"))

tab_var_perf <- tab_var %>%
  # > lab for model and type of data 
  mutate(data_type_abb = recode(data_type, "Monthly climatic predictors"="m", "Daily climatic predictors"="d"),
         model_abb = recode(model, "PCA"="pca", "FPCA"="fpca", "MFPCA"="mfpca", "PLS"="pls", "FPLS"="fpls")) %>%
  # > nb of scores 
  group_by(Country, data_type, clim.var, model) %>% 
  mutate(max_scores = max(PC)) %>% 
  mutate(PC = as.character(PC)) %>% 
  mutate(PC = if_else(PC==max_scores, if_else(max_scores %in% c("1","2","3"), PC, "all"), PC)) %>% 
  filter(PC %in% c("1", "2", "3", "all")) %>% 
  unite(col = "Model", c(model_abb, data_type_abb, PC), sep=".", remove = T) %>%
  # > merge with models' performance
  left_join(., tab_perf_labelled, by = c("data_type"="data_type", 
                                         "Country"="Country", 
                                         "Model" = "Model"))
  
tab_var_perf %>% 
  filter(Outcome == "01_Ya", 
         Country == "WORLD", 
         Type_cv == "01_YEARS", 
         pred_perf == "NSE") %>% 
  ggplot(., aes(x = nb_scores, y = pred_perf_value, 
                color = gpe_model, shape = data_type, linetype=data_type, 
                group = paste(clim.var_lab, model, gpe_model, data_type))) + 
  geom_line() +
  geom_point() + 
  facet_grid(clim.var_lab~model, switch = "y") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        strip.background = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_text(vjust=-70),
        axis.text.x = element_text(angle = 35, hjust=0.9)) +
  scale_color_manual(values = rev(c("#984EA3", "#4DAF4A")), name = "Model family") +
  scale_shape_manual(values = c(3, 20), name = "Type of predictors") + 
  scale_linetype_manual(values = c(2,1), name = "Type of predictors") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1),
         linetype = guide_legend(order = 1, override.aes = list(size=2), ncol = 1),
         shape = guide_legend(order = 1, override.aes = list(size=2), ncol = 1)) + 
  labs(x = "Scores included", y = "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)") 

#-------------------
# > Reconstruct data after data transformation
gridcode_i <- "125.25_50.75"

#-------------------
# PCA - MONTHLY DATA 

# > Objects for reconstruction
list_pca_per_variable <- analyses_world$`Monthly climatic predictors`$pca$list_pca_per_variable
init_data             <- analyses_world$`Monthly climatic predictors`$data_init %>% filter(gridcode == gridcode_i)

# > Reconstruct
tab_month_pca_reconstructed <- list_pca_per_variable %>% 
  map_dfr(., ~{ 
    # > store pca 
    pca_var <- .x$pca
    
    # > inital variables 
    init_vars <- rownames(pca_var$rotation)
    
    # > get means and sds of init data
    X <- init_data %>% 
      dplyr::select(site_year, all_of(init_vars))
    rownames(X) <- X$site_year
    mu = colMeans(X %>% dplyr::select(-site_year))
    
    # > shape table 
    X_dat <- as.data.frame(X) %>% 
      tidyr::gather(key="reconstructed.clim.var", value="init.clim.value", -site_year)
    
    # > recompute data from PCA and loadings
    Xhat_dat <- NULL
    for(i in 1:ncol(pca_var$rotation))
    {
      
      # recompute data 
      Xhat = pca_var$x[which(rownames(pca_var$x) %in% unique(X$site_year)),1:i] %*% t(pca_var$rotation[,1:i])
      Xhat = scale(Xhat, center = -mu, scale = FALSE)
      Xhat0 <- Xhat %>% 
        as.data.frame(.) %>% 
        mutate(site_year = rownames(X)) %>% 
        tidyr::gather(key   = "reconstructed.clim.var", 
                      value = "reconstructed.clim.value", 
                      -site_year) %>% 
        mutate(nComp = i)
      
      # merge with initial data 
      diff <- left_join(X_dat, 
                        Xhat0, 
                        by=c("site_year", "reconstructed.clim.var"))
      
      Xhat_dat <- rbind(Xhat_dat, diff)
      
    }
    
    Xhat_dat
  }, .id = "clim.var")  %>% 
  group_by(site_year, reconstructed.clim.var, nComp) %>% 
  mutate(rel_error_perc = ((abs(reconstructed.clim.value-init.clim.value))/init.clim.value)*100) %>% 
  mutate(rel_error_perc = ifelse(rel_error_perc > 100, 100, rel_error_perc)) %>% 
  mutate(month = substr(reconstructed.clim.var, nchar(reconstructed.clim.var), nchar(reconstructed.clim.var)))

# Reconstruction 
p1 <- tab_month_pca_reconstructed %>% 
  filter(site_year == "125.25_50.75_2000") %>% 
  ggplot(., aes(x=month)) + 
  geom_point(aes(y = init.clim.value)) + 
  geom_line(aes(y = reconstructed.clim.value, color = as.factor(nComp), group=as.factor(nComp))) + 
  facet_wrap(.~clim.var, scales="free", nrow=1) + 
  theme_bw() + 
  theme(legend.position = "bottom")

tab_month_pca_reconstructed %>% 
  filter(site_year == "125.25_50.75_2000") %>% 
  ggplot(., aes(x=month)) + 
  geom_line(aes(y = rel_error_perc, color = as.factor(nComp), group=as.factor(nComp))) +
  facet_wrap(.~clim.var, scales="free", nrow=1) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  lims(y = c(0,105))

#-------------------
# FPCA - MONTHLY DATA 
fpca_var <- analyses_bra$data_month$fpca$list_fpca_per_variable$max_temp

library(fdaACF)
fd_rec_1 <- reconstruct_fd_from_PCA(pca_struct = fpca_var, scores = fpca_var$scores[1,])
fd_rec_2 <- reconstruct_fd_from_PCA(pca_struct = fpca_var, scores = fpca_var$scores[1:2,])
fd_rec_3 <- reconstruct_fd_from_PCA(pca_struct = fpca_var, scores = fpca_var$scores[1:3,])
fd_rec_4 <- reconstruct_fd_from_PCA(pca_struct = fpca_var, scores = fpca_var$scores)

# Reconstruction of the 1rst curve
plot(fd_rec_1, ylim=c(20,35), ylab="Maximum temperature (°C)", xlab="Month") # corresponding points: points(x = 1:7, y=as.vector(predict(fd_rec_1, newdata=1:7)))
lines(fd_rec_2, col="red")
lines(fd_rec_3, col="blue")
lines(fd_rec_4, col="green")
points(x = 1:7, y = as.vector(t(dplyr::select(init_data, starts_with("monthly_max"))[1,])))
legend("topright", 
       legend = c("1 score", "2 scores", "3 scores", "4 scores"),
       col = c("black","red", "blue", "green"),
       cex = 0.75,
       lty = 1)

# Reconstruction of the 1rst curve (ggplot version)
data.frame(Month = 1:7,
           predict(fd_rec_1, newdata=seq(1,7,1))) %>% 
  dplyr::rename("Temp"=2) %>% 
  ggplot(., aes(x=Month, y=Temp)) +
  geom_line() +
  theme_bw()

#-------------------
# explained variance

analyses_bra$data_month$fpca$list_fpca_per_variable %>% 
  map_dfr(., ~{
    
    data.frame(PC = 1:length(.x$values),
               eigenvalue = .x$values) %>% 
      # > total eigen value
      mutate(sum_eig=sum(eigenvalue)) %>% 
      group_by(PC) %>%
      # > compute explained variance for each component
      mutate(variance.percent=eigenvalue/sum_eig) %>% 
      ungroup() %>% 
      # > compute cumulated explained variance
      mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
      dplyr::select(-sum_eig)
        
    
  }, .id="clim.var")

#-------------------
# MFPCA - MONTHLY DATA 

mfpca_var <- analyses_bra$data_month$mfpca

# Using ALL SCORES
par(mfrow=c(1,2))
plot(mfpca_var$list_mfpca$fit[[1]], obs = 1, lty = 2, col = 2, ylim=c(15,35))
points(x = 1:7, y = t(tab_month_bra[which(tab_month_bra$site_year==unique(tab_month_bra$site_year)[1]),
                                    c("monthly_max_2m_temperature_1","monthly_max_2m_temperature_2","monthly_max_2m_temperature_3","monthly_max_2m_temperature_4","monthly_max_2m_temperature_5","monthly_max_2m_temperature_6","monthly_max_2m_temperature_7")]), 
       col = 2)
plot(mfpca_var$list_mfpca$fit[[2]], obs = 1, lty = 2, col = 1, ylim=c(15,35))
points(x = 1:7, y = t(tab_month_bra[which(tab_month_bra$site_year==unique(tab_month_bra$site_year)[1]),
                                    c("monthly_min_2m_temperature_1","monthly_min_2m_temperature_2","monthly_min_2m_temperature_3","monthly_min_2m_temperature_4","monthly_min_2m_temperature_5","monthly_min_2m_temperature_6","monthly_min_2m_temperature_7")]), 
       col = 1)

pred <- predict(object = mfpca_var, scores=mfpca_var$list_mfpca$scores[,1])

multivExpansion <- function(multiFuns, scores)
{
  if(nObs(multiFuns) != NCOL(scores))
    stop("Number of scores does not match number of eigenfunctions.")
  
  # calculate linear combination of multivariate basis functions
  univExp <- foreach::foreach(j = seq_len(length(multiFuns))) %do% { # %do% might require extra loading
    univExpansion(type = "default", 
                  scores = scores,
                  functions = multiFuns[[j]])
  }
  
  # return as multiFunData object
  return(multiFunData(univExp))
}


pred_hm <- mfpca_var$list_mfpca$meanFunction  + 
  multivExpansion(multiFuns = mfpca_var$list_mfpca$functions, scores = mfpca_var$list_mfpca$scores)

as.data.frame(pred[[1]][1])
as.data.frame(pred_hm[[1]][1])


MFPCA::expandBasisFunction(mfpca_var$list_mfpca$scores, argvals = list(1:7), mfpca_var$list_mfpca$functions)

# Only 1 score
mfpca_var$list_mfpca$functions[[1]][1]


# Reconstruction plots
plot(mfpca_var$list_mfpca$fit[[1]], obs = 1, lty = 2, col = "green", ylim=c(20,35), ylab="Maximum temperature (°C)", xlab="Month")

for(i in 1:3)
{
  
  test <- function_mfpca(type_data = "M", 
                         vars_names = vars_names[1:2,],
                         data = tab_month_bra[1:500,],
                         nb_comp = i,
                         argvals = 1:7, 
                         univExpansion = list(list(type = "fda"),
                                              list(type = "fda")))

  plot(test$list_mfpca$fit[[1]], obs = 1, lty = 1, col = c("black","red", "blue", "green")[i], add=T)
  
  
}


points(x = 1:7, y = t(tab_month_bra[which(tab_month_bra$site_year==unique(tab_month_bra$site_year)[1]),
                                    c("monthly_max_2m_temperature_1","monthly_max_2m_temperature_2","monthly_max_2m_temperature_3","monthly_max_2m_temperature_4","monthly_max_2m_temperature_5","monthly_max_2m_temperature_6","monthly_max_2m_temperature_7")]), 
       col = 1)

legend("topright", 
       legend = c("1 score", "2 scores", "3 scores", "4 scores"),
       col = c("black","red", "blue", "green"),
       cex = 0.75,
       lty = 1)

#-------------------
# explained variance 

data.frame(PC = 1:length(analyses_bra$data_month$mfpca$list_mfpca$values),
           eigenvalue = analyses_bra$data_month$mfpca$list_mfpca$values) %>% 
  # > total eigen value
  mutate(sum_eig=sum(eigenvalue)) %>% 
  group_by(PC) %>%
  # > compute explained variance for each component
  mutate(variance.percent=eigenvalue/sum_eig) %>% 
  ungroup() %>% 
  # > compute cumulated explained variance
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
  dplyr::select(-sum_eig)

data.frame(PC = 1:length(analyses_bra$data_day$mfpca$list_mfpca$values),
           eigenvalue = analyses_bra$data_day$mfpca$list_mfpca$values) %>% 
  # > total eigen value
  mutate(sum_eig=sum(eigenvalue)) %>% 
  group_by(PC) %>%
  # > compute explained variance for each component
  mutate(variance.percent=eigenvalue/sum_eig) %>% 
  ungroup() %>% 
  # > compute cumulated explained variance
  mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
  dplyr::select(-sum_eig)

#-------------------
# PLS - MONTHLY DATA 

pls_var <- analyses_bra$data_month$plsr$list_pls_per_variable$max_temp

# ----------------
# explained variance

analyses_bra$data_month$plsr$list_pls_per_variable %>% 
  map_dfr(., ~{
    
    data.frame(PC = 1:length(.x$Xvar),
               eigenvalue = .x$Xvar) %>% 
      # > total eigen value
      mutate(sum_eig=.x$Xtotvar) %>% 
      group_by(PC) %>%
      # > compute explained variance for each component
      mutate(variance.percent=eigenvalue/sum_eig) %>% 
      ungroup() %>% 
      # > compute cumulated explained variance
      mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
      dplyr::select(-sum_eig)
    
    
  }, .id = "clim.var")


#-------------------
# FPLS - MONTHLY DATA 

# ----------------
# explained variance

analyses_bra$data_month$fplsr$list_pls_per_variable %>% 
  map_dfr(., ~{
    
    data.frame(PC = 1:length(.x$Xvar),
               eigenvalue = .x$Xvar) %>% 
      # > total eigen value
      mutate(sum_eig=.x$Xtotvar) %>% 
      group_by(PC) %>%
      # > compute explained variance for each component
      mutate(variance.percent=eigenvalue/sum_eig) %>% 
      ungroup() %>% 
      # > compute cumulated explained variance
      mutate(cumulative.variance.percent = cumsum(variance.percent)) %>% 
      dplyr::select(-sum_eig)
    
    
  }, .id = "clim.var")
