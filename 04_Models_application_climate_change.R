# -------------------------------------------------------------------------
# 
#       4. Using worst and best models to estimate the yield in USA
#       Author: M. Chen, Inrae, 2023  
# 
# -------------------------------------------------------------------------

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
library(fda) ; library(MFPCA)
# > others
library(hydroGOF)

# -------------------------------------------------------------------------
# Home-made functions performing the dimension reductions
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# ----------------------------------
# Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > monthly climatic predictors
load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))

# > daily climatic predictors
load(paste0(data_path, "01_days/list_data_day_usa.rda"))
load(paste0(data_path, "01_days/list_data_day_bra.rda"))

# > Merge datasets for both countries
tab <- list(USA   = list(data = tab_usa,   name="01_USA"),
            BRA   = list(data = tab_bra,   name="02_BRA"))

# > Load FPCA to get the scores
list_fpca_usa <- list(
  "fpca_m" = loadRDa(paste0(data_path, "00_dim_red/fpca_m_usa.rda")),
  "fpca_d" = loadRDa(paste0(data_path, "00_dim_red/fpca_d_usa.rda")))

list_fpca_bra <- list(
  "fpca_m" = loadRDa(paste0(data_path, "00_dim_red/fpca_m_bra.rda")),
  "fpca_d" = loadRDa(paste0(data_path, "00_dim_red/fpca_d_bra.rda")))

# ----------------------------------
# MODELS FORMULA

# > USA models
form.pca.m.3  <- paste0(names(tab_usa %>% dplyr::select(starts_with("PC1_"),  starts_with("PC2_"),  starts_with("PC3_"))  %>% dplyr::select(contains("month"))), collapse = " + ")
form.pca.d.3  <- paste0(names(tab_usa %>% dplyr::select(starts_with("PC1_"),  starts_with("PC2_"),  starts_with("PC3_"))  %>% dplyr::select(contains("day"))),   collapse = " + ")
form.pls.m.3  <- paste0(names(tab_usa %>% dplyr::select(starts_with("PLS1_"), starts_with("PLS2_"), starts_with("PLS3_")) %>% dplyr::select(contains("month"))), collapse = " + ")
form.fpca.m.2 <- paste0(names(tab_usa %>% dplyr::select(starts_with("FPC1_"), starts_with("FPC2_"))                       %>% dplyr::select(contains("month"))), collapse = " + ")
form.pls.d.1  <- paste0(names(tab_usa %>% dplyr::select(starts_with("PLS1_"))                                             %>% dplyr::select(contains("day"))),   collapse = " + ")

# > Means by month
form.avg.m   <- "monthly_max_2m_temperature_1 + monthly_min_2m_temperature_1 + monthly_et0_1 + monthly_surface_net_solar_radiation_1 + monthly_total_precipitation_1 + monthly_vapor_pressure_deficit_1 + 
                 monthly_max_2m_temperature_2 + monthly_min_2m_temperature_2 + monthly_et0_2 + monthly_surface_net_solar_radiation_2 + monthly_total_precipitation_2 + monthly_vapor_pressure_deficit_2 + 
                 monthly_max_2m_temperature_3 + monthly_min_2m_temperature_3 + monthly_et0_3 + monthly_surface_net_solar_radiation_3 + monthly_total_precipitation_3 + monthly_vapor_pressure_deficit_3 + 
                 monthly_max_2m_temperature_4 + monthly_min_2m_temperature_4 + monthly_et0_4 + monthly_surface_net_solar_radiation_4 + monthly_total_precipitation_4 + monthly_vapor_pressure_deficit_4 + 
                 monthly_max_2m_temperature_5 + monthly_min_2m_temperature_5 + monthly_et0_5 + monthly_surface_net_solar_radiation_5 + monthly_total_precipitation_5 + monthly_vapor_pressure_deficit_5 + 
                 monthly_max_2m_temperature_6 + monthly_min_2m_temperature_6 + monthly_et0_6 + monthly_surface_net_solar_radiation_6 + monthly_total_precipitation_6 + monthly_vapor_pressure_deficit_6 + 
                 monthly_max_2m_temperature_7 + monthly_min_2m_temperature_7 + monthly_et0_7 + monthly_surface_net_solar_radiation_7 + monthly_total_precipitation_7 + monthly_vapor_pressure_deficit_7"

# > List of models
list_models_usa <- list(
  #"pca.m.3"  = list(formula = form.pca.m.3),
  #"pca.d.3"  = list(formula = form.pca.d.3),
  "pls.m.3"  = list(formula = form.pls.m.3),
  "pls.d.1"  = list(formula = form.pls.d.1),
  "fpca.m.2" = list(formula = form.fpca.m.2))#,
  #"avg.m"    = list(formula = form.avg.m))

list_models_usa 


# ----------------------------------
# FIT THE MODELS TO COMPARE 

# > USA
list_rf_models_usa <- list_models_usa %>% 
  map(., ~{
    set.seed(101)
    ranger(as.formula(paste0("Ya ~ irrigated_portion + ", .x$formula)),
           data = tab_usa, num.tree=500, importance = "impurity")
  })

list_lm_models_usa <- list_models_usa %>% 
  map(., ~{
    
    lm(as.formula(paste0("Ya ~ irrigated_portion + ", .x$formula)),
       data = tab_usa)
    
  })

# ----------------------------------
# ANALYSE DE SENSIBILITE CHANGEMENT CLIMATIQUE - USA

tab_init <- tab_usa
dim(tab_init)
tab_day_init <- list_data_day_usa

# > Different climate change scenarios 
scenarios <- expand.grid(dT = c(0, 1, 2, 3, 4))

# > Simulate new data 
preds.sensi.usa.cc <- list()

for(i in 1:nrow(scenarios))
{
  
  # > Select the delta for each variable
  dT_i    <- scenarios[i,1]
  
  # > Apply to each dataset
  tab_cc <- tab_init %>%
    mutate(
      # > monthly averages
      monthly_max_2m_temperature_1 = monthly_max_2m_temperature_1+dT_i,
      monthly_max_2m_temperature_2 = monthly_max_2m_temperature_2+dT_i,
      monthly_max_2m_temperature_3 = monthly_max_2m_temperature_3+dT_i,
      monthly_max_2m_temperature_4 = monthly_max_2m_temperature_4+dT_i,
      monthly_max_2m_temperature_5 = monthly_max_2m_temperature_5+dT_i,
      monthly_max_2m_temperature_6 = monthly_max_2m_temperature_6+dT_i,
      monthly_max_2m_temperature_7 = monthly_max_2m_temperature_7+dT_i,
      monthly_min_2m_temperature_1 = monthly_min_2m_temperature_1+dT_i,
      monthly_min_2m_temperature_2 = monthly_min_2m_temperature_2+dT_i,
      monthly_min_2m_temperature_3 = monthly_min_2m_temperature_3+dT_i,
      monthly_min_2m_temperature_4 = monthly_min_2m_temperature_4+dT_i,
      monthly_min_2m_temperature_5 = monthly_min_2m_temperature_5+dT_i,
      monthly_min_2m_temperature_6 = monthly_min_2m_temperature_6+dT_i,
      monthly_min_2m_temperature_7 = monthly_min_2m_temperature_7+dT_i,
      # > annual averages
      year_max_2m_temperature = year_max_2m_temperature+dT_i,
      year_min_2m_temperature = year_min_2m_temperature+dT_i)
  
  # > 
  tab_day_cc <- tab_day_init
  tab_day_cc$max_2m_temperature$cum_clim.value <- tab_day_cc$max_2m_temperature$cum_clim.value + dT_i
  tab_day_cc$min_2m_temperature$cum_clim.value <- tab_day_cc$min_2m_temperature$cum_clim.value + dT_i
    
  # > compute new scores based on simulated data 
  # PCA
  #scores_pca_m <- newscores_pca(type_data  = "M", 
  #                              vars_names = vars_names, 
  #                              init_data  = tab_init %>% dplyr::select(-starts_with("PC")),
  #                              new_data   = tab_cc %>% dplyr::select(-starts_with("PC")))
  #
  #scores_pca_d <- newscores_pca(type_data  = "D", 
  #                              vars_names = vars_names, 
  #                              init_data  = tab_init %>% dplyr::select(-starts_with("PC")),
  #                              new_data   = tab_cc %>% dplyr::select(-starts_with("PC")), 
  #                              init_data_day = tab_day_cc,
  #                              new_data_day  = tab_day_cc)
  # FPCA 
  scores_fpca_m <- newscores_fpca(type_data  = "M", 
                                  vars_names = vars_names, 
                                  init_data  = tab_init %>% dplyr::select(-starts_with("FPC")),
                                  init_fpca  = list_fpca_usa$fpca_m,
                                  new_data   = tab_cc %>% dplyr::select(-starts_with("FPC")), 
                                  data_day   = NULL) 
  # PLSR
  scores_plsr_m <- newscores_plsr2(type_data  = "M", 
                                   vars_names = vars_names, 
                                   init_data  = tab_init %>% dplyr::select(-starts_with("PLS")),
                                   new_data   = tab_cc %>% dplyr::select(-starts_with("PLS")),
                                   outcome = "Ya")
  
  scores_plsr_d <- newscores_plsr2(type_data     = "D", 
                                   vars_names    = vars_names, 
                                   init_data     = tab_init %>% dplyr::select(-starts_with("PLS")),
                                   new_data      = tab_cc %>% dplyr::select(-starts_with("PLS")), 
                                   init_data_day = tab_day_init,
                                   new_data_day  = tab_day_cc,
                                   outcome       = "Ya")
  
  # > Add to initial datasets
  tab_cc_new <- tab_cc %>% 
    ungroup(.) %>% 
    dplyr::select(site_year, Ya, irrigated_portion, starts_with("monthly_"), starts_with("year")) %>%
    dplyr::select(-starts_with("monthly_cum")) %>%
    cbind(., 
          #scores_pca_m$new_data  %>% dplyr::select(starts_with("PC")), 
          #scores_pca_d$new_data  %>% dplyr::select(starts_with("PC")), 
          scores_fpca_m$new_data %>% dplyr::select(starts_with("FPC")), 
          scores_plsr_m$new_data %>% dplyr::select(starts_with("PLS")), 
          scores_plsr_d$new_data %>% dplyr::select(starts_with("PLS"))
          ) 
  
  # > Use the model on the current data to predict yield using the simulated data
  list_rf_preds_usa <- list_rf_models_usa %>% 
    map_dfr(., ~{
      
      seed <-101
      preds_rf <- predict(.x, 
                          data = tab_cc_new, 
                          type = "response", 
                          seed = seed, 
                          num.trees = 500)
              
      data.frame(site_year  = tab_cc_new$site_year,
                 obs        = tab_cc_new$Ya,
                 pred       = preds_rf$predictions,
                 dT         = dT_i,
                 country    = "01_USA", 
                 gpe_model  = "Random forest")
      
    }, .id="model")
  
  list_lm_preds_usa <- list_lm_models_usa %>% 
    map_dfr(., ~{
      
      preds_lm <- predict(.x,
                          newdata = tab_cc_new, 
                          type = "response")
      
      data.frame(site_year  = tab_cc_new$site_year,
                 obs        = tab_cc_new$Ya,
                 pred       = as.numeric(as.character(preds_lm)),
                 dT         = dT_i,
                 country    = "01_USA", 
                 gpe_model  = "Multiple linear regression")
      
    }, .id="model")
  
  # > Store
  preds.sensi.usa.cc[[paste0("scenario_", i)]] <- list(
    tab_cc_new = tab_cc_new,
    preds = rbind(list_rf_preds_usa, list_lm_preds_usa))
  
}

# > Save
save(preds.sensi.usa.cc, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/preds_cc_usa_new.rda")


rm(tab_init, tab_day_init, tab_cc, tab_cc_new, preds.sensi.usa.cc)
#----------------------------------------
# BRAZIL

# > Models formula for Brazil
form.pca.m.2   <- paste0(names(tab_bra %>% dplyr::select(starts_with("PC1_"), starts_with("PC2_"))                       %>% dplyr::select(contains("month"))), collapse = " + ")
form.pca.d.all <- paste0(names(tab_bra %>% dplyr::select(starts_with("PC"))                                              %>% dplyr::select(contains("day"))),   collapse = " + ")
form.fpca.m.2  <- paste0(names(tab_bra %>% dplyr::select(starts_with("FPC1_"), starts_with("FPC2_"))                     %>% dplyr::select(contains("month"))), collapse = " + ")
form.pls.m.3   <- paste0(names(tab_bra %>% dplyr::select(starts_with("PLS1"),  starts_with("PLS2"), starts_with("PLS3")) %>% dplyr::select(contains("month"))), collapse = " + ")
form.pls.d.all <- paste0(names(tab_bra %>% dplyr::select(starts_with("PLS"))                                             %>% dplyr::select(contains("day"))),   collapse = " + ")

# > Means by month
form.avg.m   <- "monthly_max_2m_temperature_1 + monthly_min_2m_temperature_1 + monthly_et0_1 + monthly_surface_net_solar_radiation_1 + monthly_total_precipitation_1 + monthly_vapor_pressure_deficit_1 + 
                 monthly_max_2m_temperature_2 + monthly_min_2m_temperature_2 + monthly_et0_2 + monthly_surface_net_solar_radiation_2 + monthly_total_precipitation_2 + monthly_vapor_pressure_deficit_2 + 
                 monthly_max_2m_temperature_3 + monthly_min_2m_temperature_3 + monthly_et0_3 + monthly_surface_net_solar_radiation_3 + monthly_total_precipitation_3 + monthly_vapor_pressure_deficit_3 + 
                 monthly_max_2m_temperature_4 + monthly_min_2m_temperature_4 + monthly_et0_4 + monthly_surface_net_solar_radiation_4 + monthly_total_precipitation_4 + monthly_vapor_pressure_deficit_4 + 
                 monthly_max_2m_temperature_5 + monthly_min_2m_temperature_5 + monthly_et0_5 + monthly_surface_net_solar_radiation_5 + monthly_total_precipitation_5 + monthly_vapor_pressure_deficit_5 + 
                 monthly_max_2m_temperature_6 + monthly_min_2m_temperature_6 + monthly_et0_6 + monthly_surface_net_solar_radiation_6 + monthly_total_precipitation_6 + monthly_vapor_pressure_deficit_6 + 
                 monthly_max_2m_temperature_7 + monthly_min_2m_temperature_7 + monthly_et0_7 + monthly_surface_net_solar_radiation_7 + monthly_total_precipitation_7 + monthly_vapor_pressure_deficit_7"

list_models_bra <- list(
  "pca.m.2"    = list(formula = form.pca.m.2),
  "pca.d.all"  = list(formula = form.pca.d.all),
  "pls.m.3"  = list(formula = form.pls.m.all),
  "fpca.m.2"   = list(formula = form.fpca.m.2),
  "pls.d.all"  = list(formula = form.pls.d.all),
  "avg.m"      = list(formula = form.avg.m))

list_models_bra

# > Fit models on the full Brazilian data 
list_rf_models_bra <- list_models_bra %>% 
  map(., ~{
    set.seed(101)
    ranger(as.formula(paste0("Ya ~ irrigated_portion + ", .x$formula)),
           data = tab_bra, num.tree=500, importance = "impurity")
  })

list_lm_models_bra <- list_models_bra %>% 
  map(., ~{
    
    lm(as.formula(paste0("Ya ~ irrigated_portion + ", .x$formula)),
       data = tab_bra)
    
  })

# > Simulate yield in different temperature increasing scenarios  
tab_init <- tab_bra
dim(tab_init)
tab_day_init <- list_data_day_bra %>% 
  # > necessary to have 212 non na values 
  # (bissextile years have 213 observations with 29th February)
  map(., ~{ 
    .x %>% filter(day_of_year < 213)
  })

# > Different climate change scenarios 
scenarios <- expand.grid(dT = c(0, 1, 2, 3, 4))

# > Simulate new data 
preds.sensi.bra.cc <- list()

for(i in 1:nrow(scenarios))
{
  
  # > Select the delta for each variable
  dT_i    <- scenarios[i,1]
  
  # > Apply to each dataset
  tab_cc <- tab_init %>%
    mutate(
      # > monthly averages
      monthly_max_2m_temperature_1 = monthly_max_2m_temperature_1+dT_i,
      monthly_max_2m_temperature_2 = monthly_max_2m_temperature_2+dT_i,
      monthly_max_2m_temperature_3 = monthly_max_2m_temperature_3+dT_i,
      monthly_max_2m_temperature_4 = monthly_max_2m_temperature_4+dT_i,
      monthly_max_2m_temperature_5 = monthly_max_2m_temperature_5+dT_i,
      monthly_max_2m_temperature_6 = monthly_max_2m_temperature_6+dT_i,
      monthly_max_2m_temperature_7 = monthly_max_2m_temperature_7+dT_i,
      monthly_min_2m_temperature_1 = monthly_min_2m_temperature_1+dT_i,
      monthly_min_2m_temperature_2 = monthly_min_2m_temperature_2+dT_i,
      monthly_min_2m_temperature_3 = monthly_min_2m_temperature_3+dT_i,
      monthly_min_2m_temperature_4 = monthly_min_2m_temperature_4+dT_i,
      monthly_min_2m_temperature_5 = monthly_min_2m_temperature_5+dT_i,
      monthly_min_2m_temperature_6 = monthly_min_2m_temperature_6+dT_i,
      monthly_min_2m_temperature_7 = monthly_min_2m_temperature_7+dT_i,
      # > annual averages
      year_max_2m_temperature = year_max_2m_temperature+dT_i,
      year_min_2m_temperature = year_min_2m_temperature+dT_i)
  
  # > Add degrees to daily data 
  tab_day_cc <- tab_day_init
  tab_day_cc$max_2m_temperature$cum_clim.value <- tab_day_cc$max_2m_temperature$cum_clim.value + dT_i
  tab_day_cc$min_2m_temperature$cum_clim.value <- tab_day_cc$min_2m_temperature$cum_clim.value + dT_i
  
  # > compute new scores based on simulated data 
  # PCA
  scores_pca_m <- newscores_pca(type_data  = "M", 
                                vars_names = vars_names, 
                                init_data  = tab_init %>% dplyr::select(-starts_with("PC")),
                                new_data   = tab_cc %>% dplyr::select(-starts_with("PC")))
  
  scores_pca_d <- newscores_pca(type_data  = "D", 
                                vars_names = vars_names, 
                                init_data  = tab_init %>% dplyr::select(-starts_with("PC")),
                                new_data   = tab_cc %>% dplyr::select(-starts_with("PC")), 
                                init_data_day = tab_day_cc,
                                new_data_day  = tab_day_cc)
  # FPCA 
  scores_fpca_m <- newscores_fpca(type_data  = "M", 
                                  vars_names = vars_names, 
                                  init_data  = tab_init %>% dplyr::select(-starts_with("FPC")),
                                  init_fpca  = list_fpca_usa$fpca_m,
                                  new_data   = tab_cc %>% dplyr::select(-starts_with("FPC")), 
                                  data_day   = NULL) 
  # PLSR
  scores_plsr_m <- newscores_plsr2(type_data  = "M", 
                                   vars_names = vars_names, 
                                   init_data  = tab_init %>% dplyr::select(-starts_with("PLS")),
                                   new_data   = tab_cc %>% dplyr::select(-starts_with("PLS")),
                                   outcome = "Ya")
  
  scores_plsr_d <- newscores_plsr2(type_data     = "D", 
                                   vars_names    = vars_names, 
                                   init_data     = tab_init %>% dplyr::select(-starts_with("PLS")),
                                   new_data      = tab_cc %>% dplyr::select(-starts_with("PLS")), 
                                   init_data_day = tab_day_init,
                                   new_data_day  = tab_day_cc,
                                   outcome       = "Ya")
  
  # > Add to initial datasets
  tab_cc_new <- tab_cc %>% 
    ungroup(.) %>% 
    dplyr::select(site_year, Ya, irrigated_portion, starts_with("monthly_"), starts_with("year")) %>%
    dplyr::select(-starts_with("monthly_cum")) %>%
    cbind(., 
          scores_pca_m$new_data  %>% dplyr::select(starts_with("PC")), 
          scores_pca_d$new_data  %>% dplyr::select(starts_with("PC")), 
          scores_fpca_m$new_data %>% dplyr::select(starts_with("FPC")), 
          scores_plsr_m$new_data %>% dplyr::select(starts_with("PLS")), 
          scores_plsr_d$new_data %>% dplyr::select(starts_with("PLS"))
    ) 
  
  # > Use the model on the current data to predict yield using the simulated data
  list_rf_preds_bra <- list_rf_models_bra %>% 
    map_dfr(., ~{
      
      seed <-101
      preds_rf <- predict(.x, 
                          data = tab_cc_new, 
                          type = "response", 
                          seed = seed, 
                          num.trees = 500)
      
      data.frame(site_year  = tab_cc_new$site_year,
                 obs        = tab_cc_new$Ya,
                 pred       = preds_rf$predictions,
                 dT         = dT_i,
                 country    = "02_BRA", 
                 gpe_model  = "Random forest")
      
    }, .id="model")
  
  list_lm_preds_bra <- list_lm_models_bra %>% 
    map_dfr(., ~{
      
      preds_lm <- predict(.x,
                          newdata = tab_cc_new, 
                          type = "response")
      
      data.frame(site_year  = tab_cc_new$site_year,
                 obs        = tab_cc_new$Ya,
                 pred       = as.numeric(as.character(preds_lm)),
                 dT         = dT_i,
                 country    = "02_BRA", 
                 gpe_model  = "Multiple linear regression")
      
    }, .id="model")
  
  # > Store
  preds.sensi.bra.cc[[paste0("scenario_", i)]] <- list(
    tab_cc_new = tab_cc_new,
    preds = rbind(list_rf_preds_bra, list_lm_preds_bra))
  
}

# > Save
save(preds.sensi.bra.cc, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/preds_cc_bra.rda")

rm(tab_init, tab_day_init, tab_cc, tab_cc_new, preds.sensi.usa.cc)


stop()

# Loading predictions 
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/preds_cc_usa.rda")
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/preds_cc_bra.rda")

# Loading data with country
load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))

# Merge and remove the Desert points
pred.usa <- preds.sensi.usa.cc %>% 
  map_dfr(., ~{.x$preds}) %>% 
  left_join(., tab_usa %>% dplyr::select(site_year, country_name)) %>%
  filter(country_name != "Desert")

pred.bra <- preds.sensi.bra.cc %>% 
  map_dfr(., ~{.x$preds}) %>% 
  left_join(., tab_bra %>% dplyr::select(site_year, country_name)) %>%
  filter(country_name != "Desert")

preds.cc <- rbind(pred.usa, 
                  pred.bra) %>% 
  # cc scenario labels
  mutate(dT = recode(dT, "0"="+0°C", "1"="+1°C", "2"="+2°C", "3"="+3°C", "4"="+4°C"),
         dT = factor(dT, levels = c("+0°C", "+1°C", "+2°C", "+3°C", "+4°C"))) %>% 
  separate(site_year, into = c("x", "y", "year"), sep = "_") %>% 
  # compute median prediction
  group_by(country, country_name, gpe_model, model, dT, x, y) %>% 
  summarise(Median_Prediction = median(pred)) %>% 
  pivot_wider(names_from = "dT", values_from = Median_Prediction) %>% 
  rename("Ref"="+0°C") %>% 
  pivot_longer(cols = starts_with("+"), names_to = "dT", values_to = "Median_Prediction") %>% 
  # relative difference between scenario with +0°C (initial) and scenarios with increase in temperatures
  group_by(country, country_name, gpe_model, model, dT, x, y) %>% 
  mutate(Difference_Prediction_percentage = 100*((Median_Prediction-Ref)/Ref)) %>%
  # models labels
  mutate(gpe_model = factor(gpe_model, levels = c("Random forest", "Multiple linear regression"))) %>%
  #mutate(model     = factor(model, levels = c("pca.m.3",  "pca.m.2",  "pca.d.3", "pca.d.all", "pls.m.3", "fpca.m.2","pls.d.1",  "pls.d.all", "avg.m"))) %>%
  mutate(model     = factor(model, levels = c('pca.m.3','pca.m.2',
                                              'pca.d.3','pca.d.all', 
                                              'fpca.m.3','fpca.m.2',
                                              'pls.m.2','pls.m.all', 
                                              'pls.d.3','pls.d.all', 'avg.m'))) %>%
  mutate(selection_step = case_when(
    model %in% c("pca.m.3",  "pca.m.2") ~ 'Reference approach', 
    model %in% c("pca.d.3", "pca.d.all") ~ 'Variant varying in\ntemporal resolution\nof climate data',
    model %in% c("pls.m.3", "pls.m.2", "pls.m.all",
                 "fpca.m.2", "fpca.m.3") ~ "Variants varying in\ndimension reduction\ntechniques",
    model %in% c("pls.d.1",  "pls.d.3" , "pls.d.all") ~ "Variant varying in\ndimension reduction\ntechniques and\ntemporal resolution",
    model == "avg.m" ~ 'Model based on\nmonthly averages'),
    selection_step = factor(selection_step, levels = c('Reference approach', 
                                                       'Variant varying in\ntemporal resolution\nof climate data',
                                                       "Variants varying in\ndimension reduction\ntechniques",
                                                       "Variant varying in\ndimension reduction\ntechniques and\ntemporal resolution",
                                                       "Model based on\nmonthly averages"))) %>%
  mutate(x = as.numeric(as.character(x)), 
         y = as.numeric(as.character(y)))

save(preds.cc, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/tab_preds_cc.rda")

stop()

#-----------------------------------------
# STOP STOP STOP STOP 
# OLD VERSION OF THE CODE
# > BRAZIL
set.seed(101)
rf_pca.m2_bra <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.2)),
                        data = tab_bra, num.tree=500, importance = "impurity")

set.seed(101)
rf_pls.mall_bra <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.pls.m.all)),
                        data = tab_bra, num.tree=500, importance = "impurity")



# MLR
# > BRAZIL
lm_pca.m2_bra <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.2)),
                data = tab_bra)

lm_pls.mall_bra <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.pls.m.all)),
                    data = tab_bra)

# > USA
lm_pca.m3_usa <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.3)),
                data = tab_usa)

lm_pls.m2_usa <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.pls.m.2)),
                    data = tab_usa)

# MODEL INCLUDING THE MONTHLY AVERAGES 
# RF
# > BRAZIL
set.seed(101)
rf_avg.m_bra <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.avg.m)),
                   data = tab_bra, num.tree=500, importance = "impurity")
# > USA
set.seed(101)
rf_avg.m_usa <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.avg.m)),
                   data = tab_usa, num.tree=500, importance = "impurity")
# MLR
# > BRAZIL 
lm_avg.m_bra <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.avg.m)),
                data = tab_bra)
# > USA
lm_avg.m_usa <- lm(as.formula(paste0("Ya ~ irrigated_portion + ", form.avg.m)),
                data = tab_usa)

# ----------------------------------

# > Load PCA to get the scores
list_pca <- list(
  "01_USA" = loadRDa(paste0(data_path, "00_dim_red/pca_m_usa.rda")),
  "02_BRA" = loadRDa(paste0(data_path, "00_dim_red/pca_m_bra.rda")))

list_plsr <- list(
  "01_USA" = loadRDa(paste0(data_path, "00_dim_red/plsr_m_usa.rda")),
  "02_BRA" = loadRDa(paste0(data_path, "00_dim_red/plsr_m_bra.rda")))

# ----------------------------------
# ANALYSE DE SENSIBILITE CHANGEMENT CLIMATIQUE 

# > Different climate change scenarios 
scenarios <- expand.grid(dT = c(0, 1, 2, 3, 4))

# > Simulate new data 
preds.sensi.usa.cc <- list()
preds.sensi.bra.cc <- list()
 
for(i in 1:nrow(scenarios))
{
  
  # > Select the delta for each variable
  dT_i    <- scenarios[i,1]
  
  # > Apply to each dataset
  # USA
  tab_usa_cc <- tab_usa %>%
    #filter(country_name=="United States of America") %>% 
    mutate(
      # > monthly averages
      monthly_max_2m_temperature_1 = monthly_max_2m_temperature_1+dT_i,
      monthly_max_2m_temperature_2 = monthly_max_2m_temperature_2+dT_i,
      monthly_max_2m_temperature_3 = monthly_max_2m_temperature_3+dT_i,
      monthly_max_2m_temperature_4 = monthly_max_2m_temperature_4+dT_i,
      monthly_max_2m_temperature_5 = monthly_max_2m_temperature_5+dT_i,
      monthly_max_2m_temperature_6 = monthly_max_2m_temperature_6+dT_i,
      monthly_max_2m_temperature_7 = monthly_max_2m_temperature_7+dT_i,
      monthly_min_2m_temperature_1 = monthly_min_2m_temperature_1+dT_i,
      monthly_min_2m_temperature_2 = monthly_min_2m_temperature_2+dT_i,
      monthly_min_2m_temperature_3 = monthly_min_2m_temperature_3+dT_i,
      monthly_min_2m_temperature_4 = monthly_min_2m_temperature_4+dT_i,
      monthly_min_2m_temperature_5 = monthly_min_2m_temperature_5+dT_i,
      monthly_min_2m_temperature_6 = monthly_min_2m_temperature_6+dT_i,
      monthly_min_2m_temperature_7 = monthly_min_2m_temperature_7+dT_i,
      # > annual averages
      year_max_2m_temperature = year_max_2m_temperature+dT_i,
      year_min_2m_temperature = year_min_2m_temperature+dT_i)
  
  # BRA
  tab_bra_cc <- tab_bra %>% 
    mutate(
      # > monthly averages
      monthly_max_2m_temperature_1 = monthly_max_2m_temperature_1+dT_i,
      monthly_max_2m_temperature_2 = monthly_max_2m_temperature_2+dT_i,
      monthly_max_2m_temperature_3 = monthly_max_2m_temperature_3+dT_i,
      monthly_max_2m_temperature_4 = monthly_max_2m_temperature_4+dT_i,
      monthly_max_2m_temperature_5 = monthly_max_2m_temperature_5+dT_i,
      monthly_max_2m_temperature_6 = monthly_max_2m_temperature_6+dT_i,
      monthly_max_2m_temperature_7 = monthly_max_2m_temperature_7+dT_i,
      monthly_min_2m_temperature_1 = monthly_min_2m_temperature_1+dT_i,
      monthly_min_2m_temperature_2 = monthly_min_2m_temperature_2+dT_i,
      monthly_min_2m_temperature_3 = monthly_min_2m_temperature_3+dT_i,
      monthly_min_2m_temperature_4 = monthly_min_2m_temperature_4+dT_i,
      monthly_min_2m_temperature_5 = monthly_min_2m_temperature_5+dT_i,
      monthly_min_2m_temperature_6 = monthly_min_2m_temperature_6+dT_i,
      monthly_min_2m_temperature_7 = monthly_min_2m_temperature_7+dT_i,
      # > annual averages
      year_max_2m_temperature = year_max_2m_temperature+dT_i,
      year_min_2m_temperature = year_min_2m_temperature+dT_i)
  
  # > list to store the results
  list_scores_usa <- list()
  list_scores_bra <- list()
  list_scores_usa_pls <- list()
  list_scores_bra_pls <- list()
  
  # > Apply PCA loads on new data 
  for(var_j in unique(vars_names$clim.var)) 
  {
    
    clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
    
    # > PCA loads
    loads_usa <- list_pca[["01_USA"]]$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
    loads_bra <- list_pca[["02_BRA"]]$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
    
    # > PLSR 
    plsr_usa <- list_plsr[["01_USA"]]$list_pls_per_variable[[paste0(clim.var_abb_j)]]
    plsr_bra <- list_plsr[["02_BRA"]]$list_pls_per_variable[[paste0(clim.var_abb_j)]]
    
    # > Mean and sd of the original data to standardize 
    mu_usa <- tab_usa %>% ungroup(.) %>%
      dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
      colMeans(.)
    mu_bra <- tab_bra %>% ungroup(.) %>%
      dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
      colMeans(.)
    
    sd_usa <- tab_usa %>% ungroup(.) %>%
      dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
      apply(., 2, sd)
    sd_bra <- tab_bra %>% ungroup(.) %>%
      dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
      apply(., 2, sd)
    
    # > New data 
    tab_bra_cc_j <- tab_bra_cc %>% ungroup(.) %>%
      dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
    tab_usa_cc_j <- tab_usa_cc %>% ungroup(.) %>%
      dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
    
    # > Scale new data based on mean and sd of the original data 
    z_tab_bra_cc_j <- scale(tab_bra_cc_j[,-1], center = mu_bra, scale=sd_bra)
    z_tab_usa_cc_j <- scale(tab_usa_cc_j[,-1], center = mu_usa, scale=sd_usa)
    
    scores_tab_bra_cc_j <- list()
    scores_tab_usa_cc_j <- list()
    
    # > Compute PCA scores bra
    scores_tab_bra_cc_j <- list()
    for(k in 1:ncol(loads_bra))
    {
      
      scores_tab_bra_cc_j[[paste0("PC", k, "_month_")]] <- 
        data.frame(site_year = tab_bra_cc_j$site_year,
                   score = sapply(1:ncol(z_tab_bra_cc_j), function(x) z_tab_bra_cc_j[,x] * loads_bra[x,k] ) %>% apply(., 1, sum))
      
    }
    
    # > Compute PCA scores usa
    for(l in 1:ncol(loads_usa))
    {
      
      scores_tab_usa_cc_j[[paste0("PC", l, "_month_")]] <- 
        data.frame(site_year = tab_usa_cc_j$site_year,
                   score = sapply(1:ncol(z_tab_usa_cc_j), function(x) z_tab_usa_cc_j[,x] * loads_usa[x,l] ) %>% apply(., 1, sum))
      
    }
    
    # > PLSR score from newdata
    tab_scores_usa_plsr_j <- predict(plsr_usa, type="scores", newdata=z_tab_usa_cc_j)
    tab_scores_bra_plsr_j <- predict(plsr_bra, type="scores", newdata=z_tab_bra_cc_j)
    
    colnames(tab_scores_usa_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_usa_plsr_j),"_month_",clim.var_abb_j)
    colnames(tab_scores_bra_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_bra_plsr_j),"_month_",clim.var_abb_j)
    
    # > Store
    list_scores_usa[[paste0(clim.var_abb_j)]] <- plyr::ldply(scores_tab_usa_cc_j, data.frame, .id = "PC_id") 
    list_scores_bra[[paste0(clim.var_abb_j)]] <- plyr::ldply(scores_tab_bra_cc_j, data.frame, .id = "PC_id") 

    list_scores_usa_pls[[paste0(clim.var_abb_j)]] <- tab_scores_usa_plsr_j
    list_scores_bra_pls[[paste0(clim.var_abb_j)]] <- tab_scores_bra_plsr_j
    
  }
  
  # > PC scores in table
  tab_scores_usa <- plyr::ldply(list_scores_usa, data.frame, .id="var_name") %>% 
    unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
    pivot_wider(names_from = c("PC_name"), values_from = "score")
  tab_scores_bra <- plyr::ldply(list_scores_bra, data.frame, .id="var_name") %>% 
    unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
    pivot_wider(names_from = c("PC_name"), values_from = "score")
  
  tab_scores_usa_plsr <- map_dfc(list_scores_usa_pls, data.frame)
  tab_scores_bra_plsr <- map_dfc(list_scores_bra_pls, data.frame)
  
  # > Add to initial datasets
  tab_usa_cc_new <- tab_usa_cc %>% 
    ungroup(.) %>% 
    dplyr::select(site_year, Ya, irrigated_portion, starts_with("monthly_"), starts_with("year")) %>%
    dplyr::select(-starts_with("monthly_cum")) %>%
    left_join(tab_scores_usa, by = "site_year") %>% 
    cbind(., tab_scores_usa_plsr)
  
  tab_bra_cc_new <- tab_bra_cc %>% 
    ungroup(.) %>% 
    dplyr::select(site_year, Ya, irrigated_portion, starts_with("monthly_"), starts_with("year")) %>%
    dplyr::select(-starts_with("monthly_cum")) %>%
    left_join(tab_scores_bra, by = "site_year") %>% 
    cbind(., tab_scores_bra_plsr)
  
  # > Use the model on the data to predict the simulated data
  seed <-101
  pred.rf_pca.m3_usa <- predict(rf_pca.m3_usa, 
                               data = tab_usa_cc_new, 
                               type = "response", 
                               seed = seed, 
                               num.trees = 500)
  
  seed <-101
  pred.rf_pls.m2_usa <- predict(rf_pls.m2_usa, 
                                data = tab_usa_cc_new, 
                                type = "response", 
                                seed = seed, 
                                num.trees = 500)
  
  seed <-101
  pred.rf_pca.m2_bra <- predict(rf_pca.m2_bra, 
                           data = tab_bra_cc_new, 
                           type = "response", 
                           seed = seed, 
                           num.trees = 500)
  seed <-101
  pred.rf_pls.mall_bra <- predict(rf_pls.mall_bra, 
                                data = tab_bra_cc_new, 
                                type = "response", 
                                seed = seed, 
                                num.trees = 500)
  seed <-101
  pred.rf_avg.m_usa <- predict(rf_avg.m_usa, 
                                data = tab_usa_cc_new, 
                                type = "response", 
                                seed = seed, 
                                num.trees = 500)
  seed <-101
  pred.rf_avg.m_bra <- predict(rf_avg.m_bra, 
                                data = tab_bra_cc_new, 
                                type = "response", 
                                seed = seed, 
                                num.trees = 500)
  
  pred.lm_pca.m3_usa <- predict(lm_pca.m3_usa,
                            newdata = tab_usa_cc_new, 
                            type = "response")
  
  pred.lm_pls.m2_usa <- predict(lm_pls.m2_usa,
                                newdata = tab_usa_cc_new, 
                                type = "response")
  
  pred.lm_pca.m2_bra <- predict(lm_pca.m2_bra,
                            newdata = tab_bra_cc_new, 
                            type = "response")
  pred.lm_pls.mall_bra <- predict(lm_pls.mall_bra,
                                newdata = tab_bra_cc_new, 
                                type = "response")
  
  pred.lm_avg.m_usa <- predict(lm_avg.m_usa,
                                newdata = tab_usa_cc_new, 
                                type = "response")
  pred.lm_avg.m_bra <- predict(lm_avg.m_bra,
                                newdata = tab_bra_cc_new, 
                                type = "response")
  
  # > Store
  preds.sensi.usa.cc[[paste0("scenario_", i)]] <- data.frame(site_year            = tab_usa_cc_new$site_year,
                                                             obs                  = tab_usa_cc_new$Ya,
                                                             pred.rf_pca.m3_usa   = pred.rf_pca.m3_usa$predictions,
                                                             pred.rf_pls.m2_usa   = pred.rf_pls.m2_usa$predictions,
                                                             pred.lm_pca.m3_usa   = as.numeric(as.character(pred.lm_pca.m3_usa)),
                                                             pred.lm_pls.m2_usa   = as.numeric(as.character(pred.lm_pls.m2_usa)),
                                                             pred.rf_avg.m_usa    = pred.rf_avg.m_usa$predictions,
                                                             pred.lm_avg.m_usa    = as.numeric(as.character(pred.lm_avg.m_usa)),
                                                             dT        = dT_i,
                                                             country = "01_USA")
  
  preds.sensi.bra.cc[[paste0("scenario_", i)]] <- data.frame(site_year            = tab_bra_cc_new$site_year,
                                                             obs                  = tab_bra_cc_new$Ya,
                                                             pred.rf_pca.m2_bra   = pred.rf_pca.m2_bra$predictions,
                                                             pred.rf_pls.mall_bra = pred.rf_pls.mall_bra$predictions,
                                                             pred.lm_pca.m2_bra   = as.numeric(as.character(pred.lm_pca.m2_bra)),
                                                             pred.lm_pls.mall_bra = as.numeric(as.character(pred.lm_pls.mall_bra)),
                                                             pred.rf_avg.m_bra    = pred.rf_avg.m_bra$predictions,
                                                             pred.lm_avg.m_bra    = as.numeric(as.character(pred.lm_avg.m_bra)),
                                                             dT        = dT_i,
                                                             country = "02_BRA")


}

# > Format tables for both countries
# and compute median yields over years
tab.preds.sensi.usa.cc <- plyr::ldply(preds.sensi.usa.cc, data.frame) %>% 
  left_join(., tab_usa %>% dplyr::select(site_year, country_name)) %>% 
  filter(country_name != "Desert") %>% 
  # cc scenario labels
  mutate(dT = recode(dT, "0"="+0°C", "1"="+1°C", "2"="+2°C", "3"="+3°C", "4"="+4°C"),
         dT = factor(dT, levels = c("+0°C", "+1°C", "+2°C", "+3°C", "+4°C"))) %>% 
  separate(site_year, into = c("x", "y", "year"), sep = "_") %>% 
  # set in long format
  pivot_longer(cols = starts_with("pred"), names_to = "Model", values_to = "Prediction") %>% 
  # compute median prediction
  group_by(x, y, dT, Model) %>% 
  summarise(Median_Prediction = median(Prediction))  %>% 
  # relative difference between scenario with +0°C (initial) and scenarios with increase in temperatures
  split(.$Model) %>% 
  map_dfr(., ~{
    
    .x %>% 
      pivot_wider(names_from = "dT", values_from = Median_Prediction) %>% 
      rename("Ref"="+0°C") %>% 
      pivot_longer(cols = starts_with("+"), names_to = "dT", values_to = "Median_Prediction") %>% 
      group_by(x, y, dT) %>% 
      mutate(Difference_Prediction_percentage = 100*((Median_Prediction-Ref)/Ref))
          
        }, .id = "Dim_red_Model") %>% 
  # models labels
  mutate(Model = recode(Dim_red_Model, 
                        "pred.rf_pca.m3_usa"="Random forest", 
                        "pred.lm_pca.m3_usa"="Multiple linear regression", 
                        "pred.rf_pls.m2_usa"="Random forest", 
                        "pred.lm_pls.m2_usa"="Multiple linear regression",
                        "pred.rf_avg.m_usa"="Random forest",
                        "pred.lm_avg.m_usa"="Multiple linear regression")) %>% 
  mutate(Model = factor(Model, levels = c("Random forest", 
                                          "Multiple linear regression"))) %>%
  # models labels
  mutate(Dim_red = recode(Dim_red_Model, 
                          "pred.rf_pca.m3_usa"="pca.m.3", 
                          "pred.lm_pca.m3_usa"="pca.m.3",
                          "pred.rf_pls.m2_usa"="pls.m.2", 
                          "pred.lm_pls.m2_usa"="pls.m.2", 
                          "pred.rf_avg.m_usa"="avg.m",
                          "pred.lm_avg.m_usa"="avg.m")) %>% 
  mutate(Dim_red = factor(Dim_red, levels = c("pca.m.3", "pls.m.2",
                                              "avg.m"))) %>%
  mutate(x = as.numeric(as.character(x)), 
         y = as.numeric(as.character(y)))

tab.preds.sensi.bra.cc <- plyr::ldply(preds.sensi.bra.cc, data.frame) %>% 
  left_join(., tab_bra %>% dplyr::select(site_year, country_name)) %>% 
  filter(country_name != "Desert") %>% 
  mutate(dT = recode(dT, "0"="+0°C", "1"="+1°C", "2"="+2°C", "3"="+3°C", "4"="+4°C"),
         dT = factor(dT, levels = c("+0°C", "+1°C", "+2°C", "+3°C", "+4°C"))) %>% 
  separate(site_year, into = c("x", "y", "year"), sep = "_") %>% 
  # set in long format
  pivot_longer(cols = starts_with("pred"), names_to = "Model", values_to = "Prediction") %>% 
  # compute median prediction
  group_by(x, y, dT, Model) %>% 
  summarise(Median_Prediction = median(Prediction))  %>% 
  # relative difference between scenario with +0°C (initial) and scenarios with increase in temperatures
  split(.$Model) %>% 
  map_dfr(., ~{
    
    .x %>% 
      pivot_wider(names_from = "dT", values_from = Median_Prediction) %>% 
      rename("Ref"="+0°C") %>% 
      pivot_longer(cols = starts_with("+"), names_to = "dT", values_to = "Median_Prediction") %>% 
      group_by(x, y, dT) %>% 
      mutate(Difference_Prediction_percentage = 100*((Median_Prediction-Ref)/Ref))
    
  }, .id = "Dim_red_Model") %>% 
  # models labels
  mutate(Model = recode(Dim_red_Model, 
                        "pred.rf_pca.m2_bra"="Random forest", 
                        "pred.lm_pca.m2_bra"="Multiple linear regression",
                        "pred.lm_pls.mall_bra"="Multiple linear regression",
                        "pred.rf_pls.mall_bra"="Random forest", 
                        "pred.rf_avg.m_bra"="Random forest",
                        "pred.lm_avg.m_bra"="Multiple linear regression")) %>% 
  mutate(Model = factor(Model, levels = c("Random forest", 
                                          "Multiple linear regression"))) %>%
  # models labels
  mutate(Dim_red = recode(Dim_red_Model, 
                          "pred.rf_pca.m2_bra"="pca.m.2", 
                          "pred.lm_pca.m2_bra"="pca.m.2",
                          "pred.rf_pls.mall_bra"="pls.m.all", 
                          "pred.lm_pls.mall_bra"="pls.m.all",
                          "pred.rf_avg.m_bra"="avg.m",
                          "pred.lm_avg.m_bra"="avg.m")) %>% 
  mutate(Dim_red = factor(Dim_red, levels = c("pca.m.2", "pls.m.all",
                                              "avg.m"))) %>%
  mutate(x = as.numeric(as.character(x)), 
         y = as.numeric(as.character(y)))

# > Store in list 
list_preds_cc <- list("01_USA"=tab.preds.sensi.usa.cc,
                      "02_BRA"=tab.preds.sensi.bra.cc)
save(list_preds_cc, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/preds_cc.rda")


