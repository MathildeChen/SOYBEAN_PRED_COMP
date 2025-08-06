# -------------------------------------------------------------------------
# 
#       3. Evaluate predictive performance of models 
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
# > Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))
load(paste0(data_path, "tab_world.rda"))

# > daily climatic predictors (for PLSR)
load(paste0(data_path, "01_days/list_data_day_usa.rda"))
load(paste0(data_path, "01_days/list_data_day_bra.rda"))
load(paste0(data_path, "01_days/list_data_day_world.rda"))

# ---------------------------------- 
# Functions 
# > Load function to cross-validate the models (year-by-year cross-validation)
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")


# ---------------------------------- 
# Model predictive performances assessed by cross-validation on years and on sites 
# - years: data for each year are used as the test set
# - sites: sites are splited into k folds (here k=10) and all data (i.e. data for all years for selected sites) for each fold are used as test

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

system.time({  # estimate run time
  pred_performances <- list(
    #USA   = list(data = tab_usa, data_day = list_data_day_usa, name="01_USA"),
    #BRA   = list(data = tab_bra, data_day = list_data_day_bra, name="02_BRA")
    WORLD = list(data = tab_world, data_day = list_data_day_world, name="03_WORLD")
    ) %>% 
    map(., ~{
      
      # > data to use for models fitting, name of the analysis, and bootstrap procedure
      dat_pred <- .x$data
      data_day <- .x$data_day
      country <-  .x$name
      tab_sites <- dat_pred %>% 
        distinct(x, y, gridcode)
      
      # --------------------------------------
      # > FIT MODELS & ASSESS PREDICTIVE PERFORMANCES
      #   - BASED ON CROSS-VALIDATION ON YEARS (N years = 35)
      #   - BASED ON CROSS-VALIDATION ON SITES (N folds = 10)
      
      # Models list
      list_models <- list(
        ## DAILY DATA 
        ## > PCA on daily data 
        pca.d.1   = list(name = "pca.d.1",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"))), collapse = " + ")),
        pca.d.2   = list(name = "pca.d.2",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"))), collapse = " + ")),
        pca.d.3   = list(name = "pca.d.3",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"), starts_with("PC3_day"))), collapse = " + ")),
        pca.d.all = list(name = "pca.d.all",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        ## > FPCA on daily data 
        fpca.d.1   = list(name = "fpca.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"))), collapse = " + ")),
        fpca.d.2   = list(name = "fpca.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"))), collapse = " + ")),
        fpca.d.3   = list(name = "fpca.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"), starts_with("FPC3_day"))), collapse = " + ")),
        fpca.d.all = list(name = "fpca.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        ## > MFPCA on daily data 
        mfpca.d.1   = list(name = "mfpca.d.1"  , formula = "MFPC1_day"),
        mfpca.d.2   = list(name = "mfpca.d.2"  , formula = "MFPC1_day + MFPC2_day"),
        mfpca.d.3   = list(name = "mfpca.d.3"  , formula = "MFPC1_day + MFPC2_day + MFPC3_day"),
        mfpca.d.all = list(name = "mfpca.d.all", formula = "MFPC1_day + MFPC2_day + MFPC3_day + MFPC4_day"),
        ## > PLS on daily data 
        pls.d.1     = list(name = "pls.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"))), collapse = " + ")),
        pls.d.2     = list(name = "pls.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"))), collapse = " + ")),
        pls.d.3     = list(name = "pls.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"), starts_with("PLS3_day"))), collapse = " + ")),
        pls.d.all   = list(name = "pls.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        # > FPLS on daily data 
        fpls.d.1    = list(name = "fpls.d.1",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.2    = list(name = "fpls.d.2",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_"), starts_with("FPLS2_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.3    = list(name = "fpls.d.3",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_"), starts_with("FPLS2_"), starts_with("FPLS3_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.all  = list(name = "fpls.d.all",  formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        
        # MONTHLY DATA
        # > Monthly average
        avg.m       = list(name = "avg.m",       formula = "monthly_max_2m_temperature_1          + monthly_max_2m_temperature_2          + monthly_max_2m_temperature_3          + monthly_max_2m_temperature_4          + monthly_max_2m_temperature_5          + monthly_max_2m_temperature_6          + monthly_max_2m_temperature_7 +
                                                            monthly_min_2m_temperature_1          + monthly_min_2m_temperature_2          + monthly_min_2m_temperature_3          + monthly_min_2m_temperature_4          + monthly_min_2m_temperature_5          + monthly_min_2m_temperature_6          + monthly_min_2m_temperature_7 +
                                                            monthly_et0_1                         + monthly_et0_2                         + monthly_et0_3                         + monthly_et0_4                         + monthly_et0_5                         + monthly_et0_6                         + monthly_et0_7 +
                                                            monthly_vapor_pressure_deficit_1      + monthly_vapor_pressure_deficit_2      + monthly_vapor_pressure_deficit_3      + monthly_vapor_pressure_deficit_4      + monthly_vapor_pressure_deficit_5      + monthly_vapor_pressure_deficit_6      + monthly_vapor_pressure_deficit_7 +
                                                            monthly_surface_net_solar_radiation_1 + monthly_surface_net_solar_radiation_2 + monthly_surface_net_solar_radiation_3 + monthly_surface_net_solar_radiation_4 + monthly_surface_net_solar_radiation_5 + monthly_surface_net_solar_radiation_6 + monthly_surface_net_solar_radiation_7 + 
                                                            monthly_total_precipitation_1         + monthly_total_precipitation_2         + monthly_total_precipitation_3         + monthly_total_precipitation_4         + monthly_total_precipitation_5         + monthly_total_precipitation_6         + monthly_total_precipitation_7"),
        # > PCA on monthly data 
        pca.m.1      = list(name = "pca.m.1",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"))), collapse = " + ")),
        pca.m.2      = list(name = "pca.m.2",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"))), collapse = " + ")),
        pca.m.3      = list(name = "pca.m.3",    formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"), starts_with("PC3_month"))), collapse = " + ")),
        pca.m.all    = list(name = "pca.m.all",  formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        ## > FPCA on monthly data 
        fpca.m.1     = list(name = "fpca.m.1",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"))), collapse = " + ")),
        fpca.m.2     = list(name = "fpca.m.2",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"))), collapse = " + ")),
        fpca.m.3     = list(name = "fpca.m.3",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"), starts_with("FPC3_month"))), collapse = " + ")),
        fpca.m.all   = list(name = "fpca.m.all", formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        ## > MFPCA on monthly data 
        mfpca.1      = list(name = "mfpca.m.1",  formula = "MFPC1_month"),
        mfpca.2      = list(name = "mfpca.m.2",  formula = "MFPC1_month + MFPC2_month"),
        mfpca.3      = list(name = "mfpca.m.3",  formula = "MFPC1_month + MFPC2_month + MFPC3_month"),
        mfpca.m.all  = list(name = "mfpca.m.all", formula = "MFPC1_month + MFPC2_month + MFPC3_month + MFPC4_month"),
        # > PLS on monthly data 
        pls.m.1     = list(name = "pls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"))), collapse = " + ")),
        pls.m.2     = list(name = "pls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"))), collapse = " + ")),
        pls.m.3     = list(name = "pls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"), starts_with("PLS3_month"))), collapse = " + ")),
        pls.m.all   = list(name = "pls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS")) %>% dplyr::select(contains("month"))), collapse = " + "))
        # > FPLS on monthly data 
        fpls.m.1    = list(name = "fpls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.2    = list(name = "fpls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_"), starts_with("FPLS2_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.3    = list(name = "fpls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS1_"), starts_with("FPLS2_"), starts_with("FPLS3_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.all  = list(name = "fpls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPLS")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        #
        ## ANNUAL DATA
        avg.a        = list(name = "avg.a",        formula = "year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit"),
        avg.zscore.a = list(name = "avg.zscore.a", formula = "mean_year_z_scores")
      )
      
      # --------------------------------------
      
      # Models fitting and cross-validation (full data provided)
      list_fit_cv <- list_models %>% 
        map_dfr(., ~{
          
          # ----------------------------------------------------------------------------
          # --------------------------------------
          # PREDICTED OUTCOME = SOYBEAN YIELD 
          # --------------------------------------
          # model formula
          model_name    <- .x$name
          #model_formula <- paste0("Ya ~ irrigated_portion + ", .x$formula)
            
          
          ## --------------------------------------
          ## CROSS-VALIDATION ON YEARS
          ## --------------------------------------
          ## RANDOM FOREST
          mod_rf_cv <- function_cv_year(model      = model_formula, 
                                        outcome    = "Ya",
                                        data       = dat_pred, 
                                        data_day   = data_day,
                                        model_name = paste0("rf_", model_name), 
                                        res        = "perf",
                                        save       = TRUE, 
                                        path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", country, "/01_YEARS"))
          
          # LINEAR REGRESSION
          mod_lm_cv <- function_cv_year(model      = model_formula, 
                                        outcome    = "Ya",
                                        model_gpe  = "lm",
                                        data       = dat_pred, 
                                        data_day   = data_day,
                                        model_name = paste0("lm_", model_name), 
                                        res        = "perf",
                                        save       = TRUE, 
                                        path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", country, "/01_YEARS"))
          
          ## --------------------------------------
          ## CROSS-VALIDATION ON SITES
          ## --------------------------------------
          ## RANDOM FOREST
          mod_rf_cv_geo <- function_cv_geo(model      = model_formula, 
                                           outcome    = "Ya",
                                           data       = dat_pred, 
                                           data_day   = data_day,
                                           model_name = paste0("rf_", model_name), 
                                           res        = "perf",
                                           save       = TRUE, 
                                           path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", country, "/02_GEO"))
          
          # LINEAR REGRESSION 
          mod_lm_cv_geo <- function_cv_geo(model       = model_formula, 
                                           outcome    = "Ya",
                                           model_gpe  = "lm",
                                           data       = dat_pred, 
                                           data_day   = data_day,
                                           model_name = paste0("lm_", model_name), 
                                           res        = "perf",
                                           save       = TRUE, 
                                           path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", country, "/02_GEO"))
          
          rm(model_formula, mod_rf_cv, mod_lm_cv, mod_rf_cv_geo, mod_lm_cv_geo)
          
          # ----------------------------------------------------------------------------
          # --------------------------------------
          # PREDICTED OUTCOME = SOYBEAN YIELD ANOMALY
          # --------------------------------------
          # model formula
          model_formula <- paste0("Ya_ano ~ irrigated_portion + ", .x$formula)
          
          # --------------------------------------
          # CROSS-VALIDATION ON YEARS
          # --------------------------------------
          # RANDOM FOREST
          mod_rf_cv <- function_cv_year(model      = model_formula, 
                                        outcome    = "Ya_ano",
                                        data       = dat_pred, 
                                        data_day   = data_day,
                                        model_name = paste0("rf_", model_name), 
                                        res        = "perf",
                                        save       = TRUE, 
                                        path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", country, "/01_YEARS"))
          
          # LINEAR REGRESSION
          mod_lm_cv <- function_cv_year(model      = model_formula, 
                                        outcome    = "Ya_ano",
                                        model_gpe  = "lm",
                                        data       = dat_pred, 
                                        data_day   = data_day,
                                        model_name = paste0("lm_", model_name), 
                                        res        = "perf",
                                        save       = TRUE, 
                                        path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", country, "/01_YEARS"))
          
          # --------------------------------------
          # CROSS-VALIDATION ON SITES
          # --------------------------------------
          # RANDOM FOREST
          mod_rf_cv_geo <- function_cv_geo(model      = model_formula,
                                           outcome    = "Ya_ano", 
                                           data       = dat_pred, 
                                           data_day   = data_day,
                                           model_name = paste0("rf_", model_name), 
                                           res        = "perf",
                                           save       = TRUE, 
                                           path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", country, "/02_GEO"))
          
          # LINEAR REGRESSION 
          mod_lm_cv_geo <- function_cv_geo(model       = model_formula,
                                           outcome    = "Ya_ano", 
                                           model_gpe  = "lm",
                                           data       = dat_pred, 
                                           data_day   = data_day,
                                           model_name = paste0("lm_", model_name), 
                                           res        = "perf",
                                           save       = TRUE, 
                                           path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", country, "/02_GEO"))
            
          rm(model_formula, mod_rf_cv, mod_lm_cv, mod_rf_cv_geo, mod_lm_cv_geo)
          # ----------------------------------------------------------------------------
          
         }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
  })

beepr::beep(1)

