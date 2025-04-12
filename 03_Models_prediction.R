# -------------------------------------------------------------------------
# 
#       3. Predictive models 
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

# ---------------------------------- 
# Functions 
# > Load function to cross-validate the models (year-by-year cross-validation)
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# ---------------------------------- 
# RANDOM FOREST
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
    USA   = list(data = tab_usa, name="01_USA")
    , BRA   = list(data = tab_bra, name="02_BRA")
    #, WORLD = list(data = tab_world, name="03_WORLD")
    ) %>% 
    map(., ~{
      
      # > data to use for models fitting, name of the analysis, and bootstrap procedure
      dat_pred <- .x$data
      country <-  .x$name
      tab_sites <- dat_pred %>% 
        distinct(x, y, gridcode)
      
      # --------------------------------------
      # > FIT MODELS & ASSESS PREDICTIVE PERFORMANCES
      #   BASED ON CROSS-VALIDATION ON YEARS (N years = 35)
      
      # Models list
      list_models <- list(
        # DAILY DATA 
        # > PCA on daily data 
        pca.d.1   = list(name = "pca.d.1",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"))), collapse = " + ")),
        pca.d.2   = list(name = "pca.d.2",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"))), collapse = " + ")),
        pca.d.3   = list(name = "pca.d.3",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"), starts_with("PC3_day"))), collapse = " + ")),
        pca.d.all = list(name = "pca.d.all",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"), starts_with("PC3_day"), starts_with("PC4_day"), starts_with("PC5_day"), starts_with("PC6_day"), starts_with("PC7_day"), starts_with("PC8_day"), starts_with("PC9_day"), starts_with("PC10_day"), starts_with("PC11_day"), starts_with("PC12_day"))), collapse = " + ")),
        # > FPCA on daily data 
        fpca.d.1   = list(name = "fpca.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"))), collapse = " + ")),
        fpca.d.2   = list(name = "fpca.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"))), collapse = " + ")),
        fpca.d.3   = list(name = "fpca.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"), starts_with("FPC3_day"))), collapse = " + ")),
        fpca.d.all = list(name = "fpca.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"), starts_with("FPC3_day"), starts_with("FPC4_day"))), collapse = " + ")),
        # > MFPCA on daily data 
        mfpca.d.1   = list(name = "mfpca.d.1"  , formula = "MFPC1_day"),
        mfpca.d.2   = list(name = "mfpca.d.2"  , formula = "MFPC1_day + MFPC2_day"),
        mfpca.d.3   = list(name = "mfpca.d.3"  , formula = "MFPC1_day + MFPC2_day + MFPC3_day"),
        mfpca.d.all = list(name = "mfpca.d.all", formula = "MFPC1_day + MFPC2_day + MFPC3_day + MFPC4_day"),
        # > PLS on daily data 
        pls.d.1     = list(name = "pls.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"))), collapse = " + ")),
        pls.d.2     = list(name = "pls.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"))), collapse = " + ")),
        pls.d.3     = list(name = "pls.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"), starts_with("PLS3_day"))), collapse = " + ")),
        pls.d.all   = list(name = "pls.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("PLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        # > FPLS on daily data 
        fpls.d.1    = list(name = "fpls.d.1",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.2    = list(name = "fpls.d.2",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.3    = list(name = "fpls.d.3",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_"), contains("FPLS3_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.all  = list(name = "fpls.d.all",  formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        
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
        pca.m.all    = list(name = "pca.m.all",  formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"), starts_with("PC3_month"), starts_with("PC4_month"), starts_with("PC5_month"), starts_with("PC6_month"), starts_with("PC7_month"))), collapse = " + ")),
        # > FPCA on monthly data 
        fpca.m.1     = list(name = "fpca.m.1",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"))), collapse = " + ")),
        fpca.m.2     = list(name = "fpca.m.2",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"))), collapse = " + ")),
        fpca.m.3     = list(name = "fpca.m.3",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"), starts_with("FPC3_month"))), collapse = " + ")),
        fpca.m.all   = list(name = "fpca.m.all", formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"), starts_with("FPC3_month"), starts_with("FPC4_month"))), collapse = " + ")),
        # > MFPCA on monthly data 
        mfpca.1      = list(name = "mfpca.m.1",  formula = "MFPC1_month"),
        mfpca.2      = list(name = "mfpca.m.2",  formula = "MFPC1_month + MFPC2_month"),
        mfpca.3      = list(name = "mfpca.m.3",  formula = "MFPC1_month + MFPC2_month + MFPC3_month"),
        mfpca.m.all  = list(name = "mfpca.m.all", formula = "MFPC1_month + MFPC2_month + MFPC3_month + MFPC4_month"),
        # > PLS on monthly data 
        pls.m.1     = list(name = "pls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"))), collapse = " + ")),
        pls.m.2     = list(name = "pls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"))), collapse = " + ")),
        pls.m.3     = list(name = "pls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"), starts_with("PLS3_month"))), collapse = " + ")),
        pls.m.all   = list(name = "pls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("PLS")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        # > FPLS on monthly data 
        fpls.m.1    = list(name = "fpls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.2    = list(name = "fpls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.3    = list(name = "fpls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_"), contains("FPLS3_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.all  = list(name = "fpls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        
        # ANNUAL DATA
        avg.a        = list(name = "avg.a",        formula = "year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit"),
        avg.zscore.a = list(name = "avg.zscore.a", formula = "mean_year_z_scores")
      )
      
      # --------------------------------------
      
      # Models fitting and cross-validation on years (full data provided)
      list_fit_cv <- list_models %>% 
        map_dfr(., ~{
          
          # model formula
          model_formula <- paste0("Ya ~ irrigated_portion + ", .x$formula)
          model_name    <- .x$name
          
          # > fit
          set.seed(101)
          mod  <- ranger(as.formula(model_formula),
                         data=dat_pred, 
                         num.tree=500,
                         importance="impurity") 
          
          # > cross-validation & save prediction for each model and site-year
          mod_cv <- function_cv_year(model      = model_formula, 
                                     data       = dat_pred, 
                                     model_name = model_name, 
                                     res        = "perf",
                                     save       = TRUE, 
                                     path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds/", country))
          
          # --------------------------------------
          # BOOTSTRAP PROCEDURE TO COMPUTE 95% CI of RMSEP, NSE, R2, and Bias
          
          # bootstrapping
          #mod_bs <- boot(data          = tab_sites,     # table with sites for the analysis --> bootstrap on sites
          #               tab_test      = dat_pred,      # table with full data for each site-years
          #               model_formula = model_formula, # model formula
          #               statistic     = bs,            # function to run in the bootstrap procedure
          #               R             = 500)           # number of bootstrap replications
          
          # compute CI for each indicator
          #RMSEP_ci <- getCI.bs(mod_bs, 1)
          #NSE_ci   <- getCI.bs(mod_bs, 2)
          #R2_ci    <- getCI.bs(mod_bs, 3)
          #Bias_ci  <- getCI.bs(mod_bs, 4)
          
          # --------------------------------------
          # OUTPUT
          data.frame(
            Model         = .x$name,
            N_predictors  = mod$num.independent.variables,
            model_formula = .x$formula,
            # Performance indicators + corresponding 95% CIs
            RMSEP = mod_cv$RMSEP, #RMSEP.conf.low = RMSEP_ci$conf.low, RMSEP.conf.high = RMSEP_ci$conf.high,
            NSE = mod_cv$NSE,     #NSE.conf.low   = NSE_ci$conf.low,   NSE.conf.high   = NSE_ci$conf.high,
            R2 = mod_cv$R2,       #R2.conf.low    = R2_ci$conf.low,    R2.conf.high    = R2_ci$conf.high,
            Bias = mod_cv$Bias)   #Bias.conf.low = Bias_ci$conf.low,   Bias.conf.high  = Bias_ci$conf.high)
            
          
         }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
  })

beepr::beep(1)
# ---------------------------------- 
# LINEAR MODELS (for composite variables)
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
    #USA   = list(data = tab_usa, name="01_USA"),
    #BRA   = list(data = tab_bra, name="02_BRA")
    WORLD = list(data = tab_world, name="03_WORLD")
  ) %>% 
    map(., ~{
      
      # > data to use for models fitting, name of the analysis, and bootstrap procedure
      dat_pred <- .x$data
      country <-  .x$name
      tab_sites <- dat_pred %>% 
        distinct(x, y, gridcode)
      
      # --------------------------------------
      # > FIT MODELS & ASSESS PREDICTIVE PERFORMANCES
      #   BASED ON CROSS-VALIDATION ON YEARS (N years = 35)
      
      # Models list
      # Models list
      list_models <- list(
        # DAILY DATA 
        # > PCA on daily data 
        pca.d.1   = list(name = "pca.d.1",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"))), collapse = " + ")),
        pca.d.2   = list(name = "pca.d.2",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"))), collapse = " + ")),
        pca.d.3   = list(name = "pca.d.3",       formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"), starts_with("PC3_day"))), collapse = " + ")),
        pca.d.all = list(name = "pca.d.all",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_day"), starts_with("PC2_day"), starts_with("PC3_day"), starts_with("PC4_day"), starts_with("PC5_day"), starts_with("PC6_day"), starts_with("PC7_day"), starts_with("PC8_day"), starts_with("PC9_day"), starts_with("PC10_day"), starts_with("PC11_day"), starts_with("PC12_day"))), collapse = " + ")),
        # > FPCA on daily data 
        fpca.d.1   = list(name = "fpca.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"))), collapse = " + ")),
        fpca.d.2   = list(name = "fpca.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"))), collapse = " + ")),
        fpca.d.3   = list(name = "fpca.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"), starts_with("FPC3_day"))), collapse = " + ")),
        fpca.d.all = list(name = "fpca.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_day"), starts_with("FPC2_day"), starts_with("FPC3_day"), starts_with("FPC4_day"))), collapse = " + ")),
        # > MFPCA on daily data 
        mfpca.d.1   = list(name = "mfpca.d.1"  , formula = "MFPC1_day"),
        mfpca.d.2   = list(name = "mfpca.d.2"  , formula = "MFPC1_day + MFPC2_day"),
        mfpca.d.3   = list(name = "mfpca.d.3"  , formula = "MFPC1_day + MFPC2_day + MFPC3_day"),
        mfpca.d.all = list(name = "mfpca.d.all", formula = "MFPC1_day + MFPC2_day + MFPC3_day + MFPC4_day"),
        # > PLS on daily data 
        pls.d.1     = list(name = "pls.d.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"))), collapse = " + ")),
        pls.d.2     = list(name = "pls.d.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"))), collapse = " + ")),
        pls.d.3     = list(name = "pls.d.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_day"), starts_with("PLS2_day"), starts_with("PLS3_day"))), collapse = " + ")),
        pls.d.all   = list(name = "pls.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("PLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        # > FPLS on daily data 
        fpls.d.1    = list(name = "fpls.d.1",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.2    = list(name = "fpls.d.2",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.3    = list(name = "fpls.d.3",    formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_"), contains("FPLS3_")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        fpls.d.all  = list(name = "fpls.d.all",  formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS")) %>% dplyr::select(contains("day"))), collapse = " + ")),
        
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
        pca.m.all    = list(name = "pca.m.all",  formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PC1_month"), starts_with("PC2_month"), starts_with("PC3_month"), starts_with("PC4_month"), starts_with("PC5_month"), starts_with("PC6_month"), starts_with("PC7_month"))), collapse = " + ")),
        # > FPCA on monthly data 
        fpca.m.1     = list(name = "fpca.m.1",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"))), collapse = " + ")),
        fpca.m.2     = list(name = "fpca.m.2",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"))), collapse = " + ")),
        fpca.m.3     = list(name = "fpca.m.3",   formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"), starts_with("FPC3_month"))), collapse = " + ")),
        fpca.m.all   = list(name = "fpca.m.all", formula = paste0(names(dat_pred %>% dplyr::select(starts_with("FPC1_month"), starts_with("FPC2_month"), starts_with("FPC3_month"), starts_with("FPC4_month"))), collapse = " + ")),
        # > MFPCA on monthly data 
        mfpca.1      = list(name = "mfpca.m.1",  formula = "MFPC1_month"),
        mfpca.2      = list(name = "mfpca.m.2",  formula = "MFPC1_month + MFPC2_month"),
        mfpca.3      = list(name = "mfpca.m.3",  formula = "MFPC1_month + MFPC2_month + MFPC3_month"),
        mfpca.m.all  = list(name = "mfpca.m.all", formula = "MFPC1_month + MFPC2_month + MFPC3_month + MFPC4_month"),
        # > PLS on monthly data 
        pls.m.1     = list(name = "pls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"))), collapse = " + ")),
        pls.m.2     = list(name = "pls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"))), collapse = " + ")),
        pls.m.3     = list(name = "pls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(starts_with("PLS1_month"), starts_with("PLS2_month"), starts_with("PLS3_month"))), collapse = " + ")),
        pls.m.all   = list(name = "pls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("PLS")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        # > FPLS on daily data 
        fpls.m.1    = list(name = "fpls.m.1",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.2    = list(name = "fpls.m.2",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.3    = list(name = "fpls.m.3",     formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS1_"), contains("FPLS2_"), contains("FPLS3_")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        fpls.m.all  = list(name = "fpls.m.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("FPLS")) %>% dplyr::select(contains("month"))), collapse = " + ")),
        
        # ANNUAL DATA
        avg.a        = list(name = "avg.a",        formula = "year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit"),
        avg.zscore.a = list(name = "avg.zscore.a", formula = "mean_year_z_scores")
      )
      
      # --------------------------------------
      
      # Models fitting and cross-validation on years (full data provided)
      list_fit_cv <- list_models %>% 
        map_dfr(., ~{
          
          # model formula
          model_formula <- paste0("Ya ~ irrigated_portion + ", .x$formula)
          model_name    <- .x$name
          
          # > fit
          mod  <- lm(as.formula(model_formula),
                         data=dat_pred) 
          
          # > cross-validation & save prediction for each model and site-year
          mod_cv <- function_cv_year(model      = model_formula, 
                                     model_gpe  = "lm",
                                     data       = dat_pred, 
                                     model_name = paste0("lm_", model_name), 
                                     res        = "perf",
                                     save       = TRUE, 
                                     path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds/", country))
          
          # --------------------------------------
          # BOOTSTRAP PROCEDURE TO COMPUTE 95% CI of RMSEP, NSE, R2, and Bias
          
          # bootstrapping
          #mod_bs <- boot(data          = tab_sites,     # table with sites for the analysis --> bootstrap on sites
          #               tab_test      = dat_pred,      # table with full data for each site-years
          #               model_formula = model_formula, # model formula
          #               statistic     = bs,            # function to run in the bootstrap procedure
          #               R             = 500)           # number of bootstrap replications
          
          # compute CI for each indicator
          #RMSEP_ci <- getCI.bs(mod_bs, 1)
          #NSE_ci   <- getCI.bs(mod_bs, 2)
          #R2_ci    <- getCI.bs(mod_bs, 3)
          #Bias_ci  <- getCI.bs(mod_bs, 4)
          
          # --------------------------------------
          # OUTPUT
          data.frame(
            Model         = .x$name,
            N_predictors  = length(mod$coefficients)-1,
            model_formula = .x$formula,
            # Performance indicators + corresponding 95% CIs
            RMSEP = mod_cv$RMSEP, #RMSEP.conf.low = RMSEP_ci$conf.low, RMSEP.conf.high = RMSEP_ci$conf.high,
            NSE = mod_cv$NSE,     #NSE.conf.low   = NSE_ci$conf.low,   NSE.conf.high   = NSE_ci$conf.high,
            R2 = mod_cv$R2,       #R2.conf.low    = R2_ci$conf.low,    R2.conf.high    = R2_ci$conf.high,
            Bias = mod_cv$Bias)   #Bias.conf.low = Bias_ci$conf.low,   Bias.conf.high  = Bias_ci$conf.high)
          
        }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
})



stop()

# > Test with points + lines
tab_perf_labelled %>% 
  filter(pred_perf == "RMSEP") %>% 
  mutate(nb_scores = factor(nb_scores, levels = rev(c("1 score", "2 scores", "3 scores", "All scores", ""))))%>% 
  ggplot(., aes(x = pred_perf_value, y = nb_scores, color = data_type, shape = gpe_model)) +
  geom_point(position = position_dodge2(width = 0.5)) + 
  geom_linerange(aes(xmin = 0, xmax = pred_perf_value),
                 position = position_dodge2(width = 0.5)) + 
  #geom_vline(aes(xintercept = best_pred_perf_value), color = "red", linetype=2) + 
  facet_grid(dimension_reduction ~ Country_lab, scales = "free_y", space = "free", switch = "y") + 
  theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0),
        legend.position = "bottom") + 
  labs(x = "Root mean square error of prediction") + 
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3"), name = "Climatic data type") + 
  guides(color = guide_legend(ncol = 1))

# > Test with heat maps
tab_perf_labelled %>% 
  filter(pred_perf == "RMSEP") %>% 
  #mutate(dimension_reduction = factor(dimension_reduction, levels = rev(c("Averages", "PCA", "FPCA", "MFPCA")))) %>% 
  mutate(data_type = recode(data_type, "Daily climatic predictors"="Daily", "Monthly climatic predictors"="Monthly", "Annual climatic predictors"="Annual")) %>% 
  mutate(pred_perf_value_lab = round(pred_perf_value, 3))  %>% 
  mutate(nb_scores = factor(nb_scores, levels = rev(c("1 score", "2 scores", "3 scores", "All scores", ""))))%>% 
  arrange(Country_lab, pred_perf_value) %>% 
  group_by(Country_lab) %>% 
  mutate(order = row_number()) %>%
  ggplot(., aes(x = data_type, 
                y = nb_scores, 
                fill = order, 
                label = pred_perf_value_lab)) + 
  geom_raster() +
  geom_text() + 
  theme_cowplot() + 
  theme(axis.title = element_blank(), 
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0)) +
  facet_grid(dimension_reduction~Country_lab, scales = "free", space = "free", switch = "y") + 
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = pal, name = "Order based on RMSEP\n(lower value is better)")

# > Test with ordering columns

#      Model   N_predictors R_squared Pred_error  cv_RMSEP    cv_NSE
#1     pca.d.1            7 0.6371912 0.17916810 0.4527288 0.5849429
#2     pca.d.2           13 0.6671153 0.16439052 0.4787495 0.5358607
#3     pca.d.3           19 0.6956387 0.15030465 0.4922830 0.5092489
#4   pca.d.all           48 0.7070793 0.14465484 0.5185287 0.4555259
#5    fpca.d.1            7 0.6634170 0.16621688 0.4369128 0.6134362
#6    fpca.d.2           13 0.6877426 0.15420400 0.4671610 0.5580585
#7    fpca.d.3           19 0.7234046 0.13659281 0.4817411 0.5300420
#8  fpca.d.all           25 0.7372297 0.12976549 0.4899798 0.5138302
#9     mfpca.d            5 0.7275392 0.13455100 0.4400883 0.6077967
#10         mm           43 0.8229764 0.08742065 0.4535029 0.5835223
#11    pca.m.1            7 0.7600521 0.11849497 0.3997485 0.6764024
#12    pca.m.2           13 0.8036085 0.09698522 0.4276999 0.6295668
#13    pca.m.3           19 0.8121387 0.09277271 0.4409187 0.6063152
#14  pca.m.all           43 0.8204804 0.08865328 0.4636733 0.5646328
#15   fpca.m.1            7 0.7463083 0.12528213 0.4176794 0.6467210
#16   fpca.m.2           13 0.7964701 0.10051044 0.4439341 0.6009120
#17   fpca.m.3           19 0.8117603 0.09295959 0.4461818 0.5968605
#18 fpca.m.all           25 0.8157593 0.09098474 0.4520550 0.5861774
#19    mfpca.m            5 0.7211530 0.13770475 0.4767167 0.5397938
#20         ma            7 0.7527266 0.12211254 0.4042229 0.6691178
#21  mfpca.d.1            2 0.7210687 0.13774637 0.3725334 0.7189640
#22  mfpca.d.2            3 0.6598594 0.16797374 0.4236238 0.6365937
#23  mfpca.d.3            4 0.7140171 0.14122871 0.4337152 0.6190737
#24  mfpca.m.1            2 0.7470199 0.12493075 0.3539170 0.7463502
#25  mfpca.m.2            3 0.6152375 0.19000967 0.4577338 0.5757151
#26  mfpca.m.3            4 0.6785798 0.15872893 0.4727441 0.5474319


library(dplyr)

#' @param rmse Root mean squared error on your sample
#' @param df Degrees of Freedom in your model. In this case it should be the
#'   same as the number of observations in your sample.
rmse_interval <- function(rmse, deg_free, p_lower = 0.025, p_upper = 0.975){
  tibble(.pred_lower = sqrt(deg_free / qchisq(p_upper, df = deg_free)) * rmse,
         .pred_upper = sqrt(deg_free / qchisq(p_lower, df = deg_free)) * rmse)
}

load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds/02_BRA/avg.a.rda")
preds %>% 
  summarise(rmse = caret::RMSE(Ya_obs, Ya_pred), 
            n = n()) %>% 
  mutate(rmse_interval(rmse, n))



# > Function to use in boot function
bs <- function(data_pred,      # table with full climatic data and yield data for site-years corresponding to each selected sites
               indices,       # indices to select sites
               tab_bra_sites = tab_bra_sites,  # table with unique combinaison of longitude and latitude selected for boostrap
               model_formula = "year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit", 
               seed = 101) # model formula
{
  
  # > Data for bootstrap and associated climatic data 
  sites_bs <- tab_bra_sites[indices,]
  data_pred_bs  <- data %>% 
    filter(gridcode %in% unique(sites_bs$gridcode))
  
  # > Model to fit
  model <- paste0("Ya ~ irrigated_portion + ", model_formula)
  
  # ------------------
  # > Check performance of these models compared to a model using monthly average data
  # Set seed 
  seed <- seed
  
  mod_bs_cv <- function_cv_year(model         = model, 
                                data          = data_pred_bs, 
                                model_name    = "test", 
                                res           = "perf")
  
  res.bs <- c(mod_bs_cv$RMSEP, mod_bs_cv$NSE, mod_bs_cv$R2, mod_bs_cv$Bias)
  return(res.bs)  
}

# -----------
# Function to compute 95% CI for bootstrapped estimates
# - x: a list of bootstrapped statistics, produced by implementing the bs.f() function (see above) in the boot::boot() function
# - w: the index of the bootstrapped vector for which computing 95%CI 

getCI.bs <- function(x, w) {
  
  ci_bca <- coxed::bca(x$t[,w])
  tab_bca <- data.frame("index" = w, "method"="bca", "statistic" = x$t0[w], "conf.low" = ci_bca[1], "conf.high" = ci_bca[2])
  
  return(data.frame(tab_bca, row.names=NULL))
}  


# > Number of replicates
B <- 10
tab_bra_sites <- tab_bra %>% filter(year %in% 1991:1993) %>% distinct(x, y, gridcode)
tab_test <- tab_bra %>% filter(year %in% 1991:1993)
# > Bootstrap

# > Corresponding CI
boot.ci(test.bs, conf = c(0.95), type=c("perc", "norm"))




test.bs.2 <- boot(data = tab_bra_sites,
                  tab_test=tab_test,
                  #predictors="year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit",
                  model_formula = "Ya ~ irrigated_portion + year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit",
                  statistic  = bs, 
                  R = B)#, parallel = "multicore", ncpus = n.cores)

# compute CI for each indicator
RMSEP_ci <- getCI.bs(test.bs.2, 1)
NSE_ci   <- getCI.bs(test.bs.2, 2)
R2_ci    <- getCI.bs(test.bs.2, 3)
Bias_ci  <- getCI.bs(test.bs.2, 4)

# --------------------------------------
# OUTPUT
data.frame(# Performance indicators + corresponding 95% CIs
  RMSEP.conf.low = RMSEP_ci$conf.low, RMSEP.conf.high = RMSEP_ci$conf.high,
  NSE.conf.low   = NSE_ci$conf.low,   NSE.conf.high   = NSE_ci$conf.high,
  R2.conf.low    = R2_ci$conf.low,    R2.conf.high    = R2_ci$conf.high,
  Bias.conf.low = Bias_ci$conf.low,   Bias.conf.high  = Bias_ci$conf.high)
  