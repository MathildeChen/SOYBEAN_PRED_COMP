
# -------------------------------------------------------------------------
# 
#       3. Predictive models - yield anomalies
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
# Compute yield anomalies
# USA
tab_usa_max <- tab_usa %>% 
  dplyr::select(x, y, year, Ya) %>% 
  group_by(x, y) %>% 
  mutate(year_max = max(year)) %>% 
  filter(year==year_max) %>% 
  dplyr::rename("Ya_max"="Ya") %>% 
  dplyr::select(-year, -year_max)

tab_usa_2 <- left_join(tab_usa, tab_usa_max, by = c("x", "y")) %>% 
  ungroup() %>%
  mutate(Ya_ano = Ya - Ya_max)

ggplot(data = tab_usa_2) + geom_histogram(aes(x=Ya_ano), bins = 100)
summary(tab_usa_2$Ya_ano)

# BRA
tab_bra_max <- tab_bra %>% 
  dplyr::select(x, y, year, Ya) %>% 
  group_by(x, y) %>% 
  mutate(year_max = max(year)) %>% 
  filter(year==year_max) %>% 
  dplyr::rename("Ya_max"="Ya") %>% 
  dplyr::select(-year, -year_max)

tab_bra_2 <- left_join(tab_bra, tab_bra_max, by = c("x", "y")) %>% 
  ungroup() %>%
  mutate(Ya_ano = Ya - Ya_max)
summary(tab_bra_2$Ya_ano)

# WORLD
tab_world_max <- tab_world %>% 
  dplyr::select(x, y, year, Ya) %>% 
  group_by(x, y) %>% 
  mutate(year_max = max(year)) %>% 
  filter(year==year_max) %>% 
  dplyr::rename("Ya_max"="Ya") %>% 
  dplyr::select(-year, -year_max)

tab_world_2 <- left_join(tab_world, tab_world_max, by = c("x", "y")) %>% 
  ungroup() %>%
  mutate(Ya_ano = Ya - Ya_max)
summary(tab_world_2$Ya_ano)

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
    #USA   = list(data = tab_usa_2, name="01_USA"),
    #  BRA   = list(data = tab_bra_2, name="02_BRA")
    WORLD = list(data = tab_world_2, name="03_WORLD")
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
        ## > FPLS on monthly data 
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
          model_formula <- paste0("Ya_ano ~ irrigated_portion + ", .x$formula)
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
                                     path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds_anomalies/", country))
          
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
    USA   = list(data = tab_usa_2,   name="01_USA"),
    BRA   = list(data = tab_bra_2,   name="02_BRA"),
    WORLD = list(data = tab_world_2, name="03_WORLD")
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
          model_formula <- paste0("Ya_ano ~ irrigated_portion + ", .x$formula)
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
                                     path_save  = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds_anomalies/", country))

          
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





# > Function to read rda into list
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

tab_preds <- list(
  USA   = list(name="01_USA")
  #, BRA   = list(name="02_BRA"),
  #, WORLD = list(name="03_WORLD")
) %>% 
  map_dfr(., ~{
    
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds/", .x$name)
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path")
    
  }, .id="Country") %>% 
  dplyr::select(-path) %>% 
  mutate(outcome="yield")

tab_preds_ano <- list(
  #USA   = list(name="01_USA")
  #, BRA   = list(name="02_BRA"),
  WORLD = list(name="03_WORLD")
) %>% 
  map_dfr(., ~{
    
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds_anomalies/", .x$name)
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path")
    
  }, .id="Country") %>% 
  dplyr::select(-path) %>% 
  mutate(outcome="yield anomalie") %>% 
  left_join(., tab_usa_2 %>% 
              dplyr::select(x, y, year, Ya_ano) %>% 
              unite("site_year", x:year) %>% 
              dplyr::rename("Observed"="Ya_ano"), 
            by=c("preds.site_year"="site_year")) %>% 
  dplyr::select(Country, preds.year, preds.Model, preds.site_year, Observed, preds.Ya_pred, preds.N_predictors, outcome)

tab_preds_ano_yield_max <- list(
  USA   = list(name="01_USA")
  #, BRA   = list(name="02_BRA"),
  #, WORLD = list(name="03_WORLD")
) %>% 
  map_dfr(., ~{
    
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds_anomalies/", .x$name)
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path")
    
  }, .id="Country") %>% 
  dplyr::select(-path) %>% 
  mutate(outcome="yield anomalie") %>% 
  left_join(., tab_usa_2 %>% 
              dplyr::select(x, y, year, Ya_max) %>% 
              unite("site_year", x:year), 
            by=c("preds.site_year"="site_year")) %>%
  mutate(Predicted = Ya_max+preds.Ya_pred) %>% 
  mutate(outcome="yield latest year+yield anomalie") %>% 
  dplyr::select(Country, preds.year, preds.Model, preds.site_year, preds.Ya_obs, Predicted, preds.N_predictors, outcome)

colnames(tab_preds) <- c("Country", "Year", "Model", "Site_year", "Observed", "Predicted", "N_predictors", "Outcome")
colnames(tab_preds_ano) <- c("Country", "Year", "Model", "Site_year", "Observed", "Predicted", "N_predictors", "Outcome")
colnames(tab_preds_ano_yield_max) <- c("Country", "Year", "Model", "Site_year", "Observed", "Predicted", "N_predictors", "Outcome")

tab_preds_2 <- rbind(tab_preds_ano, tab_preds, tab_preds_ano_yield_max) 

tab_preds_2 %>%
  gather(key="type", value="outcome", Observed, Predicted) %>% 
  filter(Country == "USA",
         Model %in% c("avg.zscore.a", "avg.a", "pca.m.1", "mfpca.d.1", "fpca.m.1")) %>% 
  ggplot(., aes(x=outcome, fill=type)) +
  geom_density(alpha=0.5) +
  facet_wrap(.~Outcome)

tab_preds_2_density <- tab_preds_2 %>%
  filter(Country == "USA",
         Model %in% c("avg.zscore.a", "avg.a", "pca.m.1", "mfpca.d.1", "fpca.m.1")) %>% 
  split(.$Model, .$Outcome) %>%
  map_dfr(., ~{ 
    .x$density <- get_density(x = .x$Observed, y = .x$Predicted, n=100)
    .x 
    }) 
  
py<-tab_preds_2_density %>% 
  filter(Outcome=="yield") %>%
  ggplot(.) +
  geom_point(aes(x=Predicted, y=Observed, color = density)) +
  geom_abline(slope = 1, color="blue") + 
  theme_bw() + 
  theme(legend.position="bottom") +
  facet_grid(Model~Outcome) +
  scale_color_gradientn(colours = c("grey", "red", "darkred", "black"), values = seq(0,1,by=0.25))  +
  lims(x=c(0,6), y=c(0,6))

pya<-tab_preds_2_density %>% 
  filter(Outcome=="yield anomalie") %>%
  ggplot(.) +
  geom_point(aes(x=Predicted, y=Observed, color = density)) +
  geom_abline(slope = 1, color="blue") + 
  theme_bw() + 
  theme(legend.position="bottom") +
  facet_grid(Model~Outcome) +
  scale_color_gradientn(colours = c("grey", "red", "darkred", "black"), values = seq(0,1,by=0.25))  +
  lims(x=c(-2,2), y=c(-2,2))

pyay<-tab_preds_2_density %>% 
  filter(Outcome=="yield latest year+yield anomalie") %>%
  ggplot(.) +
  geom_point(aes(x=Predicted, y=Observed, color = density)) +
  geom_abline(slope = 1, color="blue") + 
  theme_bw() + 
  theme(legend.position="bottom") +
  facet_grid(Model~Outcome) +
  scale_color_gradientn(colours = c("grey", "red", "darkred", "black"), values = seq(0,1,by=0.25))  +
  lims(x=c(0,6), y=c(0,6))

p2 <- plot_grid(py, pya, pyay, nrow=1)
ggsave(p2, filename = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/pred.obs.png",
       width=10, height=10)


tab_preds_2_lab <-tab_preds_2 %>% 
  group_by(Outcome, Country, Model) %>% 
  summarise(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted),
    R2    = caret::R2(obs = Observed,        pred = Predicted),
    Bias  = Metrics::bias(actual = Observed, predicted = Predicted),
    N_pred = max(N_predictors)
  ) %>% 
  # set long format 
  pivot_longer(cols=c(RMSEP, NSE, R2, Bias), names_to = "pred_perf", values_to = "pred_perf_value") %>% 
  # groupe of models (rf or lm)
  mutate(gpe_model = ifelse(sub("\\_.*", "", Model) == "lm", "lm", "rf")) %>% 
  mutate(Model = sub(".*_", "", Model)) %>% 
  # > change order in Models
  mutate(Model = factor(Model, levels = rev(c("avg.zscore.a", 
                                              "avg.a", "avg.m", 
                                              "pca.m.1",   "pca.m.2",   "pca.m.3",   "pca.m.all",
                                              "pca.d.1",   "pca.d.2",   "pca.d.3",   "pca.d.all",
                                              "fpca.m.1",  "fpca.m.2",  "fpca.m.3",  "fpca.m.all",
                                              "fpca.d.1",  "fpca.d.2",  "fpca.d.3",  "fpca.d.all",
                                              "mfpca.m.1", "mfpca.m.2", "mfpca.m.3", "mfpca.m.all",
                                              "mfpca.d.1", "mfpca.d.2", "mfpca.d.3", "mfpca.d.all",
                                              "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all",
                                              "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                                              "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all",
                                              "fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all")))) %>% 
  # type of data (daily, month, annual)
  mutate(data_type=case_when(
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",
                 "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                 'mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 
                 'fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all', 
                 'pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all') ~ "Daily climatic predictors", 
    Model %in% c("avg.a", "avg.zscores.a") ~ "Annual climatic predictors", 
    TRUE ~ "Monthly climatic predictors")) %>% 
  mutate(data_type = factor(data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Annual climatic predictors"))) %>% 
  # dimension reduction techniques
  mutate(dimension_reduction=case_when(
    Model %in% c('pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all',   'pca.m.1',   'pca.m.2',   'pca.m.3',   'pca.m.all')  ~"PCA",
    Model %in% c('fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all',  'fpca.m.1',  'fpca.m.2',  'fpca.m.3',  'fpca.m.all') ~"FPCA",
    Model %in% c('mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 'mfpca.m.1', 'mfpca.m.2', 'mfpca.m.3', 'mfpca.m.all')~"MFPCA",
    Model %in% c("pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",   "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all")  ~"PLS",
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",  "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all") ~"FPLS",
    Model == "avg.zscore.a" ~ "Standardised\nmean",
    TRUE ~ "Averages")) %>% 
  mutate(dimension_reduction = factor(dimension_reduction, levels = c("Standardised\nmean", "Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS"))) %>% 
  # best model per country and indicator
  group_by(Country, pred_perf) %>% 
  mutate(best_pred_perf_value = if_else(pred_perf %in% c("RMSEP", "Bias"), min(abs(pred_perf_value)), max(pred_perf_value))) %>% 
  filter(pred_perf %in% c("RMSEP", "NSE")) %>%
  mutate(pred_perf_lab = case_when(
    pred_perf == "RMSEP" ~ "RMSEP\n(lower is better)", 
    pred_perf == "Bias"  ~ "Bias\n(lower is better)", 
    pred_perf == "NSE"   ~ "NSE\n(higher is better)",
    pred_perf == "R2"    ~ "R squared\n(higher is better)"),
    pred_perf_lab = factor(pred_perf_lab, levels = c("RMSEP\n(lower is better)", 
                                                     "NSE\n(higher is better)",
                                                     "R squared\n(higher is better)",
                                                     "Bias\n(lower is better)"))) %>% 
  # country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) %>% 
  # nb of scores 
  mutate(nb_scores = case_when(
    Model %in% c('pca.d.1',   'fpca.d.1',   'mfpca.d.1',   'pca.m.1',   'fpca.m.1',   'mfpca.m.1',   'pls.d.1',   'pls.m.1',   'fpls.d.1',   'fpls.m.1')  ~"1 score",
    Model %in% c('pca.d.2',   'fpca.d.2',   'mfpca.d.2',   'pca.m.2',   'fpca.m.2',   'mfpca.m.2',   'pls.d.2',   'pls.m.2',   'fpls.d.2',   'fpls.m.2')  ~"2 scores",
    Model %in% c('pca.d.3',   'fpca.d.3',   'mfpca.d.3',   'pca.m.3',   'fpca.m.3',   'mfpca.m.3',   'pls.d.3',   'pls.m.3',   'fpls.d.3',   'fpls.m.3')  ~"3 scores",
    Model %in% c('pca.d.all', 'fpca.d.all', 'mfpca.d.all', 'pca.m.all', 'fpca.m.all', 'mfpca.m.all', 'pls.d.all', 'pls.m.all', 'fpls.d.all', 'fpls.m.all')~"All scores",
    TRUE ~ "")) %>% 
  mutate(nb_scores = factor(nb_scores, levels = c("1 score", "2 scores", "3 scores", "All scores", "")))


tab_preds_2_lab %>%
  filter(#pred_perf == "NSE", 
         Country_lab == "United-States of America\n(N=29803)",
         Model %in% c("avg.zscore.a", "avg.a", "pca.m.1", "mfpca.d.1", "fpca.m.1"),
         gpe_model =="rf") %>%
  # > plot
  ggplot(., aes(x = pred_perf_value, y=Model)) +
  geom_point(aes(color = Outcome),
             size = 2, shape = 21) + 
  facet_grid(dimension_reduction ~ Country_lab, scales = "free_y", space = "free", switch = "y") + 
  theme_bw() +
  theme(legend.position = 'bottom', 
        axis.title.y = element_blank(),
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0, size=9),
        strip.text.x = element_text(size=10), 
        strip.background = element_rect(fill="grey", color = "transparent"),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10),
        #panel.grid.major.x = element_line(color = "grey", linewidth = 0.01),
        #panel.grid.minor.x = element_line(color = "grey", linewidth = 0.01),
        panel.grid = element_blank()) + 
  scale_color_manual(values = c("black", "red", "blue"), name = "Outcome") + 
  scale_x_continuous(breaks = seq(0.2, 1.2, by=0.2), labels = seq(0.2, 1.2, by=0.2)) + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1), 
         fill = guide_none(), 
         shape = guide_none())   + 
  labs(x = "Nash-Sutcliffe model efficiency (higher value indicates a better predictive performance)")
