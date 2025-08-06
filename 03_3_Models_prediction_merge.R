
# -------------------------------------------------------------------------
# 
#       3. Predictive models - load and merge predictions for each model
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

# > Function to read rda into list
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

# -------------------------------------------------------------------------
# > Load predictions for each model, each country, and both outcomes 
tab_preds <- list(
  USA   = list(name="01_USA"),
  BRA   = list(name="02_BRA"),
  WORLD = list(name="03_WORLD")) %>% 
  map_dfr(., ~{
    
    # -----------------------
    # -----------------------
    # Predictions of yield 
    # > Cross-validated on years
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", .x$name, "/01_YEARS")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="01_Ya", 
             type_cv = "01_YEARS") %>% 
      # rename for consistency among cv procedures
      dplyr::rename("preds.fold_cv"="preds.year")
    # -----------------------
    # > Cross-validated on sites
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds/", .x$name, "/02_GEO")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds_geo <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="01_Ya", 
             type_cv = "02_GEO")
    
    # -----------------------
    # -----------------------
    # Predictions of yield anomalies
    # > Cross-validated on years
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", .x$name, "/01_YEARS")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds_ano <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="02_Ya_ano", 
             type_cv = "01_YEARS") %>% 
      # rename for consistency among cv procedures
      dplyr::rename("preds.fold_cv"="preds.year")
    # -----------------------
    # > Cross-validated on sites
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_preds_anomalies/", .x$name, "/02_GEO")
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    tab_preds_ano_geo <- Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path") %>% 
      mutate(outcome="02_Ya_ano", 
             type_cv = "02_GEO")
    # -----------------------
    rbind(tab_preds, tab_preds_geo, tab_preds_ano, tab_preds_ano_geo)
    
  }, .id="Country") %>% 
  dplyr::select(-path) 

# > Rename columns' name
head(tab_preds)
colnames(tab_preds) <- c("Country", "Fold_cv", "Model", "Site_year", "Observed", "Predicted", "N_predictors", "Outcome", "Type_cv")
names(tab_preds)

# > Compute RMSEP, NSE, R2 and Bias
# for each model in each country

tab_perf <- tab_preds %>% 
  #mutate(Observed = if_else(is.na(Observed) == T, 0, Observed)) %>% 
  group_by(Outcome, Country, Model, Type_cv) %>% 
  summarise(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted),
    R2    = caret::R2(obs = Observed,        pred = Predicted),
    Bias  = Metrics::bias(actual = Observed, predicted = Predicted),
    N_pred = max(N_predictors)
  )

# > Label models and country for plots
tab_perf_labelled <- tab_perf %>% 
  # set long format 
  pivot_longer(cols=c(RMSEP, NSE, R2, Bias), names_to = "pred_perf", values_to = "pred_perf_value") %>% 
  # > group of models (rf or lm)
  mutate(gpe_model = ifelse(sub("\\_.*", "", Model) == "lm", "lm", "rf")) %>%  
  mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
  mutate(Model = sub(".*_", "", Model)) %>% 
  # > change annual by seasonal
  mutate(Model = recode(Model, "avg.zscore.a"="avg.zscore.s", "avg.a"="avg.s")) %>%
  # > change order in Models
  mutate(Model = factor(Model, levels = rev(c("avg.zscore.s", 
                                              "avg.s", "avg.m", 
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
  # > type of data (daily, month, annual)
  mutate(data_type=case_when(
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",
                 "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                 'mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 
                 'fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all', 
                 'pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all') ~ "Daily climatic predictors", 
    Model %in% c("avg.s", "avg.zscores.s") ~ "Seasonal climatic predictors", 
    TRUE ~ "Monthly climatic predictors")) %>% 
  mutate(data_type = factor(data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Seasonal climatic predictors"))) %>% 
  # > dimension reduction techniques
  mutate(dimension_reduction=case_when(
    Model %in% c('pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all',   'pca.m.1',   'pca.m.2',   'pca.m.3',   'pca.m.all')  ~"PCA",
    Model %in% c('fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all',  'fpca.m.1',  'fpca.m.2',  'fpca.m.3',  'fpca.m.all') ~"FPCA",
    Model %in% c('mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 'mfpca.m.1', 'mfpca.m.2', 'mfpca.m.3', 'mfpca.m.all')~"MFPCA",
    Model %in% c("pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",   "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all")  ~"PLS",
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",  "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all") ~"FPLS",
    Model == "avg.zscore.s" ~ "Standardised\nmean",
    TRUE ~ "Averages")) %>% 
  mutate(dimension_reduction = factor(dimension_reduction, levels = c("Standardised\nmean", "Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS"))) %>% 
  # > country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) %>% 
  # > nb of scores 
  mutate(nb_scores = case_when(
    Model %in% c('pca.d.1',   'fpca.d.1',   'mfpca.d.1',   'pca.m.1',   'fpca.m.1',   'mfpca.m.1',   'pls.d.1',   'pls.m.1',   'fpls.d.1',   'fpls.m.1')  ~"1 score",
    Model %in% c('pca.d.2',   'fpca.d.2',   'mfpca.d.2',   'pca.m.2',   'fpca.m.2',   'mfpca.m.2',   'pls.d.2',   'pls.m.2',   'fpls.d.2',   'fpls.m.2')  ~"2 scores",
    Model %in% c('pca.d.3',   'fpca.d.3',   'mfpca.d.3',   'pca.m.3',   'fpca.m.3',   'mfpca.m.3',   'pls.d.3',   'pls.m.3',   'fpls.d.3',   'fpls.m.3')  ~"3 scores",
    Model %in% c('pca.d.all', 'fpca.d.all', 'mfpca.d.all', 'pca.m.all', 'fpca.m.all', 'mfpca.m.all', 'pls.d.all', 'pls.m.all', 'fpls.d.all', 'fpls.m.all')~"All scores",
    TRUE ~ "")) %>% 
  mutate(nb_scores = factor(nb_scores, levels = c("1 score", "2 scores", "3 scores", "All scores", ""))) %>% 
  # > outcome
  mutate(Outcome_lab = recode(Outcome, "01_Ya" = "Soybean yield", "02_Ya_ano"="Soybean yield anomaly"), 
         Outcome_lab = factor(Outcome_lab, levels = c("Soybean yield", "Soybean yield anomaly"))) %>% 
  # > type cv 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) %>% 
  # > best model per country and indicator
  group_by(Type_cv_lab, Outcome_lab, Country, pred_perf) %>% 
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
  # > negative values
  mutate(
    negative_pred_perf_value = if_else(pred_perf_value < 0, 1, 0),
    pred_perf_value_lab = if_else(pred_perf_value < 0, 0, pred_perf_value))

# > Save predictions and model performances
save(tab_preds,         file = paste0(data_path, "tab_preds.rda"))
save(tab_perf_labelled, file = paste0(data_path, "tab_perf_models.rda"))

