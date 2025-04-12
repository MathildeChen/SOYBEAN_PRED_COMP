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

# ----------------------------------
# Data 
data_path <- "..."
save_path <- "..."

# USA 
load(paste0(data_path, "tab_usa.rda"))

# BRAZIL
load(paste0(data_path, "tab_bra.rda"))

# WORLD 
load(paste0(data_path, "tab_world.rda"))

# ----------------------------------
# CROSS VALIDATION ON YEARS
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
function_cv_year <- function(model, 
                             model_gpe  = "rf",
                             data, 
                             seed=101, 
                             model_name = "Give_me_a_name",
                             save = F,
                             path_save = NULL, 
                             res = "preds"){
  
  # Data to store predicted values
  data_pred <- list()
  
  # CROSS VALIDATION ON YEARS
  for(y in unique(data$year))
  {
    
    # ---------------------
    # Define test and train datasets
    data_test <- data[which(data$year == y),]
    #dim(data_test)
    
    data_train <- data[which(data$year != y),]
    #dim(data_train)
    
    # ---------------------
    # Fit model on the train dataset
    # random forest
    if(model_gpe == "rf")
    {
      # > fit on train
      set.seed(seed)
      mod_train <- ranger(as.formula(model),
                          data=data_train, 
                          num.tree=500,
                          importance="impurity")   
      
      # > predict on test
      pred.test <- predict(mod_train, 
                           data = data_test, 
                           type = "response", 
                           seed = seed, 
                           num.trees = 500)
      # > retrieve predictions & nb of predictors
      pred.test_vec <- pred.test$predictions
      N_predictors <- mod_train$num.independent.variables
      
    }
    # linear regression (benchmark)
    if(model_gpe == "lm")
    {
      # > fit on train
      mod_train <- lm(as.formula(model),
                      data=data_train)   
      
      # > predict on test
      pred.test <- predict(mod_train,
                           newdata = data_test, 
                           type = "response")
      # > retrieve predictions
      pred.test_vec <- as.numeric(as.character(pred.test))
      N_predictors <- length(mod_train$coefficients) - 1
    }
    
    
    # ---------------------
    # Store the predictions
    data_pred[[paste0(y)]] <- data.frame(Model     = model_name,
                                         site_year = data_test$site_year,
                                         Ya_obs    = data_test$Ya,
                                         Ya_pred   = pred.test_vec,
                                         N_predictors = N_predictors)
    
  }
  
  # > predictions in a table
  preds <- plyr::ldply(data_pred, data.frame, .id = "year")
  
  # > save result
  if(save==TRUE)
  {
    save(preds, file = paste0(path_save, "/", model_name, ".rda"))
  }
  
  # RESULTS RETURNED
  # > return predictions
  if(res == "preds")
  {
    res <- preds
  }
  # > return prediction performance indicators
  if(res == "perf")
  {
    # > performance indicators
    RMSEP <- caret::RMSE(obs = preds$Ya_obs,      pred = preds$Ya_pred)
    NSE   <- hydroGOF::NSE(obs = preds$Ya_obs,    sim=preds$Ya_pred)
    R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    res <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, "R2"=R2, "Bias"=Bias, "N_predictors"=N_predictors)
  }
  
  return(res)
  
}


# ---------------------------------- 
# RANDOM FOREST
system.time({  # estimate run time
  pred_performances <- list(
    USA   = list(data = tab_usa, name="01_USA"), 
    BRA   = list(data = tab_bra, name="02_BRA"), 
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
      list_models <- list(
        pls.d.all   = list(name = "pls.d.all",   formula = paste0(names(dat_pred %>% dplyr::select(contains("PLS")) %>% dplyr::select(contains("day"))), collapse = " + ")))
      
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
                                     path_save  = paste0(save_path, country))
          
          # --------------------------------------
          # OUTPUT
          data.frame(
            Model         = .x$name,
            N_predictors  = mod$num.independent.variables,
            model_formula = .x$formula,
            # Performance indicators + corresponding 95% CIs
            RMSEP = mod_cv$RMSEP, 
            NSE =   mod_cv$NSE,     
            R2 =    mod_cv$R2,       
            Bias =  mod_cv$Bias)   
          
          
        }) 
      
    })
  
  # >>> stop cluster//
  stopCluster(my.cluster)
  
})

beepr::beep(1)



