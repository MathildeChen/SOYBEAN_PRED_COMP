# functions for predictions
library(pls)

function_plsr2 <- function(type_data, 
                           data,
                           load_data    = F,
                           vars_names,
                           cum_clim.var = T,
                           scale        = T,
                           nb_comp){
  
  # -----------------
  # -----------------
  
  # > Object to store the PLS and outputs
  list_pls <- list()
  list_pls_scores <- list()
  
  # -----------------
  # -----------------
  # > Monthly data 
  if(type_data == "M")
  {
    
    # > Apply PLSR on each variable 
    for(var_i in unique(vars_names$clim.var)){
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("plsr on ", var_i_abb))
      
      # -----------------
      # > DATA SELECTION & PREPARATION FOR PLSR
      # > select data for the considered variable
      data_var_i <- data %>% 
        ungroup(.) %>% 
        dplyr::select(site_year, Ya, starts_with(paste0("monthly_", var_i)))
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
      
      #message("- ok data")
      
      # > SCALE DATA 
      if(scale == TRUE)
      {
        # > select data for the considered variable and scale it
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya) %>%
          scale(.) %>% 
          as.data.frame(.)
      }
      
      if(scale == FALSE)
      {
        # > do not scale data (not recommanded)
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya) 
        warning("Warning: data used for PLS is not scaled")
        
      }
      # > rownames are sites-year
      rownames(z_tab_pls_var_i) <- data_var_i$site_year
      z_tab_pls_var_i$site_year <- rownames(z_tab_pls_var_i)
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(z_tab_pls_var_i$site_year)), sort(unique(data$site_year)))
      
      #message("- ok scale data")
      
      # > join with yield data 
      tab_pls_var_i <- data_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        left_join(., as.data.frame(z_tab_pls_var_i), by = "site_year") 
      
      # -----------------
      # > PLSR
      plsr_var_i_cv <- plsr(Ya ~ ., 
                            ncomp = nb_comp, 
                            data = tab_pls_var_i %>% 
                              dplyr::select(-site_year), 
                            validation = "CV")
      
      #message("plsr ok")
      
      # -----------------
      # > store scores from PLSR
      var_i_plsr_scores <- plsr_var_i_cv$scores 
      colnames(var_i_plsr_scores) <- paste0("PLS", 1:ncol(var_i_plsr_scores), "_month_", var_i_abb)
      rownames(var_i_plsr_scores) <- tab_pls_var_i$site_year
      
      # -------------------------
      # > STORE PLS AND SCORES 
      list_pls_scores[[paste0(var_i_abb)]] <- var_i_plsr_scores
      list_pls[[paste0(var_i_abb)]]        <- plsr_var_i_cv
      
    }
    
    # > Extract scores for all variables 
    tab_PLS_scores <- do.call(cbind, list_pls_scores) %>% as.data.frame(.)
    message("score tab ok")
    
    # > Return the FPCA for each variable and the scores 
    res <- list("tab_PLS_scores"        = tab_PLS_scores,
                "list_pls_per_variable" = list_pls)
    
  }
  
  
  # -----------------
  # -----------------
  # > Daily data 
  if(type_data == "D")
  {
    
    # > For each variable 
    for(var_i in unique(vars_names$clim.var))
    {
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("plsr on ", var_i_abb))
      
      # -------------------------
      # > LOAD DATA (if load_data=T)
      if(load_data == T)
      {
        
        # > Load climatic data
        era5daily_var_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_data_", var_i, ".rda"))
        
        # > Check if each site-year has data
        testthat::expect_equal(sort(unique(era5daily_var_i$site_year)), sort(unique(data$site_year)))
        
      }
      # else use provided data (provided as list of the dataframe for each climatic variable)
      if(load_data == F)
      {
        era5daily_var_i <- data[[paste0(var_i)]]
      }
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR PLSR
      # > select data for the considered variable
      # > if accumulated data over growing season is used
      if(cum_clim.var == TRUE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.) 
      }
      # > if raw data is used
      if(cum_clim.var == FALSE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, Ya, day_of_year, clim.var, clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.) 
        warning("Warning: data used for PCA is not accumulated over growing season")
      }
      
      #message("- data ok")
      
      # > SCALE DATA 
      if(scale == TRUE)
      {
        # > select data for the considered variable and scale it
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya) %>%
          scale(.) %>% 
          as.data.frame(.)
      }
      
      if(scale == FALSE)
      {
        # > do not scale data (not recommanded)
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya)  %>% 
          as.data.frame(.)
        warning("Warning: data used for PLS is not scaled")
        
      }
      # > rownames are sites-year
      rownames(z_tab_pls_var_i) <- data_var_i$site_year
      z_tab_pls_var_i$site_year <- rownames(z_tab_pls_var_i)
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(z_tab_pls_var_i$site_year)), sort(unique(era5daily_var_i$site_year)))
      
      #message("- ok scale data")
      
      # > join with yield data 
      tab_pls_var_i <- data_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        left_join(., as.data.frame(z_tab_pls_var_i), by = "site_year")
      
      # -----------------
      # > PLSR
      plsr_var_i_cv <- plsr(Ya ~ ., 
                            ncomp = nb_comp, 
                            data = tab_pls_var_i %>% 
                              dplyr::select(-site_year), 
                            validation = "CV", 
                            scale = scale, center = scale)
      
      #message("plsr ok")
      
      # -----------------
      # > store scores from PLSR
      var_i_plsr_scores <- plsr_var_i_cv$scores 
      colnames(var_i_plsr_scores) <- paste0("PLS", 1:ncol(var_i_plsr_scores), "_day_", var_i_abb)
      rownames(var_i_plsr_scores) <- tab_pls_var_i$site_year
      
      # -------------------------
      # > STORE PLS AND SCORES 
      list_pls_scores[[paste0(var_i_abb)]] <- var_i_plsr_scores
      list_pls[[paste0(var_i_abb)]]        <- plsr_var_i_cv
      
    }
    
    # > Extract scores for all variables 
    tab_PLS_scores <- do.call(cbind, list_pls_scores) %>% as.data.frame(.)
    #message("score tab ok")
    
    # > Return the FPCA for each variable and the scores 
    res <- list("tab_PLS_scores"        = tab_PLS_scores,
                "list_pls_per_variable" = list_pls)
    
  }
  
  return(res)
  
}

# ----------------------------------
# CROSS VALIDATION ON YEARS
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"
function_cv_year <- function(model, 
                             outcome,
                             model_gpe  = "rf",
                             data, 
                             data_day,
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
    
    # --------------------------------------
    # PLSR: compute scores from train dataset and predict 
    # > monthly data 
    if(model %in% c("pls.m.1", "pls.m.2", "pls.m.3", "pls.m.all"))
    {
      
      message("refit PLSR scores for monthly data")
      # Fit PLSR model 
      plsr_train <- function_plsr2(type_data  = "M", 
                                   vars_names = vars_names,
                                   data       = data_train,
                                   nb_comp    = 7)
      
      # > list to store the results
      list_scores_pls <- list()
      
      # > Apply PLSR loads on new data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # > PLSR 
        plsr_j <- plsr_train$list_pls_per_variable[[paste0(clim.var_abb_j)]]
        
        # > Mean and sd of the original data to standardize 
        mu_train <- data_train %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
          colMeans(.)
        
        sd_train <- data_train %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
          apply(., 2, sd)
        
        # > New data 
        tab_test_j <- data_test %>% ungroup(.) %>%
          dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
        
        # > Scale new data based on mean and sd of the original data 
        z_tab_test_j <- scale(tab_test_j[,-1], center = mu_train, scale=sd_train)
        
        # > PLSR score from newdata
        tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_tab_test_j)
        
        colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_month_",clim.var_abb_j)
        
        # > Store
        list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
        
      }
      
      tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
      
      
      # Merge with data train
      data_train <- data_train %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., plsr_train$tab_PLS_scores)
      
      # Merge with data test 
      data_test <- data_test %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., tab_scores_plsr)
      
    }
    # > daily data 
    if(model %in% c("pls.d.1", "pls.d.2", "pls.d.3", "pls.d.all"))
    {
      
      message("refit PLSR scores for daily data")
      # Daily train and test data 
      data_day_train <- data_day %>% 
        map(., ~ {
          .x %>% filter(site_year %in% unique(data_train$site_year))  
        })
      
      data_day_test <- data_day %>% 
        map(., ~ {
          .x %>% filter(site_year %in% unique(data_test$site_year))  
        })
      
      
      # Fit PLSR model  
      plsr_train <- function_plsr2(type_data  = "D", 
                                   vars_names = vars_names,
                                   data       = data_day_train,
                                   nb_comp    = 214)
      
      
      # > list to store the results
      list_scores_pls <- list()
      
      # > Apply PCA loads on new data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # > PLSR 
        plsr_j <- plsr_train$list_pls_per_variable[[paste0(clim.var_abb_j)]]
        
        # > Select data for the var_j 
        data_day_train_j <- data_day_train[[paste0(var_j)]]
        data_day_test_j <- data_day_test[[paste0(var_j)]]
        
        # > Prepare for scaling 
        data_var_j <- data_day_train_j %>% 
          dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.)
        
        # > Mean and sd of the original data to standardize 
        mu_train <- data_var_j %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("day_", var_j))) %>% 
          colMeans(.)
        
        sd_train <- data_var_j %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("day_", var_j))) %>% 
          apply(., 2, sd)
        
        # > New data 
        tab_test_j <- data_day_test_j %>% 
          dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.)
        
        # > Scale new data based on mean and sd of the original data 
        z_tab_test_j <- scale(tab_test_j[,-1], center = mu_train, scale=sd_train)
        
        # > PLSR score from newdata
        tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_tab_test_j)
        
        colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_day_",clim.var_abb_j)
        
        # > Store
        list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
        
      }
      
      tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
      
      # Merge with data train
      data_train <- data_train %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., plsr_train$tab_PLS_scores)
      
      # Merge with data test 
      data_test <- data_test %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., tab_scores_plsr)
      
    }  
    
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
                                         Ya_obs    = data_test[paste0(outcome)],
                                         Ya_pred   = pred.test_vec,
                                         N_predictors = N_predictors) %>% 
      dplyr::rename("Ya_obs" = 3)
    
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
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    res <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(res)
  
}


# ----------------------------------
# CROSS VALIDATION ON SITES
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"

function_cv_geo <- function(model, 
                            outcome,
                            model_gpe  = "rf",
                            k_nb_folds=10,
                            data, 
                            seed=101, 
                            model_name = "Give_me_a_name",
                            save = F,
                            path_save = NULL, 
                            res = "preds"){
  
  # Data to store predicted values
  data_pred <- list()
  
  # CROSS VALIDATION ON LOCALISATION
  # Splits sites in k_nb_folds groups
  data_sites <- data %>% 
    distinct(gridcode, x, y)
  
  set.seed(seed); folds.cv <- caret::createFolds(y = data_sites$gridcode, 
                                                 k = k_nb_folds, 
                                                 list = F)
  data_sites$fold_for_cv <- folds.cv
  
  # Merge with initial data 
  data_for_cv <- data %>% 
    left_join(., data_sites %>% dplyr::select(-gridcode), 
              by = c("x", "y"))
  
  # Remove the sites from the dataset
  for(i in 1:length(unique(data_for_cv$fold_for_cv)))
  {
    
    # ---------------------
    # Define test and train datasets
    data_test <- data_for_cv[which(data_for_cv$fold_for_cv == i),]
    #dim(data_test)
    
    data_train <- data_for_cv[which(data_for_cv$fold_for_cv != i),]
    #dim(data_train)
    
    # --------------------------------------
    # PLSR: compute scores from train dataset and predict 
    # > monthly data 
    if(model %in% c("pls.m.1", "pls.m.2", "pls.m.3", "pls.m.all"))
    {
      
      message("refit PLSR scores for monthly data")
      # Fit PLSR model 
      plsr_train <- function_plsr2(type_data  = "M", 
                                   vars_names = vars_names,
                                   data       = data_train,
                                   nb_comp    = 7)
      
      # > list to store the results
      list_scores_pls <- list()
      
      # > Apply PLSR loads on new data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # > PLSR 
        plsr_j <- plsr_train$list_pls_per_variable[[paste0(clim.var_abb_j)]]
        
        # > Mean and sd of the original data to standardize 
        mu_train <- data_train %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
          colMeans(.)
        
        sd_train <- data_train %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
          apply(., 2, sd)
        
        # > New data 
        tab_test_j <- data_test %>% ungroup(.) %>%
          dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
        
        # > Scale new data based on mean and sd of the original data 
        z_tab_test_j <- scale(tab_test_j[,-1], center = mu_train, scale=sd_train)
        
        # > PLSR score from newdata
        tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_tab_test_j)
        
        colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_month_",clim.var_abb_j)
        
        # > Store
        list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
        
      }
      
      tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
      
      
      # Merge with data train
      data_train <- data_train %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., plsr_train$tab_PLS_scores)
      
      # Merge with data test 
      data_test <- data_test %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., tab_scores_plsr)
      
    }
    # > daily data 
    if(model %in% c("pls.d.1", "pls.d.2", "pls.d.3", "pls.d.all"))
    {
      
      message("refit PLSR scores for daily data")
      # Daily train and test data 
      data_day_train <- data_day %>% 
        map(., ~ {
          .x %>% filter(site_year %in% unique(data_train$site_year))  
        })
      
      data_day_test <- data_day %>% 
        map(., ~ {
          .x %>% filter(site_year %in% unique(data_test$site_year))  
        })
      
      
      # Fit PLSR model  
      plsr_train <- function_plsr2(type_data  = "D", 
                                   vars_names = vars_names,
                                   data       = data_day_train,
                                   nb_comp    = 214)
      
      
      # > list to store the results
      list_scores_pls <- list()
      
      # > Apply PCA loads on new data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # > PLSR 
        plsr_j <- plsr_train$list_pls_per_variable[[paste0(clim.var_abb_j)]]
        
        # > Select data for the var_j 
        data_day_train_j <- data_day_train[[paste0(var_j)]]
        data_day_test_j <- data_day_test[[paste0(var_j)]]
        
        # > Prepare for scaling 
        data_var_j <- data_day_train_j %>% 
          dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.)
        
        # > Mean and sd of the original data to standardize 
        mu_train <- data_var_j %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("day_", var_j))) %>% 
          colMeans(.)
        
        sd_train <- data_var_j %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("day_", var_j))) %>% 
          apply(., 2, sd)
        
        # > New data 
        tab_test_j <- data_day_test_j %>% 
          dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
          pivot_wider(names_from = c("clim.var", "day_of_year"), 
                      values_from = "cum_clim.value", 
                      names_prefix = "day_", 
                      names_sep = ".") %>% 
          ungroup(.)
        
        # > Scale new data based on mean and sd of the original data 
        z_tab_test_j <- scale(tab_test_j[,-1], center = mu_train, scale=sd_train)
        
        # > PLSR score from newdata
        tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_tab_test_j)
        
        colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_day_",clim.var_abb_j)
        
        # > Store
        list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
        
      }
      
      tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
      
      # Merge with data train
      data_train <- data_train %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., plsr_train$tab_PLS_scores)
      
      # Merge with data test 
      data_test <- data_test %>% 
        dplyr::select(-starts_with("PLS")) %>% 
        cbind(., tab_scores_plsr)
      
    }
    
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
    data_pred[[paste0(i)]] <- data.frame(Model     = model_name,
                                         site_year = data_test$site_year,
                                         Ya_obs    = data_test[,paste0(outcome)],
                                         Ya_pred   = pred.test_vec,
                                         N_predictors = N_predictors) %>% 
      dplyr::rename("Ya_obs" = 3)
    
  }
  
  # > predictions in a table
  preds <- plyr::ldply(data_pred, data.frame, .id = "fold_cv")
  
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
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    res <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(res)
  
}