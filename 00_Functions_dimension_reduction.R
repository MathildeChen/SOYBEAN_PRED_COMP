# -------------------------------------------------------------------------
# 
#       FUNCTIONS 
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
library(caret) ; library(ranger) ; library(fastshap)
# > functional analysis & partial least square
library(fda) ; library(MFPCA) ; library(pls) ; library(fda.usc)
# > others
library(hydroGOF) # for NSE computation

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
                               "surface_net_solar_radiation"="Solar radiation",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))

# ----------------------------------
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# ----------------------------------------
# Compute monthly averages from daily averages

monthly_average <- function(var_i, 
                            load = T, 
                            data_var_i = NULL, 
                            crop = NULL){ # crop = "maize" or "soybean"
  
  # > load data 
  if(load==TRUE)
  {
    if(is.null(crop) == T){ message("No crop is provided in the arguments") } 
    data <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/", crop, "/era5daily_data_", var_i, ".rda"))
  }
  
  # > or use the data provided
  if(load==FALSE)
  {
    data <- data_var_i
    
  }
  
  # > compute means for raw and accumulated variables
  data_m <- data %>% 
    # store year in another variable 
    #rename("year_Ya" = "year") %>% 
    # > retrieve month from date 
    mutate(month = month(date),
           year = year(date)) %>% 
    # > compute monthly averages for the variable (raw and accumulated)
    group_by(site_year, x, y, gridcode, country_name, year, month, clim.var, Ya, Ya_ano) %>% 
    summarise(mean_clim.value = mean(clim.value),
              mean_cum_clim.value = mean(cum_clim.value))
  
  return(data_m)
  
}

# ----------------------------------------
# Compute annual averages from daily averages

annual_average <- function(var_i, 
                           load = T, 
                           data_var_i = NULL, 
                           crop = NULL){ # crop = "maize" or "soybean"
  
  # > load data 
  if(load==TRUE)
  {
    if(is.null(crop) == T){ message("No crop is provided in the arguments") } 
    data <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/", crop, "/era5daily_data_", var_i, ".rda"))
  }
  # > or use the data provided
  if(load==FALSE)
  {
    data <- data_var_i
  }
  
  # > compute means for raw and accumulated variables
  data_a <- data %>% 
    # store year in another variable 
    rename("year_Ya" = "year") %>% 
    # > retrieve year from date 
    mutate(year = year(date)) %>% 
    # > compute annual averages for the variable (raw and accumulated)
    group_by(site_year, x, y, gridcode, country_name, year_Ya, year, clim.var, Ya, Ya_ano) %>% 
    summarise(mean_clim.value = mean(clim.value),
              mean_cum_clim.value = mean(cum_clim.value))
  
  return(data_a)
  
}

# ----------------------------------
# Functions to run PCA
# Argument: 
# - data_for_pca: 

# > Function to run PCA and store results 
do.pca <- function(data_for_pca){
  
  # > Check if data are standardized
  #   Mean of means should be = 0, mean of sd should be = 1
  mean_data <- mean(round(apply(data_for_pca, MARGIN = 2, FUN = mean)))
  sd_data <- mean(round(apply(data_for_pca, MARGIN = 2, FUN = sd)))
  testthat::expect_equal(mean_data, 0)
  testthat::expect_equal(sd_data, 1)
  
  # > Check if any NA
  #   There should be no na
  na_data <- sum(apply(is.na(data_for_pca),2,sum))
  testthat::expect_equal(na_data, 0)
  
  # > Compute PCA
  pca <- prcomp(x = data_for_pca, tol = 0.05)
  
  # > Get results from PCA
  # Eigen values
  pca.eig <- get_eigenvalue(pca)
  
  # Variables results
  res.var <- suppressWarnings(get_pca_var(pca))
  pca.var <- data.frame(
    coord = res.var$coord,        # Coordinates
    contrib = res.var$contrib,    # Axis contribution
    cos2 = res.var$cos2,          # Representation quality 
    cor = res.var$cor             # Correlation
  )
  
  # Individuals results
  res.ind <- get_pca_ind(pca)
  pca.ind <- data.frame(
    ind.coord = res.ind$coord,        # Coordinates
    ind.contrib = res.ind$contrib,    # Axis contribution
    ind.cos2 = res.ind$cos2           # Representation quality
  )
  
  # > Store results
  pca.res <- list(pca = pca,
                  pca.eig = pca.eig,
                  pca.var = pca.var,
                  pca.ind = pca.ind)
  
  return(pca.res)
  
}

# ----------------------------------
# PRINCIPAL COMPONENT ANALYSIS 
# - type_data   : "D" for daily averages, "M" for monthly averages
# - data        : table including climatic features for each site-year 
#                 if type_data = "M", a data.frame with 1 line per site-year and 1 column per climatic monthly averages 
#                 if type_data = "D", a list containing 1 data.frame per climatic feature, in which 1 line corresponds to 1 site-year*date (long-format table)
# - load_data   : should data need to be loaded? if load_data=TRUE, provide to data argument a data.frame containing the site-year ID for checking whether all site-year have data
# - vars_names  : table including the name, full name, and abbreviations of climatic features
# - cum_clim.var: do climatic data need to be accumulated over the growing season (typically for daily data)
# - scale       : do data need to be scaled (default=TRUE)

function_pca <- function(type_data, 
                         data,
                         load_data    = F,
                         vars_names   = vars_names,
                         cum_clim.var = T,
                         scale        = T){
   
  # -------------------------
  # ------------------------- 
  
  # > Object to store the PCA and outputs
   list_pca        <- list()
   list_scores_pca <- list()
   
   # -------------------------
   # -------------------------
   # > Monthly data 
   if(type_data == "M")
   {
     
     # > For each variable 
     for(var_i in unique(vars_names$clim.var))
     {
       
       # > variable abbreviation
       var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
       message(paste0("pca on ", var_i_abb))
       
       # -------------------------
       # > DATA SELECTION & PREPARATION FOR PCA
       # > select data for the considered variable
       data_var_i <- data %>%  
         dplyr::select(site_year, starts_with(paste0("monthly_", var_i)))
       
       # > check if each site-year has data
       testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
       
       #message("ok data")
       
       # > SCALE DATA 
       if(scale == TRUE)
       {
         # > select data for the considered variable and scale it
         z_tab_pca_var_i <- data_var_i %>% 
           dplyr::select(-site_year) %>%
           scale(.) 
       }
       
       if(scale == FALSE)
       {
         # > do not scale data (not recommanded)
         z_tab_pca_var_i <- data_var_i %>% 
           dplyr::select(-site_year) 
         warning("Warning: data used for PCA is not scaled")
         
       }
       # > rownames are sites-year
       rownames(z_tab_pca_var_i) <- data_var_i$site_year
       
       #message("ok scale data")
       
       # -------------------------
       # > PCA
       # > do pca on data 
       pca_var_i <- do.pca(z_tab_pca_var_i)
       
       #message("ok pca")
       
       # > extract scores derived from PCA
       var_i_pca_scores           <- pca_var_i$pca$x
       colnames(var_i_pca_scores) <- paste0("PC", 1:ncol(var_i_pca_scores), "_month_", var_i_abb)
       rownames(var_i_pca_scores) <- data_var_i$site_year
       
       #message("ok scores")
       
       # > STORE PCA AND SCORES  
       list_pca[[paste0(var_i_abb)]]        <- pca_var_i
       list_scores_pca[[paste0(var_i_abb)]] <- var_i_pca_scores
       
     }
     
     # > Extract scores for all variables 
     tab_PCA_scores <- do.call(cbind, list_scores_pca) %>% as.data.frame(.)
     
     # > Return the PCA for each variable and the scores 
     res <- list("tab_PCA_scores"        = tab_PCA_scores,
                 "list_pca_per_variable" = list_pca)
     
   }
   
   # -------------------------
   # -------------------------
   # > Daily data
   if(type_data == "D")
   {
     
     # > For each variable 
     for(var_i in unique(vars_names$clim.var)){
       
       # > variable abbreviation
       var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
       message(paste0("pca on ", var_i_abb))
       
       # -------------------------
       # > LOAD DATA (if load_data=T)
       if(load_data == T)
       {
         
         # > Load climatic data
         era5daily_var_i <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_daily/soybean/era5daily_data_", var_i, ".rda"))
         
         # > Check if each site-year has data
         testthat::expect_equal(sort(unique(era5daily_var_i$site_year)), sort(unique(data$site_year)))
         
       }
       if(load_data == F)
       {
         era5daily_var_i <- data[[paste0(var_i)]]
       }
        
       # -------------------------
       # > DATA SELECTION & PREPARATION FOR PCA
       # > select data for the considered variable
       # > if accumulated data over growing season is used
       if(cum_clim.var == TRUE)
       {
         data_var_i <- era5daily_var_i %>% 
           dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
           pivot_wider(names_from = c("clim.var", "day_of_year"), 
                       values_from = "cum_clim.value", 
                       names_prefix = "day_", 
                       names_sep = ".")
       }
       
       # > if raw data is used
       if(cum_clim.var == FALSE)
       {
         data_var_i <- era5daily_var_i %>% 
           dplyr::select(site_year, day_of_year, clim.var, clim.value) %>% 
           pivot_wider(names_from = c("clim.var", "day_of_year"), 
                       values_from = "clim.value", 
                       names_prefix = "day_", 
                       names_sep = ".")
         warning("Warning: data used for PCA is not accumulated over growing season")
       }
     
       # > SCALE DATA 
       if(scale == TRUE)
       {
         # > select data for the considered variable and scale it
         z_tab_pca_var_i <- data_var_i %>% 
           ungroup() %>% 
           dplyr::select(starts_with(paste0("day_"))) %>% 
           scale(.) 
       }
       
       if(scale == FALSE)
       {
         # > not scaled data (not recommended)
         z_tab_pca_var_i <- data_var_i %>% 
           ungroup() %>% 
           dplyr::select(starts_with(paste0("day_")))
         warning("Warning: data used for PCA is not scaled")
         
       }
       
       # > rownames are sites-year
       rownames(z_tab_pca_var_i) <- data_var_i$site_year
       
       # -------------------------
       # > PCA
       # > do pca on data 
       pca_var_i <- do.pca(z_tab_pca_var_i)
       
       # > extract scores derived from PCA
       var_i_pca_scores           <- pca_var_i$pca$x
       colnames(var_i_pca_scores) <- paste0("PC", 1:ncol(var_i_pca_scores), "_day_", var_i_abb)
       rownames(var_i_pca_scores) <- data_var_i$site_year
       
       # > STORE PCA AND SCORES  
       list_pca[[paste0(var_i_abb)]]        <- pca_var_i
       list_scores_pca[[paste0(var_i_abb)]] <- var_i_pca_scores
     
     }
     
     # > Extract scores for all variables 
     tab_PCA_scores <- do.call(cbind, list_scores_pca) %>% as.data.frame(.)
     
     # > Return the PCA for each variable and the scores 
     res <- list("tab_PCA_scores"        = tab_PCA_scores,
                 "list_pca_per_variable" = list_pca)
   }
   
   # -------------------------
   # -------------------------

   return(res)
   
}

# Function to compute new score 
newscores_pca <- function(type_data, 
                          vars_names, 
                          init_data,
                          new_data, 
                          init_data_day = NULL,
                          new_data_day = NULL){
  
  # > Check if data has PLS scores 
  test_name_init_data <- init_data %>% 
    dplyr::select(starts_with("PC"))
  test_name_new_data <- new_data %>% 
    dplyr::select(starts_with("PC"))
  
  if(ncol(test_name_init_data) != 0){ stop("Error: the initial data already contains PC scores.") }
  if(ncol(test_name_new_data) != 0){ stop("Error: the new data already contains PC scores.") }
  
  # ----- Monthly data ------
  if(type_data=="M")
  {
    
    # > PCA on initial data 
    init_pca <- function_pca(type_data  = "M", 
                             vars_names = vars_names,
                             data       = init_data)
    
    # > List to store the future scores 
    list_scores_pca <- list()
    
    # > Apply PCA loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > Loads from PCA from initial data for variable j 
      loads_pca_j <- init_pca$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
      
      # > Mean and sd of the original data to standardize new data from these values
      mu_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      new_data_j <- new_data %>% ungroup(.) %>%
        dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
      z_new_data_j <- scale(new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PCA scores from newdata
      list_scores_pca_j <- list()
      for(k in 1:ncol(loads_pca_j))
      {
        
        list_scores_pca_j[[paste0("PC", k, "_month_")]] <- 
          data.frame(site_year = new_data_j$site_year,
                     score = sapply(1:ncol(z_new_data_j), function(x) z_new_data_j[,x] * loads_pca_j[x,k] ) %>% apply(., 1, sum))
        
      }
      
      # > Store
      list_scores_pca[[paste0(clim.var_abb_j)]] <- plyr::ldply(list_scores_pca_j, data.frame, .id = "PC_id") 
      #colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_month_",clim.var_abb_j)
      
    }
    
    # > Scores for all variables 
    tab_scores_pca <- plyr::ldply(list_scores_pca, data.frame, .id="var_name") %>% 
      unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
      pivot_wider(names_from = c("PC_name"), values_from = "score")
    
    # > Merge with initial data 
    init_data_pca <- init_data  %>% 
      cbind(., init_pca$tab_PCA_scores)
    
    # > Merge with new data
    new_data_pca <- new_data %>% 
      cbind(., tab_scores_pca)
    
    # > Output 
    res <- list(init_data = init_data_pca, 
                new_data  = new_data_pca)
    
  }
  # ----- Daily data ------
  if(type_data=="D")
  {
    
    # > Detect data for days 
    if(is.null(init_data_day)==T){ stop("Error: No initial daily data provided (init_data_day = NULL)")}
    if(is.null(new_data_day)==T){  stop("Error: No new daily data provided (new_data_day = NULL)")}
    
    # > Daily train and test data 
    data_day_init <- init_data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(init_data$site_year))  
      })
    
    data_day_new <- new_data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(new_data$site_year))  
      })
    
    message("refit PCA scores for daily data")
    
    # > Apply PCA loads on new data 
    init_pca <- function_pca(type_data  = "D",
                             vars_names = vars_names,
                             data       = data_day_init)
    
    # > List to store the future scores 
    list_scores_pca <- list()
    
    # > Apply PLSR loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > Loads from PCA from initial data for variable j 
      loads_pca_j <- init_pca$list_pca_per_variable[[paste0(clim.var_abb_j)]]$pca$rotation
      
      # > Select daily data for the var_j 
      init_data_j <- data_day_init[[paste0(var_j)]]
      new_data_j  <- data_day_new[[paste0(var_j)]]
      
      # > Prepare data for scaling 
      tab_init_data_j <- init_data_j %>% 
        dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from   = c("clim.var", "day_of_year"), 
                    values_from  = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep    = ".") %>% 
        ungroup(.)
      
      # > Mean and sd of the original data to standardize new data 
      mu_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      tab_new_data_j <- new_data_j %>% 
        dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from   = c("clim.var", "day_of_year"), 
                    values_from  = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep    = ".") %>% 
        ungroup(.)
      z_tab_new_data_j <- scale(tab_new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PCA score from scaled newdata
      list_scores_pca_j <- list()
      for(k in 1:ncol(loads_pca_j))
      {
        
        list_scores_pca_j[[paste0("PC", k, "_day_")]] <- 
          data.frame(site_year = tab_new_data_j$site_year,
                     score = sapply(1:ncol(z_tab_new_data_j), function(x) z_tab_new_data_j[,x] * loads_pca_j[x,k] ) %>% apply(., 1, sum))
        
      }
      
      # > Store
      list_scores_pca[[paste0(clim.var_abb_j)]] <- plyr::ldply(list_scores_pca_j, data.frame, .id = "PC_id")
      
    }
    
    # > Scores for all variables
    tab_scores_pca <- plyr::ldply(list_scores_pca, data.frame, .id="var_name") %>% 
      unite(col = "PC_name", c("PC_id", "var_name"), sep = "", remove=T) %>% 
      pivot_wider(names_from = c("PC_name"), values_from = "score")
    
    # > Merge with initial data 
    init_data_pca <- init_data  %>% 
      cbind(., init_pca$tab_PCA_scores)
    
    # > Merge with new data
    new_data_pca <- new_data %>% 
      cbind(., tab_scores_pca)
    
    # > Output 
    res <- list(init_data = init_data_pca, 
                new_data  = new_data_pca)
    
  }
  
  # > Out
  return(res)
  
}


# ----------------------------------
# FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS 
# - type_data: "D" for daily averages, "M" for monthly averages
# - data        : table including climatic features for each site-year 
#                 if type_data = "M", a data.frame with 1 line per site-year and 1 column per climatic monthly averages 
#                 if type_data = "D", a list containing 1 data.frame per climatic feature, in which 1 line corresponds to 1 site-year*date (long-format table)
# - load_data   : should data need to be loaded? if load_data=TRUE, provide to data argument a data.frame containing the site-year ID for checking whether all site-year have data
# - vars_names  : table including the name, full name, and abbreviations of climatic features
# - cum_clim.var: do climatic data need to be accumulated over the growing season (typically for daily data)
# - scale       : do data need to be scaled (default=TRUE)
# - nb_comp     : number of harmonics to compute (4 is generally enought to capture 100% of the dataset variability)
# - basis       : family basis used to smooth discrete observation from functional observations
#                 if type_data = "M", basis = create.bspline.basis(rangeval = c(1,7),   nbasis = 7,   norder = 4)
#                 if type_data = "D", basis = create.bspline.basis(rangeval = c(1,214), nbasis = 150, norder = 4)

function_fpca <- function(type_data, 
                          data,
                          load_data    = F,
                          vars_names   = vars_names,
                          cum_clim.var = T,
                          scale        = T,
                          nb_comp      = 4,
                          basis){
   
   # -----------------
   # -----------------
   
   # > Object to store the FPCA and outputs
   list_fpca        <- list()
   list_scores_fpca <- list()
   
   # -----------------
   # -----------------
   # > Monthly data 
   if(type_data == "M")
   {
     
     # > For each variable 
     for(var_i in unique(vars_names$clim.var)){
       
       # > variable abbreviation
       var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
       message(paste0("fpca on ", var_i_abb))
       
       # -------------------------
       # > DATA SELECTION & PREPARATION FOR FPCA
       # > select data for the considered variable
       data_var_i <- data %>% 
         dplyr::select(site_year, starts_with(paste0("monthly_", var_i)))
       
       # > check if each site-year has data
       testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
       
       # > shaping data for FPCA
       tab_fpca_var_i <- data_var_i %>% 
         pivot_longer(cols = starts_with("monthly_"), 
                      names_to = "clim.var",
                      values_to = "clim.value") %>% 
         pivot_wider(names_from = "site_year", 
                     values_from = "clim.value") %>% 
         dplyr::select(-clim.var)
       
       # -------------------------
       # > SMOOTHED CLIMATIC OBSERVATIONS
       month_fd_var_i <- smooth.basis(y = as.matrix(tab_fpca_var_i), 
                                      fdParobj = basis,
                                      returnMatrix=TRUE)$fd
       
       # -------------------------
       # > FPCA ON SMOOTHED CLIMATIC OBSERVATIONS
       monthpcaobj <- pca.fd(month_fd_var_i, 
                             nharm=nb_comp, 
                             centerfns = TRUE)
       
       # > EXTRACT SCORES
       var_i_fpca_scores <- monthpcaobj$scores
       colnames(var_i_fpca_scores) <- paste0("FPC", 1:ncol(var_i_fpca_scores), "_month_", var_i_abb)
       rownames(var_i_fpca_scores) <- unique(data_var_i$site_year)
       
       # -------------------------
       # > STORE FPCA AND SCORES 
       list_fpca[[paste0(var_i_abb)]] <- monthpcaobj
       list_scores_fpca[[paste0(var_i_abb)]] <- var_i_fpca_scores
   
     }
     
     # > Extract scores for all variables 
     tab_FPCA_scores <- do.call(cbind, list_scores_fpca) %>% as.data.frame(.)
     
     # > Return the PCA for each variable and the scores 
     res <- list("tab_FPCA_scores"        = tab_FPCA_scores,
                 "list_fpca_per_variable" = list_fpca)
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
       message(paste0("fpca on ", var_i_abb))
       
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
       
       #message("data ok")
       
       # -------------------------
       # > DATA SELECTION & PREPARATION FOR FPCA
       # > select data for the considered variable
       # > if accumulated data over growing season is used
       if(cum_clim.var == TRUE)
       {
         data_var_i <- era5daily_var_i %>% 
           dplyr::select(site_year, cum_clim.value, day_of_year) %>% 
           pivot_wider(names_from = "site_year", 
                       values_from = "cum_clim.value") %>% 
           dplyr::select(-day_of_year)
       }
       
       # > if raw data is used
       if(cum_clim.var == FALSE)
       {
         data_var_i <- era5daily_var_i %>% 
           dplyr::select(site_year, clim.value, day_of_year) %>% 
           pivot_wider(names_from = "site_year", 
                       values_from = "clim.value") %>% 
           dplyr::select(-day_of_year)
       }
       
       #message("data scale ok")
       
       # -------------------------
       # > SMOOTHED CLIMATIC OBSERVATIONS
       day_fd_var_i <- smooth.basis(y = as.matrix(data_var_i), 
                                    fdParobj = basis,
                                    returnMatrix=TRUE)$fd
       
       #message("data smooth ok")
       
       # -------------------------
       # > FPCA
       daypcaobj <- pca.fd(day_fd_var_i, 
                           nharm=nb_comp, 
                           centerfns = TRUE)
       
       #message("FPCA ok")
       
       # > EXTRACT SCORES 
       var_i_fpca_scores <- daypcaobj$scores
       colnames(var_i_fpca_scores) <- paste0("FPC", 1:ncol(var_i_fpca_scores), "_day_", var_i_abb)
       rownames(var_i_fpca_scores) <- names(data_var_i)
       
       # -------------------------
       # > STORE FPCA AND SCORES 
       list_fpca[[paste0(var_i_abb)]]        <- daypcaobj
       list_scores_fpca[[paste0(var_i_abb)]] <- var_i_fpca_scores
       
       #message(var_i, paste0(" ok"))
       
     }
     
     # > Extract scores for all variables 
     tab_FPCA_scores <- do.call(cbind, list_scores_fpca) %>% as.data.frame(.)
     message("score tab ok")
     # > Return the FPCA for each variable and the scores 
     res <- list("tab_FPCA_scores"        = tab_FPCA_scores,
                 "list_fpca_per_variable" = list_fpca)
     
   }
   
   # -----------------
   # -----------------
   
   return(res)
   
 }


# Function to compute new score 
newscores_fpca <- function(type_data, 
                           vars_names, 
                           init_data,
                           init_fpca,
                           new_data, 
                           data_day=NULL){
  
  # > Check if data has FPCA scores 
    test_name_init_data <- init_data %>% 
      dplyr::select(starts_with("FPC"))
    test_name_new_data <- new_data %>% 
      ungroup(.) %>% 
      dplyr::select(starts_with("FPC"))
    
    if(ncol(test_name_init_data) != 0){ stop("Error: the initial data already contains FPC scores.") }
    if(ncol(test_name_new_data) != 0){ stop("Error: the new data already contains FPC scores.") }
    
    # ----- Monthly data ------
    if(type_data=="M")
    {
      
      
      message("refit FPCA scores for monthly data")
      
      # > FPCA on initial data 
      init_fpca <- init_fpca
      
      # > List to store the future scores 
      list_scores_fpca <- list()
      
      # > Project new data into splines based on the initial data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        # > Variable j 
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # Get the FPC
        harms = init_fpca$list_fpca_per_variable[[paste0(clim.var_abb_j)]]$harmonics
        # Get the mean curve
        meanfd = init_fpca$list_fpca_per_variable[[paste0(clim.var_abb_j)]]$meanfd
        
        # evaluate FPC 
        Stempvals = eval.fd(1:7,harms)
        # evaluate the mean curve 
        mtempvals = eval.fd(1:7,meanfd)
        
        Mdat = new_data %>% ungroup(.) %>%
          dplyr::select(starts_with(paste0("monthly_", var_j)))
        #length(Mdat); nrow(Mdat)
        
        list_scores_sites <- list()
        for(i in 1:nrow(Mdat))
        {
          # Remove Mean curve from the initial data
          Mdat2 = Mdat[i,]-mtempvals
          
          # Obtain the FPC scores of the new curves
          coef = lm(t(Mdat2)~Stempvals-1)$coef
          list_scores_sites[[paste0(i)]] <- coef
          
        }
        
        # > EXTRACT SCORES
        var_i_fpca_scores <- plyr::ldply(list_scores_sites, .id = "ID") %>% dplyr::select(-ID)
        colnames(var_i_fpca_scores) <- paste0("FPC", 1:ncol(var_i_fpca_scores), "_month_", clim.var_abb_j)
        
        # -------------------------
        # > STORE SCORES 
        list_scores_fpca[[paste0(clim.var_abb_j)]] <- var_i_fpca_scores
        
      }
      
      # > Extract scores for all variables 
      tab_scores_fpca <- map_dfc(list_scores_fpca, data.frame)
      
      # > Merge with initial data 
      init_data_fpca <- init_data #%>% cbind(init_fpca$tab_FPCA_scores)
      
      # > Merge with new data
      new_data_fpca <- new_data %>% 
        cbind(., tab_scores_fpca)
      
      # > Output 
      res <- list(init_data = init_data_fpca, 
                  new_data  = new_data_fpca)
    }
    
    # ----- Daily data ------ does not work properly :(
    if(type_data=="D")
    {
      
      message("refit FPCA scores for daily data")
      
      # > FPCA on initial data 
      init_fpca <- init_fpca
      
      # > List to store the future scores 
      list_scores_fpca <- list()
      
      # > Project new data into splines based on the initial data 
      for(var_j in unique(vars_names$clim.var)) 
      {
        
        # > Variable j 
        clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
        
        # Get the FPC
        harms = init_fpca$list_fpca_per_variable[[paste0(clim.var_abb_j)]]$harmonics
        # Get the mean curve
        meanfd = init_fpca$list_fpca_per_variable[[paste0(clim.var_abb_j)]]$meanfd
        
        # evaluate FPC 
        Stempvals = predict.fd(harms)
        # evaluate the mean curve 
        mtempvals = predict.fd(meanfd)
        
        # shape daily data 
        Mdat = data_day[[paste0(var_j)]] %>% ungroup(.) %>%
          dplyr::select(site_year, cum_clim.value, day_of_year) %>% 
          pivot_wider(names_from = "day_of_year", 
                      values_from = "cum_clim.value") %>% 
          dplyr::select(-site_year)
        #length(Mdat); nrow(Mdat)
        
        # predict scores for each curve
        list_scores_sites <- list()
        for(i in 1:nrow(Mdat))
        {
          # Remove Mean curve from the initial data
          Mdat2 = Mdat[i,]-mtempvals
          
          # Obtain the FPC scores of the new curves
          coef = lm(t(Mdat2)~Stempvals-1)$coef
          list_scores_sites[[paste0(i)]] <- coef
          
        }
        
        # > EXTRACT SCORES
        var_i_fpca_scores <- plyr::ldply(list_scores_sites, .id = "ID") %>% dplyr::select(-ID)
        colnames(var_i_fpca_scores) <- paste0("FPC", 1:ncol(var_i_fpca_scores), "_day_", clim.var_abb_j)
        
        # -------------------------
        # > STORE SCORES 
        list_scores_fpca[[paste0(clim.var_abb_j)]] <- var_i_fpca_scores
        
      }
      
      # > Extract scores for all variables 
      tab_scores_fpca <- map_dfc(list_scores_fpca, data.frame)
      
      # > Merge with initial data 
      init_data_fpca <- init_data %>%
        cbind(init_fpca$tab_FPCA_scores)
      
      # > Merge with new data
      new_data_fpca <- new_data %>% 
        cbind(., tab_scores_fpca)
      
      # > Output 
      res <- list(init_data = init_data_fpca, 
                  new_data  = new_data_fpca)
    }
    
    return(res)
}




#fpca_d_bra <- loadRDa(paste0(path_to_dim_red, "fpca_D_bra.rda"))
#
#list_test <- list_data_day_bra %>% 
#  map(., ~{ .x %>% filter(day_of_year < 213)})
#
#td <- newscores_fpca(type_data = "D", 
#                     vars_names = vars_names[1,],
#                     init_data = tab_bra %>% dplyr::select(-starts_with("FPC")),
#                     new_data  = tab_bra %>% dplyr::select(-starts_with("FPC")), 
#                     init_fpca = fpca_d_bra$list_fpca_per_variable, 
#                     data_day = list_test)
#
#summary(td$new_data_pca$FPC1_day_max_temp)
#summary(td$init_data_pca %>% 
#          filter(site_year %in% unique(tab_usa$site_year)[1:10]) %>% 
#          pull(FPC1_day_max_temp))
#
#summary(td$new_data_pca$FPC1_day_rad)
#summary(td$init_data_pca %>% 
#          filter(site_year %in% unique(tab_usa$site_year)[1:10]) %>% 
#          pull(FPC1_day_rad))
#
#
# ----------------------------------
# MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS 
# - type_data    : "D" for daily averages, "M" for monthly averages
# - data         : table including climatic features for each site-year 
#                  if type_data = "M", a data.frame with 1 line per site-year and 1 column per climatic monthly averages 
#                  if type_data = "D", a list containing 1 data.frame per climatic feature, in which 1 line corresponds to 1 site-year*date (long-format table)
# - load_data    : should data need to be loaded? if load_data=TRUE, provide to data argument a data.frame containing the site-year ID for checking whether all site-year have data
# - vars_names   : table including the name, full name, and abbreviations of climatic features
# - cum_clim.var : do climatic data need to be accumulated over the growing season (typically for daily data)
# - scale        : do data need to be scaled (default=TRUE)
# - nb_comp      : number of components computed
# - argvals      : timing of observations, 
#                  if type_data = "M", argsvals <- 1:7 
#                  if type_data = "D", argsvals <- 1:214 
#                  can also be entered as: argvals <- unique(data$day_of_year)  
# - univExpansion: see ?MFPCA::MFPCA

function_mfpca <- function(type_data, 
                           data,
                           load_data    = F,
                           vars_names,
                           cum_clim.var = T,
                           scale        = T,
                           nb_comp      = 4,
                           argvals,
                           univExpansion = list(list(type = "fda"),
                                                list(type = "fda"),
                                                list(type = "fda"),
                                                list(type = "fda"),
                                                list(type = "fda"),
                                                list(type = "fda"))){
  
  # -----------------
  # -----------------
  
  # > Monthly data 
  if(type_data == "M")
  {
    
    # > Object to store multivariate functional observations 
    list_fun_vars <-list()
    
    for(var_i in unique(vars_names$clim.var))
    {
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR MFPCA
      # > select data for the considered variable
      X <- data %>% 
        ungroup(.) %>% 
        dplyr::select(starts_with(paste0("monthly_", var_i))) %>%
        as.matrix(.)
      
      # > STORE MULTIVARIATE FUNCTIONAL OBSERVATIONS IN A LIST 
      list_fun_vars[[paste0(var_i_abb)]] <- funData(argvals = argvals, X = X)
      
    }
    
    # -------------------------
    # > MERGE ALL FUNCTIONAL OBSERVATIONS
    multi_fun_vars <- multiFunData(list_fun_vars)
    
    # > MFPCA BASED ON UNIVARTIATE FPCA
    MFPCA_m <- MFPCA(multi_fun_vars, 
                     M             = nb_comp, 
                     uniExpansions = univExpansion, 
                     fit           = T, 
                     approx.eigen  = T)
    
    # -------------------------
    # > SCORES 
    MFPCA_scores <- MFPCA_m$scores
    colnames(MFPCA_scores) <- paste0("MFPC", 1:nb_comp, "_month")
    tab_MFPCA_scores <- data.frame(site_year = unique(data$site_year), 
                                   MFPCA_scores) 
    
    # -------------------------
    # > Return the PCA for each variable and the scores 
    res <- list("tab_MFPCA_scores" = tab_MFPCA_scores,
                "list_mfpca"       = MFPCA_m)
  }
  
  # -----------------
  # -----------------
  # > Daily data 
  if(type_data == "D")
  {
    
    # -------------------------
    # > STORE MULTIVARIATE FUNCTIONAL OBSERVATIONS IN A LIST 
    list_fun_vars <- data %>% 
      map(., ~{
        
        # > Shape data
        tab_for_X <- .x %>%
          arrange(site_year, date) %>% 
          ungroup() %>%
          dplyr::select(site_year, day_of_year, cum_clim.value) 
        
        # > Observations
        X <- tab_for_X %>% 
          pivot_wider(names_from = "day_of_year", values_from = "cum_clim.value") %>% 
          dplyr::select(-site_year) %>% 
          as.matrix(.)
        
        # > funData
        funData(argvals = argvals,
                X = X)
        
      })
    
    # -------------------------
    # > MERGE ALL FUNCTIONAL OBSERVATIONS
    multi_fun_vars <- multiFunData(list_fun_vars)
    
    # -------------------------
    # > MFPCA BASED ON UNIVARTIATE FPCA
    MFPCA_d <- MFPCA(multi_fun_vars, 
                     M             = nb_comp, 
                     uniExpansions = univExpansion, 
                     fit           = T, 
                     approx.eigen  = T)
    
    # -------------------------
    # > SCORES 
    MFPCA_scores <- MFPCA_d$scores
    colnames(MFPCA_scores) <- paste0("MFPC", 1:nb_comp, "_day")
    tab_MFPCA_scores <- data.frame(site_year = unique(data[[1]]$site_year), 
                                   MFPCA_scores) 
    
    # -------------------------
    # > Return the PCA for each variable and the scores 
    res <- list("tab_MFPCA_scores" = tab_MFPCA_scores,
                "list_mfpca"       = MFPCA_d)
    
  }
  
  # -----------------
  # -----------------
  
  return(res)
}


# ----------------------------------
# PARTIAL LEAST SQUARE REGRESSION

function_plsr <- function(type_data, 
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
      
      plsr_var_i <- plsr(Ya ~ ., 
                         ncomp = nb_comp, 
                         data = tab_pls_var_i %>% 
                           dplyr::select(-site_year), 
                         validation = "CV")
      
      # > selection nb_comp based on cross-validation
      # The approach "onesigma" simply returns the first model where the optimal CV 
      # is within one standard error of the absolute optimum (Hastie, Tibshirani and Friedman, 2009).
      
      ncomp.onesigma <- selectNcomp(plsr_var_i, 
                                    method = "onesigma")
      
      plsr_var_i_cv <- plsr(Ya ~ ., 
                            ncomp = ncomp.onesigma, 
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
      
      plsr_var_i <- plsr(Ya ~ ., 
                         ncomp = nb_comp, 
                         data = tab_pls_var_i %>% 
                           dplyr::select(-site_year), 
                         validation = "CV") # > if scale==T, data is centered and scaled
      
      # > selection nb_comp based on cross-validation
      # The approach "onesigma" simply returns the first model where the optimal CV 
      # is within one standard error of the absolute optimum (Hastie, Tibshirani and Friedman, 2009).
      
      ncomp.onesigma <- selectNcomp(plsr_var_i, 
                                    method = "onesigma")
      
      plsr_var_i_cv <- plsr(Ya ~ ., 
                            ncomp = ncomp.onesigma, 
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

# New function without choosing nb of component by cv 
function_plsr2 <- function(type_data, 
                           data,
                           load_data    = F,
                           vars_names,
                           cum_clim.var = T,
                           scale        = T,
                           outcome      = NULL,
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
        dplyr::select(site_year, Ya, Ya_ano, starts_with(paste0("monthly_", var_i)))
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
      
      #message("- ok data")
      
      # > SCALE DATA 
      if(scale == TRUE)
      {
        # > select data for the considered variable and scale it
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya, -Ya_ano) %>%
          scale(.) %>% 
          as.data.frame(.)
      }
      
      if(scale == FALSE)
      {
        # > do not scale data (not recommanded)
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya, -Ya_ano) 
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
        dplyr::select(site_year, Ya, Ya_ano) %>% 
        left_join(., as.data.frame(z_tab_pls_var_i), by = "site_year") 
      
      # -----------------
      # > PLSR
      # on yield
      if(outcome == "Ya")
      {
        
        plsr_var_i_cv <- plsr(Ya ~ ., 
                              ncomp = nb_comp, 
                              data = tab_pls_var_i %>% 
                                dplyr::select(-site_year, -Ya_ano), 
                              validation = "CV")
        
      }
      # or on yield anomaly
      if(outcome == "Ya_ano")
      {
        
        plsr_var_i_cv <- plsr(Ya_ano ~ ., 
                              ncomp = nb_comp, 
                              data = tab_pls_var_i %>% 
                                dplyr::select(-site_year, -Ya), 
                              validation = "CV")
        
      }
      if(outcome != "Ya_ano" & outcome != "Ya")
      {
        stop("No outcome provided")
      }
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
    #message("score tab ok")
    
    # > Return the PLS for each variable and the scores 
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
          dplyr::select(site_year, Ya, Ya_ano, day_of_year, clim.var, cum_clim.value) %>% 
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
          dplyr::select(site_year, Ya, Ya_ano, day_of_year, clim.var, clim.value) %>% 
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
          dplyr::select(-site_year, -Ya, -Ya_ano) %>%
          scale(.) %>% 
          as.data.frame(.)
      }
      
      if(scale == FALSE)
      {
        # > do not scale data (not recommanded)
        z_tab_pls_var_i <- data_var_i %>% 
          dplyr::select(-site_year, -Ya, -Ya_ano)  %>% 
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
        dplyr::select(site_year, Ya, Ya_ano) %>% 
        left_join(., as.data.frame(z_tab_pls_var_i), by = "site_year")
      
      # -----------------
      # > PLSR
      # on yield
      if(outcome == "Ya")
      {
        
        plsr_var_i_cv <- plsr(Ya ~ ., 
                              ncomp = nb_comp, 
                              data = tab_pls_var_i %>% 
                                dplyr::select(-site_year, -Ya_ano), 
                              validation = "CV", 
                              scale = scale, center = scale)
        
      }
      # or on yield anomaly
      if(outcome == "Ya_ano")
      {
        
        plsr_var_i_cv <- plsr(Ya_ano ~ ., 
                              ncomp = nb_comp, 
                              data = tab_pls_var_i %>% 
                                dplyr::select(-site_year, -Ya), 
                              validation = "CV", 
                              scale = scale, center = scale)
        
      }
      if(outcome != "Ya_ano" & outcome != "Ya")
      {
        stop("No outcome provided")
      }
      #message("plsr ok")
      
      # -----------------
      # > store scores from PLSR
      var_i_plsr_scores <- plsr_var_i_cv$scores 
      colnames(var_i_plsr_scores) <- paste0("PLS", 1:ncol(var_i_plsr_scores), "_day_", var_i_abb)
      #rownames(var_i_plsr_scores) <- tab_pls_var_i$site_year
      
      # -------------------------
      # > STORE PLS AND SCORES 
      list_pls_scores[[paste0(var_i_abb)]] <- var_i_plsr_scores
      list_pls[[paste0(var_i_abb)]]        <- plsr_var_i_cv
      
    }
    
    # > Extract scores for all variables 
    tab_PLS_scores <- do.call(cbind, list_pls_scores) %>% as.data.frame(.)
    #message("score tab ok")
    
    # > Return the PLS for each variable and the scores 
    res <- list("tab_PLS_scores"        = tab_PLS_scores,
                "list_pls_per_variable" = list_pls)
    
  }
  
  return(res)
  
}

# Function to compute new score 
newscores_plsr2 <- function(type_data, 
                            vars_names, 
                            init_data,
                            new_data, 
                            init_data_day = NULL, 
                            new_data_day = NULL, 
                            outcome = NULL){
  
  # > Check if data has PLS scores 
  test_name_init_data <- init_data %>% 
    dplyr::select(starts_with("PLS"))
  test_name_new_data <- new_data %>% 
    dplyr::select(starts_with("PLS"))
  
  if(ncol(test_name_init_data) != 0){ stop("Error: the initial data already contains PLS scores.") }
  if(ncol(test_name_new_data) != 0){ stop("Error: the new data already contains PLS scores.") }
  
  # ----- Monthly data ------
  if(type_data=="M")
  {
    # > Fit PLSR model on initial data 
    init_plsr <- function_plsr2(type_data  = "M", 
                                vars_names = vars_names,
                                data       = init_data,
                                nb_comp    = 7, 
                                outcome    = outcome)
    
    # > List to store the future scores 
    list_scores_pls <- list()
    
    # > Apply PLSR loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > PLSR from initial data for variable j 
      plsr_j <- init_plsr$list_pls_per_variable[[paste0(clim.var_abb_j)]]
      
      # > Mean and sd of the original data to standardize new data from these values
      mu_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- init_data %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("monthly_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      new_data_j <- new_data %>% ungroup(.) %>%
        dplyr::select(site_year, starts_with(paste0("monthly_", var_j)))
      z_new_data_j <- scale(new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PLSR score from newdata
      tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_new_data_j)
      colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_month_",clim.var_abb_j)
      
      # > Store
      list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
      
    }
    
    # > Scores for all variables 
    tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
    
    # > Merge with initial data 
    init_data_pls <- init_data  %>% 
      cbind(., init_plsr$tab_PLS_scores)
    
    # > Merge with new data
    new_data_pls <- new_data %>% 
      cbind(., tab_scores_plsr)
    
    # > Output 
    res <- list(init_data = init_data_pls, 
                new_data  = new_data_pls)
    
  }
  # ----- Daily data ------
  if(type_data=="D")
  {
    
    # > Detect data for days 
    if(is.null(init_data_day)==T){ stop("Error: No init daily data provided (init_data_day = NULL)")}
    if(is.null(new_data_day)==T){ stop("Error: No new daily data provided (new_data_day = NULL)")}
    
    # > Daily train and test data 
    data_day_init <- init_data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(init_data$site_year))  
      })
    
    data_day_new <- new_data_day %>% 
      map(., ~ {
        .x %>% filter(site_year %in% unique(new_data$site_year))  
      })
    
    #message("refit PLSR scores for daily data")
    
    # > Fit PLSR model on initial data 
    init_plsr <- function_plsr2(type_data  = "D",
                                vars_names = vars_names,
                                data       = data_day_init,
                                nb_comp    = 20, 
                                outcome    = outcome)
    
    # > List to store the future scores 
    list_scores_pls <- list()
    
    # > Apply PLSR loads on new data 
    for(var_j in unique(vars_names$clim.var)) 
    {
      
      # > Variable j 
      clim.var_abb_j <- vars_names[which(vars_names$clim.var==var_j),]$clim.var_abb
      
      # > PLSR from initial data for variable j 
      plsr_j <- init_plsr$list_pls_per_variable[[paste0(clim.var_abb_j)]]
      
      # > Select daily data for the var_j 
      init_data_j <- data_day_init[[paste0(var_j)]]
      new_data_j  <- data_day_new[[paste0(var_j)]]
      
      # > Prepare data for scaling 
      tab_init_data_j <- init_data_j %>% 
        dplyr::select(site_year, Ya, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from = c("clim.var", "day_of_year"), 
                    values_from = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep = ".") %>% 
        ungroup(.)
      
      # > Mean and sd of the original data to standardize new data 
      mu_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        colMeans(.)
      sd_init_data <- tab_init_data_j %>% ungroup(.) %>%
        dplyr::select(starts_with(paste0("day_", var_j))) %>% 
        apply(., 2, sd)
      
      # > Scale new data based on mean and sd of the original data 
      tab_new_data_j <- new_data_j %>% 
        dplyr::select(site_year, day_of_year, clim.var, cum_clim.value) %>% 
        pivot_wider(names_from = c("clim.var", "day_of_year"), 
                    values_from = "cum_clim.value", 
                    names_prefix = "day_", 
                    names_sep = ".") %>% 
        ungroup(.)
      z_tab_new_data_j <- scale(tab_new_data_j[,-1], center = mu_init_data, scale=sd_init_data)
      
      # > PLSR score from scaled newdata
      tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=z_tab_new_data_j)
      colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_day_", clim.var_abb_j)
      
      # > Store
      list_scores_pls[[paste0(clim.var_abb_j)]] <- tab_scores_plsr_j
      
    }
    
    # > Scores for all variables
    tab_scores_plsr <- map_dfc(list_scores_pls, data.frame)
    
    # > Merge with initial data 
    init_data_pls <- init_data  %>% 
      cbind(., init_plsr$tab_PLS_scores)
    
    # > Merge with new data
    new_data_pls <- new_data %>% 
      cbind(., tab_scores_plsr)
    
    # > Output 
    res <- list(init_data = init_data_pls, 
                new_data  = new_data_pls)
    
  }
  
  # > Out
  return(res)
  
}

# ----------------------------------
# FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION

function_fplsr <- function(type_data, 
                           data,
                           load_data    = F,
                           vars_names,
                           cum_clim.var = T,
                           scale        = T,
                           nb_comp,
                           basis){
  
  # -----------------
  # -----------------
  
  # > Object to store the PLS and outputs
  list_fpls <- list()
  list_fpls_scores <- list()
  
  # -----------------
  # -----------------
  # > Monthly data 
  if(type_data == "M")
  {
    
    # > For each variable 
    for(var_i in unique(vars_names$clim.var)){
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("fplsr on ", var_i_abb))
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR FPLSR
      # > select data for the considered variable
      data_var_i <- data %>% 
        dplyr::select(site_year, Ya, starts_with(paste0("monthly_", var_i)))
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
      
      # > select Y (yield)
      y <- data_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        distinct() %>% 
        pull(Ya)
      
      # > shaping data for FPLSR
      tab_fplsr_var_i <- data_var_i %>% 
        dplyr::select(-Ya) %>% 
        pivot_longer(cols = starts_with("monthly_"), 
                     names_to = "clim.var",
                     values_to = "clim.value") %>% 
        pivot_wider(names_from = "site_year", 
                    values_from = "clim.value") %>% 
        dplyr::select(-clim.var)
      
      # -------------------------
      # > SMOOTHED CLIMATIC OBSERVATIONS
      fd_var_i <- smooth.basis(y = as.matrix(tab_fplsr_var_i), 
                               fdParobj = basis,
                               returnMatrix=TRUE)$fd
      
      # > transform fd into fdata for fdata2pls function
      x <- fdata(fd_var_i)
      
      # -------------------------
      # > SCALE SMOOTHED ESTIMATED (asked in the fdata2pls function help)
      if(scale == TRUE){
        
        # > center data 
        z_x <- fdata.cen(x)$Xcen
        
      }
      if(scale == FALSE){
        
        # > not center data 
        z_x <- x
        warning("Warning: data used for FPLSR is not scaled")
      }
      
      # -------------------------
      # > FPLSR ON SMOOTHED CLIMATIC OBSERVATIONS
      
      # > optimal hyperparameters (ncomp + penalization)
      ncomp.cv <- fregre.pls.cv(fdataobj = z_x, y = y, 
                                kmax = nb_comp, 
                                lambda=TRUE, 
                                P=c(0,0,1))
      
      # > Fit FPLSR model using optimal set of hyperparameters
      fplsobj <- fdata2pls(fdataobj = z_x, y = y, 
                           ncomp = max(ncomp.cv$pls.opt), 
                           lambda = ncomp.cv$lambda.opt)
      
      # > EXTRACT SCORES
      var_i_fpls_scores <- as.data.frame(fplsobj$x)
      colnames(var_i_fpls_scores) <- paste0("FPLS", 1:ncol(var_i_fpls_scores), "_month_", var_i_abb)
      rownames(var_i_fpls_scores) <- unique(data_var_i$site_year)
      
      # -------------------------
      # > STORE FPLS AND SCORES 
      list_fpls[[paste0(var_i_abb)]] <- fplsobj
      list_fpls_scores[[paste0(var_i_abb)]] <- var_i_fpls_scores
      
    }
    
    # > Extract scores for all variables 
    tab_FPLS_scores <- do.call(cbind, list_fpls_scores) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_FPLS_scores"        = tab_FPLS_scores,
                "list_fpls_per_variable" = list_fpls)
    
    
    
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
      message(paste0("fpls on ", var_i_abb))
      
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
      
      #message("data ok")
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR FPLS
      
      # > select Y (yield)
      y <- era5daily_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        distinct() %>% 
        pull(Ya)
      
      # > shaping data for FPLSR
      
      # > select data for the considered variable
      # > if accumulated data over growing season is used
      if(cum_clim.var == TRUE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, cum_clim.value, day_of_year) %>% 
          pivot_wider(names_from = "site_year", 
                      values_from = "cum_clim.value") %>% 
          dplyr::select(-day_of_year)
      }
      
      # > if raw data is used
      if(cum_clim.var == FALSE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, clim.value, day_of_year) %>% 
          pivot_wider(names_from = "site_year", 
                      values_from = "clim.value") %>% 
          dplyr::select(-day_of_year)
      }
      
      #message("data scale ok")
      
      # -------------------------
      # > SMOOTHED CLIMATIC OBSERVATIONS
      day_fd_var_i <- smooth.basis(y = as.matrix(data_var_i), 
                                   fdParobj = basis,
                                   returnMatrix=TRUE)$fd
      
      #message("data smooth ok")
      
      # > transform fd into fdata for fdata2pls function
      x <- fdata(day_fd_var_i)
      
      # -------------------------
      # > SCALE SMOOTHED ESTIMATED (asked in the fdata2pls function help)
      if(scale == TRUE){
        
        # > center data 
        z_x <- fdata.cen(x)$Xcen
        
      }
      if(scale == FALSE){
        
        # > not center data 
        z_x <- x
        warning("Warning: data used for FPLSR is not scaled")
      }
      
      # -------------------------
      # > FPLSR ON SMOOTHED CLIMATIC OBSERVATIONS
      
      # > optimal hyperparameters (ncomp + penalization)
      ncomp.cv <- fregre.pls.cv(fdataobj = z_x, y = y, 
                                kmax = nb_comp, 
                                lambda=TRUE, 
                                P=c(0,0,1))
      
      # > Fit FPLSR model using optimal set of hyperparameters
      fplsobj <- fdata2pls(fdataobj = z_x, y = y, 
                           ncomp = max(ncomp.cv$pls.opt), 
                           lambda = ncomp.cv$lambda.opt)
      
      # > EXTRACT SCORES
      var_i_fpls_scores <- as.data.frame(fplsobj$x)
      colnames(var_i_fpls_scores) <- paste0("FPLS", 1:ncol(var_i_fpls_scores), "_day_", var_i_abb)
      rownames(var_i_fpls_scores) <- names(data_var_i)
      
      # -------------------------
      # > STORE FPLS AND SCORES 
      list_fpls[[paste0(var_i_abb)]] <- fplsobj
      list_fpls_scores[[paste0(var_i_abb)]] <- var_i_fpls_scores
      
      
    }
    
    # > Extract scores for all variables 
    tab_FPLS_scores <- do.call(cbind, list_fpls_scores) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_FPLS_scores"        = tab_FPLS_scores,
                "list_fpls_per_variable" = list_fpls)
    
  }

  # -----------------
  # -----------------
  
  return(res)
  
}


# ----------------------------------
# FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION


function_fplsr2 <- function(type_data, 
                           data,
                           load_data    = F,
                           vars_names,
                           cum_clim.var = T,
                           scale        = T,
                           nb_comp,
                           basis){
  
  # -----------------
  # -----------------
  
  # > Object to store the PLS and outputs
  list_fpls <- list()
  list_fpls_scores <- list()
  
  # -----------------
  # -----------------
  # > Monthly data 
  if(type_data == "M")
  {
    
    # > For each variable 
    for(var_i in unique(vars_names$clim.var)){
      
      # > variable abbreviation
      var_i_abb <- vars_names[which(vars_names$clim.var==var_i), "clim.var_abb"]
      message(paste0("fplsr on ", var_i_abb))
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR FPLSR
      # > select data for the considered variable
      data_var_i <- data %>% 
        dplyr::select(site_year, Ya, starts_with(paste0("monthly_", var_i)))
      
      # > check if each site-year has data
      testthat::expect_equal(sort(unique(data_var_i$site_year)), sort(unique(data$site_year)))
      
      # > select Y (yield)
      y <- data_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        distinct() %>% 
        pull(Ya)
      
      # > shaping data for FPLSR
      tab_fplsr_var_i <- data_var_i %>% 
        dplyr::select(-Ya) %>% 
        pivot_longer(cols = starts_with("monthly_"), 
                     names_to = "clim.var",
                     values_to = "clim.value") %>% 
        pivot_wider(names_from = "site_year", 
                    values_from = "clim.value") %>% 
        dplyr::select(-clim.var)
      
      # -------------------------
      # > SMOOTHED CLIMATIC OBSERVATIONS
      fd_var_i <- smooth.basis(y = as.matrix(tab_fplsr_var_i), 
                               fdParobj = basis,
                               returnMatrix=TRUE)$fd
      
      # > transform fd into fdata for fdata2pls function
      x <- fdata(fd_var_i)
      
      # -------------------------
      # > SCALE SMOOTHED ESTIMATED (asked in the fdata2pls function help)
      if(scale == TRUE){
        
        # > center data 
        z_x <- fdata.cen(x)$Xcen
        
      }
      if(scale == FALSE){
        
        # > not center data 
        z_x <- x
        warning("Warning: data used for FPLSR is not scaled")
      }
      
      # -------------------------
      # > FPLSR ON SMOOTHED CLIMATIC OBSERVATIONS
      
      # > Fit FPLSR model
      fplsobj <- fdata2pls(fdataobj = z_x, y = y, 
                           ncomp = nb_comp, 
                           lambda = 0)
      
      # > EXTRACT SCORES
      var_i_fpls_scores <- as.data.frame(fplsobj$x)
      colnames(var_i_fpls_scores) <- paste0("FPLS", 1:ncol(var_i_fpls_scores), "_month_", var_i_abb)
      rownames(var_i_fpls_scores) <- unique(data_var_i$site_year)
      
      # -------------------------
      # > STORE FPLS AND SCORES 
      list_fpls[[paste0(var_i_abb)]] <- fplsobj
      list_fpls_scores[[paste0(var_i_abb)]] <- var_i_fpls_scores
      
    }
    
    # > Extract scores for all variables 
    tab_FPLS_scores <- do.call(cbind, list_fpls_scores) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_FPLS_scores"        = tab_FPLS_scores,
                "list_fpls_per_variable" = list_fpls)
    
    
    
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
      message(paste0("fpls on ", var_i_abb))
      
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
      
      #message("data ok")
      
      # -------------------------
      # > DATA SELECTION & PREPARATION FOR FPLS
      
      # > select Y (yield)
      y <- era5daily_var_i %>% 
        dplyr::select(site_year, Ya) %>% 
        distinct() %>% 
        pull(Ya)
      
      # > shaping data for FPLSR
      
      # > select data for the considered variable
      # > if accumulated data over growing season is used
      if(cum_clim.var == TRUE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, cum_clim.value, day_of_year) %>% 
          pivot_wider(names_from = "site_year", 
                      values_from = "cum_clim.value") %>% 
          dplyr::select(-day_of_year)
      }
      
      # > if raw data is used
      if(cum_clim.var == FALSE)
      {
        data_var_i <- era5daily_var_i %>% 
          dplyr::select(site_year, clim.value, day_of_year) %>% 
          pivot_wider(names_from = "site_year", 
                      values_from = "clim.value") %>% 
          dplyr::select(-day_of_year)
      }
      
      #message("data scale ok")
      
      # -------------------------
      # > SMOOTHED CLIMATIC OBSERVATIONS
      day_fd_var_i <- smooth.basis(y = as.matrix(data_var_i), 
                                   fdParobj = basis,
                                   returnMatrix=TRUE)$fd
      
      #message("data smooth ok")
      
      # > transform fd into fdata for fdata2pls function
      x <- fdata(day_fd_var_i)
      
      # -------------------------
      # > SCALE SMOOTHED ESTIMATED (asked in the fdata2pls function help)
      if(scale == TRUE){
        
        # > center data 
        z_x <- fdata.cen(x)$Xcen
        
      }
      if(scale == FALSE){
        
        # > not center data 
        z_x <- x
        warning("Warning: data used for FPLSR is not scaled")
      }
      
      # -------------------------
      # > FPLSR ON SMOOTHED CLIMATIC OBSERVATIONS
      
      # > Fit FPLSR model using optimal set of hyperparameters
      fplsobj <- fdata2pls(fdataobj = z_x, y = y, 
                           ncomp = nb_comp, 
                           lambda = 0)
      
      # > EXTRACT SCORES
      var_i_fpls_scores <- as.data.frame(fplsobj$x)
      colnames(var_i_fpls_scores) <- paste0("FPLS", 1:ncol(var_i_fpls_scores), "_day_", var_i_abb)
      rownames(var_i_fpls_scores) <- names(data_var_i)
      
      # -------------------------
      # > STORE FPLS AND SCORES 
      list_fpls[[paste0(var_i_abb)]] <- fplsobj
      list_fpls_scores[[paste0(var_i_abb)]] <- var_i_fpls_scores
      
      
    }
    
    # > Extract scores for all variables 
    tab_FPLS_scores <- do.call(cbind, list_fpls_scores) %>% as.data.frame(.)
    
    # > Return the PCA for each variable and the scores 
    res <- list("tab_FPLS_scores"        = tab_FPLS_scores,
                "list_fpls_per_variable" = list_fpls)
    
  }
  
  # -----------------
  # -----------------
  
  return(res)
  
}


# ----------------------------------
# WRAPPER FOR COMPUTING NEW SCORES 
newscores_wraper <- function(model_init    = NULL, 
                             init_data     = NULL, 
                             pred_data     = NULL, 
                             init_data_day = NULL, 
                             pred_data_day = NULL,
                             outcome       = NULL, 
                             vars_names    = NULL){
  
  # Checks
  if(is.null(model_init)) { stop("No model provided, fill 'model_init' argument") }
  if(is.null(init_data) | is.null(pred_data)){ stop("No initial or predicted data provided, fill arguments 'init_data' or 'pred_data'") }
  if(is.null(vars_names)){
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
                                   "surface_net_solar_radiation"="Solar radiation",
                                   "total_precipitation"        ="Precipitation",
                                   "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))
  }
  
  # --------------------------------------
  # Compute scores from train dataset and predict for test data set 
  # --------------------------------------
  # PCA: 
  # > monthly data
  if(model_init %in% c("pca.m.1", "pca.m.2", "pca.m.3", "pca.m.all"))
  {
    
    message("refit PCA scores for monthly data")
    
    # > Compute new scores 
    new_scores_pca_M <- newscores_pca(type_data = "M", 
                                      vars_names = vars_names, 
                                      init_data = init_data %>% dplyr::select(-starts_with("PC")),
                                      new_data =  pred_data  %>% dplyr::select(-starts_with("PC")), 
                                      data_day = NULL)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_pca_M$init_data
    new_pred_data  <- new_scores_pca_M$new_data
    
  }
  # > daily data
  if(model_init %in% c("pca.d.1", "pca.d.2", "pca.d.3", "pca.d.all"))
  {
    
    message("refit PCA scores for daily data")
    
    # > Compute new scores 
    new_scores_pca_D <- newscores_pca(type_data  = "M", 
                                      vars_names = vars_names, 
                                      init_data  = init_data %>% dplyr::select(-starts_with("PC")),
                                      new_data   = pred_data %>% dplyr::select(-starts_with("PC")),
                                      data_day   = init_data_day)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_pca_D$init_data
    new_pred_data  <- new_scores_pca_D$new_data
    
  }
  
  # --------------------------------------
  # PLSR: compute scores from train dataset and predict for test data set 
  # > monthly data 
  if(model_init %in% c("pls.m.1", "pls.m.2", "pls.m.3", "pls.m.all"))
  {
    
    if(is.null(outcome)) { stop("No outcome provided, fill 'outcome' argument (either Ya or Ya_ano)") }
    
    message("refit PLSR scores for monthly data")
    
    # > Compute new scores 
    new_scores_plsr2_M <- newscores_plsr2(type_data     = "M", 
                                          vars_names    = vars_names, 
                                          init_data     = init_data %>% dplyr::select(-starts_with("PLS")),
                                          new_data      = pred_data %>% dplyr::select(-starts_with("PLS")), 
                                          init_data_day = NULL, 
                                          new_data_day  = NULL, 
                                          outcome       = outcome)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_plsr2_M$init_data
    new_pred_data  <- new_scores_plsr2_M$new_data
    
  }
  # > daily data 
  if(model_init %in% c("pls.d.1", "pls.d.2", "pls.d.3", "pls.d.all"))
  {
    
    if(is.null(outcome)) { stop("No outcome provided, fill 'outcome' argument (either Ya or Ya_ano)") }
    
    message("refit PLSR scores for daily data")
    
    # > Compute new scores
    new_scores_plsr2_D <- newscores_plsr2(type_data     = "D", 
                                          vars_names    = vars_names, 
                                          init_data     = init_data %>% dplyr::select(-starts_with("PLS")),
                                          new_data      = pred_data %>% dplyr::select(-starts_with("PLS")),
                                          init_data_day = init_data_day, 
                                          new_data_day  = pred_data_day, 
                                          outcome       = outcome)
    
    # > Data train and test with new scores
    new_init_data  <- new_scores_plsr2_D$init_data
    new_pred_data  <- new_scores_plsr2_D$new_data
    
  }
  
  # Output
  res <- list(new_init_data = new_init_data, 
              new_pred_data = new_pred_data)
  return(res)
  
  
}


# ----------------------------------
# CROSS VALIDATION ON YEARS
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"
function_cv_year <- function(model      = NULL, 
                             outcome    = NULL,
                             model_gpe  = "rf",
                             recompute_scores = FALSE,
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
    
    # Name of the model 
    model_init <- substr(model_name, 4, nchar(model_name))
    
    # ---------------------
    # Define test and train datasets
    data_train_0 <- data[which(data$year != y),]
    data_test_0 <- data[which(data$year == y),]
    
    # --------------------------------------
    # If necessary, compute scores from train dataset and predict for test data set 
    if(recompute_scores == TRUE)
    {
      
      new_scores <- newscores_wraper(model_init = model_init, 
                                     init_data  = data_train_0, 
                                     pred_data  = data_test_0, 
                                     data_day   = data_day, 
                                     outcome    = outcome)
      
      # Test and train datasets with new scores
      data_train <- new_scores$new_init_data
      data_test  <- new_scores$new_pred_data
    }
    if(recompute_scores == FALSE)
    {
      
      # Test and train datasets provided 
      data_train <- data_train_0
      data_test  <- data_test_0
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
    out <- preds
  }
  # > return prediction performance indicators
  if(res == "perf")
  {
    # > performance indicators
    RMSEP <- caret::RMSE(obs = preds$Ya_obs,      pred = preds$Ya_pred)
    NSE   <- hydroGOF::NSE(obs = preds$Ya_obs,    sim=preds$Ya_pred)
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    out <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(out)
  
}


# ----------------------------------
# CROSS VALIDATION ON SITES
# res: "preds" or "perf"
# model_gpe: "rf or "lm",
# outcome: name of the predicted variable "Ya" or "Ya_ano"

function_cv_geo <- function(model, 
                            outcome,
                            model_gpe  = "rf",
                            recompute_scores = FALSE,
                            k_nb_folds=10,
                            data, 
                            data_day,
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
    
    # --------------------------------------
    # Name of the model 
    model_init <- substr(model_name, 4, nchar(model_name))
    
    # ---------------------
    # Define test and train datasets
    data_train_0 <- data_for_cv[which(data_for_cv$fold_for_cv != i),]
    data_test_0  <- data_for_cv[which(data_for_cv$fold_for_cv == i),]
    
    # --------------------------------------
    # If necessary, compute scores from train dataset and predict for test data set 
    if(recompute_scores == TRUE)
    {
      
      new_scores <- newscores_wraper(model_init = model_init, 
                                     init_data  = data_train_0, 
                                     pred_data  = data_test_0, 
                                     data_day   = data_day, 
                                     outcome    = outcome)
      
      # Test and train datasets with new scores
      data_train <- new_scores$new_init_data
      data_test  <- new_scores$new_pred_data
    }
    if(recompute_scores == FALSE)
    {
      
      # Test and train datasets provided 
      data_train <- data_train_0
      data_test  <- data_test_0
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
    out <- preds
  }
  # > return prediction performance indicators
  if(res == "perf")
  {
    # > performance indicators
    RMSEP <- caret::RMSE(obs = preds$Ya_obs,      pred = preds$Ya_pred)
    NSE   <- hydroGOF::NSE(obs = preds$Ya_obs,    sim=preds$Ya_pred)
    #R2    <- caret::R2(obs = preds$Ya_obs,        pred = preds$Ya_pred)
    #Bias  <- Metrics::bias(actual = preds$Ya_obs, predicted = preds$Ya_pred)
    
    out <- data.frame("RMSEP"=RMSEP, "NSE"=NSE, #"R2"=R2, "Bias"=Bias, 
                      "N_predictors"=N_predictors)
  }
  
  return(out)
  
}

# -----------------------------------------------------------------------------
# In case of bootstrap procedure
# > Function to use in boot function
bs <- function(data,      # table with unique combinaison of longitude and latitude
               indices,   # indices to select sites
               tab_test,  # table with full climatic data for each site-year
               #seed=102, num.trees = 500,
               model_formula){
  
  # > Data for bootstrap
  sites_b <- data[indices,]
  data_b  <- tab_test %>% 
    filter(gridcode %in% unique(sites_b$gridcode))
  
  # ------------------
  # > Year-by-year cross-validation
  mod_bs_cv <- function_cv_year(model         = model_formula, 
                                data          = data_b, 
                                res           = "perf", 
                                save = F)
  
  # > Output
  res.bs <- c(mod_bs_cv$RMSEP, mod_bs_cv$NSE, mod_bs_cv$R2, mod_bs_cv$Bias)
  return(res.bs)  
}

