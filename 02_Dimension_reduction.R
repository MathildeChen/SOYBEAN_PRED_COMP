# -------------------------------------------------------------------------
# 
#           DIMENSION REDUCTION IN CLIMATIC PREDICTORS 
#           UNITED-STATES, BRAZIL, FULL DATASET
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
library(fda) ; library(MFPCA)
# > others
library(hydroGOF)

# -------------------------------------------------------------------------
# Home-made functions performing the dimension reductions
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# Functions features 
# > basis monthly
basis_month <- create.bspline.basis(rangeval = c(1,7), nbasis = 7, norder = 4)
# > basis day
basis_day  <- create.bspline.basis(rangeval = c(1,214), nbasis = 150, norder = 4)
basis_day2 <- create.bspline.basis(rangeval = c(1,212), nbasis = 150, norder = 4)
# >
univExpansion_6_variables = list(list(type = "fda"),
                                 list(type = "fda"),
                                 list(type = "fda"),
                                 list(type = "fda"),
                                 list(type = "fda"),
                                 list(type = "fda"))

# ----------------------------------
# Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > daily climatic predictors 
load(paste0(data_path, "01_days/list_data_day_usa.rda"))
load(paste0(data_path, "01_days/list_data_day_bra.rda"))
load(paste0(data_path, "01_days/list_data_day_world.rda"))

# > monthly climatic predictors 
load(paste0(data_path, "02_month/tab_month_usa.rda"))
load(paste0(data_path, "02_month/tab_month_bra.rda"))
load(paste0(data_path, "02_month/tab_month_world.rda"))

# ----------------------------------
# ANALYSES USA 

n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

scores_usa <- list(
  # > DAILY PREDICTORS
  "data_day" = list(type_data     = "D", 
                    data          = list_data_day_usa, 
                    load_data     = F, 
                    vars_names    = vars_names, 
                    basis         = basis_day,
                    nb_comp       = 20,
                    argvals       = 1:214,
                    univExpansion = univExpansion_6_variables,
                    analysis_set  = "usa"),
  # > MONTHLY PREDICTORS 
  "data_month" = list(type_data     = "M", 
                      data          = tab_month_usa,     
                      vars_names    = vars_names, 
                      basis         = basis_month, 
                      nb_comp       = 7,
                      argvals       = 1:7,
                      univExpansion = univExpansion_6_variables,
                      analysis_set  = "usa")) %>% 
  map(., ~{
    
    # > Apply each technique on the dataset 
    # TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS")
    pca_usa <- function_pca(type_data = .x$type_data, 
                            vars_names = .x$vars_names,
                            data = .x$data, 
                            scale = T, 
                            cum_clim.var = T)
    save(pca_usa, file = paste0(data_path, "00_dim_red/pca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    fpca_usa <- function_fpca(type_data = .x$type_data, 
                              vars_names = .x$vars_names,
                              data = .x$data,
                              nb_comp = .x$nb_comp,
                              basis = .x$basis)
    save(fpca_usa, file = paste0(data_path, "00_dim_red/fpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    mfpca_usa <- function_mfpca(type_data = .x$type_data, 
                                vars_names = .x$vars_names,
                                data = .x$data,
                                nb_comp = 4,
                                argvals = .x$argvals, 
                                univExpansion = .x$univExpansion)
    save(mfpca_usa, file = paste0(data_path, "00_dim_red/mfpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION")
    plsr_usa <- function_plsr2(type_data = .x$type_data, 
                              vars_names = .x$vars_names,
                              data = .x$data,
                              nb_comp = .x$nb_comp)
    save(plsr_usa, file = paste0(data_path, "00_dim_red/plsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION")
    fplsr_usa <- function_fplsr(type_data = .x$type_data, 
                                vars_names = .x$vars_names,
                                data = .x$data,
                                nb_comp = .x$nb_comp,
                                basis = .x$basis)
    save(fplsr_usa, file = paste0(data_path, "00_dim_red/fplsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # > Extract scores derived from each dimension reduction technique
    pca_scores <- pca_usa$tab_PCA_scores %>% 
      mutate(site_year = rownames(.))
    
    fpca_scores <- fpca_usa$tab_FPCA_scores %>% 
      mutate(site_year = rownames(.))
    
    mfpca_scores <- mfpca_usa$tab_MFPCA_scores
    
    plsr_scores <- plsr_usa$tab_PLS_scores %>% 
      mutate(site_year = rownames(.))
             
    fplsr_scores <- fplsr_usa$tab_FPLS_scores %>% 
      mutate(site_year = rownames(.))
    
    # > Table with all scores  
    tab_scores <- pca_scores %>% 
      left_join(., fpca_scores,  by = "site_year") %>% 
      left_join(., mfpca_scores, by = "site_year") %>% 
      left_join(., plsr_scores,  by = "site_year") %>% 
      left_join(., fplsr_scores, by = "site_year")
    
    # > Check that all site-years have scores from each techniques
    testthat::expect_equal(nrow(tab_scores %>% drop_na()), nrow(tab_scores), info = "Some site-years are missing in the produced data")
    
    # > Merge 
    list(
      "data"  = .x$data,
      "scores"= tab_scores
      #"pca"   = pca_usa,
      #"fpca"  = fpca_usa,
      #"mfpca" = mfpca_usa,
      #"plsr"  = plsr_usa,
      #"fplsr" = fplsr_usa
    )
    
  })

# >>> stop cluster//
stopCluster(my.cluster)
beepr::beep(1)

# > Save scores
tab_scores_day_usa <- scores_usa$data_day$scores
save(tab_scores_day_usa,   file = paste0(data_path, "01_days/tab_scores_day_usa.rda"))
tab_scores_month_usa <- scores_usa$data_month$scores
save(tab_scores_month_usa, file = paste0(data_path, "02_month/tab_scores_month_usa.rda"))

# ----------------------------------
# ANALYSES BRAZIL 

scores_bra <- list(
  # > DAILY PREDICTORS
  "data_day" = list(type_data     = "D", 
                    data          = list_data_day_bra %>% 
                      # > necessary to have 212 non na values 
                      # (bissextile years have 213 observations with 29th February)
                      map(., ~{ 
                        .x %>% filter(day_of_year < 213)
                        }), 
                    load_data     = F, 
                    vars_names    = vars_names, 
                    basis         = basis_day2,
                    nb_comp       = 20,
                    argvals       = 1:212,
                    univExpansion = univExpansion_6_variables,
                    analysis_set  = "bra"),
  # > MONTHLY PREDICTORS 
  "data_month" = list(type_data     = "M", 
                      data          = tab_month_bra,     
                      vars_names    = vars_names, 
                      basis         = basis_month, 
                      nb_comp       = 7,
                      argvals       = 1:7,
                      univExpansion = univExpansion_6_variables,
                      analysis_set  = "bra")) %>% 
  map(., ~{
    
    # > Apply each technique on the dataset 
    # TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS")
    pca_bra <- function_pca(type_data = .x$type_data, 
                            vars_names = .x$vars_names,
                            data = .x$data, 
                            scale = T, 
                            cum_clim.var = T)
    save(pca_bra, file = paste0(data_path, "00_dim_red/pca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    fpca_bra <- function_fpca(type_data = .x$type_data, 
                              vars_names = .x$vars_names,
                              data = .x$data, 
                              nb_comp = 4,
                              basis = .x$basis)
    save(fpca_bra, file = paste0(data_path, "00_dim_red/fpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    mfpca_bra <- function_mfpca(type_data = .x$type_data, 
                                vars_names = .x$vars_names,
                                data = .x$data,
                                nb_comp = 4,
                                argvals = .x$argvals, 
                                univExpansion = .x$univExpansion)
    save(mfpca_bra, file = paste0(data_path, "00_dim_red/mfpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION")
    plsr_bra <- function_plsr2(type_data = .x$type_data, 
                              vars_names = .x$vars_names,
                              data = .x$data,
                              nb_comp = .x$nb_comp)
    save(plsr_bra, file = paste0(data_path, "00_dim_red/plsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION")
    fplsr_bra <- function_plsr(type_data = .x$type_data, 
                               vars_names = .x$vars_names,
                               data = .x$data,
                               nb_comp = .x$nb_comp,
                               basis = .x$basis)
    save(fplsr_bra, file = paste0(data_path, "00_dim_red/plsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # > Extract scores derived from each dimension reduction technique
    pca_scores <- pca_bra$tab_PCA_scores %>% 
      mutate(site_year = rownames(.))
    
    fpca_scores <- fpca_bra$tab_FPCA_scores %>% 
      mutate(site_year = rownames(.))
    
    plsr_scores <- plsr_bra$tab_PLS_scores %>% 
      mutate(site_year = rownames(.))
    
    fplsr_scores <- fplsr_bra$tab_FPLS_scores %>% 
      mutate(site_year = rownames(.))
    
    # > Table with all scores  
    tab_scores <- pca_scores %>% 
      left_join(., fpca_scores,  by = "site_year") %>% 
      left_join(., mfpca_scores, by = "site_year") %>% 
      left_join(., plsr_scores,  by = "site_year") %>% 
      left_join(., fplsr_scores, by = "site_year")
    
    # > Check that all site-years have scores from each techniques
    testthat::expect_equal(nrow(tab_scores %>% drop_na()), nrow(tab_scores), info = "Some site-years are missing in the produced data")
    
    # > Merge 
    list(
      "data"  = .x$data,
      "scores"= tab_scores
      #"pca"   = pca_bra,
      #"fpca"  = fpca_bra,
      #"mfpca" = mfpca_bra,
      #"plsr"  = plsr_bra,
      #"fplsr" = fplsr_bra
    )
    
    output
    
  })


# > Save
tab_scores_day_bra <- scores_bra$data_day$scores
save(tab_scores_day_bra,   file = paste0(data_path, "01_days/tab_scores_day_bra.rda"))
tab_scores_month_bra <- scores_bra$data_month$scores
save(tab_scores_month_bra, file = paste0(data_path, "02_month/tab_scores_month_bra.rda"))

# ----------------------------------
# GLOBAL ANALYSES

# As the dataset is quite big, compute first the monthly scores, then the daily scores
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# ----------
# > MONTHLY DATA 
list(
  ## > MONTHLY PREDICTORS 
  #"data_month" = list(type_data     = "M", 
  #                    data          = tab_month_world,     
  #                    vars_names    = vars_names, 
  #                    basis         = basis_month, 
  #                    nb_comp       = 7,
  #                    argvals       = 1:7,
  #                    univExpansion = univExpansion_6_variables, 
  #                    analysis_set  = "world"),
  "data_day" = list(type_data     = "D", 
                    data          = list_daily_data_world %>% 
                      # > necessary to have 212 non na values 
                      # (bisextile years have 213 observations with 29th February)
                      map(., ~{ 
                        .x %>% filter(day_of_year < 213)
                      }), 
                    load_data     = F, 
                    vars_names    = vars_names, 
                    basis         = basis_day2,
                    nb_comp       = 4,
                    argvals       = 1:212,
                    univExpansion = univExpansion_6_variables, 
                    analysis_set  = "world")) %>% 
  map(., ~{
    
    # > Apply each technique on the dataset 
    # TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS")
    pca_world <- function_pca(type_data  = .x$type_data, 
                            vars_names = .x$vars_names,
                            data  = .x$data, 
                            scale = T, 
                            cum_clim.var = T)
    save(pca_world, file = paste0(data_path, "00_dim_red/pca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    fpca_world <- function_fpca(type_data  = .x$type_data, 
                              vars_names = .x$vars_names,
                              data  = .x$data,  
                              nb_comp = 4,
                              basis = .x$basis)
    save(fpca_world, file = paste0(data_path, "00_dim_red/fpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
    message("TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
    mfpca_world <- function_mfpca(type_data = .x$type_data, 
                                vars_names = .x$vars_names,
                                data = .x$data, 
                                nb_comp = 4,
                                argvals = .x$argvals, 
                                univExpansion = .x$univExpansion)
    save(mfpca_world, file = paste0(data_path, "00_dim_red/mfpca_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 4 - PARTIAL LEAST SQUARE REGRESSION")
    plsr_world <- function_plsr(type_data = .x$type_data, 
                                vars_names = .x$vars_names,
                                data = .x$data,
                                nb_comp = .x$nb_comp)
    save(plsr_world, file = paste0(data_path, "00_dim_red/plsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
    # TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION
    message("TECHNIQUE 5 - FUNCTIONAL PARTIAL LEAST SQUARE REGRESSION")
    fplsr_world <- function_plsr(type_data = .x$type_data, 
                                 vars_names = .x$vars_names,
                                 data = .x$data,
                                 nb_comp = .x$nb_comp)
    save(fplsr_world, file = paste0(data_path, "00_dim_red/fplsr_", .x$type_data, "_", .x$analysis_set, ".rda"))
    
   ## > Extract scores derived from each dimension reduction technique
   #pca_scores <- pca_world$tab_PCA_scores %>% 
   #  mutate(site_year = rownames(.))
   #
   #fpca_scores <- fpca_world$tab_FPCA_scores %>% 
   #  mutate(site_year = rownames(.))
   #
   #mfpca_scores <- mfpca_world$tab_MFPCA_scores
   #
   #plsr_scores <- plsr_world$tab_PLS_scores %>% 
   #  mutate(site_year = rownames(.))
   #
   #fplsr_scores <- fplsr_world$tab_FPLS_scores %>% 
   #  mutate(site_year = rownames(.))
   #
   ## > Table with all scores  
   #tab_scores <- pca_scores %>% 
   #  left_join(., fpca_scores,  by = "site_year") %>% 
   #  left_join(., mfpca_scores, by = "site_year") %>% 
   #  left_join(., plsr_scores,  by = "site_year") %>% 
   #  left_join(., fplsr_scores, by = "site_year")
   #
   ## > Check that all site-years have scores from each techniques
   #testthat::expect_equal(nrow(tab_scores %>% drop_na()), nrow(tab_scores), info = "Some site-years are missing in the produced data")
   #
   ## > Merge 
   #list(
   #  "data"  = .x$data,
   #  "scores"= tab_scores
   #  #"pca"   = pca_world,
   #  #"fpca"  = fpca_world,
   #  #"mfpca" = mfpca_world,
   #  #"plsr"  = plsr_world,
   #  #"fplsr" = fplsr_world
   #)
    
    #output
    
  })

# >>> stop cluster//
stopCluster(my.cluster)
beepr::beep(1)

# > Save 
tab_scores_month_world <- scores_world_month$data_month$scores
save(tab_scores_month_world, file = paste0(data_path, "02_month/tab_scores_month_world.rda"))

# ----------
# DAILY DATA

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
  
  scores_world_day <- list(
    "data_day" = list(type_data     = "D", 
                      data          = list_daily_data_world %>% 
                        # > necessary to have 212 non na values 
                        # (bisextile years have 213 observations with 29th February)
                        map(., ~{ 
                          .x %>% filter(day_of_year < 213)
                        }), 
                      load_data     = F, 
                      vars_names    = vars_names, 
                      basis         = basis_day2,
                      nb_comp       = 4,
                      argvals       = 1:212,
                      univExpansion = univExpansion_6_variables)) %>% 
    map(., ~{
      
      # > Apply each technique on the dataset 
      # TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS
      message("TECHNIQUE 1 - PRINCIPAL COMPONENT ANALYSIS")
      pca_world <- function_pca(type_data  = .x$type_data, 
                                vars_names = .x$vars_names,
                                data  = .x$data, 
                                scale = T, 
                                cum_clim.var = T)
      
      # TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
      message("TECHNIQUE 2 - FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
      fpca_world <- function_fpca(type_data  = .x$type_data, 
                                  vars_names = .x$vars_names,
                                  data  = .x$data, 
                                  basis = .x$basis)
      
      # > Extract scores derived from each dimension reduction technique
      pca_scores <- pca_world$tab_PCA_scores %>% 
        mutate(site_year = rownames(.))
      
      fpca_scores <- fpca_world$tab_FPCA_scores %>% 
        mutate(site_year = rownames(.))
      
      # > Merge 
      output <- pca_scores %>% 
        left_join(., fpca_scores,  by = "site_year") 
      
      # > Check that all site-years have scores from each techniques
      testthat::expect_equal(nrow(output %>% drop_na()), nrow(output), info = "Some site-years are missing in the produced data")
      
      output
      
    })
  
  
  scores_world_day2 <- list(
    "data_day" = list(type_data     = "D", 
                      data          = list_daily_data_world %>% 
                        # > necessary to have 212 non na values 
                        # (bisextile years have 213 observations with 29th February)
                        map(., ~{ 
                          .x %>% filter(day_of_year < 213)
                        }), 
                      load_data     = F, 
                      vars_names    = vars_names, 
                      basis         = basis_day2,
                      nb_comp       = 4,
                      argvals       = 1:212,
                      univExpansion = univExpansion_6_variables)) %>% 
    map(., ~{
      
      # # TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS
      message("TECHNIQUE 3 - MULTIVARIATE FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS")
      mfpca_world <- function_mfpca(type_data = .x$type_data, 
                                    vars_names = .x$vars_names,
                                    data = .x$data,
                                    nb_comp = .x$nb_comp,
                                    argvals = .x$argvals, 
                                    univExpansion = .x$univExpansion)
      
      # > Extract scores derived from each dimension reduction technique
      mfpca_scores <- mfpca_world$tab_MFPCA_scores
      
      # > Check that all site-years have scores from each techniques
      testthat::expect_equal(nrow(output %>% drop_na()), nrow(output), info = "Some site-years are missing in the produced data")
      
      mfpca_scores
      
    })
  
  
})

# >>> stop cluster//
stopCluster(my.cluster)

# > Save
tab_scores_day_world <- scores_world_day$data_day %>%
  left_join(., scores_world_day2$data_day, by = "site_year")
save(tab_scores_day_world,   file = paste0(data_path, "01_days/tab_scores_day_world.rda"))


# ----------------------------------
# MERGE ALL CLIMATIC DATA 
# ----------------------------------
# Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# USA 
# > daily climatic predictors 
# scores
load(paste0(data_path, "01_days/tab_scores_day_usa.rda"))
tab_scores_day_usa_fpls <- loadRDa(paste0(data_path, "01_days/tab_scores_day_usa_fpls.rda"))
# > monthly climatic predictors 
# averages
load(paste0(data_path, "02_month/tab_month_usa.rda"))
# scores
load(paste0(data_path, "02_month/tab_scores_month_usa.rda"))
tab_scores_month_usa_fpls <- loadRDa(paste0(data_path, "02_month/tab_scores_month_usa_fpls.rda"))
# > annual predictors 
load(paste0(data_path, "03_annual/tab_year_usa.rda"))

# BRAZIL
# > daily climatic predictors 
# scores
load(paste0(data_path, "01_days/tab_scores_day_bra.rda"))
tab_scores_day_bra_fpls <- loadRDa(paste0(data_path, "01_days/tab_scores_day_bra_fpls.rda"))
# > monthly climatic predictors 
# averages
load(paste0(data_path, "02_month/tab_month_bra.rda"))
# scores
load(paste0(data_path, "02_month/tab_scores_month_bra.rda"))
tab_scores_month_bra_fpls <- loadRDa(paste0(data_path, "02_month/tab_scores_month_bra_fpls.rda"))
# > annual predictors 
load(paste0(data_path, "03_annual/tab_year_bra.rda"))

# WORLD 
# > daily climatic predictors 
load(paste0(data_path, "01_days/tab_scores_day_world.rda"))
load(paste0(data_path, "01_days/tab_scores_mfpca_day_world.rda"))
load(paste0(data_path, "01_days/tab_scores_day_world_pls.rda"))
#tab_scores_day_world_fpls <- loadRDa(paste0(data_path, "01_days/tab_scores_day_world_fpls.rda"))
# > monthly climatic predictors 
# averages
load(paste0(data_path, "02_month/tab_month_world.rda"))
# scores
load(paste0(data_path, "02_month/tab_scores_month_world.rda"))
load(paste0(data_path, "02_month/tab_scores_month_world_pls.rda"))
#tab_scores_month_world_fpls <- loadRDa(paste0(data_path, "02_month/tab_scores_month_world_fpls.rda"))
# > annual predictors 
load(paste0(data_path, "03_annual/tab_year_world.rda"))

# ----------------------------------
# Merge data 
tab_usa <- tab_month_usa %>% 
  left_join(., tab_year_usa %>% ungroup(.) %>% dplyr::select(site_year, starts_with("year_"), mean_year_z_scores), by="site_year") %>% 
  left_join(., tab_scores_day_usa, by="site_year") %>% 
  left_join(., tab_scores_day_usa_fpls, by="site_year") %>% 
  left_join(., tab_scores_month_usa, by="site_year") %>% 
  left_join(., tab_scores_month_usa_fpls, by="site_year")
dim(tab_usa) # N=29803
tab_usa %>% drop_na() %>% dim(.) # should be the same as previous, no NA allowed

tab_bra <- tab_month_bra %>% 
  left_join(., tab_year_bra %>% ungroup(.) %>% dplyr::select(site_year, starts_with("year_"), mean_year_z_scores), by="site_year") %>% 
  left_join(., tab_scores_day_bra, by="site_year") %>% 
  left_join(., tab_scores_day_bra_fpls, by="site_year") %>% 
  left_join(., tab_scores_month_bra, by="site_year") %>% 
  left_join(., tab_scores_month_bra_fpls, by="site_year")
dim(tab_bra) # N=14757
tab_bra %>% drop_na() %>% dim(.) # should be the same as previous, no NA allowed

tab_world <- tab_month_world %>% 
  left_join(., tab_year_world %>% ungroup(.) %>% dplyr::select(site_year, starts_with("year_"), mean_year_z_scores), by="site_year") %>% 
  left_join(., tab_scores_day_world, by="site_year") %>% 
  left_join(., tab_scores_month_world, by="site_year") %>% 
  left_join(., tab_scores_day_world_pls, by = "site_year") %>% 
  left_join(., tab_scores_month_world_pls, by = "site_year") %>% 
  left_join(., world_mfpca_scores, by = "site_year")
dim(tab_world) # N=122229
tab_world %>% dplyr::select(-country_code, -region) %>% drop_na() %>% dim(.) # should be the same as previous, no NA allowed

save(tab_usa, file = paste0(data_path, "tab_usa.rda"))
#save(tab_bra, file = paste0(data_path, "tab_bra.rda"))
#save(tab_world, file = paste0(data_path, "tab_world.rda"))
