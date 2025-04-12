# -------------------------------------------------------------------------
# 
#       Test FPCA Approach to detect error
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
# > functional analysis 
library(fda)
# > others
library(hydroGOF) # for NSE computation

# -------------------------------------------------------------------------
# Home-made functions
source("E:/POSTDOC INRAE/ANALYSES/00_TOOLS/00_Functions.R")
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

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

# Colors
pal <- wesanderson::wes_palette("Zissou1", 6, type="continuous")
# -------------------------------------------------------------------------
# Addins

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
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit_1"   ="Vapor pressure deficit"))

# ----------------------------------
# Data
# - climatic variables (from ERA5 database)
# - yield data (from GHDY)

load("E:/POSTDOC INRAE/ANALYSES/00_DATA/data_soybean_month.rda")
load("E:/POSTDOC INRAE/ANALYSES/00_DATA/data_soybean_daily.rda")

tab_month_test <- tab_month %>% 
  # > USA
  filter(country_name == "United States of America") %>% 
  # > years 1990-2000
  filter(year %in% 1990:2000) 

tab_month_test2 <- tab_month %>% 
  # > USA
  filter(country_name == "United States of America") %>% 
  # > years 1990-2000
  filter(year == 2001) 

dim(tab_month_test) # 9130 sites-years

load("C:/Users/benni/Documents/Post doc/ERA5_daily/era5day_test.rda")

# daily

load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/list_data_day_usa.rda")

list_data_test_day <- list_data_day_usa %>% 
  map(., ~{
    
    .x %>% 
      # > USA
      filter(country_name != "Desert") %>% 
      # > years 1990-2000
      filter(year %in% 1990:2000) 
    
    
  })

list_data_test_day2 <- list_data_day_usa %>% 
  map(., ~{
    
    .x %>% 
      # > USA
      filter(country_name != "Desert") %>% 
      # > years 2001
      filter(year == 2001) 
    
    
  })

# ----------------------------------
# Home-made functions performing the dimension reductions
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# ----------------------------------
# PCA

# test on USA - 1990:2000 data
# daily
test_pca <- function_pca(type_data = "D", vars_names = vars_names[1,],
                         data = list_data_test_day, cum_clim.var = T)

summary(test_pca$tab_PCA_scores$PC1_day_max_temp)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-37.0762 -11.3508   0.4897   0.0000  10.7260  41.7315

# monthly
test_pca_m <- function_pca(type_data = "M", vars_names = vars_names[1,],
                           data = tab_month_test)

summary(test_pca_m$tab_PCA_scores$PC1_month_max_temp)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -5.95814 -1.87488 -0.05382  0.00000  1.84553  6.10693 

# ----------------------------------
# FPCA
# test on USA - 1990:2000 data
# monthly
M_basis <- create.bspline.basis(rangeval = c(1,7), nbasis = 7, norder = 4)

test_fpca_m <- function_fpca(type_data="M", 
                             data=tab_month_test2,
                             vars_names=vars_names,
                             scale = T,
                             basis=M_basis) 

summary(test_fpca_m$tab_FPCA_scores$FPC1_month_max_temp)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -17.556848  -5.602494  -0.003194   0.000000   5.352117  18.520998 

# daily
nbasis = 150
argvals = seq(1,214,len=214)
daybasis <- create.bspline.basis(rangeval = c(1,214), 
                                 nbasis = nbasis, 
                                 norder = 4)

test_fpca <- function_fpca(type_data="D", 
                           data=list_data_test_day,
                           vars_names=vars_names[1,],
                           cum_clim.var = T,
                           scale = T,
                           basis=daybasis) 

summary(test_fpca$tab_FPCA_scores$FPC1_day_max_temp)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -13784.9  -4092.5   -180.3      0.0   4219.7  13917.2 

# ----------------------------------
# MFPCA
# test on data in USA from 1990:2000
test_mfpca_m <- function_mfpca(type_data = "M", 
                               data = tab_month_test2, 
                               vars_names = vars_names, 
                               argvals = 1:7,
                               univExpansion = list(list(type = "fda"),
                                                    list(type = "fda"),
                                                    list(type = "fda"),
                                                    list(type = "fda"),
                                                    list(type = "fda"),
                                                    list(type = "fda")))

summary(test_mfpca_m$tab_MFPCA_scores$MFPC1_month)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -24.9765  -7.7983  -0.3604   0.0000   7.8356  27.0977 

test_mfpca <- function_mfpca(type_data = "D", 
                             data = list_data_test_day, 
                             vars_names, 
                             argvals = 1:214)



# > daily climatic predictors 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"
load(paste0(data_path, "01_days/list_data_day_world.rda"))
list_daily_data <- list_daily_data_world %>% map(., ~{ .x %>% filter(year == 1991, 
                                                                     day_of_year < 213,
                                                                     country_name %in% c("Brazil", "Italy")) })

test_mfpca <- function_mfpca(type_data = "D", 
                             data = list_daily_data, 
                             vars_names, 
                             argvals = 1:212)

# ----------------------------------

# PLS
# monthly
test_pls_m <- function_plsr(type_data = "M", 
                            data = tab_month_test, 
                            vars_names = vars_names, 
                            nb_comp = 7, scale = T)


test_pls_m$tab_PLS_scores$site_year <- rownames(test_pls_m$tab_PLS_scores)


pls.m.test <- test_pls_m$list_pls_per_variable$max_temp
plot(RMSEP(pls.m.test), legendpos = "topright", type="p")

par(mfrow = c(2,4))
for(i in 1:7){ 
  plot(x = predict(pls.m.test, comp = 1:i), 
       y = pls.m.test$model$Ya, xlab = "Ya reconstructed", ylab="Ya observed", 
       main = paste0(i, " PLS scores"))
  abline(a = 0, b=1, col="red")
  }

plot(test_pls_m$list_pls_per_variable$max_temp, ncomp = ncol(test_pls_m$tab_PLS_scores), line = TRUE)




# daily 
test_pls_d <- function_plsr(type_data = "D",
                            load_data = F, cum_clim.var = T,
                            data = list_data_test_day, 
                            vars_names = vars_names[1:3,], 
                            nb_comp = 20, scale = T)



head(test_pls_d$tab_PLS_scores)


par(mfrow=c(1,1))
ncomp.onesigma <- selectNcomp(test_pls_d$list_pls_per_variable$et0[[1]], 
                              method = "onesigma", plot = TRUE, ylim=c(0,1))

for(var in 1:3)
{
  
  pls.d.test <- test_pls_d$list_pls_per_variable[[var]][[1]]
  
  par(mfrow = c(4,5))
  
  for(i in 1:ncol(pls.d.test$scores)){ 
    plot(x = predict(pls.d.test, comp = 1:i), 
         y = pls.d.test$model$Ya, xlab = paste0("Ya reconstructed\nusing", i, " PLS scores"), ylab="Ya observed", 
         main = paste0(unique(vars_names$clim.var)[var]), cex = 0.2)
    abline(a = 0, b=1, col="red")
  }
  
}


testscores <- newscores_plsr2(type_data="D", 
                              vars_names=vars_names[c(1,5),], 
                              init_data=tab_month_test,
                              new_data=tab_month_test2, 
                              init_data_day = list_data_test_day, 
                              new_data_day = list_data_test_day2, 
                              outcome = "Ya")

init_data=tab_month_test
new_data=tab_month_test2 
init_data_day = list_data_test_day
new_data_day = list_data_test_day2


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
                            vars_names = vars_names[c(1,5),],
                            data       = data_day_init,
                            nb_comp    = 20, 
                            outcome    = "Ya")

# > List to store the future scores 
list_scores_pls <- list()

# > Apply PLSR loads on new data 
var_j <- unique(vars_names[c(1,5),]$clim.var)[1] 

  
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
  tab_scores_plsr_j <- predict(plsr_j, type="scores", newdata=as.data.frame(z_tab_new_data_j))
  colnames(tab_scores_plsr_j) <- paste0("PLS", 1:ncol(tab_scores_plsr_j),"_day_", clim.var_abb_j)
  

# ----------------------------------

# FPLS
# monthly
# > basis
monthbasis7 <- create.bspline.basis(rangeval = c(1,7), nbasis = 7, norder = 4)

test_fpls_m <- function_fplsr(type_data = "M", 
                            data = tab_month_test, 
                            vars_names = vars_names[1,], 
                            nb_comp = 7, 
                            scale = T, 
                            basis = monthbasis7)

head(test_fpls_m$tab_FPLS_scores)


# daily
basis_day  <- create.bspline.basis(rangeval = c(1,214), nbasis = 150, norder = 4)

test_fpls_d <- function_fplsr(type_data = "D", 
                              data = list_data_test_day, 
                              vars_names = vars_names[1,], 
                              nb_comp = 20, 
                              scale = T, 
                              basis = basis_day)



# test cv 
# model formula
# Models list
list_models <- list(avg.zscore.a = list(name = "avg.zscore.a", formula = "mean_year_z_scores"))

model_formula <- paste0("Ya ~ irrigated_portion + ", list_models$avg.zscore.a$formula)
model_name    <- list_models$avg.zscore.a$name

# > fit
mod  <- lm(as.formula(model_formula),
           data=tab_bra) 

# > cross-validation & save prediction for each model and site-year
mod_cv <- function_cv_year(model      = model_formula, 
                           model_gpe  = "lm",
                           data       = tab_bra, 
                           model_name = paste0("lm_", model_name), 
                           res        = "perf")



