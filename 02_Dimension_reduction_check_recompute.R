# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)

# Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

list_scores <- list()
# Load 
for(country in c("world", "usa", "bra"))
{
  
  for(data_type in c("M", "D"))
  {
    
    for(model in c("pca", "fpca", "mfpca", "plsr"))
    {
      
      mod <- loadRDa(paste0(data_path, "00_dim_red/", model, "_", data_type, "_", country, ".rda"))
      
      list_scores[[paste0(country)]][[paste0(data_type)]][[paste0(model)]] <- mod[1]
      
    }
    
    
  }
  
}

# > Compare with initial data 
load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))
load(paste0(data_path, "tab_world.rda"))

list_scores_init <- list()
# Load 
for(country in c("world", "usa", "bra"))
{
  
  if(country == "world"){
    dat <- tab_world
  }
  if(country == "usa"){
    dat <- tab_usa
  }
  if(country == "bra"){
    dat <- tab_bra
  }
  
  for(data_type in c("month", "day"))
  {
    
    for(model in c("PC", "FPC", "MFPC", "PLS"))
    {
      
      list_scores_init[[paste0(country)]][[paste0(data_type)]][[paste0(model)]] <- dat %>% 
        dplyr::select(contains(data_type)) %>% 
        dplyr::select(starts_with(model))
      
    }
    
    
  }
  
}







