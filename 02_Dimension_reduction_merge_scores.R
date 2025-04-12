# -------------------------------------------------------------------------
# 
#           DIMENSION REDUCTION IN CLIMATIC PREDICTORS 
#             MERGE COMPUTED SCORES FOR EACH ANALYSES
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

# ----------------------------------
# > Path to files
path_to_dim_red <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/"
patterns <- c("M_usa", "D_usa", "M_bra", "D_bra", "M_world", "D_world")

# > Load and read
for(p in patterns)
{
  
  # > Files to read 
  files <- list.files(path_to_dim_red, pattern = paste0(p, ".rda$"))
  
  # > Object to store the scores
  tab_scores_i <- list()
  
  # > Load and retrieve the scores from FPLS
  for(i in 1:length(unique(files)))
  {
    
    # Load 
    dim_red_i <- loadRDa(file = paste0(path_to_dim_red, files[i]))
    
    # Scores
    scores_i  <- dim_red_i[[1]] %>% 
      dplyr::select(-starts_with("site_year"))
    
    # Remove all before and up to ":":
    colnames(scores_i) <- gsub("^.*\\.","", names(scores_i))
    
    # Store
    tab_scores_i[[paste0(files[i])]] <- scores_i
    rm(scores_i, dim_red_i)
    gc()
    
  }
  
  # > Merge all scores and save
  tab_scores_p <- map_dfc(tab_scores_i, ~{ .x })
  save(tab_scores_p, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_scores/tab_scores_", p, ".rda"))
  rm(tab_scores_p)
  gc()

}

# ---------------------------------------
# EXTRACT SCORES FOR FPLSR for global dataset
# > on monthly data 
# > Path to files
path_to_fplsr_M <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/fplsr_M_world/"
files_M <- list.files(path_to_fplsr_M, pattern = ".rda$")

# > Object to store the scores
fplsr_M_scores_l <- list()

# > Load and retrieve the scores from FPLS
for(i in 1:length(unique(files_M)))
{
  
  fplsr_i <- loadRDa(file = paste0(path_to_fplsr_M, files_M[i]))
  scores_i <- fplsr_i$tab_FPLS_scores
  
  # Remove all before and up to ":":
  colnames(scores_i) <- gsub("^.*\\.","", names(scores_i))
  
  fplsr_M_scores_l[[paste0(files_M[i])]] <- scores_i
  rm(scores_i, fplsr_i)
  
}

# > Merge all scores
tab_scores_fplsr_M_world <- map_dfc(fplsr_M_scores_l, ~{ .x })
# > Save
save(tab_scores_fplsr_M_world, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_scores/tab_scores_fplsr_M_world.rda"))

# > daily data 
# > Path to files
path_to_fplsr_D <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/fplsr_D_world/"
files_D <- list.files(path_to_fplsr_D, pattern = ".rda$")

# > Object to store the scores
fplsr_D_scores_l <- list()

# > Load and retrieve the scores from FPLS
for(i in 1:length(unique(files_D)))
{
  
  fplsr_i <- loadRDa(file = paste0(path_to_fplsr_D, files_D[i]))
  scores_i <- fplsr_i$tab_FPLS_scores
  
  # Remove all before and up to ":":
  colnames(scores_i) <- gsub("^.*\\.","", names(scores_i))
  
  fplsr_D_scores_l[[paste0(files_D[i])]] <- scores_i
  rm(scores_i, fplsr_i)
  
}

# > Merge all scores
tab_scores_fplsr_D_world <- map_dfc(fplsr_D_scores_l, ~{ .x })
# > Save
save(tab_scores_fplsr_D_world, file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_scores/tab_scores_fplsr_D_world.rda"))


