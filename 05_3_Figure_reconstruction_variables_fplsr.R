# -------------------------------------------------------------------------
# 
#           DIMENSION REDUCTION IN CLIMATIC PREDICTORS 
# 
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
library(fda.usc)

# ----------------------------------
# > Path to files
path_to_dim_red <- "/work_home/machen/02_data/00_dim_red/"

# ---------------------------------------
# EXTRACT SCORES FOR FPLSR for global dataset
# > on monthly data 
# > Path to files
message("fplsr_M_world")
path_to_fplsr_M <- "fplsr_M_world/"
files_M <- list.files(paste0(path_to_dim_red, path_to_fplsr_M), pattern = ".rda$")

for(i in 1:length(files_M))
{
  
  fplsr_i <- loadRDa(file = paste0(path_to_dim_red, path_to_fplsr_M, files_M[i]))
  cat(paste0(names(fplsr_i$list_fpls_per_variable), "\n"))
  tryCatch(print(summary(fplsr_i$list_fpls_per_variable[[1]], draw=F)), error = function(e) {})
  cat("\n")
}

# > on daily data 
message("fplsr_D_world")
path_to_fplsr_D <- "fplsr_D_world/"
files_D <- list.files(paste0(path_to_dim_red, path_to_fplsr_D), pattern = ".rda$")

for(i in 1:length(files_D))
{
  
  fplsr_i <- loadRDa(file = paste0(path_to_dim_red, path_to_fplsr_D, files_D[i]))
  cat(paste0(names(fplsr_i$list_fpls_per_variable), "\n"))
  tryCatch(print(summary(fplsr_i$list_fpls_per_variable[[1]], draw=F)), error = function(e) {})
  cat("\n")
  
}

# ---------------------------------------
# EXTRACT SCORES FOR FPLSR for USA
message("fplsr_M_usa")
fplsr_M_usa <- loadRDa(file = paste0(path_to_dim_red, "fplsr_M_usa.rda"))
message("et0")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$et0)), error = function(e) {})
message("max_temp")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$max_temp)), error = function(e) {})
message("min_temp")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$min_temp)), error = function(e) {})
message("prec")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$prec)), error = function(e) {})
message("rad")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$rad)), error = function(e) {})
message("vapor_pressure_deficit")
tryCatch(print(summary(fplsr_M_usa$list_fpls_per_variable$vapor_pressure_deficit)), error = function(e) {})

message("fplsr_D_usa")
fplsr_D_usa <- loadRDa(file = paste0(path_to_dim_red, "fplsr_D_usa.rda"))
message("et0")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$et0)), error = function(e) {})
message("max_temp")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$max_temp)), error = function(e) {})
message("min_temp")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$min_temp)), error = function(e) {})
message("prec")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$prec)), error = function(e) {})
message("rad")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$rad)), error = function(e) {})
message("vapor_pressure_deficit")
tryCatch(print(summary(fplsr_D_usa$list_fpls_per_variable$vapor_pressure_deficit)), error = function(e) {})

# ---------------------------------------
# EXTRACT SCORES FOR FPLSR for BRA
message("fplsr_M_bra")
fplsr_M_bra <- loadRDa(file = paste0(path_to_dim_red, "fplsr_M_bra.rda"))
message("et0")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$et0)), error = function(e) {})
message("max_temp")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$max_temp)), error = function(e) {})
message("min_temp")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$min_temp)), error = function(e) {})
message("prec")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$prec)), error = function(e) {})
message("rad")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$rad)), error = function(e) {})
message("vapor_pressure_deficit")
tryCatch(print(summary(fplsr_M_bra$list_fpls_per_variable$vapor_pressure_deficit)), error = function(e) {})

message("fplsr_D_bra")
fplsr_D_bra <- loadRDa(file = paste0(path_to_dim_red, "fplsr_D_bra.rda"))
message("et0")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$et0)), error = function(e) {})
message("max_temp")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$max_temp)), error = function(e) {})
message("min_temp")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$min_temp)), error = function(e) {})
message("prec")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$prec)), error = function(e) {})
message("rad")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$rad)), error = function(e) {})
message("vapor_pressure_deficit")
tryCatch(print(summary(fplsr_D_bra$list_fpls_per_variable$vapor_pressure_deficit)), error = function(e) {})
