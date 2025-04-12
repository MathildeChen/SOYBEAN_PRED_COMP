# -------------------------------------------------------------------------
#
#       ERA5-Land daily daily climatic data 2000-2017 for EU
#       Resolution: 0.5*0.5°
#       Merge data from individual .nc to global dataset 
#         
# -------------------------------------------------------------------------

# ----------------------------------------
# Packages & tools
library(tidyverse)
library(stringr)
library(lubridate)
library(terra) ; library(rnaturalearth)
library(parallel) ; library(doParallel); library(foreach)
library(CCMHr)

# Homemade function
source("E:/POSTDOC INRAE/DATA/01_CLIMATE/ERA5/functions_to_read_era5.R")

# Home-made functions performing the dimension reductions
source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# ----------------------------------------
# Data 

# > path to climatic data
path <- "C:/Users/benni/Documents/Post doc/Test"

# > load 1 initial yield file to resample era5 data 
yield_ref <- rast("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/maize/yield_1981.nc4")

# > path to yield data (and coordinates of pixels for which we want to load climate data)
dat_coords_USA <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_1_yield_usa.rda") 
dat_coords_BRA <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_2_yield_bra.rda") 
#dat_coords_WORLD <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_3_yield_world.rda") 

# > count nb of cells per crop (including and excluding desert)
dim(dat_coords_USA) ; length(unique(dat_coords_USA$gridcode))
# 37147 sites-years, 1034 sites
dim(dat_coords_USA[which(dat_coords_USA$country_name != "Desert"),]) ; length(unique(dat_coords_USA[which(dat_coords_USA$country_name != "Desert"),]$gridcode))
# 29803 sites-years, 830 sites

dim(dat_coords_BRA) ; length(unique(dat_coords_BRA$gridcode))
# 18429 sites-years, 526 sites
dim(dat_coords_BRA[which(dat_coords_BRA$country_name != "Desert"),]) ; length(unique(dat_coords_USA[which(dat_coords_USA$country_name != "Desert"),]$gridcode))
# 14575 sites-years, 424 sites

# ----------------------------------------
# Individual .nc files with climatic ERA5 data 

# > extract all the names of the files
filenames <- list.files(path, pattern="*.nc", full.names = TRUE)

# > split the files among the different variables 
filetable <- data.frame(filename = filenames) %>% 
  # > add variable 
  mutate(var = case_when(
    str_detect(filename, "10m_u_component_of_wind")         == T ~ "10m_u_component_of_wind",
    str_detect(filename, "10m_v_component_of_wind")         == T ~ "10m_v_component_of_wind",
    str_detect(filename, "total_precipitation")             == T ~ "total_precipitation",
    str_detect(filename, "mean_2m_temperature")             == T ~ "2m_temperature",
    str_detect(filename, "mean_2m_dewpoint_temperature")    == T ~ "2m_dewpoint_temperature",
    str_detect(filename, "minimum_2m_temperature")          == T ~ "min_2m_temperature",
    str_detect(filename, "maximum_2m_temperature")          == T ~ "max_2m_temperature",
    str_detect(filename, "minimum_2m_dewpoint_temperature") == T ~ "min_2m_dewpoint_temperature",
    str_detect(filename, "maximum_2m_dewpoint_temperature") == T ~ "max_2m_dewpoint_temperature",
    str_detect(filename, "surface_pressure")                == T ~ "surface_pressure",
    str_detect(filename, "surface_net_solar_radiation")     == T ~ "surface_net_solar_radiation"
  )) %>% 
  # > add month and year 
  mutate(year  = substr(substr(filename, nchar(filename)-10, nchar(filename)), 2, 5),
         month = substr(substr(filename, nchar(filename)-10, nchar(filename)), 7, 8))

# > examine data 
filetable %>% 
  group_by(var, year) %>%
  summarise(n_files=n()) %>% 
  ggplot(., aes(x=as.numeric(as.character(year)), y=var, fill=as.factor(n_files))) +
  geom_tile(colour="white") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_viridis_d(direction = -1, name="Number of months available") + 
  labs(x="Years")

# For the moment, daily data are split between months
# 1 file = daily data for each pixel / month / year
# Temporal range: 1980-2017
# Spatial coverage: global, 0.5° resolution

# ----------------------------------------
# PREPARATION FOR DATA LOADING
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
                               "vapor_pressure_deficit"     ="vpd_1"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature",
                               "max_2m_temperature"         ="Maximum temperature",
                               "et0"                        ="Evapotranspiration ref",
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit"   ="Vapor pressure deficit"))

# > ERA5 variables to compute VPD, ET0, and wind speed
var_vpd_1 <- c("min_2m_temperature", "max_2m_temperature", "min_2m_dewpoint_temperature", "max_2m_dewpoint_temperature")

var_et0 <- c("10m_u_component_of_wind", "10m_v_component_of_wind", "min_2m_temperature", "max_2m_temperature", 
             "2m_dewpoint_temperature", "surface_net_solar_radiation", "surface_pressure")

# > select files to merge for each variable 
files_to_merge <- list()

# test set 
#files_to_merge[[paste0("max_temp")]] <- filetable %>% mutate(to_keep = case_when(var == "max_2m_temperature"  & year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 
#files_to_merge[[paste0("min_temp")]] <- filetable %>% mutate(to_keep = case_when(var == "min_2m_temperature"  & year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 
#files_to_merge[[paste0("prec")]]     <- filetable %>% mutate(to_keep = case_when(var == "total_precipitation"& year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 
#files_to_merge[[paste0("rad")]]      <- filetable %>% mutate(to_keep = case_when(var == "surface_net_solar_radiation" & year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 
#files_to_merge[[paste0("et0")]]      <- filetable %>% mutate(to_keep = case_when(var %in% var_et0  & year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 
#files_to_merge[[paste0("vpd_1")]]    <- filetable %>% mutate(to_keep = case_when(var %in% var_vpd_1& year == 2000 ~ 1, TRUE ~ 0)) %>% filter(to_keep == 1) 

# full set 
files_to_merge[[paste0("max_temp")]] <- filetable %>% mutate(to_keep = case_when(var == "max_2m_temperature" & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var == "max_2m_temperature" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("min_temp")]] <- filetable %>% mutate(to_keep = case_when(var == "min_2m_temperature" & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var == "min_2m_temperature" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("prec")]]     <- filetable %>% mutate(to_keep = case_when(var == "total_precipitation" & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var == "total_precipitation" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("rad")]]      <- filetable %>% mutate(to_keep = case_when(var == "surface_net_solar_radiation" & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var == "surface_net_solar_radiation" & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("et0")]]      <- filetable %>% mutate(to_keep = case_when(var %in% var_et0  & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var %in% var_et0  & year %in% 1981:2016 ~ 1, 
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 
files_to_merge[[paste0("vpd_1")]]    <- filetable %>% mutate(to_keep = case_when(var %in% var_vpd_1 & year == 1980 & month %in% c("11", "12") ~ 1,
                                                                                 var %in% var_vpd_1 & year %in% 1981:2016 ~ 1,
                                                                                 TRUE ~ 0)) %>% filter(to_keep == 1) 

# > Check selected files 
files_to_merge %>% 
  map(., ~{
    .x %>% group_by(var, year) %>% count() %>% 
      ggplot(., aes(x=as.numeric(as.character(year)), y=n, color=var)) + geom_line() + geom_point() + theme_bw() + theme(legend.position = "none") + facet_wrap(.~var)
  })

# ----------------------------------------
# USA 

# Coordinates
dat_coord_i <- dat_coords_USA %>% 
  ungroup() %>%
  distinct(x, y, gridcode, country_name) %>%
  filter(country_name == "Desert")

dim(dat_coord_i) # N=204

# Folder to save
save_i <- "02_1_clim_days_usa"

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop
# >>> create the cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
# >>> register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Load, merge and recompute variable for the set of gridcells used for test
for(v in c("min_temp", "prec", "rad"))#, "max_temp", "et0", "vpd_1"))
{
  
  if(v %in% c("max_temp", "min_temp", "prec", "rad"))
  {
    # > Split files by decade
    # 1981-1990
    files_to_merge_1981_1990 <- files_to_merge[[paste0(v)]] %>% 
      mutate(year = as.numeric(as.character(year))) %>% 
      filter(year <1991)
    # 1991-2000
    files_to_merge_1991_2000 <- files_to_merge[[paste0(v)]] %>% 
      mutate(year = as.numeric(as.character(year))) %>% 
      filter(year > 1990) %>% 
      filter(year < 2001)
    # 2001-2010
    files_to_merge_2001_2010 <- files_to_merge[[paste0(v)]] %>% 
      mutate(year = as.numeric(as.character(year))) %>% 
      filter(year >=2001) %>% 
      filter(year <2011)
    # 2011-2016
    files_to_merge_2011_2016 <- files_to_merge[[paste0(v)]] %>% 
      mutate(year = as.numeric(as.character(year))) %>% 
      filter(year >=2011)
    
    files_to_merge_v <- list("1981_1990"=files_to_merge_1981_1990,
                             "1991_2000"=files_to_merge_1991_2000,
                             "2001_2010"=files_to_merge_2001_2010,
                             "2011_2016"=files_to_merge_2011_2016)
    
    # > Load data 
    for(d in c("1981_1990", "1991_2000", "2001_2010", "2011_2016"))
    {
      
      # > Load from files 
      era5daily_init <- merge_era5_data(var = v,
                                        crop = "Soybean",
                                        files_to_merge = files_to_merge_v[[paste0(d)]]$filename, 
                                        dat.coords = dat_coord_i,
                                        yield_ref = yield_ref, 
                                        save_output = F)
      
      # > (Re)Compute the variables
      era5daily_correct<- correct_era5_data(clim.var = v, 
                                            data.clim.var = era5daily_init, 
                                            cum.value = T)
      
      # > save
      save(era5daily_correct, 
           file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_usa.rda"))
      
      # > remove
      rm(era5daily_init, era5daily_correct)
      gc()
      
      
    }
  }
  if(v %in% c("et0", "vpd_1"))
  {
    # > Split files by year
    files_to_merge_v <- files_to_merge[[paste0(v)]] %>% split(.$year)
    
    # > Load data 
    for(d in unique(files_to_merge[[paste0(v)]]$year))
    {
      
      # > Load from files 
      era5daily_init <- merge_era5_data(var = v,
                                        crop = "Soybean",
                                        files_to_merge = files_to_merge_v[[paste0(d)]]$filename, 
                                        dat.coords = dat_coord_i,
                                        yield_ref = yield_ref, 
                                        save_output = F)
      
      # > (Re)Compute the variables
      era5daily_correct<- correct_era5_data(clim.var = v, 
                                            data.clim.var = era5daily_init, 
                                            cum.value = T)
      
      # > save
      save(era5daily_correct, 
           file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_usa.rda"))
      
      # > remove
      rm(era5daily_init, era5daily_correct)
      gc()
      
      
    }
  }
}  

# >>> stop cluster//
stopCluster(my.cluster)

# > Merge data 
# Load, merge and recompute variable for the set of gridcells used for test

# > Object to store desert data 
list_data_day_usa_desert <- list()

for(v in c("max_temp", "min_temp", "prec", "rad", "et0", "vpd_1"))
{
  # > List variables
  list_v <- list()
  
  if(v %in% c("max_temp", "min_temp", "prec", "rad"))
  {
    # > Load data 
    for(d in c("1981_1990", "1991_2000", "2001_2010", "2011_2016"))
    {
      
      # > Load data
      era5daily_correct_v_d <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_usa.rda"))
      
      # > Store in a list
      list_v[[paste0(d)]] <- era5daily_correct_v_d
      
      # > remove unused files
      rm(era5daily_correct_v_d)
      
    }
  }
  if(v %in% c("et0", "vpd_1"))
  {
    # > Load data 
    for(d in unique(files_to_merge[[paste0(v)]]$year))
    {
      
      # > Load data
      era5daily_correct_v_d <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_usa.rda"))
      
      # > Store in a list
      list_v[[paste0(d)]] <- era5daily_correct_v_d
      
      # > remove unused files
      rm(era5daily_correct_v_d)
      
    }
  }
  
  # > Merge data for the variable 
  era5daily_correct_v <- map_dfr(list_v, data.frame)
  list_data_day_usa_desert[[paste0(v)]] <- era5daily_correct_v

}  

# > Check
list_data_day_usa_desert %>% map(., ~ {  
  dim(.x)
  range(.x$date)
  })

# > Merge with initial usa data 
# Load initial data 
list_data_day_usa_init <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/02_1_clim_days_usa/list_data_day_usa.rda")

# New data 
list_data_day_usa <- list()

for(i in 1:6)
{
  
  # > USA data 
  list_data_day_usa_v <- list_data_day_usa_init[[paste0(vars_names[i,"clim.var"])]] %>% 
    # Add the Ya_ano data 
    left_join(., dat_coords_USA %>% 
                mutate(year = as.numeric(as.character(year))) %>% 
                dplyr::select(x, y, year, Ya_ano), by = c("x", "y", "year")) %>% 
    dplyr::select(site_year, x, y, year, gridcode, country_name, date, day_of_year, clim.var, clim.value, cum_clim.value, Ya, Ya_ano)
  
  # > Desert data
  list_data_day_usa_desert_v <- list_data_day_usa_desert[[paste0(vars_names[i,"clim.var_abb"])]] %>% 
    mutate(year = year(ymd(date))) %>% 
    # Add the Ya_ano data 
    left_join(., dat_coords_USA %>%
                mutate(year = as.numeric(as.character(year))) %>% 
                dplyr::select(x, y, year, Ya, Ya_ano), by = c("x", "y", "year")) %>% 
    dplyr::select(site_year, x, y, year, gridcode, country_name, date, day_of_year, clim.var, clim.value, cum_clim.value, Ya, Ya_ano)
  
  # > Check
  testthat::expect_equal(unique(list_data_day_usa_v$clim.var), unique(list_data_day_usa_desert_v$clim.var))
  
  # > Add both datasets
  list_data_day_usa[[paste0(vars_names[i,"clim.var"])]] <- rbind(list_data_day_usa_v, list_data_day_usa_desert_v)
  
}
 
# Check
list_data_day_usa %>% 
  map(., ~{
    
    unique(.x %>% group_by(site_year) %>% count() %>% pull(n))
    
  })

# Save new data (USA + Desert)
save(list_data_day_usa, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/list_data_day_usa.rda")

# > Compute monthly data 
list_data_month_usa <- list_data_day_usa %>% 
  map(., ~{
    
    monthly_average(var_i = unique(.x$clim.var),
                    load = F, 
                    data_var_i = .x)
    
    
  })

# Table with monthly data
tab_month_usa <- plyr::ldply(list_data_month_usa, data.frame) %>% 
  # > rename mean columns
  rename("monthly"     = "mean_clim.value",
         "monthly_cum" = "mean_cum_clim.value") %>% 
  # > identify month id
  group_by(site_year, x, y, gridcode, country_name, clim.var) %>%
  arrange(year, month) %>%
  mutate(month_nb = row_number()) %>% 
  # > remove useless column
  ungroup() %>% 
  dplyr::select(-.id, -month) %>% 
  # > pivot 
  pivot_wider(names_from = c("clim.var", "month_nb"), values_from = c("monthly", "monthly_cum")) 

head(tab_month_usa)
dim(tab_month_usa) # N=37147 
dim(dat_coords_USA) # should have the same number of lines

# Save
save(tab_month_usa, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/02_month/tab_month_usa.rda")

# Table with annual averages
tab_year_usa <- tab_month_usa %>% 
  group_by(site_year, x, y, gridcode, country_name, year, Ya, Ya_ano) %>% 
  summarise(year_max_2m_temperature         =mean(c(monthly_max_2m_temperature_1, monthly_max_2m_temperature_2, monthly_max_2m_temperature_3, monthly_max_2m_temperature_4, monthly_max_2m_temperature_5, monthly_max_2m_temperature_6, monthly_max_2m_temperature_7)), 
            year_min_2m_temperature         =mean(c(monthly_min_2m_temperature_1, monthly_min_2m_temperature_2, monthly_min_2m_temperature_3, monthly_min_2m_temperature_4, monthly_min_2m_temperature_5, monthly_min_2m_temperature_6, monthly_min_2m_temperature_7)),
            year_et0                        =mean(c(monthly_et0_1, monthly_et0_2, monthly_et0_3, monthly_et0_4, monthly_et0_5, monthly_et0_6, monthly_et0_7)),
            year_surface_net_solar_radiation=mean(c(monthly_surface_net_solar_radiation_1, monthly_surface_net_solar_radiation_2, monthly_surface_net_solar_radiation_3, monthly_surface_net_solar_radiation_4, monthly_surface_net_solar_radiation_5, monthly_surface_net_solar_radiation_6, monthly_surface_net_solar_radiation_7)),
            year_total_precipitation        =mean(c(monthly_total_precipitation_1, monthly_total_precipitation_2, monthly_total_precipitation_3, monthly_total_precipitation_4, monthly_total_precipitation_5, monthly_total_precipitation_6, monthly_total_precipitation_7)),
            year_vapor_pressure_deficit     =mean(c(monthly_vapor_pressure_deficit_1, monthly_vapor_pressure_deficit_2, monthly_vapor_pressure_deficit_3, monthly_vapor_pressure_deficit_4, monthly_vapor_pressure_deficit_5, monthly_vapor_pressure_deficit_6, monthly_vapor_pressure_deficit_7))) 

dim(tab_year_usa) # N=37147

# Table with zscore of annual averages
tab_year_zscore_usa <- tab_year_usa %>% 
  ungroup(.) %>% 
  dplyr::select(site_year, starts_with("year_")) %>% 
  mutate_if(is.numeric, scale) %>% 
  pivot_longer(cols = starts_with("year_"), names_to = "clim.var", values_to = "clim.z_scores") %>% 
  group_by(site_year) %>%
  summarise(mean_year_z_scores = mean(clim.z_scores))

dim(tab_year_zscore_usa) # N=37147

# > Merge
tab_year_usa <- left_join(tab_year_usa, tab_year_zscore_usa, by = "site_year")

# > Save 
save(tab_year_usa, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/03_annual/tab_year_usa.rda")

# -----------------------------
# Brazil 
dat_coord_i <- dat_coords_BRA %>% 
  ungroup() %>% 
  distinct(x, y, gridcode, country_name) %>%
  filter(country_name == "Desert")

dim(dat_coord_i) # N=102

save_i <- "02_2_clim_days_bra"

# Load, merge and recompute variable for the set of gridcells used for test
for(v in c("max_temp", "min_temp", "prec", "rad", "et0", "vpd_1"))
{
  
  # > Split files by decade
  # 1981-1990
  files_to_merge_1981_1990 <- files_to_merge[[paste0(v)]] %>% 
    filter(year <1991)
  # 1991-2000
  files_to_merge_1991_2000 <- files_to_merge[[paste0(v)]] %>% 
    filter(year >=1991) %>% 
    filter(year <2001)
  # 2001-2010
  files_to_merge_2001_2010 <- files_to_merge[[paste0(v)]] %>% 
    filter(year >=2001) %>% 
    filter(year <2011)
  # 2011-2016
  files_to_merge_2011_2016 <- files_to_merge[[paste0(v)]] %>% 
    filter(year >=2011)
  
  files_to_merge_v <- list("1981_1990"=files_to_merge_1981_1990,
                           "1991_2000"=files_to_merge_1991_2000,
                           "2001_2010"=files_to_merge_2001_2010,
                           "2011_2016"=files_to_merge_2011_2016)
  
  # > Load data 
  for(d in c("1981_1990", "1991_2000", "2001_2010", "2011_2016"))
  {
    
    # > Load from files 
    era5daily_init <- merge_era5_data(var = v,
                                      crop = "Soybean",
                                      files_to_merge = files_to_merge_v[[paste0(d)]]$filename, 
                                      dat.coords = dat_coord_i,
                                      yield_ref = yield_ref, 
                                      save_output = F)
    
    # > (Re)Compute the variables
    era5daily_correct<- correct_era5_data(clim.var = v, 
                                          data.clim.var = era5daily_init, 
                                          cum.value = T)
    
    # > save
    save(era5daily_correct, 
         file = paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_bra.rda"))
    
    # > remove
    rm(era5daily_init, era5daily_correct)
    gc()
    
  }
  
}  

# > Merge data 
# Load, merge and recompute variable for the set of gridcells used for test

# > Object to store desert data 
list_data_day_bra_desert <- list()

for(v in c("max_temp", "min_temp", "prec", "rad", "et0", "vpd_1"))
{
  # > List variables
  list_v <- list()
  
  if(v %in% c("max_temp", "min_temp", "prec", "rad", "et0", "vpd_1"))
  {
    # > Load data 
    for(d in c("1981_1990", "1991_2000", "2001_2010", "2011_2016"))
    {
      
     # > Load data
      era5daily_correct_v_d <- loadRDa(paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/", save_i, "/era5daily_", v, "_", d ,"_desert_bra.rda"))
      
      # > Store in a list
      list_v[[paste0(d)]] <- era5daily_correct_v_d
      
      # > remove unused files
      rm(era5daily_correct_v_d)
      
    }
  }
  
  # > Merge data for the variable 
  era5daily_correct_v <- map_dfr(list_v, data.frame)
  list_data_day_bra_desert[[paste0(v)]] <- era5daily_correct_v
  rm(era5daily_correct_v)
}  

# > Check
list_data_day_bra_desert %>% map(., ~ {  
  dim(.x)
  range(.x$date)
})

# > Merge with initial usa data 
# Load initial data 
list_data_day_bra_init <- loadRDa("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/02_2_clim_days_bra/list_data_day_bra.rda")

# New data 
list_data_day_bra <- list()

for(i in 1:6)
{
  
  # > BRA data 
  list_data_day_bra_v <- list_data_day_bra_init[[paste0(vars_names[i,"clim.var"])]] %>% 
    # Add the Ya_ano data 
    left_join(., dat_coords_BRA %>% 
                mutate(year = as.numeric(as.character(year))) %>% 
                dplyr::select(x, y, year, Ya_ano), by = c("x", "y", "year")) %>% 
    dplyr::select(site_year, x, y, year, gridcode, country_name, date, day_of_year, clim.var, clim.value, cum_clim.value, Ya, Ya_ano)
  
  # > Desert data
  list_data_day_bra_desert_v <- list_data_day_bra_desert[[paste0(vars_names[i,"clim.var_abb"])]] %>% 
    mutate(year = year(ymd(date))) %>% 
    # Add the Ya_ano data 
    left_join(., dat_coords_BRA %>%
                mutate(year = as.numeric(as.character(year))) %>% 
                dplyr::select(x, y, year, Ya, Ya_ano), by = c("x", "y", "year")) %>% 
    dplyr::select(site_year, x, y, year, gridcode, country_name, date, day_of_year, clim.var, clim.value, cum_clim.value, Ya, Ya_ano)
  
  # > Check
  testthat::expect_equal(unique(list_data_day_bra_v$clim.var), unique(list_data_day_bra_desert_v$clim.var))
  
  # > Add both datasets
  list_data_day_bra[[paste0(vars_names[i,"clim.var"])]] <- rbind(list_data_day_bra_v, list_data_day_bra_desert_v)
  
}

# Check
list_data_day_bra %>% 
  map(., ~{
    
    unique(.x %>% group_by(site_year) %>% count() %>% pull(n))
    
  })

# Save new data (BRA + Desert)
save(list_data_day_bra, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/01_days/list_data_day_bra.rda")

# > Compute monthly data 
list_data_month_bra <- list_data_day_bra %>% 
  map(., ~{
    
    monthly_average(var_i = unique(.x$clim.var),
                    load = F, 
                    data_var_i = .x)
    
  })

# Table with monthly data
tab_month_bra <- plyr::ldply(list_data_month_bra, data.frame) %>% 
  # > rename mean columns
  rename("monthly"     = "mean_clim.value",
         "monthly_cum" = "mean_cum_clim.value") %>% 
  # > correct year 
  mutate(year = if_else(country_name == "Brazil" & month %in% c(11,12), year+1, year)) %>% 
  # > identify month id
  group_by(site_year, x, y, gridcode, country_name, year, clim.var) %>%
  arrange(year, month) %>%
  mutate(month_nb = row_number()) %>% 
  # > remove useless column
  ungroup() %>% 
  dplyr::select(-.id, -month) %>% 
  # > pivot 
  pivot_wider(names_from = c("clim.var", "month_nb"), values_from = c("monthly", "monthly_cum")) 

head(tab_month_bra)
dim(tab_month_bra)
dim(dat_coords_BRA) # should have the same number of lines 

# Save
save(tab_month_bra, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/02_month/tab_month_bra.rda")

# Table with annual averages
tab_year_bra <- tab_month_bra %>% 
  group_by(site_year, x, y, gridcode, country_name, year, Ya, Ya_ano) %>% 
  summarise(year_max_2m_temperature         =mean(c(monthly_max_2m_temperature_1, monthly_max_2m_temperature_2, monthly_max_2m_temperature_3, monthly_max_2m_temperature_4, monthly_max_2m_temperature_5, monthly_max_2m_temperature_6, monthly_max_2m_temperature_7)), 
            year_min_2m_temperature         =mean(c(monthly_min_2m_temperature_1, monthly_min_2m_temperature_2, monthly_min_2m_temperature_3, monthly_min_2m_temperature_4, monthly_min_2m_temperature_5, monthly_min_2m_temperature_6, monthly_min_2m_temperature_7)),
            year_et0                        =mean(c(monthly_et0_1, monthly_et0_2, monthly_et0_3, monthly_et0_4, monthly_et0_5, monthly_et0_6, monthly_et0_7)),
            year_surface_net_solar_radiation=mean(c(monthly_surface_net_solar_radiation_1, monthly_surface_net_solar_radiation_2, monthly_surface_net_solar_radiation_3, monthly_surface_net_solar_radiation_4, monthly_surface_net_solar_radiation_5, monthly_surface_net_solar_radiation_6, monthly_surface_net_solar_radiation_7)),
            year_total_precipitation        =mean(c(monthly_total_precipitation_1, monthly_total_precipitation_2, monthly_total_precipitation_3, monthly_total_precipitation_4, monthly_total_precipitation_5, monthly_total_precipitation_6, monthly_total_precipitation_7)),
            year_vapor_pressure_deficit     =mean(c(monthly_vapor_pressure_deficit_1, monthly_vapor_pressure_deficit_2, monthly_vapor_pressure_deficit_3, monthly_vapor_pressure_deficit_4, monthly_vapor_pressure_deficit_5, monthly_vapor_pressure_deficit_6, monthly_vapor_pressure_deficit_7))) 

dim(tab_year_bra) # N=18429

# Table with zscore of annual averages
tab_year_zscore_bra <- tab_year_bra %>% 
  ungroup(.) %>% 
  dplyr::select(site_year, starts_with("year_")) %>% 
  mutate_if(is.numeric, scale) %>% 
  pivot_longer(cols = starts_with("year_"), names_to = "clim.var", values_to = "clim.z_scores") %>% 
  group_by(site_year) %>%
  summarise(mean_year_z_scores = mean(clim.z_scores))

dim(tab_year_zscore_bra) # N=18429

# > Merge
tab_year_bra <- left_join(tab_year_bra, tab_year_zscore_bra, by = "site_year")

# > Save 
save(tab_year_bra, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/03_annual/tab_year_bra.rda")
