# -------------------------------------------------------------------------
# 
#       PREPARE DATASETS FOR THE ANALYSES 
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
# ----------------------------------

# ----------------------------------
# Load climatic predictors
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > monthly predictors 
tab_month_world <- loadRDa(paste0(data_path, "02_month/tab_month_world.rda"))
tab_month_usa   <- loadRDa(paste0(data_path, "02_month/tab_month_usa.rda"))
tab_month_bra   <- loadRDa(paste0(data_path, "02_month/tab_month_bra.rda"))

# > annual predictors 
tab_year_world <- loadRDa(paste0(data_path, "03_annual/tab_year_world.rda"))
tab_year_usa   <- loadRDa(paste0(data_path, "03_annual/tab_year_usa.rda"))
tab_year_bra   <- loadRDa(paste0(data_path, "03_annual/tab_year_bra.rda"))

# > scores 
tab_scores_M_world <- loadRDa(paste0(data_path, "04_scores/tab_scores_M_world.rda"))
tab_scores_M_usa   <- loadRDa(paste0(data_path, "04_scores/tab_scores_M_usa.rda"))
tab_scores_M_bra   <- loadRDa(paste0(data_path, "04_scores/tab_scores_M_bra.rda"))

tab_scores_D_world <- loadRDa(paste0(data_path, "04_scores/tab_scores_D_world.rda"))
tab_scores_D_usa   <- loadRDa(paste0(data_path, "04_scores/tab_scores_D_usa.rda"))
tab_scores_D_bra   <- loadRDa(paste0(data_path, "04_scores/tab_scores_D_bra.rda"))

tab_scores_fplsr_M_world <- loadRDa(paste0(data_path, "04_scores/tab_scores_fplsr_M_world.rda"))
tab_scores_fplsr_D_world <- loadRDa(paste0(data_path, "04_scores/tab_scores_fplsr_D_world.rda"))

# Irrigation fraction
fractarea_soybean_irr_high_res <- raster::raster("E:/POSTDOC INRAE/DATA/02_YIELDS/SPAM/spam2010v2r0_global_harv_area.geotiff/spam2010V2r0_global_H_SOYB_I.tif")
fractarea_soybean_irr_high_res
# class      : RasterLayer 
# dimensions : 2160, 4320, 9331200  (nrow, ncol, ncell)
# resolution : 0.083333, 0.083333  (x, y)
# extent     : -180, 179.9986, -89.99928, 90  (xmin, xmax, ymin, ymax)
# crs        : +proj=longlat +datum=WGS84 +no_defs 
# source     : spam2010V2r0_global_H_SOYB_I.tif 
# names      : spam2010V2r0_global_H_SOYB_I

# > resolution reduction (from 0.08333333 -> 0.5)
fractarea_soybean_irr_low_res <- terra::aggregate(fractarea_soybean_irr_high_res, fact=6, fun="mean")

# > add coordinates
fractarea_soybean_irr <- data.frame(irrigated_portion = values(fractarea_soybean_irr_low_res),
                                    coordinates(fractarea_soybean_irr_low_res)) %>%  
  # > round x and y to 2 digits to be consistent with yield data 
  mutate(x=round(x,2),
         y=round(y,2)) %>% 
  # proportion of area in %., we need to /100 to have %
  mutate(irrigated_portion = irrigated_portion/100) %>% 
  # transform NA into 0
  mutate(irrigated_portion =if_else(is.na(irrigated_portion)==T, 0, irrigated_portion))

summary(fractarea_soybean_irr$irrigated_portion)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000000  0.005443  0.000000 22.159389 

# ----------------------------------
# Global data set 

tab_world <- cbind(
  left_join(tab_month_world, tab_year_world, by = c("x", "y", "year","site_year", "gridcode", "country_name", "Ya", "irrigated_portion")),
  tab_scores_M_world, tab_scores_D_world, tab_scores_fplsr_M_world, tab_scores_fplsr_D_world) %>%
  left_join(., fractarea_soybean_irr, by =c("x", "y"))

dim(tab_world) # N=122229

dim(na.omit(tab_world))

save(tab_world, file = paste0(data_path, "tab_world.rda"))

# ----------------------------------
# USA

tab_usa <- cbind(
  left_join(tab_month_usa, tab_year_usa, by = c('site_year', 'x', 'y', 'gridcode', 'country_name', 'year', 'Ya', 'Ya_ano')),
  tab_scores_M_usa, tab_scores_D_usa) %>%
  left_join(., fractarea_soybean_irr, by =c("x", "y"))

dim(tab_usa) # N=37147

dim(na.omit(tab_usa)) # should be same dimension 

save(tab_usa, file = paste0(data_path, "tab_usa.rda"))
# ----------------------------------
# BRAZIL

tab_bra <- cbind(
  left_join(tab_month_bra, tab_year_bra, by = c('site_year', 'x', 'y', 'gridcode', 'country_name', 'year', 'Ya', 'Ya_ano')),
  tab_scores_M_bra, tab_scores_D_bra) %>%
  left_join(., fractarea_soybean_irr, by =c("x", "y"))

dim(tab_bra) # N=18429

dim(na.omit(tab_bra)) # should be same dimension 

save(tab_bra, file = paste0(data_path, "tab_bra.rda"))

stop()
# ----------------------------------

# > Entire dataset
tab_month_world <- tab_month
length(unique(tab_month_world$site_year)) # N=122229

# > Desert 
tab_desert <- tab_month_world %>% 
  filter(country_name == "Desert")

# > USA 
tab_month_usa <- tab_month %>% filter(country_name == "United States of America")
length(unique(tab_month_usa$site_year)) # N=29803 

# > Brazil
tab_month_bra <- tab_month %>% filter(country_name == "Brazil")
length(unique(tab_month_bra$site_year)) # N=14757

# > Add desert points in USA and BRAZIL 
n_desert_usa <- round(( round( ( length(unique(tab_month_usa$site_year))/(1-0.2) ), digits = 0) - length(unique(tab_month_usa$site_year)) ) / 35, digits = 0) ; n_desert_usa
n_desert_bra <- round(( round( ( length(unique(tab_month_bra$site_year))/(1-0.2) ), digits = 0) - length(unique(tab_month_bra$site_year)) ) / 35, digits = 0) ; n_desert_bra

# > Sample random site-years in desert 
desert_usa <- sample(unique(tab_desert$gridcode), size = n_desert_usa, replace = F)
desert_bra <- sample(unique(tab_desert$gridcode), size = n_desert_bra, replace = F)

# > Desert data 
tab_desert_usa <- tab_desert %>% 
  filter(gridcode %in% desert_usa)
tab_desert_bra <- tab_desert %>% 
  filter(gridcode %in% desert_bra)

# > Add it to the initial tables
tab_month_usa_2 <- rbind(tab_month_usa, tab_desert_usa)
tab_month_bra_2 <- rbind(tab_month_bra, tab_desert_bra)

# > Maps 
world <- ne_countries(scale = "medium", returnclass = "sf")

# World 
tab_month_world %>% 
  mutate(country_name = if_else(country_name == "Desert", "Desert", "Global")) %>% 
  mutate(country_name = factor(country_name, levels = c("Global", "Desert"))) %>% 
  distinct(x, y, country_name) %>%
  group_by(country_name) %>% 
  count()

#1 Global        2784
#2 Desert         663

ggplot(data = tab_month_world %>% 
         mutate(country_name = if_else(country_name == "Desert", "Unsuitable environmental conditions (N=663)", "Global (N=2784)")) %>% 
         distinct(x, y, country_name)) + 
  geom_sf(data = world) + 
  geom_point(aes(x = x, y = y, color = country_name), size = 0.2) + 
  scale_color_manual(values = c("darkgreen", "black"), name = "Sites") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1)) + 
  theme_map() + 
  theme(legend.position = "bottom") + ggtitle("Global dataset")

# USA
tab_month_usa_2 %>% 
  distinct(x, y, country_name) %>%
  group_by(country_name) %>% 
  count()
#1 United States of America   830
#2 Desert                     213

ggplot(data = tab_month_usa_2 %>% 
         mutate(country_name = if_else(country_name == "Desert", "Unsuitable environmental conditions (N=213)", "United States of America (N=830)")) %>% 
         distinct(x, y, country_name)) + 
  geom_sf(data = world) + 
  geom_point(aes(x = x, y = y, color = country_name), size = 0.2) + 
  scale_color_manual(values = c("red", "black"), name = "Sites") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1)) + 
  theme_map() + 
  theme(legend.position = "bottom") + ggtitle("US dataset")

# Brazil 
tab_month_bra_2 %>% 
  distinct(x, y, country_name) %>%
  group_by(country_name) %>% 
  count()
#1 Brazil         424
#2 Desert         105

ggplot(data = tab_month_bra_2 %>% 
         mutate(country_name = if_else(country_name == "Desert", "Unsuitable environmental conditions (N=105)", "Brazil (N=424)")) %>%
         distinct(x, y, country_name)) + 
  geom_sf(data = world) + 
  geom_point(aes(x = x, y = y, color = country_name), size = 0.2) + 
  scale_color_manual(values = c("blue", "black"), name = "Sites") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1)) + 
  theme_map() + 
  theme(legend.position = "bottom") + ggtitle("Brazil dataset")

