# -------------------------------------------------------------------------
#           
#       Global dataset of historical yields for major crops 1981–2016  
#       Crops: Soybean & Maize
#       Time resolution: 1981 - 2016  
#         
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Loading packages 
library(tidyverse) ; library(stringr)
library(cowplot)
library(terra) ; library(raster)
library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(osmextract)

# -------------------------------------------------------------------------
# Home-made functions
source("E:/POSTDOC INRAE/ANALYSES/00_TOOLS/00_Functions.R")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 1. Loading yield .nc files

# References of yield data
# - Iizumi, T., Sakai, T. The global dataset of historical yields for major crops 1981–2016. Sci Data 7, 97 (2020). https://doi.org/10.1038/s41597-020-0433-7
# - Iizumi, T. Global dataset of historical yields v1.2 and v1.3 aligned version. (2019) PANGAEA, https://doi.org/10.1594/PANGAEA.909132 

# > path to maize (major crop) and soybean yield data
path_soybean <- "E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/gdhy_v1.2_v1.3_20190128/soybean"

# > extract all the names of the files
filenames <- list.files(c(path_soybean), full.names = TRUE)

# > split the files among the different variables 
filetable <- data.frame(filename = filenames) %>% 
  # > add year 
  mutate(year = substr(filename, nchar(filename)-7, nchar(filename)-4)) %>% 
  # > add crop
  mutate(crop = case_when(
    str_detect(string = filename, pattern = "soybean") == T ~ "soybean",
    str_detect(string = filename, pattern = "maize")   == T ~ "maize"
  ))

# > extract data for each crop
Ya_full <- filetable %>% 
  split(.$crop) %>% 
  map_dfr(., ~{ 
    
    # > for a crop, extract the data for each year
    .x %>% 
      split(.$year) %>% 
      map_dfr(., ~{
        
        # > labels for month & year
        year_i  <- unique(.x$year)
        
        # > export the data from .nc file
        raster_i <- rast(paste0(.$filename))
        
        # > transform it into a data.frame
        tab_raster_i <- as.data.frame(raster_i, xy=T) 
        
        tab_raster_i %>% 
          # > add date
          mutate(year = year_i) %>% 
          # > rename var as yield
          rename("Ya"="var")
        
      }, .id = "year")
    
  }, .id = "crop")

# > identify individual pairs of long and lat
dat_coords <- Ya_full %>% 
  distinct(x, y) %>% 
  mutate(x=if_else(x>180, x-360, x))

# > attribute country and world region to each combination of long and lat
dat_coords_country <- coords2continent(points = dat_coords)

# > merge to the full yield dataset
Ya_full_coords <- left_join(Ya_full%>% 
                              mutate(x=if_else(x>180, x-360, x)), 
                            dat_coords_country, by = c("x", "y"))


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 2. Remove yield trend
# To avoid any confusion with technological progress, soybean yield data
# are detrended in order to remove the increasing trends of soybean yield 
# time series due to improved cultivars and technological progress 

# > Check if any NA
dim(na.omit(Ya_full_coords)) ; dim(Ya_full_coords) 

# > remove NAs (normally 0 NA) # 
Ya_full_coords <- na.omit(Ya_full_coords) 

# > detrend yield data for each crop 
Ya_full_coords.notrend <- Ya_full_coords %>% 
  # > identify pairs of long and lat 
  unite("gridcode", x:y, remove=F) %>% 
  # > split by crop
  split(.$crop) %>% 
  map_dfr(., ~ {
    
    data.yield.notrend_i <- .x
    
    # > loop on grid cell (combination of long and lat)
    gridcode <- unique(data.yield.notrend_i$gridcode)
    
    # > compute trend + detrend
    for (i in gridcode){
      # loop on gridcode starts with subset on gridcode
      tmp <- subset(data.yield.notrend_i, gridcode==i)
      
      # fitting the trend needs at least 4 values
      if(nrow(tmp) >= 4)
      {
        # fit smoothing spline
        mod <- smooth.spline(tmp$year, tmp$Ya)
        
        # add yearly residuals from smoothing splines to 
        # expected value from smoothing spline
        data.yield.notrend_i$Ya[data.yield.notrend_i$gridcode==i] <- predict(mod)$y[length(tmp$Ya)] + residuals(mod)
        # yield anomaly
        data.yield.notrend_i$Ya_ano[data.yield.notrend_i$gridcode==i] <- residuals(mod)
      }
      
      # if less than 4 values, return NA
      if(nrow(tmp) < 4)
      {
        data.yield.notrend_i$Ya[data.yield.notrend_i$gridcode==i] <-  NA
      }
      
    }
    
    data.yield.notrend_i
    
  }, .id = "crop")


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 3. Identify the grid-cells in which crops cover at least 1% of the surface

# Reference of the crop mask used 
# Monfreda, C., N. Ramankutty, and J. A. Foley (2008), Farming the planet: 2. Geographic
# distribution of crop areas, yields, physiological types, and net primary production in the year
# 2000, Global Biogeochem. Cycles, 22, GB1022, doi:10.1029/2007GB002947.

# > path to crop mask
fractarea_soybean_high_res <- raster::raster("E:/POSTDOC INRAE/DATA/02_YIELDS/MONFREDA/monfreda_soybean/soybean_HarvAreaYield_Geotiff/soybean_HarvAreaYield_Geotiff/soybean_HarvestedAreaFraction.tif") 

# > resolution reduction (from 0.08333333 -> 0.5)
fractarea_soybean_low_res <- terra::aggregate(fractarea_soybean_high_res, fact=6, fun="mean")

# > add coordinates
fractarea_soybean <- data.frame(soybean = values(fractarea_soybean_low_res),
                                coordinates(fractarea_soybean_low_res)) %>% 
  # > round x and y to 2 digits to be consistent with yield data 
  mutate(x=round(x,2),
         y=round(y,2)) %>% 
  # > in MONFREDA data set, area is expressed as proportion
  # to express proportion of area in %, we need to *100
  mutate(crop_area = soybean*100) %>% 
  dplyr::select(-soybean)

# > merge area + yields 
Ya_full_coords.notrend_area <- Ya_full_coords.notrend %>% 
  left_join(., fractarea_soybean, by=c("x", "y")) %>% 
  # > identify grid-cells with more than 1% of maize or soybean harvested
  mutate(more_than_1_percent    = if_else(crop_area > 1,   1,   0))

# > check for NA values
Ya_full_coords.notrend_area %>% filter(is.na(crop_area) == T) # 0 lines: OK 

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 4. Databases

# > Soybean: select grid-cells with 
# - area of soybean of at least 1%
# - located in Argentina, Brazil, Canada, China, India, Italy, USA
data.yield.notrend.soybean <- Ya_full_coords.notrend_area %>% 
  # > keep data for soybean
  filter(crop == "soybean") %>% 
  # > identify grid-cells producing >1% soybean from Argentina, Brazil, Canada, China, India, Italy, USA
  mutate(to_keep = if_else(more_than_1_percent == 1 & country_code %in% c("ARG", "BRA", "CAN", "CHN", "IND", "ITA", "USA"), 1, 0)) %>% 
  # > remove useless variables
  dplyr::select(-continent, -more_than_1_percent)

# > check no data in other countries
data.yield.notrend.soybean$country_name <- droplevels(data.yield.notrend.soybean$country_name)
table(data.yield.notrend.soybean$to_keep, data.yield.notrend.soybean$country_name)

# > check area in selected grid-cells
summary(data.yield.notrend.soybean$crop_area[which(data.yield.notrend.soybean$to_keep==0 & data.yield.notrend.soybean$country_code %in% c("ARG", "BRA", "CAN", "CHN", "IND", "ITA", "USA"))]) # > should be <1%
summary(data.yield.notrend.soybean$crop_area[which(data.yield.notrend.soybean$to_keep==1)]) # > should be >1%

# > Split data in 3 datasets 
#   - Global 
#   - US
#   - Brazil 

soybean_world <- data.yield.notrend.soybean %>% 
  filter(to_keep == 1)  
dim(soybean_world) # N=98361
soybean_world %>% distinct(x, y) %>% dim(.) # N = 2784

soybean_usa <- data.yield.notrend.soybean %>% 
  filter(to_keep == 1) %>% 
  filter(country_code == "USA") 
dim(soybean_usa) # N=29803
soybean_usa %>% distinct(x, y) %>% dim(.) # N = 830

soybean_bra <- data.yield.notrend.soybean %>% 
  filter(to_keep == 1) %>% 
  filter(country_code == "BRA")
dim(soybean_bra) # N=14757
soybean_bra %>% distinct(x, y) %>% dim(.) # N = 424

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 5. Add locations with yield equals zero in desert and arctic areas
# In this section we use the Koppen-Geiger climate classification
# to identify desert and arctic areas in which we will set maize and 
# soybean yields equal to zero

# > read Koppen-Geiger classification in raster format
# data are from http://koeppen-geiger.vu-wien.ac.at/present.htm
# August 2022 version
KGmap <- raster("C:/Users/benni/Documents/Post doc/Map_KG-Global/Map_KG-Global/KG_1986-2010.grd") ; KGmap

#class      : RasterLayer 
#dimensions : 2160, 4320, 9331200  (nrow, ncol, ncell)
#resolution : 0.08333333, 0.08333333  (x, y)
#extent     : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#source     : KG_1986-2010.grd 
#names      : layer 
#values     : 1, 32  (min, max)

# > resolution reduction (from 0.083333 -> 0.5)
KGmap <- terra::aggregate(KGmap, fact=6, fun="modal")

# > define climate zones which we consider as desert
desert <- c(7, 8, 20, 21, 30, 31)

# > check how many data points are in these zones (less than 1%)
soybean_world %>% 
  distinct(x, y) %>% 
  left_join(., as.data.frame(KGmap, xy=T), by=c("x", "y")) %>% 
  mutate(layer=round(layer, digits = 0)) %>% 
  group_by(layer) %>% 
  count() %>%
  ungroup() %>% 
  mutate(tot = sum(n), 
         freq = (n/tot)*100) %>% 
  filter(layer %in% desert)

#  layer     n   tot   freq
#1     8     2  2784 0.0718
#2    20     1  2784 0.0359
#3    31     1  2784 0.0359

# --------------
# > Define global desert 
KGmap.desert <- KGmap
KGmap.desert[!(KGmap.desert %in% desert)] <- 0
KGmap.desert <- crop(KGmap.desert, extent(-180, 180,-60, 90)) # > remove antarctic

# > prepare data for sampling 
list_data <- list(
  "01_USA" = list(data_sample = soybean_usa %>% 
                    dplyr::select(crop, gridcode, x, y, Ya, Ya_ano, year, region, country_name, country_code),
                  KGmap.desert.test = KGmap.desert, 
                  name = "01_USA",
                  name_save = "01_USA"),
  "02_BRA" = list(data_sample = soybean_bra %>% 
                    dplyr::select(crop, gridcode, x, y, Ya, Ya_ano, year, region, country_name, country_code),
                  KGmap.desert.test = KGmap.desert, 
                  name = "02_BRA",
                  name_save = "02_BRA"),
  "03_WORLD" = list(data_sample = soybean_world %>% dplyr::select(crop, gridcode, x, y, Ya, Ya_ano, year, region, country_name, country_code),
                    KGmap.desert.test = KGmap.desert, 
                    name = "03_WORLD", 
                    name_save = "03_WORLD"))

# > sample desert points for each zone
list_data_desert <- list_data %>% 
  map(., ~{ 
      
      # > apply function
      data_sample_desert_i <- sample.desert(perc_sample       = 0.2, 
                                            data_sample       = .x$data_sample,
                                            KGmap.desert.test = .x$KGmap.desert.test, 
                                            desert            = desert) 
      
      data_sample_desert_i_coords <- coords2continent(points = subset(data_sample_desert_i, select=c(x, y))) %>%  
        unite("gridcode", x, y, remove = F) %>%
        # > identify duplicates
        mutate(duplicates = if_else(gridcode %in% unique(soybean_world$gridcode), 1, 0))
      
      # > final set of points (excluding the duplicates)
      gridcode_desert <- data_sample_desert_i_coords %>% 
        filter(duplicates == 0) %>% 
        pull(gridcode)
      
      # > expand data to get all years needed
      data_sample_desert <- expand.grid(
        crop     = unique(.x$data_sample$crop),
        gridcode = gridcode_desert,
        year     = 1981:2016) %>% 
        left_join(., data_sample_desert_i_coords, by="gridcode") %>% 
        mutate(Ya = 0, Ya_ano = 0,
               country_name = "Desert") %>% 
        arrange(crop, gridcode) %>% 
        dplyr::select(crop, gridcode, x, y, Ya, Ya_ano, year, region, country_name, country_code)
      
      # > merge with initial dataset 
      data_fin <- rbind(.x$data_sample, data_sample_desert)
      
      # > check that there is no duplicate
      testthat::expect_equal(0, data_fin %>% 
                               distinct(gridcode, country_name) %>% 
                               group_by(gridcode, country_name) %>%
                               summarise(n=n()) %>% 
                               filter(n != 1) %>% nrow(.))
      
      # > out
      list(data_sample        = .x$data_sample,
           data_desert        = data_sample_desert_i_coords,
           data_sample_desert = data_fin,
           KGmap.desert.test  = .x$KGmap.desert.test,
           name               = .x$name,
           name_save          = .x$name_save)

  })

# > Check if the number of desert points represent ~25% of the final dataset
list_data_desert %>% map_dfr(., ~{ 
    
    n_desert <- nrow(.x$data_sample_desert[which(.x$data_sample_desert$country_name == "Desert"),])
    n_init   <- nrow(.x$data_sample)
    n_total  <- nrow(.x$data_sample_desert)
    
    data.frame(
      dataset            = .x$name,            # > name of the data set
      n_total            = n_total,            # > size of the initial data + sampled desert points 
      n_init             = n_init,             # > size of the initial data, excluding the desert points 
      n_desert           = n_desert) %>%       # > size of sampled points
      mutate(prop_desert = n_desert / n_total) # > which proportion of the dataset is in desert --> should be close to 0.25
    
})

#  dataset n_total n_init n_desert prop_desert
#1   01_USA   37147  29803     7344   0.1977010
#2   02_BRA   18429  14757     3672   0.1992512
#3 03_WORLD  122229  98361    23868   0.1952728
# --> prop_desert should be around 0.20

# > Check the number of sites in total, and in desert in each region
list_data_desert %>% map_dfr(., ~{ 
    
    n_desert <- .x$data_sample_desert[which(.x$data_sample_desert$country_name == "Desert"),] %>% 
      distinct(x, y) %>% nrow(.)
    n_init   <- .x$data_sample %>% 
      distinct(x, y) %>% nrow(.)
    n_total  <- .x$data_sample_desert %>% 
      distinct(x, y) %>% nrow(.)
    
    data.frame(
      dataset            = .x$name,            # > name of the data set
      n_total            = n_total,            # > size of the initial data + sampled desert points 
      n_init             = n_init,             # > size of the initial data, excluding the desert points 
      n_desert           = n_desert) %>%       # > size of sampled points
      mutate(prop_desert = n_desert / n_total) # > which proportion of the dataset is in desert --> should be close to 0.25

})

#  dataset n_total n_init n_desert prop_desert
#1   01_USA    1034    830      204   0.1972921
#2   02_BRA     526    424      102   0.1939163
#3 03_WORLD    3444   2784      660   0.1916376

# > save datasets
yield_usa <- list_data_desert$'01_USA'$data_sample_desert
yield_bra <- list_data_desert$'02_BRA'$data_sample_desert
yield_world <- list_data_desert$'03_WORLD'$data_sample_desert

save(yield_usa, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_1_yield_usa.rda")
save(yield_bra, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_2_yield_bra.rda")
save(yield_world, file = "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_data/01_3_yield_world.rda")
