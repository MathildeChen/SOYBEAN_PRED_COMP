# -------------------------------------------------------------------------
# 
#       SUPPLEMENTARY FIGURES FOR THE ARTICLE
#       Author: M. Chen, Inrae, 2023  
# 
# -------------------------------------------------------------------------

# ----------------------------------
# > Packages
library(tidyverse) ; library(stringr) ; library(lubridate) ; library(CCMHr)
# > plots
library(cowplot) ; library(wesanderson) ; library(RColorBrewer)
# > maps 
library(terra) ; library(raster) ; library("rnaturalearth") ; library("rnaturalearthdata") ; library(sf) ; library(sp) ; library(rworldmap) ; library(rmapshaper) ; library(tidygeocoder) 

# ----------------------------------
# Graphical features 

# > Color palettes
# > gradient blue - yellow - red
pal <- wesanderson::wes_palette("Zissou1", 6, type="continuous")

# > countries
#pal_country <- RColorBrewer::brewer.pal(7, "Set1")
pal_country <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33")# reorder to have desert in yellow

# > Map for world and continents 
world         <- ne_countries(scale = "medium", returnclass = "sf")
africa        <- ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")
asia          <- ne_countries(continent = "Asia", scale = "medium", returnclass = "sf")
north_america <- ne_countries(continent = "North America", scale = "medium", returnclass = "sf")
south_america <- ne_countries(continent = "South America", scale = "medium", returnclass = "sf")
americas      <- rbind(north_america, south_america)
europe        <- ne_countries(continent = "Europe", scale = "medium", returnclass = "sf")
brazil       <- ne_countries(country = "Brazil", scale = "medium", returnclass = "sf")

# > path to save the images
save_path <- "E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/FIGURES/"

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
                               "vapor_pressure_deficit"     ="vpd"),
         clim.var_abb2 = recode(clim.var, 
                               "min_2m_temperature"         ="Min temp.",
                               "max_2m_temperature"         ="Max temp.",
                               "et0"                        ="ET0",
                               "surface_net_solar_radiation"="Solar rad.",
                               "total_precipitation"        ="Total prec.",
                               "vapor_pressure_deficit"     ="VPD"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature (°C)",
                               "max_2m_temperature"         ="Maximum temperature (°C)",
                               "et0"                        ="Evapotransp. ref (mm/day)",
                               "surface_net_solar_radiation"="Net solar radiations (MJ/m²)",
                               "total_precipitation"        ="Total precipitations (mm)",
                               "vapor_pressure_deficit"     ="Vapor pressure deficit")) %>% 
  mutate(clim.var_lab2 = recode(clim.var, 
                                "min_2m_temperature"         ="minimum temperature\n(°C)",
                                "max_2m_temperature"         ="maximum temperature\n(°C)",
                                "et0"                        ="evapotranspiration of reference\n(mm/day)",
                                "surface_net_solar_radiation"="net solar radiations\n(MJ/m²)",
                                "total_precipitation"        ="total precipitations\n(mm)",
                                "vapor_pressure_deficit"     ="vapor pressure deficit\n"))

# ----------------------------------
# Data 
data_path <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"

# > Yield data 
load("E:/POSTDOC INRAE/DATA/02_YIELDS/GDHY_v1.3/yield_full.rda")
data_soybean_init <- data_crop_init %>% 
  filter(crop != "maize")

# > Load daily climatic data 
#load(paste0(data_path, "01_days/list_data_day_usa.rda"))
#load(paste0(data_path, "01_days/list_data_day_bra.rda"))
#load(paste0(data_path, "01_days/list_data_day_world.rda"))

# > Yield and all climatic aggregated predictors 
load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))
load(paste0(data_path, "tab_world.rda"))

# > Load cumulated explained variables (for figure 3)
load(paste0(data_path, "tab_var_dim_red_full.rda"))

# > Load table with models performances
load(paste0(data_path, "tab_perf_models.rda"))

# > Predictions 
load(paste0(data_path, "tab_preds.rda"))

# > Predictions in temperature increase
load(paste0(data_path, "05_sensi_cc/tab_preds_cc.rda"))

# ----------------------------------
# Functions 
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

# ----------------------------------
# > Supplementary maps 
# USA
tab_usa %>% 
  distinct(x, y, country_name) %>% 
  group_by(country_name) %>% count()
#1 United States of America   830
#2 Desert                     204

map_usa <- tab_usa %>% 
  distinct(x, y, country_name) %>%
  filter(country_name=="United States of America") %>% 
  ggplot(.) + 
  geom_sf(data=world) + 
  geom_point(aes(x = x, y = y), 
             size=0.25, 
             color = "#E41A1C") + 
  lims(x=c(-120, -57), y=c(25,50)) + 
  theme_map() ; map_usa

map_usa2 <- ggplot(data = tab_usa %>% 
                     mutate(country_name = if_else(country_name == "Desert", 
                                                   "Selected sites in geographical areas unsuitable for soybean production (N=204)", 
                                                   "United States of America (N=830)")) %>% 
                     mutate(country_name = factor(country_name, levels = c("United States of America (N=830)", 
                                                                           "Selected sites in geographical areas unsuitable for soybean production (N=204)"))) %>% 
                     distinct(x, y, country_name)) + 
  geom_sf(data = world) + 
  geom_point(aes(x = x, y = y, color = country_name), size = 0.2) + 
  scale_color_manual(values = c("#E41A1C", "black"), name = "") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=1), ncol = 1)) + 
  theme_map() + 
  theme(plot.margin = unit(c(0,1,0,0), "cm")) +
  theme(legend.position = "bottom", 
        legend.text = element_text(size=8)) ; map_usa2

plot_grid(plot_grid(map_usa, ggplot() + theme_void(), ncol = 1, rel_heights = c(0.75, 0.25)), 
          map_usa2, nrow = 1, rel_widths = c(0.35, 0.65), align = "t", axis = "tlbr", labels = c("a.", "b."))

ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/00_MAPS/USA_map.png"), 
       width = 11, height=4, dpi=300, bg="white")

# Brazil
tab_bra %>% 
  distinct(x, y, country_name) %>%
  group_by(country_name) %>% 
  count()
#1 Brazil         424
#2 Desert         102

map_bra <- tab_bra %>% 
  distinct(x, y, country_name) %>%
  filter(country_name=="Brazil") %>% 
  ggplot(data=.) + 
  geom_sf(data=brazil) + 
  geom_point(aes(x = x, y = y), 
             size=0.25, 
             color = "#377EB8") + 
  theme_map() ; map_bra

map_bra2 <- ggplot(data = tab_bra %>% 
                     mutate(country_name = if_else(country_name == "Desert", 
                                                   "Selected sites in geographical areas unsuitable for soybean production (N=102)", 
                                                   "Brazil (N=424)")) %>%
                     mutate(country_name = factor(country_name, levels = c("Brazil (N=424)", 
                                                                           "Selected sites in geographical areas unsuitable for soybean production (N=102)"))) %>% 
                     distinct(x, y, country_name)) + 
  geom_sf(data = world) + 
  geom_point(aes(x = x, y = y, color = country_name), size = 0.2) + 
  scale_color_manual(values = c("#377EB8", "black"), name = "") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=1), ncol = 1)) + 
  theme_map() + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size=8))

plot_grid(plot_grid(map_bra, ggplot() + theme_void(), ncol = 1, rel_heights = c(0.75, 0.25)), 
          map_bra2, nrow = 1, rel_widths = c(0.35, 0.65), align = "t", axis = "tlbr", labels = c("c.", "d."))

ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/00_MAPS/BRA_map.png"), 
       width = 11, height=4, dpi=300, bg="white")

# Full data 
map_world_full <- ggplot() + 
  geom_sf(data=world) + 
  geom_point(data = tab_world %>% 
               ungroup() %>%
               distinct(x, y, country_name) %>% 
               mutate(country_name  = as.character(country_name)) %>% 
               mutate(country_color = case_when(
                 country_name == "United States of America"~ "United States of America (N=830)",
                 country_name == "Brazil"                  ~ "Brazil (N=424)",
                 #country_name == "Desert"                  ~ "Desert",
                 TRUE ~ "World (N=3447)")),# %>% mutate(is_desert     = if_else(country_name == "Desert", "Desert, i.e. 0 yield sites", "Other")), 
             aes(x = x, y = y, 
                 #shape = as.factor(is_desert),
                 color = country_color), 
             size=0.25) + 
  theme_map() + 
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("#377EB8","#E41A1C","black"), name = "Analysis") + 
  #scale_color_manual(values = c("#377EB8","darkgrey", "#E41A1C","black"), name = "Analysis") + 
  #scale_shape_manual(values = c(0,15), "Desert points") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=1), ncol = 1), 
         shape = guide_legend(order = 2, override.aes = list(size=1), ncol = 1)) ; map_world_full

ggsave(map_world_full, filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/02_Supplementary_Figure_map.png"), 
       width = 14, height=9, dpi=300, bg="white")


# ----------------------------------
# Soybean yield and yield anomaly distribution
data_soybean_init$x <- as.numeric(as.character(data_soybean_init$x))
data_soybean_init$y <- as.numeric(as.character(data_soybean_init$y))
data_soybean_init$year <- as.numeric(as.character(data_soybean_init$year))

data_Ya <- rbind( 
  # > global dataset
  left_join(tab_world %>% dplyr::select(gridcode, x, y, year, country_name, Ya, Ya_ano) %>% mutate(year = as.numeric(as.character(year))), 
            data_soybean_init[,-1], 
            by = c("gridcode", "x", "y", "year", "country_name")) %>% 
    mutate(analyse_set = "Global (N=122229)"),
  
  # > USA
  left_join(tab_usa %>% dplyr::select(gridcode, x, y, year, country_name, Ya, Ya_ano), 
            data_soybean_init[,-1], 
            by = c("gridcode", "x", "y", "year", "country_name")) %>% 
    mutate(analyse_set = "United States of America (N=29803)"),
  
  # > BRAZIL
  left_join(tab_bra %>% dplyr::select(gridcode, x, y, year, country_name, Ya, Ya_ano), 
            data_soybean_init[,-1], 
            by = c("gridcode", "x", "y", "year", "country_name")) %>% 
    mutate(analyse_set = "Brazil (N=14757)"))  %>%
  mutate(Ya_init = if_else(is.na(Ya_init)==T, 0, Ya_init),
         Ya_detrented = if_else(is.na(Ya_detrented)==T, 0, Ya_detrented))

# check
testthat::expect_equal(data_Ya$Ya, data_Ya$Ya_detrented)

data_Ya <- data_Ya %>%
  dplyr::select(-Ya) %>%
  gather(key="Yield_type", value="Yield_value", Ya_init, Ya_detrented) %>%
  mutate(Yield_type = recode(Yield_type, 
                             "Ya_init"      = "before detrending", 
                             "Ya_detrented" = "after detrending")) %>% 
  mutate(Yield_type = factor(Yield_type, levels = c("before detrending", "after detrending"))) %>% 
  mutate(analyse_set = factor(analyse_set, levels = c("Global (N=122229)", 
                                                      "United States of America (N=29803)",
                                                      "Brazil (N=14757)"))) 

# Plot for each country
for(country in unique(data_Ya$analyse_set))
{
  
  data_Ya_country_i <- data_Ya %>% 
    filter(analyse_set==country)
  
  # > Yield per year and per country (5th, 50th, and 95th percentiles)
  p1 <- data_Ya_country_i %>% 
    mutate(year = as.numeric(as.character(year))) %>% 
    group_by(analyse_set, Yield_type, year) %>%
    summarise(perc05 = quantile(Yield_value, probs = c(0.05), na.rm=T),
              median = median(Yield_value, na.rm=T),
              perc95 = quantile(Yield_value, probs = c(0.95), na.rm=T)) %>% 
    ggplot(data = .) +
    geom_hline(aes(yintercept = median(data_Ya_country_i[which(data_Ya_country_i$Yield_type=="after detrending"),]$Yield_value, na.rm=T)),
               color = "grey", lty=2) + 
    geom_ribbon(aes(x= year, 
                    ymin = perc05, 
                    ymax = perc95, 
                    fill = Yield_type),
                alpha = 0.2) + 
    geom_line(aes(x= year, y = median, 
                  color = Yield_type)) +
    theme_cowplot() + 
    theme(legend.position = "bottom",
          title = element_text(size=10), 
          axis.text = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size=10), 
          legend.text = element_text(size=10),
          plot.caption.position = "plot", 
          plot.caption = element_text(hjust = 0)) + 
    lims(y = c(0,6)) + 
    scale_color_manual(values = c("#FF7F00", "black"), name="Soybean yield") + 
    scale_fill_manual(values = c("#FF7F00", "black"), name="Soybean yield") + 
    labs(x = "Year", y = "Yield (tons/hectare)", 
         caption = "Note: Dotted line represents the median reference yield, i.e. the yield of the latest year available for each site, used for detrending") + 
    ggtitle("a. Soybean yield (median, 5th and 95th percentiles)") ; p1
  
  # > density of yield across all years and sites 
  p2 <- data_Ya_country_i %>% 
    ggplot(data = .) +
    geom_vline(aes(xintercept = median(data_Ya_country_i[which(data_Ya_country_i$Yield_type=="after detrending"),]$Yield_value, na.rm=T)),
               color = "grey", lty=2) + 
    geom_density(aes(x = Yield_value, fill = Yield_type, color=Yield_type), 
                 alpha=0.2) + 
    theme_cowplot() + 
    theme(legend.position = "bottom",
          title = element_text(size=10), 
          axis.text  = element_text(size=10),
          axis.title = element_text(size=10),
          legend.text = element_text(size=10)) + 
    scale_color_manual(values = c("#FF7F00", "black"), name="Soybean yield") + 
    scale_fill_manual(values = c("#FF7F00", "black"), name="Soybean yield") + 
    #guides(color = guide_legend(order = 1, override.aes = list(size=1), ncol = 1, title = "Soybean yield"), 
    #       fill  = guide_none())  + 
    labs(x="Yield (tons/hectare)", y = "Density", caption = "") + 
    ggtitle("") ; p2
  
  # > Yield anomalies
  p3 <- data_Ya_country_i %>% 
    filter(Yield_type=="after detrending") %>% 
    mutate(year = as.numeric(as.character(year))) %>% 
    group_by(analyse_set, year) %>%
    summarise(perc05 = quantile(Ya_ano, probs = c(0.05), na.rm=T),
              median = median(Ya_ano, na.rm=T),
              perc95 = quantile(Ya_ano, probs = c(0.95), na.rm=T)) %>% 
    ggplot(data = .) +
    geom_hline(yintercept = 0, color = "grey", lty=2) + 
    geom_ribbon(aes(x= year, 
                    ymin = perc05, 
                    ymax = perc95), 
                fill = "red",
                alpha = 0.2) + 
    geom_line(aes(x= year, y = median), 
              color = "red") +
    theme_cowplot() + 
    theme(legend.position = "bottom",
          title = element_text(size=10), 
          axis.text = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size=10), 
          legend.text = element_text(size=10)) + 
    lims(y = c(-1.5, 1.5)) + 
    labs(x = "Year", y = "Yield anomalie (tons/hectare)") + 
    ggtitle("b. Soybean yield anomalie (median, 5th and 95th percentiles)") ; p3
  
  # > boxplot of yield across all years and sites 
  p4 <- data_Ya_country_i %>% 
    filter(Yield_type=="after detrending") %>% 
    mutate(year = as.numeric(as.character(year))) %>% 
    ggplot(data = .) +
    geom_vline(xintercept = 0, color = "grey", lty=2) + 
    geom_density(aes(x = Ya_ano), 
                 fill = "red", color="red", alpha=0.2) + 
    theme_cowplot() + 
    theme(legend.position = "bottom",
          title = element_text(size=10), 
          axis.text  = element_text(size=10),
          axis.title = element_text(size=10),
          legend.text = element_text(size=10)) + 
    labs(x = "Yield anomalie (tons/hectare)", y = "Density") + 
    ggtitle("") ; p4
  
  ya_supp <- plot_grid(p1,p2, 
                       p3,p4, 
                       nrow=2, rel_widths = c(0.6, 0.4), rel_heights = c(0.6, 0.4),
                       align = "hv", axis = "l") ; ya_supp
  
  # > save
  ggsave(ya_supp, filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/01_YIELDS/Supplementary_Figure_yields_", country,".png"), 
         width = 11, height=8, dpi=300, bg="white")
  
}

# Other (not used for the moment)
# > Annual yield distribution 
# WORLD excluding desert 
ya_des <- tab_year_world %>% 
  filter(country_name != "Desert") %>% 
  ungroup() %>% 
  distinct(x, y, year, Ya) %>% 
  ggplot(., aes(x = as.numeric(as.character(year)), y = Ya, group=year)) + 
  geom_hline(yintercept = mean(tab_year_world[which(tab_year_world$country_name!="Desert"),]$Ya), 
             color = "black", linetype = 2, linewidth = 1) + 
  geom_boxplot(#outlier.shape = NA, 
    outlier.size = 0.1, width=0.75, alpha = 0.4) + 
  theme_cowplot() + 
  theme(title = element_text(size=10), 
        axis.text = element_text(size=10)) +
  scale_x_continuous(breaks = seq(1980,2020,by=5), 
                     labels = seq(1980,2020,by=5), 
                     limits = c(1980,2017)) + 
  lims(y = c(0,7)) +
  ggtitle(paste0("a. Global, excluding sites with unsuitable climatic conditions for soybean production (N=", nrow(tab_year_world %>% filter(country_name != "Desert")), ")")) +
  labs(x = "Year", y = "Yield\n(tons/hectare)")  ; ya_des

# WORLD
ya_wor <- tab_year_world %>% 
  #filter(country_name != "Desert") %>% 
  ungroup() %>% 
  distinct(x, y, year, Ya) %>% 
  ggplot(., aes(x = as.numeric(as.character(year)), y = Ya, group=year)) + 
  geom_hline(yintercept = mean(tab_year_world$Ya), 
             color = "black", linetype = 2, linewidth = 1) + 
  geom_boxplot(#outlier.shape = NA, 
    outlier.size = 0.1, width=0.75, alpha = 0.4) + 
  theme_cowplot() + 
  theme(title = element_text(size=10), 
        axis.text = element_text(size=10),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = seq(1980,2020,by=5), 
                     labels = seq(1980,2020,by=5), 
                     limits = c(1980,2017)) + 
  lims(y = c(0,7)) +
  ggtitle(paste0("a. Global (N=", nrow(tab_year_world), ")")) +
  labs(x = "Year", y = "Yield\n(tons/hectare)") 

# USA
ya_usa <- tab_year_usa %>% 
  ungroup() %>% 
  distinct(x, y, year, Ya) %>% 
  ggplot(., aes(x = as.numeric(as.character(year)), y = Ya, group=year)) + 
  geom_hline(yintercept = mean(tab_year_usa$Ya), color = "#E41A1C", linetype = 2, linewidth = 1) + 
  #geom_hline(yintercept = mean(tab_year_world$Ya), color = "black", linetype = 2) + 
  geom_boxplot(#outlier.shape = NA, 
    outlier.size = 0.1,
    width=0.75, fill = "#E41A1C", alpha = 0.5) + 
  theme_cowplot() + 
  theme(#axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    title = element_text(size=10), 
    axis.text = element_text(size=10)) +
  scale_x_continuous(breaks = seq(1980,2020,by=5), 
                     labels = seq(1980,2020,by=5), 
                     limits = c(1980,2017)) + 
  lims(y = c(0,7)) +
  ggtitle(paste0("b. United-States of America (N=", nrow(tab_year_usa), ")")) +
  labs(x = "Year", y = "Yield\n(tons/hectare)") ; ya_usa

# BRAZIL
ya_bra <- tab_year_bra %>% 
  ungroup() %>% 
  distinct(x, y, year, Ya) %>% 
  ggplot(., aes(x = as.numeric(as.character(year)), y = Ya, group=year)) + 
  geom_hline(yintercept = mean(tab_year_bra$Ya), color = "#377EB8", linetype = 2, linewidth = 1) + 
  #geom_hline(yintercept = mean(tab_year_world$Ya), color = "black", linetype = 2) + 
  geom_boxplot(#outlier.shape = NA, 
    outlier.size = 0.1,
    width=0.75, fill = "#377EB8", alpha = 0.5) + 
  theme_cowplot() + 
  theme(#axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    title = element_text(size=10), 
    axis.text = element_text(size=10)) +
  scale_x_continuous(breaks = seq(1980,2020,by=5), 
                     labels = seq(1980,2020,by=5), 
                     limits = c(1980,2017)) + 
  lims(y = c(0,7)) +
  ggtitle(paste0("c. Brazil (N=", nrow(tab_year_bra), ")")) +
  labs(x = "Year", y = "Yield\n(tons/hectare)")  ; ya_bra

ya_dist <- plot_grid(
  ya_wor, ya_usa, ya_bra, 
  #rel_heights = c(0.33,0.33, 0.4),
  ncol = 1
)

#ggsave(ya_dist, filename=paste0(save_path, "Figure_2.png"), width = 9, height=10, dpi=300, bg="white")

# ------------------------------------------------------------
# CLIMATIC VARIABLES 

# > Daily climatic data 
#   Median, 5th and 95th percentiles of cumulative 
clim_data <- list("United-States of America (N=29803)"= list_data_day_usa,
                  "Brazil (N=14757)"                   = list_data_day_bra, 
                  "Global (N=122229)"                  = list_daily_data_world
                  ) %>% 
  map_dfr(., ~{
    
    .x %>%
      map_dfr(., ~{
        
        # > Variable 
        var_i <- unique(.x$clim.var)
        var_lab_i <- vars_names[which(vars_names$clim.var==var_i), "clim.var_lab2"]
        
        .x %>% 
          # Some years have more than 212 days 
          filter(day_of_year < 212) %>% 
          # Compute 5th, 95th percentiles and median of climatic feature 
          group_by(clim.var, day_of_year) %>% 
          summarise(perc05 = quantile(cum_clim.value, probs = c(0.05)),
                    median = median(cum_clim.value),
                    perc95 = quantile(cum_clim.value, probs = c(0.95)))
        
        
      })
    
  }, .id = "Country")

rm(list_data_day_usa, list_data_day_bra, list_daily_data_world)

# > Monthly climatic data 
clim_data_month <- rbind(
  tab_world %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit")) %>% 
    mutate(Country = "Global (N=122229)"), 
  tab_usa %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit")) %>% 
    mutate(Country = "United-States of America (N=29803)"), 
  tab_bra %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit"))  %>% 
    mutate(Country = "Brazil (N=14757)")) %>% 
  pivot_longer(cols = starts_with("monthly"), values_to = "clim.value", names_to = "clim.var.month") %>% 
  mutate(clim.var = substr(clim.var.month, 9, nchar(clim.var.month)-2)) %>%
  mutate(clim.var.month = substr(clim.var.month, nchar(clim.var.month), nchar(clim.var.month))) %>%
  left_join(., vars_names, by="clim.var") %>%
  mutate(clim.var_lab = factor(clim.var_lab, 
                               levels = c("Maximum temperature (°C)", 
                                          "Minimum temperature (°C)",
                                          "Total precipitations (mm)",
                                          "Net solar radiations (MJ/m²)",
                                          "Evapotransp. ref (mm/day)",
                                          "Vapor pressure deficit"))) %>%
  group_by(Country, clim.var_lab, clim.var.month) %>% 
  summarise(perc05 = quantile(clim.value, probs = c(0.05)),
            median = median(clim.value),
            perc95 = quantile(clim.value, probs = c(0.95))) %>% 
  mutate(Country = factor(Country, levels = c("Global (N=122229)", "United-States of America (N=29803)", "Brazil (N=14757)"))) 

# > Annual climatic data 
clim_data_year <- rbind(
  tab_world %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit")) %>% 
    mutate(Country = "Global (N=122229)"), 
  tab_usa %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit")) %>% 
    mutate(Country = "United-States of America (N=29803)"), 
  tab_bra %>% dplyr::select(
    country_name, x, y, Ya, year, site_year,
    contains("monthly_max_2m"), contains("monthly_min_2m"), contains("monthly_total"), contains("monthly_et0"), contains("monthly_surface_net_solar_radiation"), contains("monthly_vapor_pressure_deficit"))  %>% 
    mutate(Country = "Brazil (N=14757)")) %>% 
  pivot_longer(cols = starts_with("monthly"), values_to = "clim.value", names_to = "clim.var.month") %>% 
  mutate(clim.var = substr(clim.var.month, 9, nchar(clim.var.month)-2)) %>%
  mutate(clim.var.month = substr(clim.var.month, nchar(clim.var.month), nchar(clim.var.month))) %>%
  left_join(., vars_names, by="clim.var") %>%
  mutate(clim.var_lab = factor(clim.var_lab, 
                               levels = c("Maximum temperature (°C)", 
                                          "Minimum temperature (°C)",
                                          "Total precipitations (mm)",
                                          "Net solar radiations (MJ/m²)",
                                          "Evapotransp. ref (mm/day)",
                                          "Vapor pressure deficit"))) %>%
  group_by(Country, clim.var_lab) %>% 
  summarise(perc05 = quantile(clim.value, probs = c(0.05)),
            median = median(clim.value),
            perc95 = quantile(clim.value, probs = c(0.95))) %>% 
  mutate(Country = factor(Country, levels = c("Global (N=122229)", "United-States of America (N=29803)", "Brazil (N=14757)"))) 

# > Plots
for(country in unique(clim_data$Country))
{
  
  # > Color for each analysis set 
  country_color <- case_when(country == "United-States of America (N=29803)" ~ "#E41A1C",
                             country == "Brazil (N=14757)" ~ "#377EB8",
                             country == "Global (N=122229)" ~ "#4DAF4A")
    
  # > DAYS 
  p_days <- clim_data %>% 
    filter(Country==country) %>% 
    left_join(., vars_names, by="clim.var") %>%
    mutate(clim.var_lab = factor(clim.var_lab, 
                                 levels = c("Maximum temperature (°C)", 
                                            "Minimum temperature (°C)",
                                            "Total precipitations (mm)",
                                            "Net solar radiations (MJ/m²)",
                                            "Evapotransp. ref (mm/day)",
                                            "Vapor pressure deficit")))%>% 
    mutate(month_day_1 = if_else(day_of_year %in% c(1, 31, 62, 92, 123, 153, 184), day_of_year, NA)) %>% 
    # Plot 
    ggplot(data = .) + 
    geom_ribbon(aes(x= day_of_year, ymin = perc05, ymax = perc95), fill = country_color,
                alpha = 0.2) + 
    geom_point(aes(x = month_day_1, y = median), color = country_color, 
               size = 0.75, shape = 3) +  
    geom_line(aes(x= day_of_year, y = median), color = country_color) +
    facet_wrap(.~clim.var_lab, scale="free", nrow = 1)  +
    theme_cowplot() + 
    theme(legend.position = "none",
          title = element_text(size=10), 
          axis.text = element_text(size=10),
          strip.text = element_text(size=10), 
          legend.text = element_text(size=10)) + 
    labs(x = "Day of the growing season", y = paste0("Cumulated climatic data\nover soybean growing season")) + 
    ggtitle("a. Cumulative daily climate data")
  
  # > MONTHS
  p_months <- clim_data_month %>% 
    filter(Country==country) %>% 
    ggplot(.) + 
    geom_ribbon(aes(x= clim.var.month, ymin = perc05, ymax = perc95, group = Country),
                alpha = 0.2, fill = country_color) + 
    geom_point(aes(x = clim.var.month, y = median, group = Country),
               color = country_color, size = 0.75) + 
    geom_line(aes(x = clim.var.month, y = median, group = Country),
              color = country_color) + 
    facet_wrap(.~clim.var_lab, scale="free", nrow = 1)  +
    theme_cowplot() + 
    theme(legend.position = "none",
          title = element_text(size=10), 
          axis.text = element_text(size=10),
          strip.text = element_text(size=10), 
          legend.text = element_text(size=10)) + 
    labs(x = "Month of the growing season", y = "Climatic data\nover soybean growing season") +
    ggtitle("b. Monthly climate data"); p_months
  
  # > YEARS
  p_year <- clim_data_year %>% 
    filter(Country==country) %>% 
    ggplot(.) + 
    geom_linerange(aes(x = 0, ymin = perc05, ymax = perc95),
                   color = country_color, position = position_dodge2(width = 0.75)) + 
    geom_point(aes(x = 0, y = median, group = Country), 
               shape = 15, color = country_color, position = position_dodge2(width = 0.75)) + 
    facet_wrap(.~clim.var_lab, scale="free", nrow = 1)  +
    theme_cowplot() + 
    theme(legend.position = "bottom",
          title = element_text(size=10), 
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=10),
          strip.text = element_text(size=10), 
          legend.text = element_text(size=10)) + 
    labs(x = "Annual median over the growing season", y = "Climatic data\nover soybean growing season") +
    ggtitle("c. Annual climate data"); p_year
  
  # Save only plot for daily/monthy cumulative variables 
  #ggsave(p_days,   filename=paste0(save_path, "Figure_3_1.png"), width = 13, height=4, dpi=300, bg="white")
  #ggsave(p_months, filename=paste0(save_path, "Figure_3_2.png"), width = 13, height=4, dpi=300, bg="white")
  
  p_clim <- plot_grid(p_days, p_months, p_year, nrow = 3)
  ggsave(p_clim,   
         filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/03_CLIMATE/Supplementary_Figure_clim_", country, ".png"), 
         width = 15, height=9.5, dpi=300, bg="white")

  
}

# > Values corresponding to dots in the plot: 
clim_data_month <- list("United-States of America (N=29803)"= list_data_day_usa,
                        "Brazil (N=14757)"                   = list_data_day_bra, 
                        "Global (N=122229)"                  = list_data_day_world) %>% 
  map_dfr(., ~{
    
    .x %>%
      map_dfr(., ~{
        
        # > Variable 
        var_i <- unique(.x$clim.var)
        var_lab_i <- vars_names[which(vars_names$clim.var==var_i), "clim.var_lab2"]
        
        .x %>% 
          # 3 months after seeding
          #filter(day_of_year == 92) %>% 
          # 6 months after seeding
          filter(day_of_year == 184) %>% 
          # Compute 5th, 95th percentiles and median of climatic feature 
          group_by(clim.var) %>% 
          summarise(perc05 = quantile(cum_clim.value, probs = c(0.05)),
                    median = median(cum_clim.value),
                    perc95 = quantile(cum_clim.value, probs = c(0.95)))
        
      })
    
  }, .id = "Country")

clim_data_month %>% 
  mutate(lab = paste0(round(median, 2), " (", round(perc05, 2), ", ", round(perc95, 2), ")")) %>% 
  dplyr::select(-median, -perc05, -perc95) %>%
  pivot_wider(names_from = "Country", values_from = lab)



# -----------------------------------------------------------
# % of total variance explained in the data by each strategy 
# of dimension reduction techniques in USA and BRAZIL datasets
# (Plot for world in the main)

# > Function to remove useless pannels
remove_facets <- function(plot, layout) {
  layout <- strsplit(layout, split = '\n')[[1]]
  layout <- lapply(layout, trimws)
  layout <- matrix(unlist(sapply(layout, strsplit, "")),
                   nrow = length(layout), byrow = T)
  layout <- which(layout == "#", arr.ind = TRUE)
  prm <- apply(layout,1,\(x) {
    c(glue::glue("panel-{x[1]}-{x[2]}"),
      glue::glue("strip-t-{x[2]}-{x[1]}"))
  })
  # https://stackoverflow.com/a/30372692/1296582
  g <- ggplot2::ggplotGrob(plot)
  rm_grobs <- g$layout$name %in% prm
  g$grobs[rm_grobs] <- NULL
  g$layout <- g$layout[!rm_grobs, ]
  ggpubr::as_ggplot(g)
}

remove <- c("aa#aa
             aa#aa
             aa#aa
             aa#aa
             aa#aa
             aa#aa
             ##a##")

Supp_Fig_var_USA <- tab_var_full %>% 
  filter(country=="USA") %>% 
  # > keep for 1 climatic variable 
  #filter(clim.var == "max_temp" | is.na(clim.var)==T) %>%
  # > select only the 10 first cp 
  filter(PC < 11) %>% 
  # > or based on the cumulative explained variance
  dplyr::select(model, data_type, clim.var_lab, PC, variance.percent, cumulative.variance.percent) %>% 
  # > levels 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               "Total\nprecipitations (mm)"="Average\nprecipitation (mm)",
                               "Net solar\nradiations (MJ/m²)"="Net solar\nradiation (MJ/m²)"),
         clim.var_lab = factor(clim.var_lab, levels = c('Maximum\ntemperature (°C)', 
                                                        'Minimum\ntemperature (°C)', 
                                                        'Average\nprecipitation (mm)', 
                                                        'Net solar\nradiation (MJ/m²)', 
                                                        'Evapotransp.\nref (mm/day)', 
                                                        'Vapor pressure\ndeficit', 
                                                        'All variables\n(only for MFPCA)'))) %>% 
  ggplot(.) + 
  geom_hline(yintercept = 90, color = "darkgrey", lty=2) + 
  geom_line(aes(y=cumulative.variance.percent, 
                x=PC, color = data_type)) + 
  geom_point(aes(y=cumulative.variance.percent, shape = data_type,
                 x=PC, color = data_type), size=1) + 
  #geom_line(aes(y=variance.percent)) + 
  #geom_point(aes(y=variance.percent, color=data_type), shape = 21) + 
  facet_grid(clim.var_lab ~ model, switch = "y") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        strip.background = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_text(vjust=-50),
        plot.caption.position = "plot", 
        plot.caption = element_text(hjust = 0)) + 
  scale_color_manual(values = c("red", "black"), name = "Temporal resolution") +
  scale_shape_manual(values = c(16,8), name = "Temporal resolution") + 
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,105)) + 
  guides(color = guide_legend(order = 1, ncol = 2),
         shape = guide_legend(order = 1, override.aes = list(size=2), ncol = 2)) + 
  labs(x = "Scores", y = "Cumulated explained variance (%)", 
       caption = "Abbreviations: PCA: principal component analysis; FPCA: functional PCA; MFPCA: multivariate FPCA; PLSR: partial least square regression; FPLSR: functional PLSR.") ; Supp_Fig_var_USA

# > Save 
ggsave(plot = remove_facets(Supp_Fig_var_USA, remove), 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/04_EXPLAINED_VAR/Supp_Fig_var_USA.png"), 
       width = 10, height=11, dpi=300, bg="white")

# > BRA
Supp_Fig_var_BRA <- tab_var_full %>% 
  filter(country=="BRA") %>% 
  # > keep for 1 climatic variable 
  #filter(clim.var == "max_temp" | is.na(clim.var)==T) %>%
  # > select only the 10 first cp 
  filter(PC < 11) %>% 
  # > or based on the cumulative explained variance
  dplyr::select(model, data_type, clim.var_lab, PC, variance.percent, cumulative.variance.percent) %>% 
  # > levels 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               "Total\nprecipitations (mm)"="Average\nprecipitation (mm)",
                               "Net solar\nradiations (MJ/m²)"="Net solar\nradiation (MJ/m²)"),
         clim.var_lab = factor(clim.var_lab, levels = c('Maximum\ntemperature (°C)', 
                                                        'Minimum\ntemperature (°C)', 
                                                        'Average\nprecipitation (mm)', 
                                                        'Net solar\nradiation (MJ/m²)', 
                                                        'Evapotransp.\nref (mm/day)', 
                                                        'Vapor pressure\ndeficit', 
                                                        'All variables\n(only for MFPCA)'))) %>% 
  ggplot(.) + 
  geom_hline(yintercept = 90, color = "darkgrey", lty=2) + 
  geom_line(aes(y=cumulative.variance.percent, 
                x=PC, color = data_type)) + 
  geom_point(aes(y=cumulative.variance.percent, shape = data_type,
                 x=PC, color = data_type), size=1) + 
  #geom_line(aes(y=variance.percent)) + 
  #geom_point(aes(y=variance.percent, color=data_type), shape = 21) + 
  facet_grid(clim.var_lab ~ model, switch = "y") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        strip.background = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        strip.switch.pad.grid = unit(1, "cm"),
        axis.title.y = element_text(vjust=-50),
        plot.caption.position = "plot", 
        plot.caption = element_text(hjust = 0)) + 
  scale_color_manual(values = c("red", "black"), name = "Temporal resolution") +
  scale_shape_manual(values = c(16,8), name = "Temporal resolution") + 
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0,105)) + 
  guides(color = guide_legend(order = 1, ncol = 2),
         shape = guide_legend(order = 1, override.aes = list(size=2), ncol = 2)) + 
  labs(x = "Scores", y = "Cumulated explained variance (%)", 
       caption = "Abbreviations: PCA: principal component analysis; FPCA: functional PCA; MFPCA: multivariate FPCA; PLSR: partial least square regression; FPLSR: functional PLSR.") ; Supp_Fig_var_BRA

# > Save 
ggsave(plot = remove_facets(Supp_Fig_var_BRA, remove), 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/04_EXPLAINED_VAR/Supp_Fig_var_BRA.png"), 
       width = 10, height=11, dpi=300, bg="white")

# --------------------------------------------------------------------------
# Scores plots 

# > Generate a table with the different scores to keep for each analysis
country <- c("01_USA", "02_BRA", "03_WORLD")
dim_red   <- c("PC", "FPC", "PLS", "FPLS")
data_type <- c("month", "day")
PC <- 1:2
clim.var <- c(unique(vars_names$clim.var_abb)[1:5], "vapor_pressure_deficit")
tab_scores_2 <- expand.grid(country=country, 
                            dim_red=dim_red, 
                            data_type=data_type, 
                            PC=PC, 
                            clim.var=clim.var) %>% 
  mutate(dim_red_full = recode(dim_red, "PC"="PCA", "FPC"="FPCA", "MFPC"="MFPCA")) %>% 
  unite("dim_red_PC", c(dim_red, PC), remove = F, sep="") %>%
  unite("score_to_keep", c(dim_red_PC, data_type, clim.var), sep="_", remove = F) %>% 
  rbind(., expand.grid(country=country,
                       score_to_keep = c("MFPC1_day", "MFPC2_day", "MFPC1_month", "MFPC2_month"),
                      clim.var = clim.var,
                      dim_red_full = "MFPCA") %>% 
          separate(col = score_to_keep, into = c("dim_red_PC", "data_type"), sep = "_", remove = F) %>% 
          mutate(dim_red = substr(dim_red_PC, 1, nchar(dim_red_PC)-1),
                 PC = substr(dim_red_PC, nchar(dim_red_PC), nchar(dim_red_PC))) %>% 
          dplyr::select(country, score_to_keep, dim_red_PC, dim_red, data_type, PC, clim.var, dim_red_full)) 

# > WORLD
tab_full_i <- tab_world

tab_scores_2 %>% 
  filter(country == "03_WORLD") %>% 
  #filter(clim.var == "vapor_pressure_deficit") %>% 
  mutate(clim.var = droplevels(.$clim.var),
         dim_red = droplevels(.$dim_red)) %>%
  # > climatic variable 
  split(.$clim.var) %>% 
  map(., ~{
    
    country_i <- unique(.x$country)
    clim.var_i  <- unique(.x$clim.var)
    clim_var2_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor_pressure_deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var"])
    clim_var3_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor pressure deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var_lab2"])
    
    # > data_type (month or day)
    list_plot_clim.var <- .x %>% 
      split(.$data_type) %>% 
      map(., ~{ 
        
        data_type_i <- unique(.x$data_type)
        
        # > dim_red (PC, FPC, MFPC, PLS, FPLS)
        list_plot_clim_data_type <- .x %>% 
          split(.$dim_red) %>% 
          map(., ~{
            
            dim_red_i <- unique(.x$dim_red)
            dim_red_full_i <- unique(.x$dim_red_full)
            
            to_keep   <- c(unique(.x %>%
                                    filter(score_to_keep %in% names(tab_full_i)) %>% 
                                    pull(score_to_keep)), 
                           paste0("year_", clim_var2_i))
            
            # > Select the correct table 
            tab_i <- tab_full_i %>% 
              dplyr::select(all_of(to_keep)) 
            
            # > If at least 2 scores
            if(ncol(tab_i) > 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 3) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_hline(yintercept = mean(tab_i[,2])) +
                geom_point(aes(x = tab_i[,1], y = tab_i[,2], color = order),
                           size=0.5) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text = element_text(size=8),
                      legend.position = "bottom") + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="Score 2") +
                ggtitle(paste0(dim_red_full_i))
            }
            
            # > If only 1 score
            if(ncol(tab_i) == 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 2) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_jitter(aes(x = tab_i[,1], y = as.factor(0), color = order),
                            size=0.5, width = 0.75) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text.x = element_text(size=8),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "bottom",
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank()) + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="") +
                ggtitle(paste0(dim_red_full_i))
              
            } 
            
            # > Plot 
            plot_i
            
          })
        
        # > Legends and Plot 
        leg <- get_legend(list_plot_clim_data_type$PC)
        plot <- plot_grid(list_plot_clim_data_type$PC + theme(legend.position = "none"),
                          list_plot_clim_data_type$FPC + theme(legend.position = "none"),
                          list_plot_clim_data_type$MFPC + theme(legend.position = "none"), 
                          list_plot_clim_data_type$PLS + theme(legend.position = "none"), 
                          list_plot_clim_data_type$FPLS + theme(legend.position = "none"), 
                          nrow = 1)# + theme(plot.background = element_rect(fill = "white", color = "transparent")),
        list("legend"=leg, "plot"=plot)
        
      })
    
    # Plot for the clim.var
    p_save <- plot_grid(
      plot_grid(ggplot() + theme_void(),
                list_plot_clim.var$month$plot, 
                ggplot() + theme_void(),
                list_plot_clim.var$day$plot, 
                ncol = 1,
                rel_heights = c(0.1, 1, 0.1, 1),
                labels = c("", "a. Monthly averages", "", "b. Daily average"), 
                label_size = 11, label_x = -0.01, label_y = 1.1),
      list_plot_clim.var$month$legend,
      ggplot() + 
        theme_void() + 
        theme(plot.caption.position = "plot", 
              plot.caption = element_text(hjust = 0, size=8)) +
        labs(caption = "Abbreviations: PCA: Principal Component Analysis; FPCA: Functional PCA; MFPCA: Multivariate FPCA; PLS: Partial Least-Square regression; FPLS: Functional PLS"),
      ncol = 1, rel_heights = c(1, 0.1, 0.05)
    )
    # > Save
    ggsave(p_save + theme(plot.background = element_rect(fill = "white", color = "transparent")), 
           filename = paste0(save_path, "/SUPPLEMENTARY_FIGURES/05_SCORES_PLOTS/", country_i,"_", clim.var_i, ".png"),
           width=12.5, height=7, dpi=300)
})

# > USA
tab_full_i <- tab_usa 

tab_scores_2 %>% 
  filter(country == "01_USA") %>% 
  #filter(clim.var == "vapor_pressure_deficit") %>% 
  mutate(clim.var = droplevels(.$clim.var),
         dim_red = droplevels(.$dim_red)) %>%
  # > climatic variable 
  split(.$clim.var) %>% 
  map(., ~{
    
    country_i <- unique(.x$country)
    clim.var_i  <- unique(.x$clim.var)
    clim_var2_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor_pressure_deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var"])
    clim_var3_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor pressure deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var_lab2"])
    
    # > data_type (month or day)
    list_plot_clim.var <- .x %>% 
      split(.$data_type) %>% 
      map(., ~{ 
        
        data_type_i <- unique(.x$data_type)
        
        # > dim_red (PC, FPC, MFPC, PLS, FPLS)
        list_plot_clim_data_type <- .x %>% 
          split(.$dim_red) %>% 
          map(., ~{
            
            dim_red_i <- unique(.x$dim_red)
            dim_red_full_i <- unique(.x$dim_red_full)
            
            to_keep   <- c(unique(.x %>%
                                    filter(score_to_keep %in% names(tab_full_i)) %>% 
                                    pull(score_to_keep)), 
                           paste0("year_", clim_var2_i))
            
            # > Select the correct table 
            tab_i <- tab_full_i %>% 
              dplyr::select(all_of(to_keep)) 
            
            # > If at least 2 scores
            if(ncol(tab_i) > 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 3) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_hline(yintercept = mean(tab_i[,2])) +
                geom_point(aes(x = tab_i[,1], y = tab_i[,2], color = order),
                           size=0.5) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text = element_text(size=8),
                      legend.position = "bottom") + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="Score 2") +
                ggtitle(paste0(dim_red_full_i))
            }
            
            # > If only 1 score
            if(ncol(tab_i) == 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 2) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_jitter(aes(x = tab_i[,1], y = as.factor(0), color = order),
                            size=0.5, width = 0.75) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text.x = element_text(size=8),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "bottom",
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank()) + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="") +
                ggtitle(paste0(dim_red_full_i))
              
            } 
            
            # > Plot 
            plot_i
            
          })
        
        # > Legends and Plot 
        leg <- get_legend(list_plot_clim_data_type$PC)
        plot <- plot_grid(list_plot_clim_data_type$PC   + theme(legend.position = "none"),
                          list_plot_clim_data_type$FPC  + theme(legend.position = "none"),
                          list_plot_clim_data_type$MFPC + theme(legend.position = "none"), 
                          list_plot_clim_data_type$PLS  + theme(legend.position = "none"), 
                          list_plot_clim_data_type$FPLS + theme(legend.position = "none"), 
                          nrow = 1)# + theme(plot.background = element_rect(fill = "white", color = "transparent")),
        list("legend"=leg, "plot"=plot)
        
      })
    
    # Plot for the clim.var
    p_save <- plot_grid(
      plot_grid(ggplot() + theme_void(),
                list_plot_clim.var$month$plot, 
                ggplot() + theme_void(),
                list_plot_clim.var$day$plot, 
                ncol = 1,
                rel_heights = c(0.1, 1, 0.1, 1),
                labels = c("", "a. Monthly averages", "", "b. Daily average"), 
                label_size = 11, label_x = -0.01, label_y = 1.1),
      list_plot_clim.var$month$legend,
      ggplot() + 
        theme_void() + 
        theme(plot.caption.position = "plot", 
              plot.caption = element_text(hjust = 0, size=8)) +
        labs(caption = "Abbreviations: PCA: Principal Component Analysis; FPCA: Functional PCA; MFPCA: Multivariate FPCA; PLS: Partial Least-Square regression; FPLS: Functional PLS"),
      ncol = 1, rel_heights = c(1, 0.1, 0.05)
    ) 
    # > Save
    ggsave(p_save + theme(plot.background = element_rect(fill = "white", color = "transparent")), 
           filename = paste0(save_path, "/SUPPLEMENTARY_FIGURES/05_SCORES_PLOTS/", country_i,"_",clim.var_i, ".png"),
           width=12.5, height=7, dpi=300)
  })

# > BRA
tab_full_i <- tab_bra 

tab_scores_2 %>% 
  filter(country == "02_BRA") %>% 
  #filter(clim.var == "vapor_pressure_deficit") %>% 
  mutate(clim.var = droplevels(.$clim.var),
         dim_red = droplevels(.$dim_red)) %>%
  # > climatic variable 
  split(.$clim.var) %>% 
  map(., ~{
    
    country_i <- unique(.x$country)
    clim.var_i  <- unique(.x$clim.var)
    clim_var2_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor_pressure_deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var"])
    clim_var3_i <- ifelse(clim.var_i == "vapor_pressure_deficit", "vapor pressure deficit", vars_names[which(vars_names$clim.var_abb == clim.var_i),"clim.var_lab2"])
    
    # > data_type (month or day)
    list_plot_clim.var <- .x %>% 
      split(.$data_type) %>% 
      map(., ~{ 
        
        data_type_i <- unique(.x$data_type)
        
        # > dim_red (PC, FPC, MFPC, PLS, FPLS)
        list_plot_clim_data_type <- .x %>% 
          split(.$dim_red) %>% 
          map(., ~{
            
            dim_red_i <- unique(.x$dim_red)
            dim_red_full_i <- unique(.x$dim_red_full)
            
            to_keep   <- c(unique(.x %>%
                                    filter(score_to_keep %in% names(tab_full_i)) %>% 
                                    pull(score_to_keep)), 
                           paste0("year_", clim_var2_i))
            
            # > Select the correct table 
            tab_i <- tab_full_i %>% 
              dplyr::select(all_of(to_keep)) 
            
            # > If at least 2 scores
            if(ncol(tab_i) > 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 3) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_hline(yintercept = mean(tab_i[,2])) +
                geom_point(aes(x = tab_i[,1], y = tab_i[,2], color = order),
                           size=0.5) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text = element_text(size=8),
                      legend.position = "bottom") + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="Score 2") +
                ggtitle(paste0(dim_red_full_i))
            }
            
            # > If only 1 score
            if(ncol(tab_i) == 2)
            {
              
              tab_i <- tab_i %>% 
                dplyr::rename("year_var" = 2) %>%
                # order based on year average
                ungroup() %>% 
                arrange(year_var) %>% 
                mutate(order = 1:nrow(.))
              
              # > Plot 
              plot_i <- ggplot(data = tab_i) +
                geom_vline(xintercept = mean(tab_i[,1])) +
                geom_jitter(aes(x = tab_i[,1], y = as.factor(0), color = order),
                            size=0.5, width = 0.75) + 
                theme_bw() + 
                theme(title = element_text(size = 10),
                      axis.title = element_text(size=10),
                      axis.text.x = element_text(size=8),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "bottom",
                      panel.grid.major.y = element_blank(),
                      panel.grid.minor.y = element_blank()) + 
                scale_color_gradientn(colours = pal, name=paste0("Rank based on annual average of ", clim_var3_i, " in each site-year (lower value: lower rank)"),
                                      guide = guide_colorbar(barwidth = 15, barheight = 0.5, title.vjust = 1)) + 
                labs(x="Score 1", y="") +
                ggtitle(paste0(dim_red_full_i))
              
            } 
            
            # > Plot 
            plot_i
            
          })
        
        # > Legends and Plot 
        leg <- get_legend(list_plot_clim_data_type$PC)
        plot <- plot_grid(list_plot_clim_data_type$PC   + theme(legend.position = "none"),
                          list_plot_clim_data_type$FPC  + theme(legend.position = "none"),
                          list_plot_clim_data_type$MFPC + theme(legend.position = "none"), 
                          list_plot_clim_data_type$PLS  + theme(legend.position = "none"), 
                          list_plot_clim_data_type$FPLS + theme(legend.position = "none"), 
                          nrow = 1)# + theme(plot.background = element_rect(fill = "white", color = "transparent")),
        list("legend"=leg, "plot"=plot)
        
      })
    
    # Plot for the clim.var
    p_save <- plot_grid(
      plot_grid(ggplot() + theme_void(),
                list_plot_clim.var$month$plot, 
                ggplot() + theme_void(),
                list_plot_clim.var$day$plot, 
                ncol = 1,
                rel_heights = c(0.1, 1, 0.1, 1),
                labels = c("", "a. Monthly averages", "", "b. Daily average"), 
                label_size = 11, label_x = -0.01, label_y = 1.1),
      list_plot_clim.var$month$legend,
      ggplot() + 
        theme_void() + 
        theme(plot.caption.position = "plot", 
              plot.caption = element_text(hjust = 0, size=8)) +
        labs(caption = "Abbreviations: PCA: Principal Component Analysis; FPCA: Functional PCA; MFPCA: Multivariate FPCA; PLS: Partial Least-Square regression; FPLS: Functional PLS"),
      ncol = 1, rel_heights = c(1, 0.1, 0.05)
    ) 
    # > Save
    ggsave(p_save + theme(plot.background = element_rect(fill = "white", color = "transparent")) + labs(caption = "Abbreviations: PCA: Principal Component Analysis; FPCA: Functional PCA; MFPCA: Multivariate FPCA; PLS: Partial Least-Square regression; FPLS: Functional PLS"), 
           filename = paste0(save_path, "/SUPPLEMENTARY_FIGURES/05_SCORES_PLOTS/", country_i,"_",clim.var_i, ".png"),
           width=12.5, height=7, dpi=300)
  })

# --------------------------------------------------------------------------
# Correlation between scores (PCA, FPCA, MFPCA, PLS, FPLS) and monthly cumulated variables 

# > Select data for plot correlations
# - data: tab_usa, tab_bra, or tab_world
# - var_i: one of unique(vars_name$clim.var)
# - dim_red_i: one of: PC, FPC, MFPC, PLS, FPLS
# -data_type: one of: "month" or "day"
tab_for_corrplot <- function(data, var_i, dim_red_i, data_type_i){
  
  # > Remove the characters before "." --> for FPLS 
  names(data) <- gsub("^.*\\.","", names(data))
  
  # > Variable labs and abbreviations
  abb <- vars_names[which(vars_names$clim.var==var_i),]$clim.var_abb
  abb <- ifelse(abb == "vpd", "vapor_pressure_deficit", abb)
  lab <- vars_names[which(vars_names$clim.var==var_i),]$clim.var_lab
  
  # > Select data for correlation
  if(dim_red_i != "MFPC")
  {
    tab_for_cor <- data %>% 
      dplyr::select(starts_with(paste0("monthly_", var_i)),
                    ends_with(paste0(data_type_i, "_", abb)),
                    irrigated_portion) %>% 
      dplyr::select(starts_with(paste0(dim_red_i)), 
                    starts_with(paste0("monthly_", var_i)),
                    irrigated_portion)
  }
  if(dim_red_i == "MFPC")
  {
    tab_for_cor <- data %>% 
      dplyr::select(ends_with(paste0(data_type_i)),
                    starts_with(paste0("monthly_", var_i)),
                    irrigated_portion)
  }
  
  # > Rename to fit in the plot
  colnames(tab_for_cor) <- c(paste0(dim_red_i, " score ", 1:(ncol(tab_for_cor)-8)),
                             "Month 1", "Month 2", "Month 3", "Month 4", "Month 5", "Month 6", "Month 7", 
                             "Irrigation")
  
  # > Title
  tab_for_cor$lab <- lab
  
  # > Plot
  return(tab_for_cor)
  
}


# GRAPHS
for(dim_red in c("PC", "FPC", "MFPC", "PLS", "FPLS"))
{
  
  for(data_type in c("month", "day"))
  {
    
    # > USA
    png(filename = paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/FIGURES/SUPPLEMENTARY_FIGURES/06_CORRPLOTS/USA_", data_type, "_", dim_red, ".png"),
        width = 30, height = 20, units = "cm", res = 300)
    
    # > Graphics
    par(mfrow=c(2,3))
    
    for(var_i in unique(vars_names$clim.var))
    {
      
      # > select data 
      tab_i <- tab_for_corrplot(data = tab_usa, 
                                var_i = var_i,
                                dim_red_i = dim_red, 
                                data_type_i = data_type)
      # > plot
      if(ncol(tab_i) > 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             #addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      if(ncol(tab_i) <= 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      
      
    }
    dev.off()
    
    # > BRAZIL
    png(filename = paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/FIGURES/SUPPLEMENTARY_FIGURES/06_CORRPLOTS/BRA_", data_type, "_", dim_red, ".png"),
        width = 30, height = 20, units = "cm", res = 300)
    
    # > Graphics
    par(mfrow=c(2,3))
    
    for(var_i in unique(vars_names$clim.var))
    {
      
      
      tab_i <- tab_for_corrplot(data = tab_bra, 
                                var_i = var_i, 
                                dim_red_i = dim_red, 
                                data_type_i = data_type)
      
      # > plot
      if(ncol(tab_i) > 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             #addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      if(ncol(tab_i) <= 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      
    }
    dev.off()
    
    # > WORLD
    png(filename = paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/FIGURES/SUPPLEMENTARY_FIGURES/06_CORRPLOTS/WORLD_", data_type, "_", dim_red, ".png"),
        width = 30, height = 20, units = "cm", res = 300)
    
    # > Graphics
    par(mfrow=c(2,3))
    
    for(var_i in unique(vars_names$clim.var))
    {
      
      # > Select data 
      tab_i <- tab_for_corrplot(data = tab_world, 
                                var_i = var_i, 
                                dim_red_i = dim_red, 
                                data_type_i = data_type)
      
      # > plot
      if(ncol(tab_i) > 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             #addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      if(ncol(tab_i) <= 13)
      {
        
        cor(tab_i %>% dplyr::select(-lab)) %>% 
          corrplot::corrplot(., method = "square", diag = F, 
                             type = "lower", mar = c(1, 1, 1, 1), 
                             addCoef.col = 'black', 
                             tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
                             title = paste0(unique(tab_i$lab)))
        
      }
      
    }
    dev.off()
  }
  
}


# --------------------------------------------------------------------------
# Performance of all models to predict soybean yield or soybean yield anomaly in USA, Brazil, globally 

# > Compute mean performance 
tab_perf_labelled_mean <- tab_perf_labelled %>% 
  group_by(Outcome, Outcome_lab, Country, Country_lab, Model, dimension_reduction, pred_perf, pred_perf_lab, gpe_model) %>% 
  # MEAN OF BOTH CROSS-VALIDATION
  summarise(mean_pred_perf = mean(pred_perf_value)) %>% 
  ungroup() %>%
  mutate(Type_cv_lab = "c. Mean of cross-validation results") %>%
  # > negative values
  mutate(
    negative_pred_perf_value = if_else(mean_pred_perf < 0, 1, 0),
    mean_pred_perf_lab = if_else(mean_pred_perf < 0, 0, mean_pred_perf))

# > Plots for each outcome (yield or yield anomaly) and each country
for(outcome in unique(tab_perf_labelled$Outcome_lab))
{
  for(country in unique(tab_perf_labelled$Country_lab))
  {
    
    # Select the correct data 
    tab_perf_labelled_i <- tab_perf_labelled %>% 
      filter(Country_lab==country) %>% 
      filter(Outcome_lab == outcome) %>% 
      mutate(Type_cv_lab = recode(Type_cv_lab, "Cross-validation on years"="a. Cross-validation on years", "Cross-validation on sites"="b. Cross-validation on sites"),
             Type_cv_lab = factor(Type_cv_lab, levels = c("a. Cross-validation on years", "b. Cross-validation on sites"))) %>%
      # > group of models 
      mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
      # > negative values
      mutate(negative_pred_perf_value = ifelse(negative_pred_perf_value == 0, 1, 0)) 
    
    tab_perf_labelled_mean_i <- tab_perf_labelled_mean %>% 
      filter(Country_lab==country) %>% 
      filter(Outcome_lab == outcome) %>%
      # > group of models 
      mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
      # > negative values
      mutate(negative_pred_perf_value = ifelse(negative_pred_perf_value == 0, 1, 0)) 
    
    # PLOT N of predictors per models
    p.n_pred <- tab_perf_labelled_i %>% 
      filter(pred_perf == "NSE") %>% 
      # > Nb of predictors per model
      group_by(Outcome_lab, Country_lab, dimension_reduction, Model) %>%
      mutate(N_pred = max(N_pred, na.rm = T))  %>%
      # > plot
      ggplot(., aes(x = pred_perf_value, y=Model)) +
      geom_text(aes(label = N_pred, x = 1.4), 
                size = 3, check_overlap = T) + 
      facet_grid(dimension_reduction ~ ., scales = "free_y", space = "free", switch = "y") + 
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_blank(), 
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank()) ; p.n_pred
    
    # PLOT NSE 
    # > individual cross-validation results
    p.nse <- tab_perf_labelled_i %>% 
      filter(pred_perf == "NSE") %>% 
      # > order of models 
      arrange(Type_cv_lab, Outcome_lab, Country_lab, desc(pred_perf_value)) %>% 
      group_by(Type_cv_lab, Outcome_lab, Country_lab) %>% 
      mutate(order = row_number()) %>%
      mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
             best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
      # > min and max among lm and rf
      group_by(Type_cv_lab, Outcome_lab, Country, Model) %>% 
      mutate(min_x_lab = min(pred_perf_value_lab),
             max_x_lab = max(pred_perf_value_lab)) %>% 
      # > x coordinates for labels 
      mutate(x_lab = case_when(
        pred_perf_value_lab==min_x_lab ~ pred_perf_value_lab-0.075*pred_perf_value_lab,
        pred_perf_value_lab==max_x_lab ~ pred_perf_value_lab+0.075*pred_perf_value_lab,
        TRUE ~ NA)) %>% 
      mutate(x_lab = if_else(x_lab <0, -0.03, x_lab),
             x_lab = if_else(x_lab >1, 1, x_lab)) %>% 
      # > plot
      ggplot(., aes(x = pred_perf_value_lab, y = Model)) +
      geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                     color = "grey8", size=0.5) +
      geom_point(aes(color = gpe_model, shape = best, size = negative_pred_perf_value)) + 
      geom_point(aes(color = gpe_model, size = abs(negative_pred_perf_value-1)), shape="|") + 
      #geom_point(aes(color = gpe_model, shape = best)) + 
      #geom_point(aes(color = gpe_model, fill = as.factor(negative_pred_perf_value)), size=2, shape = 21) + 
      geom_text(aes(x = x_lab, label = order), 
                nudge_y = 0.4, 
                size = 2) + 
      facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", space = "free", switch = "y") + 
      theme_bw() +
      theme(legend.position = 'bottom', 
            axis.title.y = element_blank(),
            strip.placement = "outside", 
            strip.text.y.left = element_text(angle = 0, size=9),
            strip.text.x = element_text(size=10), 
            strip.background = element_rect(fill="grey75", color = "black"),
            title = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            legend.text = element_text(size=10),
            panel.grid = element_blank()) + 
      scale_shape_manual(values = c(8, 20), name = "Rank model") + 
      scale_color_viridis_d(option = "A", begin = 0.2, end = 0.75, name = "Model family", direction = -1) + 
      scale_size_continuous(range=c(0,3), guide = guide_none()) +
      scale_fill_manual(values = c("white", "transparent")) + 
      scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2), limits = c(0,1)) + 
      guides(color = guide_legend(order = 1, override.aes = list(size=3, shape=20), ncol = 2),
             shape = guide_none())   +  
      labs(x = "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)") ; p.nse
    
    # > Results based on mean of results of each CV procedure
    p.nse_mean <- tab_perf_labelled_mean_i %>% 
      filter(pred_perf == "NSE") %>% 
      # > order of models 
      arrange(Outcome_lab, Country_lab, desc(mean_pred_perf)) %>% 
      group_by(Outcome_lab, Country_lab) %>% 
      mutate(order = row_number()) %>%
      mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
             best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
      # > min and max among lm and rf
      group_by(Outcome_lab, Country_lab, Model) %>% 
      mutate(min_x_lab = min(mean_pred_perf_lab),
             max_x_lab = max(mean_pred_perf_lab)) %>% 
      # > x coordinates for labels 
      mutate(x_lab = case_when(
        mean_pred_perf_lab==min_x_lab ~ mean_pred_perf_lab-0.075*mean_pred_perf_lab,
        mean_pred_perf_lab==max_x_lab ~ mean_pred_perf_lab+0.075*mean_pred_perf_lab,
        TRUE ~ NA)) %>% 
      mutate(x_lab = if_else(x_lab <0, -0.03, x_lab),
             x_lab = if_else(x_lab >1, 1, x_lab)) %>%
      # > plot
      ggplot(., aes(x = mean_pred_perf_lab, y=Model)) +
      geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                     color = "grey8", size=0.5) +
      geom_point(aes(color = gpe_model, shape = best, size = negative_pred_perf_value)) + 
      geom_point(aes(color = gpe_model, size = abs(negative_pred_perf_value-1)), shape="|") + 
      #geom_point(aes(color = gpe_model, shape = best)) + 
      #geom_point(aes(color = gpe_model, fill = as.factor(negative_pred_perf_value)), size=2, shape = 21) + 
      geom_text(aes(x = x_lab, label = order), 
                nudge_y = 0.4, 
                size = 2) +
      facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", space = "free", switch = "y") + 
      theme_bw() +
      theme(legend.position = 'bottom', 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            strip.placement = "outside", 
            strip.text.y.left = element_blank(),
            strip.text.x = element_text(size=10), 
            strip.background.x = element_rect(fill="grey75", color = "black"),
            strip.background.y = element_blank(),
            title = element_text(size=10),
            legend.text = element_text(size=10),
            #panel.grid.major.x = element_line(color = "grey", linewidth = 0.01),
            #panel.grid.minor.x = element_line(color = "grey", linewidth = 0.01),
            panel.grid = element_blank()) + 
      scale_shape_manual(values = c(8, 20), name = "Rank model") + 
      scale_color_viridis_d(option = "A", begin = 0.2, end = 0.75, name = "Model family", direction = -1) + 
      scale_fill_manual(values = c("white", "transparent")) + 
      scale_size_continuous(range=c(0,3), guide = guide_none()) +
      scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2), limits = c(0,1)) +
      guides(color = guide_legend(order = 1, override.aes = list(size=3, shape=20), ncol = 2),
             shape = guide_none()) +
      labs(x = "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)") ; p.nse_mean
    
    # > Get legend of the 1rst plot
    leg <- get_legend(p.nse_mean)
    
    # > Merge plots and add common x label
    plot <- plot_grid(p.nse + theme(legend.position = "none", axis.title.x = element_blank()), 
                      p.nse_mean + theme(legend.position = "none", axis.title.x = element_blank()), 
                      p.n_pred,
                      nrow=1, rel_widths = c(2.75, 1, 0.15), axis = "tblr", align = "h")
    
    plot_x <- ggdraw(add_sub(plot, "Nash-Sutcliffe model efficiency", 
                             vpadding=grid::unit(0,"lines"),y=0, x=0.6, vjust=0, size = 10))
    
    # > Add common legend and save 
    FigNSE <- plot_grid(plot_x, 
                        NULL,
                        leg, 
                        NULL,
                        ncol = 1, rel_heights = c(0.9,0.025,0.05,0.025))
    
    # PLOT RMSEP
    # > individual cross-validation results
    p.rmsep <- tab_perf_labelled_i %>%
      filter(pred_perf == "RMSEP")  %>% 
      # > order of models 
      arrange(Type_cv_lab, Outcome_lab, Country_lab, pred_perf_value) %>% 
      group_by(Type_cv_lab, Outcome_lab, Country_lab) %>% 
      mutate(order = row_number()) %>%
      mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
             best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
      # > min and max among lm and rf
      group_by(Type_cv_lab, Outcome_lab, Country, Model) %>% 
      mutate(min_x_lab = min(pred_perf_value),
             max_x_lab = max(pred_perf_value)) %>% 
      # > x coordinates for labels 
      mutate(x_lab = case_when(
        pred_perf_value==min_x_lab ~ pred_perf_value-0.025*pred_perf_value,
        pred_perf_value==max_x_lab ~ pred_perf_value+0.025*pred_perf_value,
        TRUE ~ NA)) %>% 
      mutate(x_lab = if_else(x_lab <0, 0, x_lab)) %>% 
      # > plot
      ggplot(., aes(x = pred_perf_value_lab, y = Model)) +
      geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                     color = "grey8", size=0.5) +
      geom_point(aes(color = gpe_model, shape = best), size = 3) + 
      #geom_point(aes(color = gpe_model), shape="|", size=3) + 
      geom_text(aes(x = x_lab, label = order), 
                nudge_y = 0.4, 
                size = 2) + 
      facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", space = "free", switch = "y") + 
      theme_bw() +
      theme(legend.position = 'bottom', 
            axis.title.y = element_blank(),
            strip.placement = "outside", 
            strip.text.y.left = element_text(angle = 0, size=9),
            strip.text.x = element_text(size=10), 
            strip.background = element_rect(fill="grey75", color = "black"),
            title = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            legend.text = element_text(size=10),
            panel.grid = element_blank()) + 
      scale_shape_manual(values = c(8, 20), name = "Rank model") + 
      scale_color_viridis_d(option = "A", begin = 0.2, end = 0.75, name = "Model family", direction = -1) + 
      scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2)) + 
      guides(color = guide_legend(order = 1, override.aes = list(size=3, shape=20), ncol = 2),
             shape = guide_none()) + 
      labs(x = "Root mean square error of prediction\n(lower value indicates a better predictive performance)") ; p.rmsep
    
    # > Results based on mean of results of each CV procedure
    p.rmsep_mean <- tab_perf_labelled_mean_i %>% 
      filter(pred_perf == "RMSEP") %>% 
      # > order of models 
      arrange(Outcome_lab, Country_lab, mean_pred_perf) %>% 
      group_by(Outcome_lab, Country_lab) %>% 
      mutate(order = row_number()) %>%
      mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
             best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
      # > min and max among lm and rf
      group_by(Outcome_lab, Country_lab, Model) %>% 
      mutate(min_x_lab = min(mean_pred_perf),
             max_x_lab = max(mean_pred_perf)) %>% 
      # > x coordinates for labels 
      mutate(x_lab = case_when(
        mean_pred_perf==min_x_lab ~ mean_pred_perf-0.025*mean_pred_perf,
        mean_pred_perf==max_x_lab ~ mean_pred_perf+0.025*mean_pred_perf,
        TRUE ~ NA)) %>% 
      mutate(x_lab = if_else(x_lab <0, 0, x_lab)) %>% 
      # > plot
      ggplot(., aes(x = mean_pred_perf_lab, y=Model)) +
      geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                     color = "grey8", size=0.5) +
      geom_point(aes(color = gpe_model, shape = best), size = 3) + 
      #geom_point(aes(color = gpe_model), shape="|", size=3) + 
      geom_text(aes(x = x_lab, label = order), 
                nudge_y = 0.4, 
                size = 2) + 
      facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", space = "free", switch = "y") + 
      theme_bw() +
      theme(legend.position = 'bottom', 
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            strip.placement = "outside", 
            strip.text.y.left = element_blank(),
            strip.text.x = element_text(size=10), 
            strip.background.x = element_rect(fill="grey75", color = "black"),
            strip.background.y = element_blank(),
            title = element_text(size=10),
            legend.text = element_text(size=10),
            #panel.grid.major.x = element_line(color = "grey", linewidth = 0.01),
            #panel.grid.minor.x = element_line(color = "grey", linewidth = 0.01),
            panel.grid = element_blank()) + 
      scale_shape_manual(values = c(8, 20), name = "Rank model") + 
      scale_color_viridis_d(option = "A", begin = 0.2, end = 0.75, name = "Model family", direction = -1) + 
      scale_size_continuous(range=c(0,3), guide = guide_none()) +
      #scale_fill_manual(values = c("red", "white")) + 
      scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2)) + 
      guides(color = guide_legend(order = 1, override.aes = list(size=3, shape=20), ncol = 2),
             shape = guide_none())    + 
      labs(x = "Root mean square error of prediction\n(lower value indicates a better predictive performance)") ; p.rmsep_mean
    
    # > Get legend of the 1rst plot
    leg <- get_legend(p.rmsep_mean)
    
    # > Merge plots and add common x label
    plot <- plot_grid(p.rmsep + theme(legend.position = "none", axis.title.x = element_blank()), 
                      p.rmsep_mean + theme(legend.position = "none", axis.title.x = element_blank()), 
                      p.n_pred,
                      nrow=1, rel_widths = c(2.75, 1, 0.15), axis = "tblr", align = "h")
    
    plot_x <- ggdraw(add_sub(plot, "Root mean square error of prediction", 
                             vpadding=grid::unit(0,"lines"),y=0, x=0.6, vjust=0, size = 10))
    
    # > Add common legend and save 
    FigRMSEP <- plot_grid(plot_x, 
                        NULL,
                        leg, 
                        NULL,
                        ncol = 1, rel_heights = c(0.9,0.025,0.05,0.025))

    # > Save 
    country_save <- unique(tab_perf_labelled_i$Country)
    ggsave(plot = FigNSE, 
           filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/Supplementary_Figure_", country_save, "_", outcome, "_nse.png"), 
           width = 9.5, height=12, dpi=300, bg="white")
    
    ggsave(plot = FigRMSEP, 
           filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/Supplementary_Figure_", country_save, "_", outcome, "_rmsep.png"), 
           width = 9.5, height=12, dpi=300, bg="white")
    
  }
  
}

# --------------------------------------------------------------------------
# MEANS and SDs of NSE of models 

# -----------------------------
# A. Per family of model 

tab_p1 <- tab_perf_labelled %>% 
  filter(pred_perf=="NSE", 
         Country=="WORLD", 
         Outcome=="01_Ya") %>% 
  mutate(Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) %>% 
  mutate(gpe_model = recode(gpe_model, "Linear regression"="Linear\nregression", "Random forest"="Random\nForest")) %>% 
  group_by(Outcome_lab, Type_cv_lab, gpe_model) %>%
  mutate(n_models = n()) %>% 
  mutate(y_lab = 1.05) %>%
  mutate(NSE_lab = paste0(round(mean(pred_perf_value),2), " (", round(sd(pred_perf_value),2), ")")) %>%
  mutate(N_lab = paste0("n=", n_models))

p0 <- tab_p1 %>% 
  filter(Type_cv_lab == "Cross-validation on years") %>% 
  ggplot(data=.) +
  geom_text(aes(x = gpe_model, y = y_lab, label=gpe_model), 
            check_overlap = T, position = position_dodge(width = 1), size=4, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") +
  ggtitle(" ") + 
  labs(y=" ") +
  coord_flip(clip = 'off') + 
  #theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(), 
        panel.background = element_rect(fill="grey75", color="black", size=0.5),
        plot.margin = unit(c(0,0.1,0,0), "cm")) ; p0

p1_left <- ggplot(data=tab_p1) +
  geom_boxplot(aes(x = gpe_model, y = pred_perf_value),
               width=0.25, position = position_dodge(width = 1), outlier.shape = NA) + 
  geom_text(aes(x = gpe_model, y = y_lab, label=NSE_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "transparent"),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill="grey75", color="black", size=0.5),
        strip.text.y = element_blank(),
        legend.position = "none",
        plot.tag.position = c(0.9, 1.025), 
        plot.tag = element_text(size=9), 
        panel.grid = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_y_continuous(breaks = seq(0,1,by=0.25), labels = seq(0,1,by=0.25), limits = c(0,1.1)) + 
  labs(y=" ") +
  ggtitle("a.") + 
  coord_flip() ; p1_left

p1_right <- tab_p1 %>% 
  filter(Type_cv_lab == "Cross-validation on years") %>% 
  ggplot(data=.) +
  geom_text(aes(x = gpe_model, y = y_lab, label=N_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") +
  ggtitle(" ") + 
  labs(y=" ") +
  coord_flip() + 
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) ; p1_right

p1 <- cowplot::plot_grid(p0, p1_left, p1_right, nrow = 1, rel_widths = c(0.15, 0.75,0.1), axis = "tblr", align = "hv") ; p1

#ggsave(plot=p1, 
#       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/mean_sd_nse_world_ya.png"), 
#       width = 10, height=5, dpi=300, bg="white")

# -----------------------------
# B. Per Dimension reduction technique

tab_p2 <- tab_perf_labelled %>% 
  filter(pred_perf=="NSE", 
         Country=="WORLD", 
         Outcome=="01_Ya") %>% 
  mutate(Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) %>% 
  mutate(gpe_model = recode(gpe_model, "Linear regression"="Linear\nregression", "Random forest"="Random\nForest")) %>% 
  mutate(dimension_reduction = recode(dimension_reduction, 'Standardised\nmean'='Std. mean')) %>%
  mutate(dimension_reduction = factor(dimension_reduction, levels = rev(c('Std. mean', 'Averages', 'PCA', 'FPCA', 'MFPCA', 'PLS', 'FPLS')))) %>% 
  mutate(dimension_reduction_lab = if_else(Type_cv_lab == "Cross-validation on years", dimension_reduction, NA)) %>%
  group_by(Outcome_lab, Type_cv_lab, gpe_model, dimension_reduction) %>%
  mutate(n_models = n()) %>% 
  mutate(y_lab = 1.05) %>%
  mutate(NSE_lab = paste0(round(mean(pred_perf_value),2), " (", round(sd(pred_perf_value),2), ")")) %>% 
  mutate(N_lab = paste0("n=", n_models))

p2_left <- tab_p2 %>%
  ggplot(.) +
  geom_boxplot(aes(x = dimension_reduction, y = pred_perf_value, fill = dimension_reduction),
               width=0.5, position = position_dodge(width = 1), outlier.shape = NA) + 
  geom_text(aes(x = dimension_reduction, y = y_lab, label=NSE_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.placement = "outside", 
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill="grey75", color="black", size=0.5),
        strip.text.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_viridis_d(begin = 0.85,
                        end=0.2,
                        option="B") +
  scale_y_continuous(breaks = seq(0,1,by=0.25), labels = seq(0,1,by=0.25), limits = c(0,1.1)) + 
  labs(y=" ") +
  ggtitle("b.") +
  coord_flip() ; p2_left

p2_right <- tab_p2 %>% 
  filter(Type_cv_lab == "Cross-validation on years") %>% 
  ggplot(data=.) +
  geom_text(aes(x = dimension_reduction, y = y_lab, label=N_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") +
  ggtitle(" ") + 
  labs(y=" ") +
  coord_flip() + 
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"))  ; p2_right

p2 <- cowplot::plot_grid(p0, p2_left, p2_right, nrow = 1, rel_widths = c(0.15, 0.75, 0.1), axis = "tblr", align = "hv") ; p2

#ggsave(plot=p2, 
#       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/mean_sd_nse_world_ya.png"), 
#       width = 10, height=5, dpi=300, bg="white")

# -----------------------------
# C. Per temporal resolution

tab_p3 <- tab_perf_labelled %>% 
  filter(pred_perf=="NSE", 
         Country=="WORLD", 
         Outcome=="01_Ya") %>% 
  mutate(gpe_model = recode(gpe_model, "Linear regression"="Linear\nregression", "Random forest"="Random\nForest")) %>%
  mutate(data_type = recode(data_type, 
                            "Annual climatic predictors"="Seasonal",
                            "Monthly climatic predictors"="Monthly",
                            "Daily climatic predictors"="Daily"),
         data_type = factor(data_type, levels = rev(c("Seasonal", "Monthly", "Daily")))) %>%
  group_by(Outcome_lab, Type_cv_lab, gpe_model, dimension_reduction) %>%
  filter(dimension_reduction != "Standardised\nmean") %>% 
  group_by(Outcome_lab, Type_cv_lab, gpe_model, data_type) %>%
  mutate(y_lab = 1.05) %>%
  mutate(n_models = n()) %>% 
  mutate(NSE_lab = paste0(round(mean(pred_perf_value),2), " (", round(sd(pred_perf_value),2), ")")) %>% 
  mutate(N_lab = paste0("n=", n_models))

p3_left <- tab_p3 %>%
  ggplot(.) +
  geom_boxplot(aes(y = pred_perf_value, x = data_type, fill = data_type),
               width=0.5, position = position_dodge(width = 1), outlier.shape = NA) + 
  geom_text(aes(x = data_type, y = y_lab, label=NSE_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") + 
  #facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.placement = "outside", 
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill="grey75", color="black", size=0.5),
        legend.position = "none",
        panel.grid = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_viridis_d(begin = 0.75, end=0.2, option="G") +
  scale_y_continuous(breaks = seq(0,1,by=0.25), labels = seq(0,1,by=0.25), limits = c(0,1.1)) + 
  labs(y="Nash-Sutcliffe model efficiency (NSE)") +
  coord_flip() + 
  ggtitle("c."); p3_left

p3_right <- tab_p3 %>% 
  filter(Type_cv_lab == "Cross-validation on years") %>% 
  ggplot(data=.) +
  geom_text(aes(x = data_type, y = y_lab, label=N_lab), 
            check_overlap = T, position = position_dodge(width = 1), size=3, show.legend = FALSE) +
  facet_grid(gpe_model ~Type_cv_lab, scales = "free") +
  ggtitle(" ") + 
  labs(y=" ") +
  coord_flip(clip = 'off') + 
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) ; p3_right

p3 <- cowplot::plot_grid(p0, p3_left, p3_right, nrow = 1, rel_widths = c(0.15, 0.75, 0.1), axis = "tblr", align = "hv") ; p3

#ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/mean_sd_nse_world_ya3.png"), 
#       width = 14, height=4, dpi=300, bg="white")

p <- plot_grid(p1, p2, p3,
          ncol = 1, 
          rel_heights = c(0.2, 0.5, 0.23)
) ; p

save_plot(plot = p, 
          filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/07_MODELS_PERF/mean_sd_nse_world_ya2.png"), 
          base_width = 14, base_height=8, dpi=300, bg="white")

# --------------------------------------------------------------------------
# Calibration, residuals of best models 

# Functions 
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

# -------------------------
# > Best models 
best_models <- tab_perf_labelled_mean %>% 
  filter(pred_perf=="NSE") %>%
  group_by(Outcome_lab, Country_lab, pred_perf_lab) %>% 
  arrange(desc(mean_pred_perf)) %>% 
  slice(1) %>% 
  mutate(Model = if_else(gpe_model=="Random forest", paste0("rf_", Model), paste0("lm_", Model)))

best_models %>% dplyr::select(Outcome_lab, Country_lab, Model, gpe_model) %>% arrange(Country_lab)

#   pred_perf_lab             Outcome_lab           Country_lab                           Model      gpe_model    
# 1 "NSE\n(higher is better)" Soybean yield         "Global\n(N=122229)"                  pca.m.2    Random forest
# 2 "NSE\n(higher is better)" Soybean yield anomaly "Global\n(N=122229)"                  pca.m.all  Random forest
# 3 "NSE\n(higher is better)" Soybean yield         "United-States of America\n(N=29803)" pca.m.3    Random forest
# 4 "NSE\n(higher is better)" Soybean yield anomaly "United-States of America\n(N=29803)" avg.m      Random forest
# 5 "NSE\n(higher is better)" Soybean yield         "Brazil\n(N=14757)"                   pca.m.2    Random forest
# 6 "NSE\n(higher is better)" Soybean yield anomaly "Brazil\n(N=14757)"                   fpca.m.all Random forest

# -------------------------
# > Select predictions of best models 
best_models_preds <- tab_preds %>% 
  # > best model for each dataset
  mutate(best = case_when(
    Outcome == "01_Ya" & Country == "USA" & Model == best_models[which(best_models$Country == "USA"   & best_models$Outcome == "01_Ya"),]$Model ~ 1,
    Outcome == "01_Ya" & Country == "BRA"       & Model == best_models[which(best_models$Country == "BRA"   & best_models$Outcome == "01_Ya"),]$Model ~ 1,
    Outcome == "01_Ya" & Country == "WORLD"     & Model == best_models[which(best_models$Country == "WORLD" & best_models$Outcome == "01_Ya"),]$Model ~ 1,
    Outcome == "02_Ya_ano" & Country == "USA"   & Model == best_models[which(best_models$Country == "USA"   & best_models$Outcome == "02_Ya_ano"),]$Model ~ 1,
    Outcome == "02_Ya_ano" & Country == "BRA"   & Model == best_models[which(best_models$Country == "BRA"   & best_models$Outcome == "02_Ya_ano"),]$Model ~ 1,
    Outcome == "02_Ya_ano" & Country == "WORLD" & Model == best_models[which(best_models$Country == "WORLD" & best_models$Outcome == "02_Ya_ano"),]$Model ~ 1,
    TRUE ~ 0
  )) %>% 
  filter(best == 1) %>% 
  # > compute residuals
  mutate(Residuals = Predicted - Observed) %>% 
  # > compute RMSEP and NSE
  group_by(Outcome, Country, Model, Type_cv) %>% 
  mutate(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
  ) %>%
  mutate(RMSEP = round(RMSEP, 3),
         NSE = round(NSE, 3)) %>% 
  mutate(label_perf = paste0("Name model: ", Model, "\nRMSEP: ", RMSEP, " tonnes/hectare\nNSE: ", NSE)) %>% 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) %>% 
  split(list(.$Country, .$Outcome))
  
# -------------------------
# > Plots
best_models_preds[6] %>% 
  map(., ~{
    
    # > Plot details and graph settings
    country <- unique(.x$Country)
    outcome <- unique(.x$Outcome)
    if(outcome == "01_Ya"){ 
      limsx   <- c(0,7) 
      pos_lab_x <- 0
      pos_lab_y <- 6
      lim_res  <- 5
    }
    if(outcome == "02_Ya_ano"){ 
      limsx <- c(-2.2,2.2) 
      pos_lab_x <- -2.2
      pos_lab_y <- 2
      lim_res  <- 2.2
      }
    
    # > Calibration plot 
    p.pred.obs <- .x %>% split(.$Type_cv) %>% 
          map_dfr(., ~{
            
            density <- get_density(x = .x$Observed, y = .x$Predicted, n=100)
            .x %>% 
              mutate(density=density)
          }, .id = "Type_cv") %>% 
      ggplot(data = .) + 
      geom_point(aes(x = Observed, y = Predicted, color = density), 
                 size = 0.5) + 
      geom_abline(color = "red") + 
      geom_text(aes(label = label_perf),
                x = pos_lab_x, y = pos_lab_y, check_overlap = TRUE, hjust = 0, size = 3) + 
      theme_cowplot() + 
      theme(legend.position = "bottom",
            strip.text.x = element_text(size=10), 
            #strip.background.x = element_blank(),
            title = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            legend.text = element_text(size=10)
            #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
      ) + 
      scale_color_viridis_c(option = 3, name="Density per combinaison of\nobserved/predicted yield values", guide = guide_colorbar(barwidth = 7.5, barheight = 0.35)) +
      facet_wrap(Type_cv_lab ~ ., ncol=1, scales = "free") +
      labs(x = "Observed yields (tons/hectare)", y = "Predicted yields (tons/hectare)") + 
      ggtitle("a. ") + 
      lims(x = limsx, y = limsx)
    
    # > Density of observation vs prediction 
    p.dens.pred.obs <- .x %>% 
      dplyr::select(-Residuals, -RMSEP, -NSE) %>% 
      gather(key="Type", value = "Yield", Predicted, Observed) %>% 
      ggplot(.) + 
      geom_density(aes(x = Yield, color = Type, fill = Type), alpha=0.2) + 
      theme_cowplot() + 
      theme(legend.position = "bottom",
            strip.text.x = element_text(size=10), 
            #strip.text.y = element_blank(), 
            #strip.background = element_blank(),
            title = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10),
            legend.text = element_text(size=10)
            #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
      ) + 
      facet_wrap(Type_cv_lab ~ ., ncol=1, scales = "free") +
      #facet_wrap(.~Country_lab, scales = "free_y") +
      scale_color_manual(values = c("black","#4DAF4A"), name ="") +
      scale_fill_manual(values = c("black","#4DAF4A"), name ="") +
      labs(x = "Yield (tons/hectare)", y = "Density") + 
      ggtitle("b. ") +
      lims(x = limsx)
    
    # > Residuals
    p.res <- .x %>% 
      # mean residuals 
      group_by(Country) %>% 
      mutate(mean_Residuals = mean(Residuals)) %>% 
      ggplot(.) + 
      geom_histogram(aes(x = Residuals), 
                     bins = 100, alpha = 0.5, fill = "darkgrey", color = "transparent") + 
      geom_vline(aes(xintercept = mean_Residuals), color = "black", 
                 linetype = 2, linewidth = 0.5) +
      theme_cowplot() + 
      theme(legend.position = "none",
            #strip.text = element_blank(), 
            #strip.background = element_blank(),
            title = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title.x = element_text(size=10)) + 
      facet_wrap(Type_cv_lab ~ ., ncol=1, scales = "free") +
      #facet_wrap(.~Country_lab, scales = "free_y") + 
      ggtitle("c. ") +
      labs(x = "Predicted - observed yields (tons/hectare)", y = "Counts") +
      lims(x = c(-lim_res, lim_res))
    
    # > FULL FIGURE
    FIG_CALIB <- plot_grid(p.pred.obs + theme(strip.text.y = element_blank(),
                                         strip.background = element_blank()), 
                      ggplot() + theme_void(),
                      p.dens.pred.obs + theme(strip.text.y = element_blank(),
                                              strip.background = element_blank()), 
                      ggplot() + theme_void(),
                      p.res + theme(strip.text.y = element_blank(),
                                    strip.background = element_blank()), 
                      nrow = 1, align = "h", axis="tblr", rel_widths = c(1, 0.1, 1, 0.1, 1))
    
    # > save
    ggsave(plot = FIG_CALIB, 
           filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/08_CALIB_PLOTS/", country, "_", outcome, ".png"), 
           width = 13, height=8, dpi=300, bg="white")
    
  })

# --------------------------------------------------------------------------
# IMPORTANCE PLOTS 

best_models %>% dplyr::select(Outcome_lab, Country_lab, Model, gpe_model) %>% arrange(Country_lab)
#   pred_perf_lab             Outcome_lab           Country_lab                           Model      gpe_model    
# 1 "NSE\n(higher is better)" Soybean yield         "Global\n(N=122229)"                  pca.m.2    Random forest
# 2 "NSE\n(higher is better)" Soybean yield anomaly "Global\n(N=122229)"                  pca.m.all  Random forest
# 3 "NSE\n(higher is better)" Soybean yield         "United-States of America\n(N=29803)" pca.m.3    Random forest
# 4 "NSE\n(higher is better)" Soybean yield anomaly "United-States of America\n(N=29803)" avg.m      Random forest
# 5 "NSE\n(higher is better)" Soybean yield         "Brazil\n(N=14757)"                   pca.m.2    Random forest
# 6 "NSE\n(higher is better)" Soybean yield anomaly "Brazil\n(N=14757)"                   fpca.m.all Random forest

# Fit 
form.avg.m   <- paste0(names(tab_world %>% dplyr::select(starts_with(paste("monthly_", unique(vars_names$clim.var), sep="")))), collapse = " + ")
form.pca.m.2 <- paste0(names(tab_world %>% dplyr::select(starts_with("PC1"), starts_with("PC2")) %>% dplyr::select(contains("month"))), collapse = " + ")
form.pca.m.3 <- paste0(names(tab_usa %>% dplyr::select(starts_with("PC1"), starts_with("PC2"), starts_with("PC3")) %>% dplyr::select(contains("month"))), collapse = " + ")
form.pca.m.all  <- paste0(names(tab_bra %>% dplyr::select(starts_with("PC")) %>% dplyr::select(contains("month"))), collapse = " + ")
form.fpca.m.all <- paste0(names(tab_bra %>% dplyr::select(starts_with("FPC")) %>% dplyr::select(contains("month"))), collapse = " + ")

library(ranger)

# WORLD
tab_fit <- tab_world
set.seed(101)
mod_world_ya <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.2)), data=tab_fit, num.tree=500, importance="impurity")
set.seed(101)
mod_world_ya_ano <- ranger(as.formula(paste0("Ya_ano ~ irrigated_portion + ", form.pca.m.all)), data=tab_fit, num.tree=500, importance="impurity")

tab_imp <- rbind(
  data.frame(Predictor = names(mod_world_ya$variable.importance), 
             Importance = mod_world_ya$variable.importance, 
             Outcome = "Soybean yield"),
  data.frame(Predictor = names(mod_world_ya_ano$variable.importance), 
             Importance = mod_world_ya_ano$variable.importance, 
             Outcome = "Soybean yield anomaly")) %>%
  mutate(Outcome = factor(Outcome, levels = c("Soybean yield", "Soybean yield anomaly"))) %>%
  group_by(Outcome) %>%
  arrange(desc(Importance)) %>%
  mutate(order = 1:n()) %>% 
  mutate(Importance_lab = ifelse(Outcome == "Soybean yield", Importance+1000, Importance+20)) 

p1 <- tab_imp %>% 
  filter(Outcome == "Soybean yield") %>% 
  mutate(Predictor_lab = recode(Predictor, 
                                'PC1_month_prec'='Score 1 - average precipitation', 
                                'PC1_month_min_temp'='Score 1 - min. temperature', 
                                'PC1_month_max_temp'='Score 1 - max. temperature', 
                                'PC1_month_vapor_pressure_deficit'='Score 1 - vapor pressure deficit', 
                                'PC1_month_rad'='Score 1 - radiation', 
                                'PC1_month_et0'='Score 1 - ref. evapotranspiration', 
                                'irrigated_portion'='Irrigation', 
                                'PC2_month_et0'='Score 2 - ref. evapotranspiration', 
                                'PC2_month_rad'='Score 2 - radiation', 
                                'PC2_month_min_temp'='Score 2 - min. temperature', 
                                'PC2_month_vapor_pressure_deficit'='Score 2 - vapor pressure deficit', 
                                'PC2_month_max_temp'='Score 2 - max. temperature', 
                                'PC2_month_prec'='Score 2 - average precipitation')) %>% 
  ggplot(., aes(y = reorder(Predictor_lab, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance), width=0.75) +
  #geom_text(aes(x = Importance_lab, label = order),size = 3) + 
  #facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

p2 <- tab_imp %>% 
  filter(Outcome != "Soybean yield") %>% 
  ggplot(., aes(y = reorder(Predictor, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance), width=0.75) +
  geom_text(aes(x = Importance_lab, label = order),
            size = 3) + 
  facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

plot_grid(p1, p2, nrow=1)

# > save
ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/WORLD_imp.png"), 
       width = 12, height=6, dpi=300, bg="white")
ggsave(p1 + 
         theme(strip.background = element_blank(), 
               strip.text = element_blank(),
               plot.caption = element_text(hjust=0, size=6),
               plot.caption.position = "plot"), #+ labs(caption = "Abbreviations: 'PC1_month_X' and 'PC2_month_X': first and second components derived from standard principal component analysis (PCA) applied on\nmonthly averages of the climatic variable X, respectively.\nClimatic variables abbreviations are; 'max_temp': maximum temperature (°C); 'min_temp': minimum temperature (°C); 'prec': total precipitations (mm);\n'rad': solar radiations (MJ/m²); 'vapor_pressure_deficit': vapor pressure deficit; 'et0': reference evapotranspiration (mm/day).\nOther predictor is 'irrigated_portion': fractional area of irrigated soybean in the grid-cell."), 
         filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/WORLD_imp1.png"), width = 6, height=5, dpi=300, bg="white")

# USA
tab_fit <- tab_usa
set.seed(101)
mod_usa_ya <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.3)), data=tab_fit, num.tree=500, importance="impurity")
set.seed(101)
mod_usa_ya_ano <- ranger(as.formula(paste0("Ya_ano ~ irrigated_portion + ", form.avg.m)), data=tab_fit, num.tree=500, importance="impurity")

tab_imp_usa <- rbind(
  data.frame(Predictor = names(mod_usa_ya$variable.importance), 
             Importance = mod_usa_ya$variable.importance, 
             Outcome = "Soybean yield"),
  data.frame(Predictor = names(mod_usa_ya_ano$variable.importance), 
             Importance = mod_usa_ya_ano$variable.importance, 
             Outcome = "Soybean yield anomaly")) %>%
  mutate(Outcome = factor(Outcome, levels = c("Soybean yield", "Soybean yield anomaly"))) %>%
  group_by(Outcome) %>%
  arrange(desc(Importance)) %>%
  mutate(order = 1:n()) %>% 
  mutate(Importance_lab = ifelse(Outcome == "Soybean yield", Importance+1000, Importance+2)) 

p1 <- tab_imp_usa %>% 
  filter(Outcome == "Soybean yield") %>% 
  mutate(Predictor_lab = recode(Predictor, 
                                'PC1_month_prec'='Score 1 - average precipitation', 
                                'PC1_month_min_temp'='Score 1 - min. temperature', 
                                'PC1_month_max_temp'='Score 1 - max. temperature', 
                                'PC1_month_vapor_pressure_deficit'='Score 1 - vapor pressure deficit', 
                                'PC1_month_rad'='Score 1 - radiation', 
                                'PC1_month_et0'='Score 1 - ref. evapotranspiration', 
                                'PC2_month_et0'='Score 2 - ref. evapotranspiration', 
                                'PC2_month_rad'='Score 2 - radiation', 
                                'PC2_month_min_temp'='Score 2 - min. temperature', 
                                'PC2_month_vapor_pressure_deficit'='Score 2 - vapor pressure deficit', 
                                'PC2_month_max_temp'='Score 2 - max. temperature', 
                                'PC2_month_prec'='Score 2 - average precipitation',
                                'PC3_month_prec'='Score 3 - average precipitation', 
                                'PC3_month_min_temp'='Score 3 - min. temperature', 
                                'PC3_month_max_temp'='Score 3 - max. temperature', 
                                'PC3_month_vapor_pressure_deficit'='Score 3 - vapor pressure deficit', 
                                'PC3_month_rad'='Score 3 - radiation', 
                                'PC3_month_et0'='Score 3 - ref. evapotranspiration', 
                                'irrigated_portion'='Irrigation')) %>% 
  ggplot(., aes(y = reorder(Predictor_lab, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance)) +
  #geom_text(aes(x = Importance_lab, label = order), size = 3) + 
  #facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

p2 <- tab_imp_usa %>% 
  filter(Outcome != "Soybean yield") %>% 
  ggplot(., aes(y = reorder(Predictor, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance)) +
  geom_text(aes(x = Importance_lab, label = order),
            size = 3) + 
  facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

plot_grid(p1, p2, nrow=1)

# > save
ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/USA_imp.png"), 
       width = 13, height=6, dpi=300, bg="white")
p1usa <- p1 + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.caption = element_text(hjust=0, size=6),
        plot.caption.position = "plot") #+ labs(caption = "Abbreviations: 'PC1_month_X', 'PC2_month_X', and 'PC3_month_X': first, second, and third components derived from standard principal component\nanalysis (PCA) applied on monthly averages of the climatic variable X, respectively.\nClimatic variables abbreviations are; 'max_temp': maximum temperature (°C); 'min_temp': minimum temperature (°C); 'prec': total precipitations (mm);\n'rad': solar radiations (MJ/m²); 'vapor_pressure_deficit': vapor pressure deficit; 'et0': reference evapotranspiration (mm/day).\nOther predictor is 'irrigated_portion': fractional area of irrigated soybean in the grid-cell.")
ggsave(p1usa, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/USA_imp1.png"), width = 6, height=5, dpi=300, bg="white")


# BRA
tab_fit <- tab_bra
set.seed(101)
mod_bra_ya <- ranger(as.formula(paste0("Ya ~ irrigated_portion + ", form.pca.m.2)), data=tab_fit, num.tree=500, importance="impurity")
set.seed(101)
mod_bra_ya_ano <- ranger(as.formula(paste0("Ya_ano ~ irrigated_portion + ", form.fpca.m.all)), data=tab_fit, num.tree=500, importance="impurity")

tab_imp_bra <- rbind(
  data.frame(Predictor = names(mod_bra_ya$variable.importance), 
             Importance = mod_bra_ya$variable.importance, 
             Outcome = "Soybean yield"),
  data.frame(Predictor = names(mod_bra_ya_ano$variable.importance), 
             Importance = mod_bra_ya_ano$variable.importance, 
             Outcome = "Soybean yield anomaly")) %>%
  mutate(Outcome = factor(Outcome, levels = c("Soybean yield", "Soybean yield anomaly"))) %>%
  group_by(Outcome) %>%
  arrange(desc(Importance)) %>%
  mutate(order = 1:n()) %>% 
  mutate(Importance_lab = ifelse(Outcome == "Soybean yield", Importance+1000, Importance+1)) 

p1 <- tab_imp_bra %>% 
  filter(Outcome == "Soybean yield") %>% 
  mutate(Predictor_lab = recode(Predictor, 
                                'PC1_month_prec'='Score 1 - average precipitation', 
                                'PC1_month_min_temp'='Score 1 - min. temperature', 
                                'PC1_month_max_temp'='Score 1 - max. temperature', 
                                'PC1_month_vapor_pressure_deficit'='Score 1 - vapor pressure deficit', 
                                'PC1_month_rad'='Score 1 - radiation', 
                                'PC1_month_et0'='Score 1 - ref. evapotranspiration', 
                                'PC2_month_et0'='Score 2 - ref. evapotranspiration', 
                                'PC2_month_rad'='Score 2 - radiation', 
                                'PC2_month_min_temp'='Score 2 - min. temperature', 
                                'PC2_month_vapor_pressure_deficit'='Score 2 - vapor pressure deficit', 
                                'PC2_month_max_temp'='Score 2 - max. temperature', 
                                'PC2_month_prec'='Score 2 - average precipitation',
                                'PC3_month_prec'='Score 3 - average precipitation', 
                                'PC3_month_min_temp'='Score 3 - min. temperature', 
                                'PC3_month_max_temp'='Score 3 - max. temperature', 
                                'PC3_month_vapor_pressure_deficit'='Score 3 - vapor pressure deficit', 
                                'PC3_month_rad'='Score 3 - radiation', 
                                'PC3_month_et0'='Score 3 - ref. evapotranspiration', 
                                'irrigated_portion'='Irrigation')) %>% 
  ggplot(., aes(y = reorder(Predictor_lab, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance)) +
  #geom_text(aes(x = Importance_lab, label = order), size = 3) + 
  #facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank()) ; p1

p2 <- tab_imp_bra %>% 
  filter(Outcome != "Soybean yield") %>% 
  ggplot(., aes(y = reorder(Predictor, order, min, decreasing = T))) + 
  geom_col(aes(x = Importance)) +
  geom_text(aes(x = Importance_lab, label = order),
            size = 3) + 
  facet_wrap(.~Outcome, scales = "free") + 
  theme_bw() + 
  theme(axis.title.y = element_blank()) ; p2

plot_grid(p1, p2, nrow=1)

# > save
ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/BRA_imp.png"), 
       width = 13, height=6, dpi=300, bg="white")

p1bra <- p1 + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.caption = element_text(hjust=0, size=6),
        plot.caption.position = "plot")# + labs(caption = "Abbreviations: 'PC1_month_X' and 'PC2_month_X': first and second components derived from standard principal component analysis (PCA) applied on\nmonthly averages of the climatic variable X, respectively.\nClimatic variables abbreviations are; 'max_temp': maximum temperature (°C); 'min_temp': minimum temperature (°C); 'prec': total precipitations (mm);\n'rad': solar radiations (MJ/m²); 'vapor_pressure_deficit': vapor pressure deficit; 'et0': reference evapotranspiration (mm/day).\nOther predictor is 'irrigated_portion': fractional area of irrigated soybean in the grid-cell.")
ggsave(p1bra, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/BRA_imp1.png"), width = 6, height=5, dpi=300, bg="white")

plot_grid(p1usa+ggtitle("a. US - pca.m.3"), 
          p1bra+ggtitle("b. Brazil - pca.m.2") + theme(plot.caption = element_blank()), nrow=1,
          align = "hv", axis="tlbr")

ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/09_IMP_PLOTS/USA_BRA_imp1.png"), 
       width = 13, height=6, dpi=300, bg="white")

# --------------------------------------------------------------------------
# PATIAL DEPENDENCY PLOTS



partial(cforest_adjusted, 
        pred.var = c("avg_mtg_duration", "avg_mtg_attd"), 
        trim.outliers = TRUE, chull = TRUE, parallel = TRUE,
        grid.resolution = 30,  paropts = list(.packages = "ranger"))


# --------------------------------------------------------------------------
# CLIMATE CHANGE 
tab.preds.sensi.usa.cc <- preds.cc %>% filter(country=="01_USA")
tab.preds.sensi.bra.cc <- preds.cc %>% filter(country=="02_BRA")

# Maps showing the median difference (only USA, not desert)
cc.usa.1 <- tab.preds.sensi.usa.cc %>% 
  filter(gpe_model=="Random forest") %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-110, -57), y=c(25,50)) + 
  facet_grid(model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        legend.position = "bottom") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage), max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)     ") +
  ggtitle("a. Random forest") 

cc.usa.2 <- tab.preds.sensi.usa.cc %>% 
  filter(gpe_model!="Random forest") %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-110, -57), y=c(25,50)) + 
  facet_grid(model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        legend.position = "bottom") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage), max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)     ") +
  ggtitle("b. Linear regression") 

usa.cc <- plot_grid(cc.usa.1, 
                    cc.usa.2 + theme(legend.position = "none"),
                    ncol=2, align = "hv", axis = "tblr") 

# > save
ggsave(usa.cc, filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/USA_cc2.png"), 
       width = 14, height=7, dpi=300, bg="white")


# Maps showing the median difference (only USA, not desert)
cc.bra.1 <- tab.preds.sensi.bra.cc %>% 
  filter(gpe_model=="Random forest") %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-85, -35), y=c(-35,5))  + 
  facet_grid(model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        legend.position = "none") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5, title = "Difference relative to median yield\nin current climatic conditions"),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage), max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)  ") +
  ggtitle("a. Random forest") 

cc.bra.2 <- tab.preds.sensi.bra.cc %>% 
  filter(gpe_model!="Random forest") %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-85, -35), y=c(-35,5))  + 
  facet_grid(model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "transparent"),
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside",
        legend.position = "bottom") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5, title = "Difference relative to median yield\nin current climatic conditions"),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage), max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)  ") +
  ggtitle("b. Linear regression")

cc.bra <- plot_grid(cc.bra.1,
                    cc.bra.2 + theme(legend.position = "none"), 
                    ncol=2, align = "hv", axis = "tblr")

# > save
ggsave(cc.bra, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/BRA_cc2.png"), 
       width = 14, height=7, dpi=300, bg="white")


# additional boxplots 
p1 <- tab.preds.sensi.usa.cc %>% 
  mutate(model     = factor(model, levels = rev(c("pca.m.3",  "pca.m.2",  "pca.d.3", "pca.d.all", 
                                                  "fpca.m.3", "fpca.m.2", "pls.m.2", "pls.m.all", 
                                                  "pls.d.3",  "pls.d.all", "avg.m")))) %>%
  group_by(model, gpe_model, dT) %>% 
  mutate(M=round(mean(Difference_Prediction_percentage), 2)) %>% 
  ggplot(., aes(x=model)) + 
  geom_hline(yintercept=0, color="black", lty=2) +
  geom_boxplot(aes(y=Difference_Prediction_percentage, color=gpe_model), outlier.size = 0.4) + 
  #geom_point(aes(y = M), color = "red") + 
  geom_text(aes(y = 90, label = M), size=3, check_overlap = T) + 
  facet_grid(gpe_model ~ dT, switch = "y") + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside", 
        axis.title.y = element_blank(),
        #plot.margin = unit(c(0,0,0,0), units = "cm"),
        strip.background.y = element_blank()) + 
  ggtitle("a. United States of America") + 
  lims(y = c(-50,100)) +
  labs(y = "Difference in median yield prediction (%)") +
  coord_flip()

p2 <- tab.preds.sensi.bra.cc %>% 
  mutate(model     = factor(model, levels = rev(c("pca.m.3",  "pca.m.2",  "pca.d.3", "pca.d.all", 
                                              "fpca.m.3", "fpca.m.2", "pls.m.2", "pls.m.all", 
                                              "pls.d.3",  "pls.d.all", "avg.m")))) %>%
  group_by(model, gpe_model, dT) %>% 
  mutate(M=round(mean(Difference_Prediction_percentage), 2)) %>% 
  ggplot(., aes(x=model)) + 
  geom_hline(yintercept=0, color="black", lty=2) +
  geom_boxplot(aes(y=Difference_Prediction_percentage, color=gpe_model), outlier.size = 0.4) + 
  #geom_point(aes(y = M), color = "red") + 
  geom_text(aes(y = 90, label = M), size=3, check_overlap = T) + 
  facet_grid(gpe_model ~ dT, switch = "y") + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle=0), 
        strip.placement = "outside", 
        axis.title.y = element_blank(),
        #plot.margin = unit(c(0,0,0,0), units = "cm"),
        strip.background.y = element_blank()) + 
  ggtitle("b. Brazil") + 
  lims(y = c(-50,100)) +
  labs(y = "Difference in median yield prediction (%)") +
  coord_flip()

boxp.CC <- plot_grid(p1,p2,
          nrow=2, align = "hv", axis="tblr")

ggsave(boxp.CC, filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/USA_BRA_boxlot.png"), 
       width = 10.5, height=8, dpi=300, bg="white")







# test 
plot_grid(
  tab.preds.sensi.usa.cc %>% 
    filter(Dim_red != "pls.m.2") %>% 
    ggplot(.) + 
    geom_density(aes(x = Median_Prediction, color = Dim_red)) +
    geom_density(aes(x=Ref), color="black", lty=2) +
    geom_density(data = tab_usa %>% filter(country_name!="Desert"),
                 aes(x=Ya), color = "red") +
    facet_grid(Model~dT, switch = "y") + 
    scale_color_viridis_d(name = "Dimension reduction technique") + 
    theme_bw() + 
    theme(legend.position = "bottom",
          strip.text.y.left = element_text(angle=0), 
          strip.placement = "outside", 
          axis.title = element_blank(),
          strip.background.y = element_blank(), 
          #plot.margin = unit(c(0,0,0,0), units = "cm"),
          #plot.caption.position = "plot", 
          plot.caption = element_text(hjust=0, size=6)) + 
    ggtitle("a. US") + labs(x="Predicted yield")+ lims(x=c(0,6)),
  
  tab.preds.sensi.bra.cc %>% 
    filter(Dim_red != "pls.m.all")%>% 
    ggplot(.) + 
    geom_density(aes(x = Median_Prediction, color = Dim_red)) +
    geom_density(aes(x=Ref), color="black", lty=2) +
    geom_density(data = tab_bra %>% filter(country_name!="Desert"),
                 aes(x=Ya), color = "red") +
    facet_grid(Model~dT, switch = "y") + 
    scale_color_viridis_d(name = "Dimension reduction technique") + 
    theme_bw() + 
    theme(legend.position = "bottom",
          strip.text.y.left = element_text(angle=0), 
          strip.placement = "outside", 
          axis.title.y = element_blank(),
          strip.background.y = element_blank(), 
          #plot.margin = unit(c(0,0,0,0), units = "cm"),
          #plot.caption.position = "plot", 
          plot.caption = element_text(hjust=0, size=6)) + 
    ggtitle("b. Brazil") + labs(x="Predicted yield") + lims(x=c(0,6)),
  
  ncol=1, rel_heights = c(0.47, 0.53)
)

ggsave(filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/USA_BRA_density.png"), 
       width = 9, height=8.5, dpi=300, bg="white")
