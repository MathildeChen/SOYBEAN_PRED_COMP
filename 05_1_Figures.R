# -------------------------------------------------------------------------
# 
#       FIGURES FOR THE ARTICLE
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
pal <- wes_palette("Zissou1", 6, type="continuous")

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
                               "vapor_pressure_deficit"     ="vpd"))  %>% 
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

# > Full datasets 
load(paste0(data_path, "tab_world.rda"))
load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))

# > Load cumulated explained variables (for figure 3)
load(paste0(data_path, "tab_var_dim_red_full.rda"))


# > Load table with models performances (for figures 4 and 5)
load(paste0(data_path, "tab_perf_models.rda"))

# > Load table containing all predictions (for figure 5)
load(paste0(data_path, "tab_preds.rda"))

# > Load table containing predictions under increasing temperatures scenarios
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/05_sensi_cc/tab_preds_cc.rda")

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
# Table 1 

# > Without Desert
tab_world %>%
  filter(country_name != "Desert") %>% 
  mutate(Ya_ano = if_else(is.na(Ya_ano) == T, 0, Ya_ano)) %>% 
  dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% 
  mutate(analysis_set = "WORLD") %>% 
  rbind(., 
        tab_usa %>% filter(country_name != "Desert") %>% dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% mutate(analysis_set = "USA"),
        tab_bra %>% filter(country_name != "Desert") %>% dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% mutate(analysis_set = "BRA")) %>% 
  gather(key = "var", value="value", -site_year, -analysis_set) %>% 
  mutate(analysis_set = factor(analysis_set, levels = c("WORLD", "USA", "BRA"))) %>% 
  group_by(analysis_set, var) %>% 
  summarise(n = n(),
            med = round(median(value), 2), 
            q5th = round(quantile(value, probs = c(0.05)), 2),
            q95th = round(quantile(value, probs = c(0.95)), 2),
            mean = round(mean(value), 2), 
            sd = round(sd(value), 2), 
            min = round(min(value), 2), 
            max = round(max(value), 2)) %>% 
  mutate(lab1 = paste0(med, " (", q5th, ", ", q95th, ")"),
         lab2 = paste0(mean, " (", sd, "), [", min, " - ", max, "]")) %>% 
  dplyr::select(analysis_set, var, lab2) %>% 
  tidyr::spread(key=analysis_set, value =lab2)

#  var                              WORLD                         USA                          BRA          
#1 Ya                               2.55 (1.02), [0.07 - 6.7]     3.38 (0.7), [0.76 - 5.49]    2.99 (0.62), [1.11 - 5.84]    
#2 Ya_ano                           0 (0.21), [-1.62 - 2.2]       0 (0.21), [-1.3 - 0.68]      0 (0.22), [-1.28 - 2.2]       
#4 year_max_2m_temperature          25.37 (3.64), [14.26 - 34.85] 24.82 (2.9), [17 - 33.54]    28.27 (1.82), [22.64 - 33.11] 
#5 year_min_2m_temperature          15.69 (4.21), [3.34 - 25.19]  14.08 (3.02), [5.93 - 22.56] 19.14 (2.1), [13.13 - 23.8]   
#7 year_total_precipitation         0.17 (0.07), [0.02 - 0.57]    0.13 (0.03), [0.02 - 0.29]   0.23 (0.05), [0.08 - 0.44]    
#6 year_surface_net_solar_radiation 0.63 (0.06), [0.41 - 0.85]    0.66 (0.04), [0.55 - 0.85]   0.66 (0.03), [0.51 - 0.76]    
#3 year_et0                         1.73 (0.44), [0.69 - 5.07]    2.01 (0.34), [1.24 - 5.07]   1.36 (0.25), [0.81 - 2.62]    
#8 year_vapor_pressure_deficit      0.74 (0.22), [0.24 - 2.34]    0.75 (0.18), [0.34 - 2.34]   0.63 (0.15), [0.28 - 1.41] 


# > With Desert
tab_world %>%
  mutate(Ya_ano = if_else(is.na(Ya_ano) == T, 0, Ya_ano)) %>% 
  dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% 
  mutate(analysis_set = "WORLD") %>% 
  rbind(., 
        tab_usa %>% dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% mutate(analysis_set = "USA"),
        tab_bra %>% dplyr::select(site_year, Ya, Ya_ano, starts_with("year_")) %>% mutate(analysis_set = "BRA")) %>% 
  gather(key = "var", value="value", -site_year, -analysis_set) %>% 
  mutate(analysis_set = factor(analysis_set, levels = c("WORLD", "USA", "BRA"))) %>% 
  group_by(analysis_set, var) %>% 
  summarise(med = round(median(value), 2), 
            q5th = round(quantile(value, probs = c(0.05)), 2),
            q95th = round(quantile(value, probs = c(0.95)), 2),
            mean = round(mean(value), 2), 
            sd = round(sd(value), 2), 
            min = round(min(value), 2), 
            max = round(max(value), 2)) %>% 
  mutate(lab1 = paste0(med, " (", q5th, ", ", q95th, ")"),
         lab2 = paste0(mean, " (", sd, "), [", min, " - ", max, "]")) %>% 
  dplyr::select(analysis_set, var, lab2) %>% 
  tidyr::spread(key=analysis_set, value =lab2)

#  var                              WORLD                           USA                            BRA
#1 Ya                               2.05 (1.37), [0 - 6.7]          2.71 (1.49), [0 - 5.49]        2.4 (1.32), [0 - 5.84]        
#2 Ya_ano                           0 (0.19), [-1.62 - 2.2]         0 (0.19), [-1.3 - 0.68]        0 (0.2), [-1.28 - 2.2]        
#4 year_max_2m_temperature          22.16 (10.15), [-20.81 - 42.42] 21.96 (9.61), [-20.48 - 40.75] 24.7 (10.45), [-20.37 - 40.75]
#5 year_min_2m_temperature          12.68 (9.21), [-27.6 - 28.41]   11.58 (8.25), [-27.46 - 28.14] 15.6 (9.51), [-26.75 - 28.14] 
#7 year_total_precipitation         0.14 (0.08), [0 - 0.72]         0.11 (0.05), [0 - 0.29]        0.2 (0.09), [0 - 0.44]        
#6 year_surface_net_solar_radiation 0.59 (0.14), [0.1 - 1]          0.62 (0.13), [0.11 - 0.89]     0.62 (0.13), [0.1 - 0.89]     
#3 year_et0                         1.73 (0.86), [0.19 - 8.17]      1.96 (0.84), [0.2 - 7.76]      1.45 (0.85), [0.2 - 7.21]     
#8 year_vapor_pressure_deficit      0.75 (0.53), [0.02 - 5.04]      0.78 (0.56), [0.02 - 4.85]     0.68 (0.57), [0.02 - 4.85]  

tab_world %>%
  filter(country_name!="Desert") %>% 
  mutate(irrigated_portion = irrigated_portion/100) %>% 
  pull(irrigated_portion) %>% 
  summary(.)

tab_world %>%
  filter(country_name!="Desert", irrigated_portion > 0) %>% 
  pull(irrigated_portion) %>% 
  summary(.)

# ----------------------------------
# Figure 1 

# schema of analyses framework 

# ----------------------------------
# Figure 2 

tab_world %>% 
  ungroup() %>%
  distinct(x, y, country_name) %>% 
  group_by(country_name) %>% 
  count()

# Map unique combination of site-years for WORLD (main analysis)
# Full data 
Fig2 <- tab_world %>% 
  ungroup() %>%
  distinct(x, y, country_name) %>% 
  mutate(country_name  = as.character(country_name)) %>% 
  mutate(country_color = case_when(
    country_name == "Desert"                  ~ "Desert",
    TRUE ~ "World (N=3447)")) %>% 
  mutate(is_desert     = if_else(country_name == "Desert", 
                                 "Selected sites in geographical areas unsuitable for soybean production (yield = 0 tons/hectares) (N=663)", 
                                 "Selected sites in major soybean producing areas (N=2784)")) %>%
  mutate(is_desert = factor(is_desert, levels = c("Selected sites in major soybean producing areas (N=2784)",
                                                  "Selected sites in geographical areas unsuitable for soybean production (yield = 0 tons/hectares) (N=663)"))) %>%
  ggplot(.) + 
  geom_sf(data=world %>% filter(continent != "Antarctica")) + 
  geom_point(aes(x = x, y = y, color = as.factor(is_desert)), 
             size=0.15) + 
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom", 
        plot.margin = unit(c(0,0,0,0),units = "cm")) + 
  scale_color_manual(values = c("#4DAF4A", "black"), name = "") + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1)) ; Fig2

ggsave(Fig2, filename=paste0(save_path, "Figure_2.png"), 
       width = 9, height=4, dpi=300, bg="white")

# ----------------------------------
# Figure 3
# % of total variance explained in the data by each dimension reduction techniques

# function to remove useless pannels
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

Fig3 <- tab_var_full %>% 
  filter(country=="WORLD") %>% 
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
  labs(x = "Scores", y = "Cumulated explained variance (%)"#, caption = "Abbreviations: PCA: principal component analysis; FPCA: functional PCA; MFPCA: multivariate FPCA; PLSR: partial least square regression; FPLSR: functional PLSR."
       ) ; Fig3

remove <- c("aa#aa
             aa#aa
             aa#aa
             aa#aa
             aa#aa
             aa#aa
             ##a##")

# > Save 
ggsave(plot = remove_facets(Fig3, remove), 
       filename=paste0(save_path, "Figure_3.png"), 
       width = 10, height=11, dpi=300, bg="white")

# ------------------------------------------
# Figure 4 - Performance of all models to predict soybean yield or soybean yield globally 
#            Measured as NSE (RMSEP in supplementary figures)

# > Compute mean performance 
tab_perf_labelled_i_mean <- tab_perf_labelled %>% 
  group_by(Outcome, Outcome_lab, Country, Country_lab, Model, dimension_reduction, pred_perf, pred_perf_lab, gpe_model) %>% 
  # MEAN OF BOTH CROSS-VALIDATION
  summarise(mean_pred_perf = mean(pred_perf_value)) %>% 
  ungroup() %>%
  mutate(Type_cv_lab = "c. Mean of cross-validation results")

# > Plot 
#   Result based on each type of cross-validation
p.nse <- tab_perf_labelled %>% 
  filter(pred_perf == "NSE", 
         Country=="WORLD", 
         Outcome == "01_Ya") %>% 
  # > order of models 
  arrange(Type_cv_lab, Outcome_lab, Country_lab, desc(pred_perf_value)) %>% 
  group_by(Type_cv_lab, Outcome_lab, Country_lab) %>% 
  mutate(order = row_number()) %>%
  mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
         best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
  # > min and max among lm and rf
  group_by(Type_cv_lab, Outcome_lab, Country, Model) %>% 
  mutate(min_x_lab = min(pred_perf_value),
         max_x_lab = max(pred_perf_value)) %>% 
  # > group of models 
  mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
  # > change labels 
  mutate(Type_cv_lab = recode(Type_cv_lab, "Cross-validation on years"="a. Cross-validation on years", "Cross-validation on sites"="b. Cross-validation on sites")) %>% 
  # > plot
  ggplot(., aes(x = pred_perf_value, y=Model)) +
  geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                 color = "black") +
  geom_point(aes(color = gpe_model, fill = best, shape = interaction(negative_pred_perf_value, gpe_model)),
  #geom_point(aes(color = gpe_model, fill = best, shape = gpe_model),
             size = 2) + 
  geom_text(aes(label = order), 
            nudge_y = 0.4, 
            size = 2) + 
  facet_grid(dimension_reduction ~ Type_cv_lab, scales = "free_y", space = "free", switch = "y") + 
  theme_bw() +
  theme(legend.position = 'bottom', 
        axis.title.y = element_blank(),
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0, size=9),
        strip.text.x = element_text(size=10), 
        strip.background = element_rect(fill="grey", color = "transparent"),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10),
        #panel.grid.major.x = element_line(color = "grey", linewidth = 0.01),
        #panel.grid.minor.x = element_line(color = "grey", linewidth = 0.01),
        panel.grid = element_blank(),
        panel.spacing.x = unit(1.25, "lines")) + 
  scale_shape_manual(values = c(21, 4)) + 
  scale_color_manual(values = pal_country[c(3, 4)], name = "Model family") + 
  scale_fill_manual(values = c("red", "white"), name = "Rank model") + 
  scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2), limits = c(0,1)) + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2, shape=21), ncol = 1), 
         fill = guide_legend(order = 1, override.aes = list(size=2, shape=21), ncol = 1),
         shape = guide_none())   + 
  labs(x = "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)") ; p.nse

#    Results based on mean of results of each CV procedure
p.nse_mean <- tab_perf_labelled_i_mean %>% 
  filter(pred_perf == "NSE", 
         Country=="WORLD", 
         Outcome == "01_Ya") %>% 
  # > order of models 
  arrange(Outcome_lab, Country_lab, desc(mean_pred_perf)) %>% 
  group_by(Outcome_lab, Country_lab) %>% 
  mutate(order = row_number()) %>%
  mutate(best = if_else(order==1, "Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"),
         best = factor(best, levels=c("Model with Rank = 1 (best predictive accuracy)", "Model with Rank > 1"))) %>%
  # > min and max among lm and rf
  group_by(Outcome_lab, Country_lab, Model) %>% 
  mutate(min_x_lab = min(mean_pred_perf),
         max_x_lab = max(mean_pred_perf)) %>% 
  # > plot
  ggplot(., aes(x = mean_pred_perf, y=Model)) +
  geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                 color = "black") +
  geom_point(aes(color = gpe_model, fill = best), 
             shape = 21,
             size = 2) + 
  geom_text(aes(label = order), 
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
        strip.background.x = element_rect(fill="grey", color = "transparent"),
        strip.background.y = element_blank(),
        title = element_text(size=10),
        legend.text = element_text(size=10),
        #panel.grid.major.x = element_line(color = "grey", linewidth = 0.01),
        #panel.grid.minor.x = element_line(color = "grey", linewidth = 0.01),
        panel.grid = element_blank()) + 
  scale_shape_manual(values = c(21, 4)) + 
  scale_color_manual(values = pal_country[c(3, 4)], name = "Model family") + 
  scale_fill_manual(values = c("red", "white"), name = "Rank model") + 
  scale_x_continuous(breaks = seq(0, 1.2, by=0.2), labels = seq(0, 1.2, by=0.2), limits = c(0,1)) + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2, shape=21), ncol = 1), 
         fill = guide_legend(order = 1, override.aes = list(size=2, shape=21), ncol = 1),
         shape = guide_none())   + 
  labs(x = "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)") ; p.nse_mean

# N predictor
p.n_pred <- tab_perf_labelled %>% 
  filter(pred_perf == "NSE") %>% 
  # > Nb of predictors per model
  group_by(Type_cv_lab, Outcome_lab, Country_lab, dimension_reduction, Model) %>%
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

# > Merge all plots and add the legend
leg <- get_legend(p.nse)

plot <- plot_grid(p.nse + theme(legend.position = "none", axis.title.x = element_blank()), 
          p.nse_mean + theme(legend.position = "none", axis.title.x = element_blank()), 
          nrow=1, rel_widths = c(2.75, 1), axis = "tblr", align = "h")

plot_x <- ggdraw(add_sub(plot, "Nash-Sutcliffe model efficiency\n(higher value indicates a better predictive performance)", 
               vpadding=grid::unit(0,"lines"),y=0, x=0.6, vjust=0, size = 10))

Fig4 <- plot_grid(plot_x, 
                  NULL,
                  leg, 
                  ncol = 1, rel_heights = c(0.95,0.05,0.05))

# > Save 
ggsave(plot = Fig4, 
       filename=paste0(save_path, "Figure_4.png"), 
       width = 9, height=11, dpi=300, bg="white")


# ------------------------------------------
# Figure 5 
# Observation vs Prediction & Residues 
# for the best models for each country (USA, Brazil, global)

# > Best models 
tab_perf_labelled %>% 
  filter(pred_perf == "RMSEP", Outcome == "01_Ya") %>% 
  # > order of models 
  arrange(Country, pred_perf_value) %>% 
  group_by(Country, Type_cv) %>% 
  mutate(order = row_number()) %>%
  filter(order == 1)

tab_perf_labelled %>% 
  filter(pred_perf == "NSE", Outcome == "01_Ya") %>% 
  # > order of models 
  arrange(Country, desc(pred_perf_value)) %>% 
  group_by(Country, Type_cv) %>% 
  mutate(order = row_number()) %>%
  filter(order == 1)

# > Select predictions of best models 
best_models_preds <- tab_preds %>% 
  # > best model for each dataset
  mutate(best = case_when(
    Country == "USA"   & Model == "mfpca.m.1" ~ 1,
    Country == "BRA"   & Model == "avg.a" ~ 1,
    Country == "WORLD" & Model == "pca.m.2" ~ 1,
    TRUE ~ 0
  )) %>% 
  filter(best == 1) %>% 
  # > compute residuals
  mutate(Residuals = Predicted - Observed) %>% 
  # > compute RMSEP and NSE
  group_by(Country, Model) %>% 
  mutate(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
  ) %>%
  split(.$Country) 

# > Density of points 
best_models_preds$USA$density <- get_density(x = best_models_preds$USA$Observed, y = best_models_preds$USA$Predicted, n=100)
best_models_preds$BRA$density <- get_density(x = best_models_preds$BRA$Observed, y = best_models_preds$BRA$Predicted, n=100)
best_models_preds$WORLD$density <- get_density(x = best_models_preds$WORLD$Observed, y = best_models_preds$WORLD$Predicted, n=100)

tab_best_models_preds <- plyr::ldply(best_models_preds, data.frame, .id = "Country")

# > Pred vs Obs
p.pred.obs <- tab_best_models_preds %>% 
  mutate(RMSEP = round(RMSEP, 3),
         NSE = round(NSE, 3)) %>% 
  mutate(label_perf = paste0("Name model: ", Model, "\nRMSEP: ", RMSEP, " tonnes/hectare\nNSE: ", NSE)) %>% 
  # country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) %>% 
  #filter(Site_year %in% c(unique(best_models_preds$USA$Site_year)[1:1000], unique(best_models_preds$BRA$Site_year)[1:1000], unique(best_models_preds$WORLD$Site_year)[1:1000])) %>% 
  ggplot(data = .) + 
  geom_point(aes(x = Observed, y = Predicted, color = density), 
             size = 0.5) + 
  geom_abline(color = "red") + 
  geom_text(aes(label = label_perf),
            x = 0, y = 6, check_overlap = TRUE, hjust = 0, size = 3) + 
  theme_cowplot() + 
  theme(legend.position = "none",
        strip.text.x = element_text(size=10), 
        strip.background.x = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10)
        #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
        ) + 
  facet_wrap(.~Country_lab, scales = "free_y") +
  #scale_color_gradientn(colours = rev(pal)) + 
  scale_color_viridis_c(option = 3) +
  labs(x = "Observed yields (tons/hectare)", y = "Predicted yields (tons/hectare)") + 
  ggtitle("a. ") + 
  lims(x = c(0,7), y = c(0,7)) ; p.pred.obs
  
# > Residuals
p.res <- tab_best_models_preds %>% 
  # country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) %>% 
  # mean residuals 
  group_by(Country) %>% 
  mutate(mean_Residuals = mean(Residuals)) %>% 
  ggplot(.) + 
  geom_histogram(aes(x=Residuals, fill = Country), 
               bins = 100, alpha = 0.5, color = "transparent") + 
  geom_vline(aes(xintercept = mean_Residuals), color = "black", 
             linetype = 2, linewidth = 0.5) +
  theme_cowplot() + 
  theme(legend.position = "none",
        #strip.text = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10)) + 
  facet_wrap(.~Country_lab, scales = "free_y") + 
  scale_fill_manual(values = c("#377EB8","#E41A1C",  "black")) + 
  scale_color_manual(values = c( "#377EB8", "#E41A1C","black")) + 
  labs(x = "Predicted - observed yields (tons/hectare)", y = "Counts") +
  ggtitle("b. ") ; p.res

p6 <- plot_grid(p.pred.obs, p.res, 
                ncol = 1, align = "hv")

ggsave(plot = p6, 
       filename=paste0(save_path, "Figure_6.png"), 
       width = 13, height=9, dpi=300, bg="white")

# ------------------------------------------
# Calibration plot + residuals distribution for each model 

# > Select data 
tab_preds_i <- tab_preds %>% 
  filter(Country=="WORLD",
         Outcome == "01_Ya",
         Model %in% c("pca.m.2", "rf_pca.m.2")) %>% 
    # > compute residuals
    mutate(Residuals = Predicted - Observed) %>% 
    # > compute RMSEP and NSE
    group_by(Type_cv, Country, Model) %>% 
    mutate(
      RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
      NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
    ) %>% 
  mutate(RMSEP = round(RMSEP, 3),
         NSE = round(NSE, 3)) %>% 
  mutate(Model = "random forest model including\nthe two first scores derived from\nPCA on monthly data (pca.m.2)") %>% 
  mutate(label_perf = paste0("Name model: ", Model, "\nRMSEP: ", RMSEP, " tonnes/hectare\nNSE: ", NSE)) %>% 
  # > country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) 

# > Density of points 
tab_preds_i_dens <- tab_preds_i %>% 
  split(.$Type_cv) %>% 
  map_dfr(., ~{
    
    density <- get_density(x = .x$Observed, y = .x$Predicted, n=100)
    .x %>% 
      mutate(density=density)
  }, .id = "Type_cv") %>% 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) 
  
# > Pred vs Obs (calibration plot)
p.pred.obs <- tab_preds_i_dens %>% 
  ggplot(data = .) + 
  geom_point(aes(x = Observed, y = Predicted, color = density), 
             size = 0.5) + 
  geom_abline(color = "red") + 
  #geom_text(aes(label = label_perf), x = 0, y = 6, check_overlap = TRUE, hjust = 0, size = 3) + 
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
  scale_color_viridis_c(option = 3, name="Density per combination of\nobserved/predicted yield values", 
                        guide = guide_colorbar(barwidth = 7.5, barheight = 0.5)) +
  facet_wrap(Type_cv_lab ~ ., ncol=1, scales = "free") +
  labs(x = expression(paste("Observed yields (t  ", h^{-1}, ")")), y = expression(paste("Predicted yields (t  ", h^{-1}, ")"))) + 
  ggtitle("a. ") + 
  lims(x = c(0,7), y = c(0,7))

# > Density of observation vs prediction 
p.dens.pred.obs <- tab_preds_i_dens %>% 
  dplyr::select(-density, -Residuals, -RMSEP, -NSE) %>% 
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
  scale_color_viridis_d(name ="", option = "C") +
  scale_fill_viridis_d(name ="", option = "C") +
  labs(x = expression(paste("Yield (t  ", h^{-1}, ")")), y = "Density") + 
  ggtitle("b. ")+
  lims(x = c(0,7)); p.dens.pred.obs

# > Residuals
p.res <- tab_preds_i_dens %>% 
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
  labs(x = expression(paste("Predicted - observed yields (t  ", h^{-1}, ")")), y = "Counts") ; p.res
  
Fig5 <- plot_grid(p.pred.obs + theme(strip.text.y = element_blank(),
                                     strip.background = element_blank()), 
                  ggplot() + theme_void(),
                  p.dens.pred.obs + theme(strip.text.y = element_blank(),
                                          strip.background = element_blank()), 
                  ggplot() + theme_void(),
                  p.res + theme(strip.text.y = element_blank(),
                                strip.background = element_blank()), 
                  nrow = 1, align = "h", axis="tblr", rel_widths = c(1, 0.1, 1, 0.1, 1))
# > save
ggsave(plot = Fig5, 
       filename=paste0(save_path, "Figure_5.png"), 
       width = 13, height=7, dpi=300, bg="white")
  
# ------------------------------------------
# Figure 7 & 8  
# Maps of prediction in different climate change scenarios 

# CLIMATE CHANGE 
tab.preds.sensi.usa.cc <- preds.cc %>% filter(country=="01_USA")
tab.preds.sensi.bra.cc <- preds.cc %>% filter(country=="02_BRA")

# Maps showing the median difference (only USA, not desert)
cc.usa.1 <- tab.preds.sensi.usa.cc %>% 
  filter(gpe_model=="Random forest") %>%
  #filter(model %in% c("pca.m.3", "pca.d.3", "pls.m.2", "pls.d.3")) %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-110, -57), y=c(25,50)) + 
  facet_grid(selection_step + model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "grey"),
        strip.text.y.left = element_text(angle=0, size=8), 
        strip.text.x = element_text(size=8), 
        strip.placement = "outside",
        legend.position = "bottom", 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage), max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)     ") +
  ggtitle("United States", subtitle = "a. Random forest") 

cc.usa.2 <- tab.preds.sensi.usa.cc %>% 
  filter(gpe_model!="Random forest") %>%
  #filter(model %in% c("pca.m.3", "pca.d.3", "pls.m.2", "pls.d.3")) %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-110, -57), y=c(25,50)) + 
  facet_grid(selection_step + model ~ dT) +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "grey"),
        #strip.text.y.left = element_text(angle=0, size = 8), 
        strip.text.y = element_blank(), 
        strip.background.y = element_blank(),
        strip.text.x = element_text(size=8), 
        strip.placement = "outside",
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.usa.cc$Difference_Prediction_percentage), max(tab.preds.sensi.usa.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)     ") +
  ggtitle("", subtitle = "b. Multiple linear regression") 

leg_cc.usa.1 <- get_legend(cc.usa.1)

cc.usa <- plot_grid(
  plot_grid(cc.usa.1 + 
              theme(legend.position = "none"), 
            cc.usa.2 + 
              theme(legend.position = "none"),
            ncol=2, align = "hv", axis = "tblr"),
  leg_cc.usa.1, 
  ncol=1, align = "hv", axis = "tblr", rel_heights = c(0.95, 0.05)
)
                    #, rel_widths = c(0.55, 0.45)) 

# > save
#ggsave(usa.cc, filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/USA_cc2.png"), 
#       width = 14, height=7, dpi=300, bg="white")
#

ggsave(cc.usa, filename="E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/USA_cc2.png", 
       width = 14, height=6.5, dpi=300, bg="white")

# Maps showing the median difference (only brazil, not desert)
cc.bra.1 <- tab.preds.sensi.bra.cc %>% 
  filter(gpe_model=="Random forest") %>%
  #filter(model %in% c("pca.m.2", "pca.d.all",  "pls.m.all", "pls.d.all")) %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-85, -35), y=c(-35,5))  + 
  facet_grid(selection_step + model ~ dT, switch = "y") +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "grey"),
        strip.text.y.left = element_text(angle=0, size=8), 
        strip.text.x = element_text(size=8), 
        strip.placement = "outside",
        legend.position = "bottom") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage), max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)  ") +
  ggtitle("Brazil", "a. Random forest") 

cc.bra.2 <- tab.preds.sensi.bra.cc %>% 
  filter(gpe_model!="Random forest") %>%
  #filter(model %in% c("pca.m.2", "pca.d.all",  "pls.m.all", "pls.d.all")) %>%
  ggplot(.) + 
  geom_sf(data = world) + 
  geom_point(aes(x =  x, y = y, color = Difference_Prediction_percentage),
             size = 0.1) + 
  lims(x=c(-85, -35), y=c(-35,5))  + 
  facet_grid(selection_step + model ~ dT) +
  theme_map() + 
  theme(strip.background.x = element_rect(colour = "grey"),
        #strip.text.y.left = element_text(angle=0, size = 8), 
        strip.text.y = element_blank(), 
        strip.background.y = element_blank(),
        strip.text.x = element_text(size=8), 
        strip.placement = "outside",
        legend.position = "bottom") + 
  scale_color_gradientn(colors = c("purple4","royalblue", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[1:4], "white", wesanderson::wes_palette("Zissou1", n = 10, type = "continuous")[7:10], "darkred", "black"),
                        values = c(1, (0 - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage))/(max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage) - min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 0), 
                        guide = guide_colorbar(barwidth = 15, barheight = 0.5),
                        breaks = seq(-80, 80, by=20), 
                        limits = c(min(tab.preds.sensi.bra.cc$Difference_Prediction_percentage), max(tab.preds.sensi.bra.cc$Difference_Prediction_percentage)), 
                        name = "Difference in median yield prediction (%)  ") +
  ggtitle("", "b. Multiple linear regression")

leg_cc.bra.1 <- get_legend(cc.bra.1)

cc.bra <- plot_grid(
  plot_grid(cc.bra.1 + 
              theme(legend.position = "none"), 
            cc.bra.2 + 
              theme(legend.position = "none"),
            ncol=2, align = "hv", axis = "tblr"),
  leg_cc.bra.1, 
  ncol=1, align = "hv", axis = "tblr", rel_heights = c(0.95, 0.05)
)

ggsave(cc.bra, filename="E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/BRA_cc2.png", 
       width = 14, height=6.5, dpi=300, bg="white")


# > save

cc.plot <- plot_grid(cc.usa, cc.bra, 
                    ncol=1, align = "hv", axis = "tblr")

ggsave(cc.usa, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/USA_cc3.png"), 
       width = 13, height=6, dpi=300, bg="white")
ggsave(cc.bra, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/BRA_cc3.png"), 
       width = 14, height=7, dpi=300, bg="white")

ggsave(cc.plot, 
       filename=paste0(save_path, "SUPPLEMENTARY_FIGURES/11_CC_SIMULATION/cc.plot.png"), 
       width = 14, height=10, dpi=300, bg="white")


stop()

# V2 for Figure 3
# > Plot
clim_data_plot <- clim_data %>% 
  mutate(month_day_1 = if_else(day_of_year %in% c(1, 31, 62, 92, 123, 153, 184), day_of_year, NA)) %>% 
  #filter(clim.var %in% c("max_2m_temperature", "min_2m_temperature", "total_precipitation")) %>% 
  mutate(Country = recode(Country, "United-States of America (N=29803)"="United-States of America\n(N=29803)",
                          "Brazil (N=14757)"="Brazil\n(N=14757)", 
                          "Global (N=122229)"="Global\n(N=122229)")) %>% 
  mutate(Country = factor(Country, levels = c("United-States of America\n(N=29803)", "Brazil\n(N=14757)", "Global\n(N=122229)"))) %>% 
  left_join(., vars_names, by="clim.var") %>% 
  # Plot 
  ggplot(data = .) + 
  # takes many time to be produced
  #geom_line(aes(x = day_of_year, y = cum_clim.value, group = site_year), 
  #          color = "grey") + 
  geom_ribbon(aes(x= day_of_year, ymin = perc05, ymax = perc95, fill = Country),
              alpha = 0.2) +
  geom_line(aes(x = day_of_year, y = median, color = Country)) +
  geom_point(aes(x = month_day_1, y = median, color = Country), size = 0.75) + 
  lemon::facet_rep_grid(clim.var_lab~Country, scale="free", switch = "y")  +
  theme_cowplot() + 
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none", 
        legend.text = element_text(size=10),
        strip.background.x = element_blank(),
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0, size=10),
        strip.text.x = element_text(size=10), 
        title = element_text(size=10), 
        panel.grid.major = element_line(color = "grey", linewidth = 0.1)) + 
  scale_color_manual(values = c("#E41A1C", "#377EB8", "black")) + 
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "black")) + 
  labs(x = "Day of the growing season", y = paste0("Accumulated climatic data\nover the growing season")) ; clim_data_plot

ggsave(clim_data_plot, filename=paste0(save_path, "Figure_3_2.png"), 
       width = 10, height=12, dpi=300, bg="white")

# -----------------------------
# -----------------------------
# Figures asked after Revision 1


# Predicted vs Observed yields for pca.m.2, avg.s, and avg.m 
# > Select data 
tab_preds_i <- tab_preds %>% 
  filter(Country=="WORLD",
         Outcome == "01_Ya",
         #Model %in% c("rf_pca.m.2", "rf_avg.a", "rf_avg.m")) %>% 
         Model %in% c("rf_pca.m.2", "lm_pca.m.2", 
                      "rf_pca.d.2", "lm_pca.d.2")) %>% 
  # > compute residuals
  mutate(Residuals = Predicted - Observed) %>% 
  # > compute RMSEP and NSE
  group_by(Type_cv, Country, Model) %>% 
  mutate(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
  ) %>% 
  mutate(RMSEP = round(RMSEP, 3),
         NSE = round(NSE, 3)) %>% 
  #mutate(Model = "random forest model including\nthe two first scores derived from\nPCA on monthly data (pca.m.2)") %>% 
  mutate(label_perf = paste0("Name model: ", Model, "\nRMSEP: ", RMSEP, " tonnes/hectare\nNSE: ", NSE)) %>% 
  # > country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) 

# > Density of points 
tab_preds_i_dens <- tab_preds_i %>% 
  split(.$Type_cv) %>% 
  map_dfr(., ~{
    
    .x %>% 
      split(.$Model) %>% 
      map_dfr(., ~{
        
        density <- get_density(x = .x$Observed, y = .x$Predicted, n=100)
        .x %>% 
          mutate(density=density)
        
        
      }, .id = "Model")
    
  }, .id="Type_cv") %>% 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) 

# > correlation between obs and preds for each model
tab_preds_i_corr <- tab_preds_i %>% 
  split(.$Type_cv) %>% 
  map_dfr(., ~{
    
    .x %>% 
      split(.$Model) %>% 
      map_dfr(., ~{
        
        cor(.x$Observed, .x$Predicted)
        
        
      }, .id = "Model")
    
  }, .id="Type_cv") %>% 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) 


# > Pred vs Obs (calibration plot)
p.pred.obs <- tab_preds_i_dens %>% 
  #mutate(Model = recode(Model, "rf_avg.m" = "avg.m", "rf_avg.a" = "avg.s", "rf_pca.m.2"="pca.m.2"), Model = factor(Model, levels = c("pca.m.2", "avg.m", "avg.s"))) %>% 
  mutate(Model = recode(Model, 
                        "rf_pca.m.2" = "pca.m.2\n(random forest)", 
                        "lm_pca.m.2" = "pca.m.2\n(linear regression)", 
                        "rf_pca.d.2" = "pca.d.2\n(random forest)",
                        "lm_pca.d.2" = "pca.d.2\n(linear regression)"), 
         Model = factor(Model, levels = c("pca.m.2\n(random forest)", 
                                          "pca.m.2\n(linear regression)", 
                                          "pca.d.2\n(random forest)", 
                                          "pca.d.2\n(linear regression)"))) %>% 
  ggplot(data = .) + 
  geom_point(aes(x = Observed, y = Predicted, color = density), 
             size = 0.5) + 
  geom_hline(yintercept = 0, color="grey") +
  geom_abline(color = "red") + 
  #geom_text(aes(label = label_perf), x = 0, y = 6, check_overlap = TRUE, hjust = 0, size = 3) + 
  theme_cowplot() + 
  theme(legend.position = "bottom",
        strip.text.x = element_text(size=10), 
        strip.background.x = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10)
        #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
  ) + 
  #scale_color_viridis_c(option = 3, name="Density per combination of observed/predicted yield values", guide = guide_colorbar(barwidth = 17, barheight = 0.5), direction = -1) +
  scale_color_gradientn(colours = c(wes_palette("Zissou1", 10, "continuous"), "darkred", "black"), name="Density per combination of observed/predicted yield values", guide = guide_colorbar(barwidth = 17, barheight = 0.5)) +
  facet_wrap(Type_cv_lab ~ Model, scales = "free", nrow=2) +
  labs(x = expression(paste("Observed yields (t  ", h^{-1}, ")")), y = expression(paste("Predicted yields (t  ", h^{-1}, ")"))) +
  lims(x = c(0,7), y = c(-5,7))
  #ggtitle("a. ") + 

# > save
ggsave(plot = p.pred.obs, 
       filename="E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_pred_vs_obs.png", 
       width = 13, height=8, dpi=300, bg="white")

# > Density of observation vs prediction 
# split by model
p.dens.pred.obs1 <- tab_preds_i_dens %>% 
  #mutate(Model = recode(Model, "rf_avg.m" = "avg.m", "rf_avg.a" = "avg.s", "rf_pca.m.2"="pca.m.2"),Model = factor(Model, levels = c("pca.m.2", "avg.m", "avg.s"))) %>% 
  mutate(Model = recode(Model, 
                        "rf_pca.m.2" = "pca.m.2\n(random forest)", 
                        "lm_pca.m.2" = "pca.m.2\n(linear regression)", 
                        "rf_pca.d.2" = "pca.d.2\n(random forest)",
                        "lm_pca.d.2" = "pca.d.2\n(linear regression)"), 
         Model = factor(Model, levels = c("pca.m.2\n(random forest)", 
                                          "pca.m.2\n(linear regression)", 
                                          "pca.d.2\n(random forest)", 
                                          "pca.d.2\n(linear regression)"))) %>%
  #dplyr::select(-density, -Residuals, -RMSEP, -NSE) %>% 
  #gather(key="Type", value = "Yield", Predicted, Observed) %>% 
  ggplot(.) + 
  geom_vline(xintercept = 0, color = "grey") +
  geom_density(aes(x = Predicted, fill = "grey"), alpha=0.1) + 
  geom_density(aes(x = Observed), color = "darkgreen") + 
  theme_cowplot() + 
  theme(legend.position = "none",
        strip.text.x = element_text(size=10), 
        #strip.text.y = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10)
        #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
  ) + 
  facet_wrap(Type_cv_lab ~ Model, scales = "free", nrow=2) +
  #facet_wrap(.~Country_lab, scales = "free_y") +
  #scale_color_viridis_d(name ="", option = "D") +
  #scale_color_manual(colors = , name ="", option = "D") +
  #scale_fill_viridis_d(name ="", option = "D") +
  scale_fill_viridis_d(name ="", option = "D") +
  labs(x = expression(paste("Yield (t  ", h^{-1}, ")")), y = "Density") +
  #ggtitle("a. Split by model") +
  lims(x = c(-0,7), y=c(0,1)); p.dens.pred.obs1

# with all models
p.dens.pred.obs2 <- tab_preds_i_dens %>% 
  mutate(Model = recode(Model, "rf_avg.m" = "avg.m", "rf_avg.a" = "avg.s", "rf_pca.m.2"="pca.m.2"),
         Model = factor(Model, levels = c("pca.m.2", "avg.m", "avg.s"))) %>% 
  mutate(Lab = " ") %>%
  #dplyr::select(-density, -Residuals, -RMSEP, -NSE) %>% 
  #gather(key="Type", value = "Yield", Predicted, Observed) %>% 
  ggplot(.) + 
  geom_density(aes(x = Predicted, color = Model), alpha=0.1) + 
  geom_density(data = tab_preds_i_dens %>% filter(Model == "rf_avg.m"), 
               aes(x = Observed), color = "darkgreen", lty=2,linewidth=1) + 
  theme_cowplot() + 
  theme(legend.position = "none",
        strip.text.x = element_text(size=10), 
        #strip.text.y = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=10)
        #, panel.grid.major = element_line(color = "lightgrey", linewidth = 0.1)
  ) + 
  facet_wrap(Type_cv_lab ~ Lab, scales = "free", ncol=1) +
  #facet_wrap(.~Country_lab, scales = "free_y") +
  scale_color_viridis_d(name ="", option = "D") +
  labs(x = expression(paste("Yield (t  ", h^{-1}, ")")), y = "Density") +
  ggtitle("b. All models") +
  lims(x = c(0,7), y=c(0,1)); p.dens.pred.obs2

p.dens.pred.obs <- plot_grid(p.dens.pred.obs1, p.dens.pred.obs2, 
                             ncol=2, rel_widths = c(0.725,0.275),
                             align = "v", axis = "tbrl") ; p.dens.pred.obs

# > save
ggsave(plot = p.dens.pred.obs1, 
       filename="E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_pred_vs_obs_density.png", 
       width = 13, height=6, dpi=300, bg="white")

# > Residuals
p.res <- tab_preds_i_dens %>% 
  mutate(Model = recode(Model, "rf_avg.m" = "avg.m", "rf_avg.a" = "avg.s", "rf_pca.m.2"="pca.m.2"),
         Model = factor(Model, levels = c("pca.m.2", "avg.m", "avg.s"))) %>% 
  # mean residuals 
  group_by(Type_cv, Model) %>% 
  mutate(mean_Residuals = mean(Residuals)) %>% 
  ggplot(.) + 
  geom_histogram(aes(x = Residuals, fill = Model), 
                 bins = 100, alpha = 0.5, color = "transparent") + 
  geom_vline(aes(xintercept = mean_Residuals), color = "black", 
             linetype = 2, linewidth = 0.5) +
  theme_cowplot() + 
  theme(legend.position = "none",
        #strip.text = element_blank(), 
        strip.background = element_blank(),
        title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title.x = element_text(size=10)) + 
  facet_wrap(Type_cv_lab ~ Model, scales = "free") +
  scale_fill_viridis_d(name ="", option = "D") +
  #facet_wrap(.~Country_lab, scales = "free_y") + 
  #ggtitle("c. ") +
  lims(x=c(-5,5), y=c(0, 35000)) + 
  labs(x = expression(paste("Predicted - observed yields (t  ", h^{-1}, ")")), y = "Counts") ; p.res

# > save
ggsave(plot = p.res, 
       filename="E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/Supp_Figure_res.png", 
       width = 13, height=7, dpi=300, bg="white")

# Tentative figure rendement 
for(var_i in unique(vars_names$clim.var))
{
 
  # Dynamique variable name
  var_i_quo     <- quo_name(var_i) 
  abb_var_i_quo <- quo_name(vars_names[which(vars_names$clim.var==var_i),"clim.var_abb"])
  lab_vari_i    <- vars_names[which(vars_names$clim.var==var_i),"clim.var_lab"]
  lab2_vari_i    <- vars_names[which(vars_names$clim.var==var_i),"clim.var_lab2"]
  
  # Table with data
  tab_world_i <- tab_world %>%
    dplyr::select(-starts_with("PLS"), -starts_with("FPLS"),-starts_with("FPC"),-starts_with("MFP"),-starts_with("PCA_day"), -starts_with("monthly_cum")) %>% 
    rename(year_average  := !! paste0("year_", var_i_quo),
           PC1_var       := !! paste0("PC1_month_", abb_var_i_quo),
           PC2_var       := !! paste0("PC2_month_", abb_var_i_quo),
           PC3_var       := !! paste0("PC3_month_", abb_var_i_quo),
           monthly_var_1 := !! paste0("monthly_", var_i, "_1"),
           monthly_var_2 := !! paste0("monthly_", var_i, "_2"),
           monthly_var_3 := !! paste0("monthly_", var_i, "_3"),
           monthly_var_4 := !! paste0("monthly_", var_i, "_4"),
           monthly_var_5 := !! paste0("monthly_", var_i, "_5"),
           monthly_var_6 := !! paste0("monthly_", var_i, "_6"),
           monthly_var_7 := !! paste0("monthly_", var_i, "_7"))
  
  # Site selection
  # > 5 sites with lower seasonal values per country
  sites_min <- tab_world_i %>%
    group_by(country_name, gridcode) %>%
    summarise(mean_per_site = mean(year_average)) %>%
    group_by(country_name) %>% 
    slice_min(mean_per_site, n=1) 
  
  # > 5 sites with lower seasonal values per country
  sites_max <- tab_world_i %>%
    group_by(country_name, gridcode) %>%
    summarise(mean_per_site = mean(year_average)) %>%
    group_by(country_name) %>% 
    slice_max(mean_per_site, n=1) 
  
  # Data and yield predictions for each sites and each years 
  sites_var_i <- tab_world_i %>%
    filter(gridcode %in% c(unique(sites_min$gridcode), 
                           unique(sites_max$gridcode))) %>% 
    dplyr::select(gridcode, site_year, Ya, country_code, 
                  year_average, starts_with('monthly_var'), 
                  PC1_var, PC2_var) %>%
    left_join(tab_preds %>% 
                filter(Country=="WORLD", 
                       Model == "rf_pca.m.2", 
                       Type_cv=="01_YEARS", 
                       Outcome=="01_Ya"), 
              by=c("site_year"="Site_year"))
  
  # Scale factors for dual axis
  scaleFactor_1 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_1)
  scaleFactor_2 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_2)
  scaleFactor_3 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_3)
  scaleFactor_4 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_4)
  scaleFactor_5 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_5)
  scaleFactor_6 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_6)
  scaleFactor_7 <- max(sites_var_i$Predicted) / max(sites_var_i$monthly_var_7)
  scaleFactor   <- max(sites_var_i$Predicted) / max(sites_var_i$year_average)
  
  # Plots
  # > Seasonal averages
  ps <- sites_var_i %>% 
    gather(key = "PC", value = "score", starts_with("PC")) %>% 
    #mutate(PC = substr(PC,1,3)) %>% 
    ggplot(data = .) + 
    geom_line(aes(x=score, y=year_average*scaleFactor), color="blue") + 
    geom_point(aes(x = score, 
                   y = Predicted), size=0.3) + 
    facet_wrap(.~PC, scales="free_x") +
    scale_y_continuous(name="Predicted yield", sec.axis=sec_axis(~./scaleFactor, name=paste0("Seasonal\n", lab2_vari_i))) +
    labs(x = "PC score")
  
  # > Monthly averages
  pm <- sites_var_i %>% 
    gather(key = "PC", value = "score", starts_with("PC")) %>% 
    mutate(monthly_var_1_rec = monthly_var_1*scaleFactor_1,
           monthly_var_2_rec = monthly_var_2*scaleFactor_2,
           monthly_var_3_rec = monthly_var_3*scaleFactor_3,
           monthly_var_4_rec = monthly_var_4*scaleFactor_4,
           monthly_var_5_rec = monthly_var_5*scaleFactor_5,
           monthly_var_6_rec = monthly_var_6*scaleFactor_6,
           monthly_var_7_rec = monthly_var_7*scaleFactor_7) %>% 
    gather(key = "monthly_averages", value = "averages_recalculated", ends_with("rec")) %>%
    mutate(month_lab = paste0("Month ", substr(monthly_averages, nchar(monthly_averages)-4, nchar(monthly_averages)-4))) %>%
    ggplot(data = .) + 
    geom_line(aes(x=score, y=averages_recalculated, color=month_lab), linetype = 1) + 
    geom_point(aes(x = score, 
                   y = Predicted), size=0.3) + 
    facet_grid(month_lab~PC, scales="free_x") +
    scale_color_manual(values=c(pal_country), name="") + 
    scale_y_continuous(name="Predicted yield", sec.axis=sec_axis(~./scaleFactor, name=paste0(lab_vari_i))) +
    theme(legend.position = "none", 
          strip.placement = "outside", 
          strip.text.y = element_text(angle=0)) +
    labs(x = "PC score")
  
  psm <- plot_grid(
     ps + ggtitle("A. Predicted yields and average min temperature over the growing season"),
     pm + ggtitle("B. Predicted yields and average min temperature by month"), 
     ncol=1, rel_heights = c(0.2,0.8), align = "hv", axis = "tbrl"
    
  )
  
  # Save
  # > both
  ggsave(plot = psm, 
         filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_expl_", abb_var_i_quo, ".png"), 
         width = 7, height=12, dpi=300, bg="white")
  
  # > seasonal only 
  ggsave(plot = ps, 
         filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_expl_seasonal_", abb_var_i_quo, ".png"), 
         width = 7, height=2, dpi=300, bg="white")
  
  
}


# Home-made functions performing the dimension reductions
#source("E:/POSTDOC INRAE/ANALYSES/A_MODEL_COMP/00_Functions_dimension_reduction.R")

# PCA on monthly data
load("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/00_dim_red/pca_M_world.rda")

# Eigenvalues of the scores associated 
# with each principal component
variable.eig.tab <- pca_world$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Eigenvalues
    data.frame(CP = paste0("Dim.", 1:length(.x$pca.eig[,2])), 
               variance.percent = .x$pca.eig[,2])
    
  }, .id='var') %>%
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb")) 

# Contribution and correlation of averages climate data
# with the scores associated with principal components
variable.contrib.tab <- pca_world$list_pca_per_variable %>% 
  map_dfr(., ~{
    
    # > Contribution + correlation of each variable to the first dimensions
    .x$pca.var %>% 
      mutate(metric = rownames(.)) %>%
      gather(key = "pca_feature_full", value = "value", -metric) %>%
      separate(col = "pca_feature_full", into = c("pca_feature", "B", "Dim"), remove = T) %>% 
      unite("CP", B:Dim, sep = ".") %>%
      spread(key = "pca_feature", value = "value") 
    
  }, .id='var') %>% 
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  mutate(month_lab = paste0("Month ", substr(metric, nchar(metric), nchar(metric)))) %>% 
  mutate(CP_lab = paste0("Score ", substr(CP, nchar(CP), nchar(CP)))) %>% 
  left_join(., vars_names, by = c("var"="clim.var_abb"))

# First tentative of graphic
contrib_var_2 <- variable.contrib.tab %>%
  filter(CP %in% c("Dim.1", "Dim.2")) %>%
  # Add explained variance into the graph
  left_join(variable.eig.tab, by = c("CP", "CP_lab", "var", "clim.var",  "clim.var_lab", "clim.var_lab2")) %>%
  mutate(CP_lab2 = paste0(CP_lab, " (", round(variance.percent, 1), "%)"),
         contrib_lab = round(contrib, 1)) %>%
  # Plot
  ggplot(., aes(x = month_lab, y = cor, fill = contrib, label = contrib_lab)) + 
  geom_col(color="black") + 
  #geom_point(size = 2) +
  #geom_text(aes(y = cor), color = "black", size = 2, nudge_y = 0.1) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
  scale_fill_gradientn(colours = c("white", pal[2], "blue", "darkblue", "black"),
                        #breaks = c(-1, -0.5, 0, 0.5, 1),
                        name = "Contribution of variables to dimensions (%)") +
  coord_flip() + 
  facet_grid(CP_lab ~ clim.var_lab, scales = "free", switch = "y", space = "free") +
  labs(y = "Correlation between monthly averages and scores associated with each principal component") + 
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        axis.title.y = element_blank()) ; contrib_var_2

# Partial dependency 

# > fit
library(ranger)
library(pdp)
set.seed(101)
mod  <- ranger(Ya ~ irrigated_portion +
                    PC1_month_max_temp + PC1_month_min_temp + PC1_month_et0 + PC1_month_rad + PC1_month_prec + PC1_month_vapor_pressure_deficit + 
                    PC2_month_max_temp + PC2_month_min_temp + PC2_month_et0 + PC2_month_rad + PC2_month_prec + PC2_month_vapor_pressure_deficit,
               data=tab_world, 
               num.tree=500,
               importance="impurity") 

# > variables importance 
variable.imp.tab <- data.frame(Predictor  = names(mod$variable.importance), 
                               Importance = mod$variable.importance) %>%
  mutate(Score_lab = case_when(Predictor %in% c("PC1_month_prec", "PC1_month_min_temp", "PC1_month_max_temp", "PC1_month_vapor_pressure_deficit", "PC1_month_rad", "PC1_month_et0") ~ "Score 1",
                                   Predictor %in% c("PC2_month_prec", "PC2_month_min_temp", "PC2_month_max_temp", "PC2_month_vapor_pressure_deficit", "PC2_month_rad", "PC2_month_et0") ~ "Score 2", 
                                   Predictor == "irrigated_portion" ~ "Irrigation fraction"),
         var = case_when(Predictor %in% c("PC1_month_prec", "PC2_month_prec") ~ "prec",
                         Predictor %in% c("PC1_month_min_temp", "PC2_month_min_temp") ~ "min_temp", 
                         Predictor %in% c("PC1_month_max_temp", "PC2_month_max_temp") ~ "max_temp", 
                         Predictor %in% c("PC1_month_vapor_pressure_deficit", "PC2_month_vapor_pressure_deficit") ~ "vpd", 
                         Predictor %in% c("PC1_month_rad", "PC2_month_rad") ~ "rad", 
                         Predictor %in% c("PC1_month_et0", "PC2_month_et0") ~ "et0",
                         Predictor == "irrigated_portion" ~ "Irrigation fraction")) %>% 
  arrange(desc(Importance)) %>%
  mutate(Rank = 1:length(names(mod$variable.importance))) 
rownames(variable.imp.tab) <- NULL
variable.imp.tab %>% dplyr::select(Predictor, Importance, Rank)

#                           Predictor Importance Rank
# 1                    PC1_month_prec  45809.835    1
# 2                PC1_month_min_temp  29666.185    2
# 3                PC1_month_max_temp  27543.866    3
# 4  PC1_month_vapor_pressure_deficit  26385.484    4
# 5                     PC1_month_rad  22174.020    5
# 6                     PC1_month_et0  15337.147    6
# 7                 irrigated_portion  11941.108    7
# 8                     PC2_month_et0  10906.681    8
# 9                     PC2_month_rad  10851.362    9
# 10               PC2_month_min_temp   8371.250   10
# 11 PC2_month_vapor_pressure_deficit   8140.138   11
# 12               PC2_month_max_temp   5250.271   12
# 13                   PC2_month_prec   4265.243   13

# > variables pdp
variable.pdp <- list()

# > setting for parallelization
n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop

# >>> create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

for(var in names(mod$variable.importance))
{ 
  
  # > compute pdp for each variable 
  # (takes long time!)
  variable.pdp[[paste0(var)]] <- mod %>% 
    partial(pred.var = var) %>% 
    dplyr::rename("value"=1) %>%
    mutate(variable = var)
  
}

# > transform list into dataframe
variable.pdp.tab <- plyr::ldply(variable.pdp, data.frame, .id = "variable") %>% 
  mutate(clim.var = if_else(variable != "irrigated_portion", substr(variable, 11, nchar(variable)), variable)) %>%
  mutate(var = if_else(var == "vapor_pressure_deficit", "vpd", var)) %>% 
  left_join(., vars_names, by = "clim.var")

# >>> stop cluster//
stopCluster(my.cluster)

save(variable.pdp.tab, file = paste0(data_path, "06_pdp_imp/best_03_WORLD.rda"))

# First tentative graphic pdp
load(paste0(data_path, "06_pdp_imp/best_03_WORLD.rda"))

variable.pdp.tab2 <- variable.pdp.tab %>% 
  # importance values 
  left_join(variable.imp.tab, by = c("variable" = "Predictor")) %>% 
  # variable names
  left_join(vars_names, by = c("var" = "clim.var_abb")) %>% 
  # label for irrigation 
  mutate(clim.var = if_else(variable == "irrigated_portion", "Irrigation fraction", clim.var),
         clim.var_lab = if_else(variable == "irrigated_portion", "Irrigation fraction", clim.var_lab),
         clim.var_lab2 = if_else(variable == "irrigated_portion", "irrigation fraction", clim.var_lab2)) 


# Graphics 
explain_plots <- variable.pdp.tab2 %>% 
  mutate(clim.var_lab = recode(clim.var_lab, 
                               "Maximum temperature (°C)"="Maximum temperature",
                               "Minimum temperature (°C)" = "Minimum temperature",
                               "Evapotransp. ref (mm/day)" = "Evapotransp. ref",
                               "Net solar radiations (MJ/m²)"="Net solar radiation",
                               "Total precipitations (mm)" = "Average precipitation")) %>% 
  split(.$var) %>% 
  map(., ~{
    
    .x %>% 
      split(.$Score_lab) %>% 
      map(., ~{
        
        var_i <- unique(.x$var)
        lab_var_i <- unique(.x$clim.var_lab)
        score_i <- unique(.x$Score_lab)
        pred_i <- unique(.x$Predictor)
        imp_i <- round(unique(.x$Importance), 1)
        rank_i <- unique(.x$Rank)
        
        # Contribution and correlation of averages climate data
        # with the scores associated with principal components
        p1 <- variable.contrib.tab %>%
          filter(var == var_i, 
                 CP_lab == score_i) %>%
          # Add explained variance into the graph
          left_join(variable.eig.tab, by = c("CP", "CP_lab", "var", "clim.var",  "clim.var_lab", "clim.var_lab2")) %>%
          mutate(CP_lab2 = paste0(CP_lab, " (", round(variance.percent, 1), "%)"),
                 contrib_lab = round(contrib, 1),
                 cor_lab = round(cor, 2)) %>%
          # Lab with contribution 
          mutate(nudge_x_lab = if_else(cor < 0, -1, 1),
                 x_lab = cor + 0.1*nudge_x_lab) %>%
          # reverse y lab order
          mutate(month_lab = factor(month_lab, levels = rev(c('Month 1', 'Month 2', 'Month 3', 'Month 4', 'Month 5', 'Month 6', 'Month 7')))) %>% 
          # recode
          mutate(clim.var_lab = recode(clim.var_lab, 
                                       "Maximum temperature (°C)"="Maximum\ntemperature",
                                       "Minimum temperature (°C)" = "Minimum\ntemperature",
                                       "Evapotransp. ref (mm/day)" = "Evapotransp.\nref",
                                       "Net solar radiations (MJ/m²)"="Net solar\nradiation",
                                       "Total precipitations (mm)" = "Average\nprecipitation",
                                       "Vapor pressure deficit"="Vapor pressure\ndeficit"), 
                 clim.var_lab = paste0(clim.var_lab, "\n\nImportance score:\n", imp_i)) %>% 
          # Plot
          ggplot(., aes(x = month_lab, y = cor, 
                        #fill = contrib, label = contrib_lab)) + 
                        label = cor_lab)) + 
          geom_hline(yintercept = 0, color="darkgrey", linetype = 2) + 
          geom_col(fill = "lightgrey", color="black", width = 0.75) + 
          #geom_point(size = 2) +
          geom_text(aes(y = x_lab), color = "black", size=2.5) + 
          theme_bw() +
          theme(panel.grid = element_blank(),
                #strip.background = element_blank(),
                #plot.title = element_text(face = "bold"),
                legend.position = "bottom",
                strip.placement = "outside", 
                strip.text.y.left = element_text(angle = 0, size=10),
                axis.title.y = element_blank()) +
          scale_fill_gradientn(colours = c("white", pal[2],  "blue", "darkblue", "black"),
                               limits = c(0,100),
                               name = "Contribution of variables to dimensions (%)") +
          coord_flip() + 
          facet_grid(clim.var_lab~., switch = "y") +
          #ggtitle(paste0(score_i, " - ", lab_var_i, "\nImportance: ", imp_i, " (rank: ", rank_i, ")")) + 
          labs(y = paste0("Correlation between monthly averages and ", score_i)) + 
          lims(y = c(-1.1, 1.1))
        
        # plot PDP
        p2 <- .x %>% 
          ggplot(.) +
          geom_line(aes(x = value, y = yhat)) + 
          #facet_wrap(Score_lab ~ clim.var_lab, scales="free") + 
          theme_bw() + 
          labs(x = paste0(score_i), y = "Predicted yield\n(t/ha)") + 
          lims(y = c(1, 2.5))
        
        list("contrib" = p1, 
             "pdp" = p2)
        
        
      })
    
  })



# Figure in the main

plot_grid(
  plot_grid(
    # contrib
    plot_grid(explain_plots$prec$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$min_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$max_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              ncol=1),
    # pdp
    plot_grid(explain_plots$prec$`Score 1`$pdp, 
              explain_plots$min_temp$`Score 1`$pdp, 
              explain_plots$max_temp$`Score 1`$pdp, 
              ncol=1), 
    # settings
    align="hv", axis="trlb", labels = c('a.', 'b.')),
  # legend
  get_legend(explain_plots$prec$`Score 1`$contrib), 
  ncol=1, rel_heights = c(0.95, 0.05)
)


ggsave(filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Figure_expl_3_best.png"), 
       width = 10, height=10, dpi=300, bg="white")

# Supplementary Figure 
plot_grid(
  #plot_grid(
    # contrib
    plot_grid(explain_plots$prec$`Score 1`$contrib + 
                theme(legend.position = "none") +
                ggtitle("a."), 
              explain_plots$min_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$max_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$vpd$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$rad$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$et0$`Score 1`$contrib + theme(legend.position = "none"), 
              ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
    # pdp
    plot_grid(explain_plots$prec$`Score 1`$pdp +
                ggtitle("b."), 
              explain_plots$min_temp$`Score 1`$pdp, 
              explain_plots$max_temp$`Score 1`$pdp, 
              explain_plots$vpd$`Score 1`$pdp, 
              explain_plots$rad$`Score 1`$pdp, 
              explain_plots$et0$`Score 1`$pdp,
              ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
    # settings
    align="hv", axis="trlb", 
    rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 1`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_expl_1.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

plot_grid(
  #plot_grid(
    # contrib
    plot_grid(
              explain_plots$et0$`Score 2`$contrib + 
                theme(legend.position = "none") +
                ggtitle("a."), 
              explain_plots$rad$`Score 2`$contrib + theme(legend.position = "none"), 
              explain_plots$min_temp$`Score 2`$contrib + theme(legend.position = "none"), 
              explain_plots$vpd$`Score 2`$contrib + theme(legend.position = "none"), 
              explain_plots$max_temp$`Score 2`$contrib + theme(legend.position = "none"), 
              explain_plots$prec$`Score 2`$contrib + theme(legend.position = "none"), 
              ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
    # pdp
    plot_grid(explain_plots$et0$`Score 2`$pdp +
                ggtitle("b."),
              explain_plots$rad$`Score 2`$pdp, 
              explain_plots$min_temp$`Score 2`$pdp, 
              explain_plots$vpd$`Score 2`$pdp, 
              explain_plots$max_temp$`Score 2`$pdp, 
              explain_plots$prec$`Score 2`$pdp, 
              ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
    # settings
    align="hv", axis="trlb", 
    rel_widths = c(0.55,0.45)
  # legend
  #, plot_grid(get_legend(explain_plots$prec$`Score 2`$contrib), ggplot() + theme_void(), nrow=1),
  #ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_expl_2.png"), 
       width = 11.5, height=13, dpi=300, bg="white")


# Supplementary Figure 
plot_grid(
  plot_grid(
    # contrib scores 1
    plot_grid(explain_plots$prec$`Score 1`$contrib + 
                theme(legend.position = "none") +
                ggtitle("a."), 
              explain_plots$min_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$max_temp$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$vpd$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$rad$`Score 1`$contrib + theme(legend.position = "none"), 
              explain_plots$et0$`Score 1`$contrib + theme(legend.position = "none"), 
              ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)),
    # contrib scores 2
    plot_grid(
      explain_plots$et0$`Score 2`$contrib + 
        theme(legend.position = "none") +
        ggtitle("b."), 
      explain_plots$rad$`Score 2`$contrib + theme(legend.position = "none"), 
      explain_plots$min_temp$`Score 2`$contrib + theme(legend.position = "none"), 
      explain_plots$vpd$`Score 2`$contrib + theme(legend.position = "none"), 
      explain_plots$max_temp$`Score 2`$contrib + theme(legend.position = "none"), 
      explain_plots$prec$`Score 2`$contrib + theme(legend.position = "none"), 
      ncol=1, rel_heights = c(0.17, 0.15, 0.15, 0.15, 0.15, 0.15)), 
    # settings
    align="hv", axis="trlb"),
  # legend
  get_legend(explain_plots$prec$`Score 1`$contrib),
  ncol=1, rel_heights = c(0.95, 0.05)
)

ggsave(filename=paste0("E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/01_ERL/R1/new_figs/Supp_Figure_expl_3_for_reviewers_response.png"), 
       width = 11.5, height=13, dpi=300, bg="white")

# -----------------------------
# -----------------------------
# Figures asked after Revision 2

# Statistical differences in model performances 
# Anova to test the difference between groups of models 

tab_preds_folds <- tab_preds %>% 
  # > only World data 
  filter(Country=="WORLD", Outcome == "01_Ya") %>%
  # > compute performance metric per fold
  group_by(Type_cv, Model, Fold_cv) %>% 
  summarise(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted)
  ) %>% 
  # Labels 
  # > group of models (rf or lm)
  mutate(gpe_model = ifelse(sub("\\_.*", "", Model) == "lm", "lm", "rf")) %>%  
  mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
  mutate(Model = sub(".*_", "", Model)) %>% 
  # > change annual by seasonal
  mutate(Model = recode(Model, "avg.zscore.a"="avg.zscore.s", "avg.a"="avg.s")) %>%
  # > change order in Models
  mutate(Model = factor(Model, levels = rev(c("avg.zscore.s", 
                                              "avg.s", "avg.m", 
                                              "pca.m.1",   "pca.m.2",   "pca.m.3",   "pca.m.all",
                                              "pca.d.1",   "pca.d.2",   "pca.d.3",   "pca.d.all",
                                              "fpca.m.1",  "fpca.m.2",  "fpca.m.3",  "fpca.m.all",
                                              "fpca.d.1",  "fpca.d.2",  "fpca.d.3",  "fpca.d.all",
                                              "mfpca.m.1", "mfpca.m.2", "mfpca.m.3", "mfpca.m.all",
                                              "mfpca.d.1", "mfpca.d.2", "mfpca.d.3", "mfpca.d.all",
                                              "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all",
                                              "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                                              "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all",
                                              "fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all")))) %>% 
  # > type of data (daily, month, annual)
  mutate(data_type=case_when(
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",
                 "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                 'mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 
                 'fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all', 
                 'pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all') ~ "Daily climatic predictors", 
    Model %in% c("avg.s", "avg.zscores.s") ~ "Seasonal climatic predictors", 
    TRUE ~ "Monthly climatic predictors")) %>% 
  mutate(data_type = factor(data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Seasonal climatic predictors"))) %>% 
  # > dimension reduction techniques
  mutate(dimension_reduction=case_when(
    Model %in% c('pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all',   'pca.m.1',   'pca.m.2',   'pca.m.3',   'pca.m.all')  ~"PCA",
    Model %in% c('fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all',  'fpca.m.1',  'fpca.m.2',  'fpca.m.3',  'fpca.m.all') ~"FPCA",
    Model %in% c('mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 'mfpca.m.1', 'mfpca.m.2', 'mfpca.m.3', 'mfpca.m.all')~"MFPCA",
    Model %in% c("pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",   "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all")  ~"PLS",
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",  "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all") ~"FPLS",
    Model == "avg.zscore.s" ~ "Standardised\nmean",
    TRUE ~ "Averages")) %>% 
  mutate(dimension_reduction = factor(dimension_reduction, levels = c("Standardised\nmean", "Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS"))) %>% 
  # > nb of scores 
  mutate(nb_scores = case_when(
    Model %in% c('pca.d.1',   'fpca.d.1',   'mfpca.d.1',   'pca.m.1',   'fpca.m.1',   'mfpca.m.1',   'pls.d.1',   'pls.m.1',   'fpls.d.1',   'fpls.m.1')  ~"1 score",
    Model %in% c('pca.d.2',   'fpca.d.2',   'mfpca.d.2',   'pca.m.2',   'fpca.m.2',   'mfpca.m.2',   'pls.d.2',   'pls.m.2',   'fpls.d.2',   'fpls.m.2')  ~"2 scores",
    Model %in% c('pca.d.3',   'fpca.d.3',   'mfpca.d.3',   'pca.m.3',   'fpca.m.3',   'mfpca.m.3',   'pls.d.3',   'pls.m.3',   'fpls.d.3',   'fpls.m.3')  ~"3 scores",
    Model %in% c('pca.d.all', 'fpca.d.all', 'mfpca.d.all', 'pca.m.all', 'fpca.m.all', 'mfpca.m.all', 'pls.d.all', 'pls.m.all', 'fpls.d.all', 'fpls.m.all')~"All scores",
    TRUE ~ "")) %>% 
  mutate(nb_scores = factor(nb_scores, levels = c("1 score", "2 scores", "3 scores", "All scores", ""))) %>% 
  # > type cv 
  mutate(Type_cv_lab = recode(Type_cv, "01_YEARS" = "Cross-validation on years", "02_GEO" = "Cross-validation on sites"),
         Type_cv_lab = factor(Type_cv_lab, levels = c("Cross-validation on years", "Cross-validation on sites"))) 

# Test significance of differences between groups 
tab_preds_folds$gpe_model<-as.factor(tab_preds_folds$gpe_model)
tab_preds_folds$dimension_reduction<-as.factor(tab_preds_folds$dimension_reduction)
tab_preds_folds$data_type<-as.factor(tab_preds_folds$data_type)

# ANOVA
list_aov <- list()

for(cv_i in c("01_YEARS", "02_GEO"))
{
  
  # Select data for this CV 
  # and compute mean performance by model
  data_cv_i <- tab_preds_folds %>% filter(Type_cv==cv_i) %>%
    group_by(data_type, gpe_model, dimension_reduction, Model) %>%
    summarise(NSE = mean(NSE),
              RMSEP = mean(RMSEP))
  
  data_cv_lm <- data_cv_i %>% 
    filter(gpe_model=="Linear regression")
  
  data_cv_rf <- data_cv_i %>% 
    filter(gpe_model!="Linear regression")
  
  # Define factors
  data_cv_i$gpe_model           <-as.factor(data_cv_i$gpe_model)
  data_cv_i$dimension_reduction <-as.factor(data_cv_i$dimension_reduction)
  data_cv_i$data_type           <-as.factor(data_cv_i$data_type)
  
  # Fit models 
  # > NSE
  aov_nse_1  <- aov(NSE ~ gpe_model, data = data_cv_i)
  aov_nse_2  <- aov(NSE ~ dimension_reduction, data = data_cv_i)
  aov_nse_2a <- aov(NSE ~ dimension_reduction, data = data_cv_lm)
  aov_nse_2b <- aov(NSE ~ dimension_reduction, data = data_cv_rf)
  aov_nse_3  <- aov(NSE ~ data_type, data = data_cv_i)
  aov_nse_3a <- aov(NSE ~ data_type, data = data_cv_lm)
  aov_nse_3b <- aov(NSE ~ data_type, data = data_cv_rf)
  
  # > RMSE
  aov_rmse_1  <- aov(RMSEP ~ gpe_model, data = data_cv_i)
  aov_rmse_2a <- aov(RMSEP ~ dimension_reduction, data = data_cv_lm)
  aov_rmse_2b <- aov(RMSEP ~ dimension_reduction, data = data_cv_rf)
  aov_rmse_3a <- aov(RMSEP ~ data_type, data = data_cv_lm)
  aov_rmse_3b <- aov(RMSEP ~ data_type, data = data_cv_rf)
  
  # > Outputs
  list_aov[[paste0(cv_i)]] <- list(aov_nse_GM   =aov_nse_1, 
                                   aov_nse_DR   =aov_nse_2,
                                   aov_nse_DR_LM=aov_nse_2a,
                                   aov_nse_DR_RF=aov_nse_2b,
                                   aov_nse_DT   =aov_nse_3,
                                   aov_nse_DT_LM=aov_nse_3a,
                                   aov_nse_DT_RF=aov_nse_3b,
                                   
                                   aov_rmse_GM   =aov_rmse_1, 
                                   aov_rmse_DR_LM=aov_rmse_2a, 
                                   aov_rmse_DR_RF=aov_rmse_2b, 
                                   aov_rmse_DT_LM=aov_rmse_3a,
                                   aov_rmse_DT_RF=aov_rmse_3b)
  
}


library(broom)

list_aov %>% 
  map_dfr(., ~{
    
    .x %>% 
      map_dfr(., ~{
        
        broom::tidy(.x)
        
        
      }, .id="model")
    
    
  }, .id="cv") %>% 
  # remove unused terms 
  filter(term %in% c("gpe_model", 
                     "dimension_reduction",
                     "data_type")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  dplyr::select(-df, -sumsq, -meansq, -statistic) %>% 
  spread(key="cv", value="p.value") 

#   model          term                `01_YEARS` `02_GEO` 
# 7 aov_nse_GM     gpe_model           <0.001     <0.001  
# 1 aov_nse_DR     dimension_reduction 0.065      0.010   
# 2 aov_nse_DR_LM  dimension_reduction 0.006      0.006   
# 3 aov_nse_DR_RF  dimension_reduction 0.019      0.004   
# 4 aov_nse_DT     data_type           0.005      0.006   
# 5 aov_nse_DT_LM  data_type           0.002      0.002   
# 6 aov_nse_DT_RF  data_type           0.001      0.012   

# 12 aov_rmse_GM    gpe_model           <0.001     <0.001   
# 8 aov_rmse_DR_LM dimension_reduction 0.007      0.008   
# 9 aov_rmse_DR_RF dimension_reduction 0.018      0.005   
# 10 aov_rmse_DT_LM data_type           0.001      0.001   
# 11 aov_rmse_DT_RF data_type           <0.001     0.005   

# Mean and sd of NSE
tab_preds_folds %>% 
  group_by(Type_cv,dimension_reduction) %>% 
  summarise(mean_NSE=round(mean(NSE),2),
            sd_NSE=round(sd(NSE),2)) %>% 
  unite(lab, mean_NSE:sd_NSE, sep=" (", remove=T) %>% 
  mutate(lab=paste0(lab, ")")) %>% 
  spread(key="Type_cv", "lab") 
  
tab_preds_folds %>% 
  group_by(Type_cv,data_type) %>% 
  summarise(mean_NSE=round(mean(NSE),2),
            sd_NSE=round(sd(NSE),2)) %>% 
  unite(lab, mean_NSE:sd_NSE, sep=" (", remove=T) %>% 
  mutate(lab=paste0(lab, ")")) %>% 
  spread(key="Type_cv", "lab") 


# Define factors
tab_preds_folds <- tab_perf_labelled %>% 
  filter(Outcome=="01_Ya", Country=="WORLD", pred_perf=="NSE")

tab_preds_folds$gpe_model           <-factor(tab_preds_folds$gpe_model, levels = c("Linear regression", "Random forest"))
tab_preds_folds$dimension_reduction <-factor(tab_preds_folds$dimension_reduction, levels=c("Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS", "Standardised\nmean"))
tab_preds_folds$data_type           <-factor(tab_preds_folds$data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Seasonal climatic predictors"))

mod_nse <- lm(pred_perf_value ~ dimension_reduction + gpe_model + data_type 
              + N_pred + Type_cv, 
              data = tab_preds_folds)

broom::tidy(mod_nse, conf.int=T) %>%
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(term, lab, p.value) 

AIC(mod_nse)
BIC(mod_nse)
summary(mod_nse)$r.squared

# Linear models 
list_mods <- list()

for(cv_i in c("01_YEARS", "02_GEO"))
{
  
  # Select data for this CV 
  # and compute mean performance by model
  data_cv_i <- tab_preds_folds %>% filter(Type_cv==cv_i) %>%
    group_by(data_type, gpe_model, dimension_reduction, Model) %>%
    summarise(NSE = mean(NSE),
              RMSEP = mean(RMSEP))
  
  # Define factors
  data_cv_i$gpe_model           <-factor(data_cv_i$gpe_model, levels = c("Linear regression", "Random forest"))
  data_cv_i$dimension_reduction <-factor(data_cv_i$dimension_reduction, levels=c("Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS", "Standardised\nmean"))
  data_cv_i$data_type           <-factor(data_cv_i$data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Seasonal climatic predictors"))
  
  # Fit models 
  # > NSE
  # > without model equation
  mod_nse_0 <- lm(NSE ~ gpe_model + Model, data = data_cv_i)
  mod_nse_01a <- lm(NSE ~ dimension_reduction, data = data_cv_i)
  mod_nse_01b <- lm(NSE ~ dimension_reduction + gpe_model, data = data_cv_i)
  mod_nse_02a <- lm(NSE ~ data_type, data = data_cv_i)
  mod_nse_02b <- lm(NSE ~ data_type + gpe_model, data = data_cv_i)
  mod_nse_03 <- lm(NSE ~ dimension_reduction + gpe_model + data_type, data = data_cv_i)
  
  # > including model's equation
  mod_nse_1 <- lm(NSE ~ gpe_model + Model, data = data_cv_i)
  
  mod_nse_2a <- lm(NSE ~ dimension_reduction + Model, data = data_cv_i)
  mod_nse_2b <- lm(NSE ~ dimension_reduction + gpe_model + Model, data = data_cv_i)
  
  mod_nse_3a <- lm(NSE ~ data_type + Model, data = data_cv_i)
  mod_nse_3b <- lm(NSE ~ data_type + gpe_model + Model, data = data_cv_i)
  
  mod_nse_4 <- lm(NSE ~ dimension_reduction + gpe_model + data_type + Model, data = data_cv_i)
  
  mod_rmse_4 <- lm(RMSEP ~ dimension_reduction + gpe_model + data_type + Model, data = data_cv_i)
  # > Outputs
  list_mods[[paste0(cv_i)]] <- list(mod_nse_1=mod_nse_1,  
                                    mod_nse_2a=mod_nse_2a,
                                    mod_nse_2b=mod_nse_2b,
                                    mod_nse_3a=mod_nse_3a,
                                    mod_nse_3b=mod_nse_3b,
                                    mod_nse_4=mod_nse_4,
                                    mod_nse_0=mod_nse_0,
                                    mod_nse_01a=mod_nse_01a,
                                    mod_nse_01b=mod_nse_01b,
                                    mod_nse_02a=mod_nse_02a,
                                    mod_nse_02b=mod_nse_02b,
                                    mod_nse_03=mod_nse_03,
                                    mod_rmse_4=mod_rmse_4)
  
}

# AIC BIC of the models 
list_mods %>% 
  map_dfr(., ~{
    
    .x %>% 
      map_dfr(., ~{
        
        data.frame(AIC = AIC(.x),
                   BIC = BIC(.x),
                   R2 = summary(.x)$r.squared,
                   N_predictors = nrow(summary(.x)$coefficient)-1)
        
        
      }, .id="model")
    
    
  }, .id="cv") %>% 
  arrange(cv, AIC, BIC, desc(R2), N_predictors)

#         cv       model         AIC        BIC        R2 N_predictors
#1  01_YEARS  mod_nse_2b -192.098899  -81.65327 0.9672321           43 *
#2  01_YEARS   mod_nse_1 -192.098899  -81.65327 0.9672321           43 *
#3  01_YEARS   mod_nse_4 -192.098899  -81.65327 0.9672321           43 *
#4  01_YEARS   mod_nse_0 -192.098899  -81.65327 0.9672321           43 *
#5  01_YEARS  mod_nse_3b -192.098899  -81.65327 0.9672321           43 *
#6  01_YEARS  mod_nse_03 -111.173302  -84.17548 0.8148543            9
#7  01_YEARS mod_nse_02b  -75.200616  -62.92888 0.6765757            3
#8  01_YEARS mod_nse_01b  -71.482676  -49.39355 0.6922851            7
#9  01_YEARS mod_nse_02a    8.811324   18.62871 0.1207195            2
#10 01_YEARS mod_nse_01a   15.260939   34.89572 0.1364289            6
#11 01_YEARS  mod_nse_3a   54.298135  162.28942 0.4113759           42
#12 01_YEARS  mod_nse_2a   54.298135  162.28942 0.4113759           42

#13   02_GEO   mod_nse_1 -234.114866 -123.66924 0.9812001           43 *
#14   02_GEO   mod_nse_0 -234.114866 -123.66924 0.9812001           43 *
#15   02_GEO  mod_nse_3b -234.114866 -123.66924 0.9812001           43 *
#16   02_GEO  mod_nse_2b -234.114866 -123.66924 0.9812001           43 *
#17   02_GEO   mod_nse_4 -234.114866 -123.66924 0.9812001           43 *
#18   02_GEO  mod_nse_03 -105.271562  -78.27374 0.8145611            9
#19   02_GEO mod_nse_01b  -65.318107  -43.22898 0.6908545            7
#20   02_GEO mod_nse_02b  -54.887669  -42.61593 0.6169677            3
#21   02_GEO mod_nse_02a   15.200523   25.01791 0.1143213            2
#22   02_GEO mod_nse_01a   15.709038   35.34382 0.1882082            6
#23   02_GEO  mod_nse_2a   49.642185  157.63347 0.4785538           42
#24   02_GEO  mod_nse_3a   49.642185  157.63347 0.4785538           42

# * models with similar performances

# Coefs of the full model 
broom::tidy(list_mods$"01_YEARS"$mod_nse_4, conf.int=T) %>%
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(term, lab, p.value) %>%
  # remove unused terms 
  filter(term %in% c("gpe_modelRandom forest", 
                     "dimension_reductionAverages", "dimension_reductionPCA", "dimension_reductionFPCA", "dimension_reductionMFPCA", "dimension_reductionPLS", "dimension_reductionFPLS", "dimension_reductionStandardised\nmean",
                     "data_typeDaily climatic predictors", "data_typeMonthly climatic predictors", "data_typeSeasonal climatic predictors")) 

#  term                                    lab                  p.value
#1 "dimension_reductionPCA"                -0.06 (-0.19, 0.08)  0.387  
#2 "dimension_reductionFPCA"               -0.07 (-0.2, 0.07)   0.315  
#3 "dimension_reductionMFPCA"              -0.3 (-0.44, -0.17)  <0.001 
#4 "dimension_reductionPLS"                -0.06 (-0.19, 0.08)  0.404  
#5 "dimension_reductionFPLS"               -0.38 (-0.57, -0.19) <0.001 
#6 "dimension_reductionStandardised\nmean" -0.34 (-0.48, -0.21) <0.001 
#7 "gpe_modelRandom forest"                0.39 (0.36, 0.42)    <0.001 
#8 "data_typeMonthly climatic predictors"  0.05 (-0.08, 0.19)   0.438  
#9 "data_typeSeasonal climatic predictors" -0.01 (-0.2, 0.18)   0.939  

broom::tidy(list_mods$"02_GEO"$mod_nse_4, conf.int=T) %>%
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(term, lab, p.value) %>%
  # remove unused terms 
  filter(term %in% c("gpe_modelRandom forest", 
                     "dimension_reductionAverages", "dimension_reductionPCA", "dimension_reductionFPCA", "dimension_reductionMFPCA", "dimension_reductionPLS", "dimension_reductionFPLS", "dimension_reductionStandardised\nmean",
                     "data_typeDaily climatic predictors", "data_typeMonthly climatic predictors", "data_typeSeasonal climatic predictors")) 

#  term                                    lab                  p.value
#1 "dimension_reductionPCA"                -0.09 (-0.2, 0.01)   0.088  
#2 "dimension_reductionFPCA"               -0.1 (-0.21, 0.01)   0.064  
#3 "dimension_reductionMFPCA"              -0.47 (-0.58, -0.37) <0.001 
#4 "dimension_reductionPLS"                -0.08 (-0.19, 0.02)  0.115  
#5 "dimension_reductionFPLS"               -0.35 (-0.5, -0.2)   <0.001 
#6 "dimension_reductionStandardised\nmean" -0.52 (-0.62, -0.41) <0.001 
#7 "gpe_modelRandom forest"                0.38 (0.36, 0.4)     <0.001 
#8 "data_typeMonthly climatic predictors"  0.06 (-0.04, 0.17)   0.234  
#9 "data_typeSeasonal climatic predictors" -0.03 (-0.18, 0.12)  0.694  

# > RMSEP
# Coefs of the full model 
broom::tidy(list_mods$"01_YEARS"$mod_rmse_4, conf.int=T) %>%
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(term, lab, p.value) %>%
  # remove unused terms 
  filter(term %in% c("gpe_modelRandom forest", 
                     "dimension_reductionAverages", "dimension_reductionPCA", "dimension_reductionFPCA", "dimension_reductionMFPCA", "dimension_reductionPLS", "dimension_reductionFPLS", "dimension_reductionStandardised\nmean",
                     "data_typeDaily climatic predictors", "data_typeMonthly climatic predictors", "data_typeSeasonal climatic predictors")) 

broom::tidy(list_mods$"02_GEO"$mod_rmse_4, conf.int=T) %>%
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(term, lab, p.value) %>%
  # remove unused terms 
  filter(term %in% c("gpe_modelRandom forest", 
                     "dimension_reductionAverages", "dimension_reductionPCA", "dimension_reductionFPCA", "dimension_reductionMFPCA", "dimension_reductionPLS", "dimension_reductionFPLS", "dimension_reductionStandardised\nmean",
                     "data_typeDaily climatic predictors", "data_typeMonthly climatic predictors", "data_typeSeasonal climatic predictors")) 



# All coefs
library(broom)

tab_coefs <- list_mods %>% 
  map_dfr(., ~{
    
    .x %>% 
      map_dfr(., ~{
        
        broom::tidy(.x, conf.int=T)
        
        
      }, .id="model")
    
    
  }, .id="cv") %>% 
  # remove unused terms 
  filter(term %in% c("gpe_modelRandom forest", 
                     "dimension_reductionAverages", "dimension_reductionPCA", "dimension_reductionFPCA", "dimension_reductionMFPCA", "dimension_reductionPLS", "dimension_reductionFPLS", "dimension_reductionStandardised\nmean",
                     "data_typeDaily climatic predictors", "data_typeMonthly climatic predictors", "data_typeSeasonal climatic predictors")) %>% 
  # label with mean estimate (95% CI)
  mutate(lab = paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")) %>% 
  mutate(p.value = if_else(p.value < 0.001, "<0.001", format(round(p.value, 3), scientific=F, digits = 3))) %>%
  # remove old columns
  dplyr::select(cv, model, term, lab, p.value)


