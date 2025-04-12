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






# ----------------------------------
# Model predictive performance based on RMSEP
# > Load predictions for each model

# > Function to read rda into list
rda2list <- function(file) {
  e <- new.env()
  load(file, envir = e)
  as.list(e)
}

# > Apply on all preds stored 
tab_preds <- list(
  USA   = list(name="01_USA"),
  BRA   = list(name="02_BRA"),
  WORLD = list(name="03_WORLD")
) %>% 
  map_dfr(., ~{
    
    # > folder where preds are stored
    folder <- paste0("C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/04_preds/", .x$name)
    files <- list.files(folder, pattern = ".rda$")
    
    # > merge all preds for 1 country
    Map(rda2list, file.path(folder, files)) %>% 
      plyr::ldply(., data.frame, .id = "path")
    
  }, .id="Country") %>% 
  dplyr::select(-path)

# > Rename columns' name
head(tab_preds)
colnames(tab_preds) <- c("Country", "Year", "Model", "Site_year", "Observed", "Predicted", "N_predictors")

# > Compute RMSEP, NSE, R2 and Bias
# for each model in each country

tab_perf <- tab_preds %>% 
  group_by(Country, Model) %>% 
  summarise(
    RMSEP = caret::RMSE(obs = Observed,      pred = Predicted),
    NSE   = hydroGOF::NSE(obs = Observed,    sim = Predicted),
    R2    = caret::R2(obs = Observed,        pred = Predicted),
    Bias  = Metrics::bias(actual = Observed, predicted = Predicted),
    N_pred = max(N_predictors)
  )

# > Label models and country for plots
tab_perf_labelled <- tab_perf %>% 
  # set long format 
  pivot_longer(cols=c(RMSEP, NSE, R2, Bias), names_to = "pred_perf", values_to = "pred_perf_value") %>% 
  # groupe of models (rf or lm)
  mutate(gpe_model = ifelse(sub("\\_.*", "", Model) == "lm", "lm", "rf")) %>% 
  mutate(Model = sub(".*_", "", Model)) %>% 
  # > change order in Models
  mutate(Model = factor(Model, levels = rev(c("avg.zscore.a", 
                                              "avg.a", "avg.m", 
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
  # type of data (daily, month, annual)
  mutate(data_type=case_when(
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",
                 "pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",
                 'mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 
                 'fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all', 
                 'pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all') ~ "Daily climatic predictors", 
    Model %in% c("avg.a", "avg.zscores.a") ~ "Annual climatic predictors", 
    TRUE ~ "Monthly climatic predictors")) %>% 
  mutate(data_type = factor(data_type, levels = c("Daily climatic predictors", "Monthly climatic predictors", "Annual climatic predictors"))) %>% 
  # dimension reduction techniques
  mutate(dimension_reduction=case_when(
    Model %in% c('pca.d.1',   'pca.d.2',   'pca.d.3',   'pca.d.all',   'pca.m.1',   'pca.m.2',   'pca.m.3',   'pca.m.all')  ~"PCA",
    Model %in% c('fpca.d.1',  'fpca.d.2',  'fpca.d.3',  'fpca.d.all',  'fpca.m.1',  'fpca.m.2',  'fpca.m.3',  'fpca.m.all') ~"FPCA",
    Model %in% c('mfpca.d.1', 'mfpca.d.2', 'mfpca.d.3', 'mfpca.d.all', 'mfpca.m.1', 'mfpca.m.2', 'mfpca.m.3', 'mfpca.m.all')~"MFPCA",
    Model %in% c("pls.d.1",   "pls.d.2",   "pls.d.3",   "pls.d.all",   "pls.m.1",   "pls.m.2",   "pls.m.3",   "pls.m.all")  ~"PLS",
    Model %in% c("fpls.d.1",  "fpls.d.2",  "fpls.d.3",  "fpls.d.all",  "fpls.m.1",  "fpls.m.2",  "fpls.m.3",  "fpls.m.all") ~"FPLS",
    Model == "avg.zscore.a" ~ "Standardised\nmean",
    TRUE ~ "Averages")) %>% 
  mutate(dimension_reduction = factor(dimension_reduction, levels = c("Standardised\nmean", "Averages", "PCA", "FPCA", "MFPCA", "PLS", "FPLS"))) %>% 
  # best model per country and indicator
  group_by(Country, pred_perf) %>% 
  mutate(best_pred_perf_value = if_else(pred_perf %in% c("RMSEP", "Bias"), min(abs(pred_perf_value)), max(pred_perf_value))) %>% 
  filter(pred_perf %in% c("RMSEP", "NSE")) %>%
  mutate(pred_perf_lab = case_when(
    pred_perf == "RMSEP" ~ "RMSEP\n(lower is better)", 
    pred_perf == "Bias"  ~ "Bias\n(lower is better)", 
    pred_perf == "NSE"   ~ "NSE\n(higher is better)",
    pred_perf == "R2"    ~ "R squared\n(higher is better)"),
    pred_perf_lab = factor(pred_perf_lab, levels = c("RMSEP\n(lower is better)", 
                                                     "NSE\n(higher is better)",
                                                     "R squared\n(higher is better)",
                                                     "Bias\n(lower is better)"))) %>% 
  # country label 
  mutate(Country_lab = case_when(
    Country == "USA" ~ "United-States of America\n(N=29803)", 
    Country == "BRA" ~ "Brazil\n(N=14757)", 
    TRUE ~ "Global\n(N=122229)"
  )) %>% 
  mutate(Country_lab = factor(Country_lab, levels = c("Global\n(N=122229)",
                                                      "United-States of America\n(N=29803)", 
                                                      "Brazil\n(N=14757)"))) %>% 
  # nb of scores 
  mutate(nb_scores = case_when(
    Model %in% c('pca.d.1',   'fpca.d.1',   'mfpca.d.1',   'pca.m.1',   'fpca.m.1',   'mfpca.m.1',   'pls.d.1',   'pls.m.1',   'fpls.d.1',   'fpls.m.1')  ~"1 score",
    Model %in% c('pca.d.2',   'fpca.d.2',   'mfpca.d.2',   'pca.m.2',   'fpca.m.2',   'mfpca.m.2',   'pls.d.2',   'pls.m.2',   'fpls.d.2',   'fpls.m.2')  ~"2 scores",
    Model %in% c('pca.d.3',   'fpca.d.3',   'mfpca.d.3',   'pca.m.3',   'fpca.m.3',   'mfpca.m.3',   'pls.d.3',   'pls.m.3',   'fpls.d.3',   'fpls.m.3')  ~"3 scores",
    Model %in% c('pca.d.all', 'fpca.d.all', 'mfpca.d.all', 'pca.m.all', 'fpca.m.all', 'mfpca.m.all', 'pls.d.all', 'pls.m.all', 'fpls.d.all', 'fpls.m.all')~"All scores",
    TRUE ~ "")) %>% 
  mutate(nb_scores = factor(nb_scores, levels = c("1 score", "2 scores", "3 scores", "All scores", "")))


# RMSEP
p.rmsep <- tab_perf_labelled %>% 
  # > For 1 indicator
  filter(pred_perf == "RMSEP")  %>% 
  # > change order in data type and in dimension reduction
  #mutate(dimension_reduction = factor(dimension_reduction, levels = rev(c("Averages", "PCA", "FPCA", "MFPCA")))) %>% 
  #mutate(data_type = recode(data_type, "Daily climatic predictors"="Daily", "Monthly climatic predictors"="Monthly", "Annual climatic predictors"="Annual")) %>% 
  # > round indicators
  #mutate(pred_perf_value_lab = round(pred_perf_value, 3))  %>% 
  # > scores
  #mutate(nb_scores = factor(nb_scores, levels = rev(c("1 score", "2 scores", "3 scores", "All scores", ""))))%>% 
  # > order of models 
  arrange(Country_lab, pred_perf_value) %>% 
  group_by(Country_lab) %>% 
  mutate(order = row_number()) %>%
  mutate(best = if_else(order==1, "Best model", "")) %>%
  # > min and max among lm and rf
  group_by(Country, Model) %>% 
  mutate(min_x_lab = min(pred_perf_value),
         max_x_lab = max(pred_perf_value)) %>% 
  # > group of models 
  mutate(gpe_model = recode(gpe_model, "lm"="Linear regression", "rf"="Random forest")) %>% 
  # > Nb of predictors per model
  group_by(Country_lab, dimension_reduction, Model) %>%
  mutate(N_pred = max(N_pred, na.rm = T)) %>% 
  # > plot
  ggplot(., aes(x = pred_perf_value, y=Model)) +
  geom_linerange(aes(xmin = min_x_lab, xmax = max_x_lab), 
                 color = "black") +
  geom_point(aes(color = gpe_model, fill = best),
             size = 2, shape = 21) + 
  geom_text(aes(label = order), 
            nudge_y = 0.4, 
            size = 2) + 
  geom_text(aes(label = N_pred, x = 1.4), 
            size = 3, check_overlap = T) + 
  #lemon::facet_rep_grid(dimension_reduction ~ Country_lab, scales = "free_y", space = "free", switch = "y", repeat.tick.labels = F) + 
  #theme_cowplot() + 
  facet_grid(dimension_reduction ~ Country_lab, scales = "free_y", space = "free", switch = "y") + 
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
        panel.grid = element_blank()) + 
  scale_color_manual(values = pal_country[c(3, 4)], name = "Model family") + 
  scale_fill_manual(values = c("grey", "red"), name = "Rank model") + 
  #scale_shape_manual(values = c(21, 19), name = "Rank model") + 
  scale_x_continuous(breaks = seq(0.2, 1.2, by=0.2), labels = seq(0.2, 1.2, by=0.2)) + 
  guides(color = guide_legend(order = 1, override.aes = list(size=2), ncol = 1), 
         fill = guide_none(), 
         shape = guide_none())   + 
  labs(x = "Root mean square error of prediction (lower value indicates a better predictive performance)") ; p.rmsep

ggsave(plot = p.rmsep, 
       filename=paste0(save_path, "Figure_4.png"), 
       width = 10, height=11.5, dpi=300, bg="white")