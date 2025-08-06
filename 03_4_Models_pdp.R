# -------------------------------------------------------------------------
# 
#           PARTIAL DEPENDENCY PLOTS OF BEST MODELS  
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
library(caret) ; library(ranger) ; library(fastshap) ; library(boot) ; library(coxed) 
library(pdp)
# > functional analysis 
library(fda) ; library(MFPCA)
# > others
library(hydroGOF)

# -------------------------------------------------------------------------

# ----------------------------------
# Data 
data_path     <- "C:/Users/benni/Documents/Post doc/ERA5_data_comp_models/"
supp_fig_path <- "E:/POSTDOC INRAE/PAPERS/01_MODEL_COMP/FIGURES/SUPPLEMENTARY_FIGURES/"

load(paste0(data_path, "tab_usa.rda"))
load(paste0(data_path, "tab_bra.rda"))
load(paste0(data_path, "tab_world.rda"))

# ----------------------------------
# > CORRELATION BETWEEN SCORES AND CLIMATIC VARIABLES 

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
                               "vapor_pressure_deficit"   ="vapor_pressure_deficit"))  %>% 
  mutate(clim.var_lab = recode(clim.var, 
                               "min_2m_temperature"         ="Minimum temperature",
                               "max_2m_temperature"         ="Maximum temperature",
                               "et0"                        ="Evapotranspiration ref",
                               "surface_net_solar_radiation"="Solar radiations",
                               "total_precipitation"        ="Precipitation",
                               "vapor_pressure_deficit"   ="Vapor pressure deficit"))

# WORLD
# > Correlation between all PCA scores and monthly averages 
png(filename = paste0(supp_fig_path, "Correlation plots/WORLD_all_PCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_world %>% 
    dplyr::select(starts_with(paste0("PC1_month_", clim.var_abb_i)), 
                  starts_with(paste0("PC2_month_", clim.var_abb_i)), 
                  starts_with(paste0("PC3_month_", clim.var_abb_i)),
                  starts_with(paste0("PC4_month_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i)), 
                  irrigated_portion) %>% 
    dplyr::rename("PCA score 1"= 1,
                  "PCA score 2"= 2,
                  "PCA score 3"= 3,
                  "PCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11,
                  "Irrigation"=irrigated_portion) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# --------------------------------------------
# USA

# > Correlation between all MFPCA scores and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "Correlation plots/USA_all_MFPCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  
  tab_usa %>% 
    dplyr::select(MFPC1_month, MFPC2_month, MFPC3_month, MFPC4_month, starts_with(paste0("monthly_", var_i)), irrigated_portion) %>% 
    dplyr::rename("MFPCA score 1"=MFPC1_month,
                  "MFPCA score 2"=MFPC2_month,
                  "MFPCA score 3"=MFPC3_month,
                  "MFPCA score 4"=MFPC4_month,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11,
                  "Irrigation"=irrigated_portion) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# --------------------------------------------
# BRAZIL
# > Correlation between year and monthly averages in Brazil
png(filename = paste0(supp_fig_path, "Correlation plots/BRA_year_vs_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  
  tab_bra %>% 
    dplyr::select(starts_with(paste0("year_", var_i)), starts_with(paste0("monthly_", var_i)), irrigated_portion) %>% 
    dplyr::rename("Mean over 7 months"=1,
                  "Month 1"=2, 
                  "Month 2"=3, 
                  "Month 3"=4, 
                  "Month 4"=5, 
                  "Month 5"=6, 
                  "Month 6"=7, 
                  "Month 7"=8,
                  "Irrigation"=irrigated_portion) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# --------------------------------------------
# --------------------------------------------
# > IMPORTANCE AND PDP 
# for the best models  
imp.pdp <- list(
  USA   = list(data = tab_usa, 
               formula = "MFPC1_month",
               name="01_USA"),
  BRA   = list(data = tab_bra, 
               formula = "year_max_2m_temperature + year_min_2m_temperature + year_et0 + year_surface_net_solar_radiation + year_total_precipitation + year_vapor_pressure_deficit",
               name="02_BRA"),
  WORLD = list(data = tab_world, 
               formula = "PC1_month_max_temp + PC1_month_min_temp + PC1_month_et0 + PC1_month_rad + PC1_month_prec + PC1_month_vapor_pressure_deficit + PC2_month_max_temp + PC2_month_min_temp + PC2_month_et0 + PC2_month_rad + PC2_month_prec + PC2_month_vapor_pressure_deficit",
               name="03_WORLD")) %>%
  map(., ~{ 
    
    # > data to use for models fitting, name of the analysis, and bootstrap procedure
    dat_pred      <- .x$data
    model_name    <- .x$name
    tab_sites     <- dat_pred %>% 
      distinct(x, y, gridcode)
    
    # > model formula
    model_formula <- paste0("Ya ~ irrigated_portion + ", .x$formula)
    
    # > fit
    set.seed(101)
    mod  <- ranger(as.formula(model_formula),
                   data=dat_pred, 
                   num.tree=500,
                   importance="impurity") 
    
    # > variables importance 
    variable.imp.tab <- data.frame(predictor = names(mod$variable.importance), imp = mod$variable.importance)
    
    # > variables pdp
    variable.pdp <- list()
    for(var in names(mod$variable.importance))
    { 
      
      # > setting for parallelization
      n.cores <- parallel::detectCores() - 3 ; n.cores # 7 on the computer, 11 on the laptop
      
      # >>> create the cluster
      my.cluster <- parallel::makeCluster(
        n.cores, 
        type = "PSOCK"
      )
      
      # register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      # > compute pdp for each variable 
      # (takes long time!)
      variable.pdp[[paste0(var)]] <- mod %>% 
        partial(pred.var = var) %>% 
        dplyr::rename("variable"=1)
      
    }
    
    # > transform list into dataframe
    variable.pdp.tab <- plyr::ldply(variable.pdp, data.frame, .id = "variable")
    
    # >>> stop cluster//
    stopCluster(my.cluster)
    
    
    # > output 
    out <- list(
      "name" = model_name,
      "mod" = mod,
      "imp" = variable.imp.tab,
      "pdp" = variable.pdp.tab
    )
    
    # > save
    save(out, file = paste0(data_path, "06_pdp_imp/best_", model_name, ".rda"))
    
    out
})

# LOAD 
library(CCMHr)
imp.pdp <- list()
imp.pdp[[paste0("USA")]] <- loadRDa(file = paste0(data_path, "06_pdp_imp/best_01_USA.rda"))
imp.pdp[[paste0("BRA")]] <- loadRDa(file = paste0(data_path, "06_pdp_imp/best_02_BRA.rda"))
imp.pdp[[paste0("WORLD")]] <- loadRDa(file = paste0(data_path, "06_pdp_imp/best_03_WORLD.rda"))

# IMPORTANCE
library(tidytext)
rbind(
  data.frame(imp.pdp$WORLD$imp, country="Global"),
  data.frame(imp.pdp$USA$imp, country="United States"),
  data.frame(imp.pdp$BRA$imp, country="Brazil")) %>% 
  mutate(country = factor(country, levels = c("Global", "United States", "Brazil"))) %>% 
  mutate(country_lab = recode(country, 
                              "Global"="Global (N=122229)\nModel: pca.m.2", 
                              "United States"="United States (N=29803)\nModel: mfpca.m.1", 
                              "Brazil"="Brazil (N=14757)\nModel: avg.m")) %>% 
  mutate(country_lab = factor(country_lab, 
                              levels = c("Global (N=122229)\nModel: pca.m.2", 
                              "United States (N=29803)\nModel: mfpca.m.1", 
                              "Brazil (N=14757)\nModel: avg.m"))) %>% 
  arrange(desc(country_lab), desc(imp)) %>%
  group_by(country_lab) %>% 
  mutate(order = 1:n()) %>% 
  ggplot(., aes(x = reorder_within(predictor, order, country_lab, min),
                y = imp)) + 
  geom_col(fill="grey40",  width = 0.75) +
  scale_x_reordered() +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle = 0), 
        strip.placement = 'outside',
        panel.grid = element_blank())+
  labs(y="Variable importance") +
  facet_grid(country_lab ~ ., 
             scales="free", 
             space = "free", switch = "y")

ggsave(filename = paste0(supp_fig_path, "imp_plots.png"), 
      dpi=300, width =8, height = 6)

# Separate plots
data.frame(imp.pdp$WORLD$imp, country_lab="Global (N=122229)\nModel: pca.m.2") %>% 
  ggplot(., aes(y = reorder(predictor, imp, min),
                x = imp)) + 
  geom_col(fill="grey40",  width = 0.75) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())+
  labs(x="Variable importance") +
  facet_grid(. ~ country_lab)

ggsave(filename = paste0(supp_fig_path, "Importance plots/imp_plots_world.png"), 
       dpi=300, width =8, height = 4)

data.frame(imp.pdp$USA$imp, country_lab="United States (N=29803)\nModel: mfpca.m.1") %>% 
  ggplot(., aes(y = reorder(predictor, imp, min),
                x = imp)) + 
  geom_col(fill="grey40",  width = 0.75) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())+
  labs(x="Variable importance") +
  facet_grid(. ~ country_lab)

ggsave(filename = paste0(supp_fig_path, "Importance plots/imp_plots_usa.png"), 
       dpi=300, width =8, height = 2)

data.frame(imp.pdp$BRA$imp, country_lab="Brazil (N=14757)\nModel: avg.m") %>% 
  ggplot(., aes(y = reorder(predictor, imp, min),
                x = imp)) + 
  geom_col(fill="grey40",  width = 0.75) +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank())+
  labs(x="Variable importance") +
  facet_grid(. ~ country_lab)

ggsave(filename = paste0(supp_fig_path, "Importance plots/imp_plots_bra.png"), 
       dpi=300, width =8, height = 5)

# --------------------------------------------
# --------------------------------------------
# PARTIAL DEPENDENCY PLOTS of the best models 
# Global 
names(imp.pdp$WORLD$mod$variable.importance)
#[1] "irrigated_portion"               
#[2] "PC1_month_max_temp"              
#[3] "PC1_month_min_temp"              
#[4] "PC1_month_et0"                   
#[5] "PC1_month_rad"                   
#[6] "PC1_month_prec"                  
#[7] "PC1_month_vapor_pressure_deficit"
#[8] "PC2_month_max_temp"              
#[9] "PC2_month_min_temp"              
#[10] "PC2_month_et0"                   
#[11] "PC2_month_rad"                   
#[12] "PC2_month_prec"                  
#[13] "PC2_month_vapor_pressure_deficit"

p1 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_max_temp", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly maximum temperature (°C)", 
              ylab = "Predicted yield (tons/hectare)") ; p1

p2 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_min_temp", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly minimum temperature (°C)", 
              ylab = "Predicted yield (tons/hectare)") ; p2

p3 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_prec", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly total precipitations (mm)", 
              ylab = "Predicted yield (tons/hectare)") ; p3

p4 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_rad", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly solar radiations (M2/m2)", 
              ylab = "Predicted yield (tons/hectare)") ; p4

p5 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_et0", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly evapotranspiration of reference (mm/day)", 
              ylab = "Predicted yield (tons/hectare)") ; p5

p6 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC1_month_vapor_pressure_deficit", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "First score derived from principal component analysis\napplied on monthly vapor pressure deficit", 
              ylab = "Predicted yield (tons/hectare)") ; p6

p7 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_max_temp", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly maximum temperature (°C)", 
              ylab = "Predicted yield (tons/hectare)") ; p7

p8 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_min_temp", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly minimum temperature (°C)", 
              ylab = "Predicted yield (tons/hectare)") ; p8

p9 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_prec", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly total precipitations (mm)", 
              ylab = "Predicted yield (tons/hectare)") ; p9

p10 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_rad", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly solar radiations (M2/m2)", 
              ylab = "Predicted yield (tons/hectare)") ; p10

p11 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_et0", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly evapotranspiration of reference (mm/day)", 
              ylab = "Predicted yield (tons/hectare)") ; p11

p12 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "PC2_month_vapor_pressure_deficit", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Second score derived from principal component analysis\napplied on monthly vapor pressure deficit", 
              ylab = "Predicted yield (tons/hectare)") ; p12

p13 <- imp.pdp$WORLD$mod %>% 
  partial(pred.var = "irrigated_portion", train = tab_world) %>% 
  plotPartial(rug=T, 
              train = tab_world,
              lwd = 2, 
              #ylim = c(2,5),
              xlab = "Irrigated portion of soybean in the grid cell ()", 
              ylab = "Predicted yield (tons/hectare)") ; p13

p_scores1 <- plot_grid(p1, p2, p3, p4, p5, p6,
                       nrow=2)

p_scores2 <- plot_grid(p7, p8, p9, p10, p11, p12, 
                       nrow=2)

ggsave(p_scores1, 
       filename = paste0(supp_fig_path, "pdp_world1.png"), 
       dpi=300, width = 15, height = 8)

ggsave(p_scores2, 
       filename = paste0(supp_fig_path, "pdp_world2.png"), 
       dpi=300, width = 15, height = 8)

ggsave(plot_grid(p13), 
       filename = paste0(supp_fig_path, "pdp_world3.png"), 
       dpi=300, width = 6, height = 4)

# USA 
names(imp.pdp$USA$mod$variable.importance)

p1 <- imp.pdp$USA$mod %>% 
   partial(pred.var = "MFPC1_month", train = tab_usa) %>% 
   plotPartial(rug=T, 
               train = tab_usa,
               lwd = 2, 
               #ylim = c(2,5),
               xlab = "First score derived from multivariate functional\nprincipal component analysis applied on monthly averages", 
               ylab = "Predicted yield (tons/hectare)") ; p1

p2 <- imp.pdp$USA$mod %>% 
 partial(pred.var = "irrigated_portion", train = tab_usa) %>% 
 plotPartial(rug=T, 
             train = tab_usa,
             lwd = 2, 
             #ylim = c(2,5),
             xlab = "Irrigated portion of soybean in the grid cell (%)", 
             ylab = "Predicted yield (tons/hectare)") ; p2
 
plot_grid(p1, p2, ncol = 2)
ggsave(filename = paste0(supp_fig_path, "pdp_usa.png"), 
       dpi=300, width = 10, height = 4)

# Brazil 
names(imp.pdp$BRA$mod$variable.importance)
# [1] "irrigated_portion"                "year_max_2m_temperature"         
# [3] "year_min_2m_temperature"          "year_et0"                        
# [5] "year_surface_net_solar_radiation" "year_total_precipitation"        
# [7] "year_vapor_pressure_deficit"     

p1 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_max_2m_temperature", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Maximum temperature (°C)\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p1 

p2 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_min_2m_temperature", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Minimum temperature (°C)\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p2 

p3 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_total_precipitation", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Total precipitations (mm)\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p3

p4 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_surface_net_solar_radiation", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Net surface solar radiations (MJ/m2)\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p4

p5 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_et0", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Evapotranspiration of reference (mm/day)\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p5

p6 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "year_vapor_pressure_deficit", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Vapor pressure deficit\naverage over the growing season", 
              ylab = "Predicted yield (tons/hectare)") ; p6

p7 <- imp.pdp$BRA$mod %>% 
  partial(pred.var = "irrigated_portion", train = tab_bra) %>% 
  plotPartial(rug=T, 
              train = tab_bra,
              lwd = 2, 
              cex=3,
              xlab = "Irrigated portion of soybean in the grid cell (%)", 
              ylab = "Predicted yield (tons/hectare)") ; p7

plot_grid(p1, p2, p3, p4, p5, p6, p7, ncol = 4)
ggsave(filename = paste0(supp_fig_path, "pdp_bra.png"), 
       dpi=300, width =15, height = 8)


# OTHER CORRELATION PLOTS 
# > Correlation between all PCA scores (day) and monthly averages
png(filename = paste0(supp_fig_path, "WORLD_all_PCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_world %>% 
    dplyr::select(starts_with(paste0("PC1_day_", clim.var_abb_i)), 
                  starts_with(paste0("PC2_day_", clim.var_abb_i)), 
                  starts_with(paste0("PC3_day_", clim.var_abb_i)),
                  starts_with(paste0("PC4_day_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("Month 1"=paste0("monthly_", var_i, "_1"), 
                  "Month 2"=paste0("monthly_", var_i, "_2"), 
                  "Month 3"=paste0("monthly_", var_i, "_3"), 
                  "Month 4"=paste0("monthly_", var_i, "_4"), 
                  "Month 5"=paste0("monthly_", var_i, "_5"), 
                  "Month 6"=paste0("monthly_", var_i, "_6"), 
                  "Month 7"=paste0("monthly_", var_i, "_7")) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all FPCA scores and monthly averages
png(filename = paste0(supp_fig_path, "WORLD_all_FPCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_world %>% 
    dplyr::select(starts_with(paste0("FPC1_month_", clim.var_abb_i)), 
                  starts_with(paste0("FPC2_month_", clim.var_abb_i)), 
                  starts_with(paste0("FPC3_month_", clim.var_abb_i)),
                  starts_with(paste0("FPC4_month_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    
    dplyr::rename("FPCA score 1"= 1,
                  "FPCA score 2"= 2,
                  "FPCA score 3"= 3,
                  "FPCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all FPCA scores (day) and monthly averages
png(filename = paste0(supp_fig_path, "WORLD_all_FPCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_world %>% 
    dplyr::select(starts_with(paste0("FPC1_day_", clim.var_abb_i)), 
                  starts_with(paste0("FPC2_day_", clim.var_abb_i)), 
                  starts_with(paste0("FPC3_day_", clim.var_abb_i)),
                  starts_with(paste0("FPC4_day_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    
    dplyr::rename("FPCA score 1"= 1,
                  "FPCA score 2"= 2,
                  "FPCA score 3"= 3,
                  "FPCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all MFPCA scores and monthly averages
png(filename = paste0(supp_fig_path, "WORLD_all_MFPCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  
  tab_world %>% 
    dplyr::select(MFPC1_month, MFPC2_month, MFPC3_month, MFPC4_month, starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("MFPCA score 1"=MFPC1_month,
                  "MFPCA score 2"=MFPC2_month,
                  "MFPCA score 3"=MFPC3_month,
                  "MFPCA score 4"=MFPC4_month,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all MFPCA scores (day) and monthly averages
png(filename = paste0(supp_fig_path, "WORLD_all_MFPCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  
  tab_world %>% 
    dplyr::select(MFPC1_day, MFPC2_day, MFPC3_day, MFPC4_day, 
                  starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("MFPCA score 1"=MFPC1_day,
                  "MFPCA score 2"=MFPC2_day,
                  "MFPCA score 3"=MFPC3_day,
                  "MFPCA score 4"=MFPC4_day,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()
#-----------------------------------------
# USA
# > Correlation between all PCA scores and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "USA_all_PCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_usa %>% 
    dplyr::select(starts_with(paste0("PC1_month_", clim.var_abb_i)), 
                  starts_with(paste0("PC2_month_", clim.var_abb_i)), 
                  starts_with(paste0("PC3_month_", clim.var_abb_i)),
                  starts_with(paste0("PC4_month_", clim.var_abb_i)),
                  starts_with(paste0("PC5_month_", clim.var_abb_i)),
                  starts_with(paste0("PC6_month_", clim.var_abb_i)),
                  starts_with(paste0("PC7_month_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("PCA score 1"= 1,
                  "PCA score 2"= 2,
                  "PCA score 3"= 3,
                  "PCA score 4"= 4,
                  "PCA score 5"= 5,
                  "PCA score 6"= 6,
                  "PCA score 7"= 7,
                  "Month 1"=8, 
                  "Month 2"=9, 
                  "Month 3"=10, 
                  "Month 4"=11, 
                  "Month 5"=12, 
                  "Month 6"=13, 
                  "Month 7"=14) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all PCA scores (day) and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "USA_all_PCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_usa %>% 
    dplyr::select(starts_with(paste0("PC1_day_", clim.var_abb_i)), 
                  starts_with(paste0("PC2_day_", clim.var_abb_i)), 
                  starts_with(paste0("PC3_day_", clim.var_abb_i)),
                  starts_with(paste0("PC4_day_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("PCA score 1"= 1,
                  "PCA score 2"= 2,
                  "PCA score 3"= 3,
                  "PCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all FPCA scores and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "USA_all_FPCA_monthly.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_usa %>% 
    dplyr::select(starts_with(paste0("FPC1_month_", clim.var_abb_i)), 
                  starts_with(paste0("FPC2_month_", clim.var_abb_i)), 
                  starts_with(paste0("FPC3_month_", clim.var_abb_i)),
                  starts_with(paste0("FPC4_month_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    
    dplyr::rename("FPCA score 1"= 1,
                  "FPCA score 2"= 2,
                  "FPCA score 3"= 3,
                  "FPCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all FPCA scores (day) and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "USA_all_FPCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  clim.var_abb_i <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_abb"])
  
  tab_usa %>% 
    dplyr::select(starts_with(paste0("FPC1_day_", clim.var_abb_i)), 
                  starts_with(paste0("FPC2_day_", clim.var_abb_i)), 
                  starts_with(paste0("FPC3_day_", clim.var_abb_i)),
                  starts_with(paste0("FPC4_day_", clim.var_abb_i)),
                  starts_with(paste0("monthly_", var_i))) %>% 
    
    dplyr::rename("FPCA score 1"= 1,
                  "FPCA score 2"= 2,
                  "FPCA score 3"= 3,
                  "FPCA score 4"= 4,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

# > Correlation between all MFPCA scores (day) and monthly averages in the USA 
png(filename = paste0(supp_fig_path, "USA_all_MFPCA_daily.png"), 
    width = 12, height = 8, 
    units = "in", 
    res = 500)

par(mfrow=c(2,3))
for(var_i in unique(vars_names$clim.var))
{
  
  # > Label
  lab <- unique(vars_names[which(vars_names$clim.var == var_i), "clim.var_lab"])
  
  tab_usa %>% 
    dplyr::select(MFPC1_day, MFPC2_day, MFPC3_day, MFPC4_day, 
                  starts_with(paste0("monthly_", var_i))) %>% 
    dplyr::rename("MFPCA score 1"=MFPC1_day,
                  "MFPCA score 2"=MFPC2_day,
                  "MFPCA score 3"=MFPC3_day,
                  "MFPCA score 4"=MFPC4_day,
                  "Month 1"=5, 
                  "Month 2"=6, 
                  "Month 3"=7, 
                  "Month 4"=8, 
                  "Month 5"=9, 
                  "Month 6"=10, 
                  "Month 7"=11) %>% 
    cor(.) %>% 
    corrplot(., method = "square", diag = F, 
             type = "lower", 
             addCoef.col = 'black', tl.cex = 0.7, number.cex = 0.7, tl.col = "black", 
             title = paste0(lab), mar = c(1,0,1,0))
  
}

dev.off()

