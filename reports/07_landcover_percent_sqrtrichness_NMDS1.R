#### Load required libraries ################################################################
library(tidyverse)
library(cowplot) # for visuals
library(ggplot2) # for visuals
library(ggpubr) # for visuals
library(scales) # for visuals
library(ggpmisc) # ggplot extensions, for visuals
library(grid) # for visuals
library(gridExtra) # for visuals
library(car) # Companion to Applied Regression
library(lme4) # Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(effects) # Graphical and tabular effect displays, e.g., of interactions, for various statistical models with linear predictors
library(MuMIn)#AIC, R2
library(dplyr)
library(tidyr)
# load package to visualise table in html format
library(sjPlot) # plot nice tables of model output
library(sjmisc)
library(sjlabelled)
library(performance) # the best r2 for model

library(insight) # Get variance components from random effects models

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified") # if including unspecified/other category

#### load data ################################################################

# run script 05 so data is loaded

taxonomy_Insecta <- read.delim("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/cleaned-data/DK_taxonomy_Insecta.txt")

### Covariation check #########################################

allInsects <- allInsects_totsample # check format for covers, should be 0-100 not 0-1
asvs <- totsample_asvs

#check whether explanatory variables are strongly correlated
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_1000",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_500",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_250",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_50",names(allInsects))])])
#correlations between stops and urban cover...

# change land covers to be 0-100 instead of 0-1
allInsects[, c(30:141,143)] <- allInsects[, c(30:141,143)]*100
allInsects$NMDS1 <- allInsects$NMDS1*100 

### Figure XX: the data - community comp ##########################
# richness
allInsects.long <- allInsects %>% 
  select(NMDS1, Urban_1000, Agriculture_500, Open.uncultivated.land_1000, Wetland_50, Forest_1000) %>% pivot_longer(-c(NMDS1), names_to = "landcover", values_to = "cover")

head(allInsects.long)

allInsects.long$cover_int <- floor((allInsects.long$cover - min(allInsects.long$cover)) / 10) + 1

# calculating the mean community variation across 5 land cover intervals
meancover <- allInsects.long %>%
  dplyr::group_by(cover_int, landcover) %>%
  dplyr::summarise(mean = mean(NMDS1)) %>% ungroup(landcover)

withsd <- meancover %>% dplyr::group_by(landcover) %>%
  dplyr::mutate(Sp_SD = mapply(function(x)sd(mean[-x]), 1:dplyr::n()), se = mapply(function(x)sd(mean[-x]) / sqrt(length(mean[-x])), 1:dplyr::n()))

# calculate the mean variation for each land cover, disregarding each proportional interval
meancover <- withsd %>% dplyr::group_by(landcover) %>% dplyr::summarise(moremean = mean(mean), meansd = mean(Sp_SD), meanse = mean(se)) %>% ungroup(landcover)

meancover$landcover_names <- plyr::mapvalues(meancover$landcover, from = c("Urban_1000", "Agriculture_500", "Open.uncultivated.land_1000", "Wetland_50", "Forest_1000"), to = c("Urban", "Farmland", "Grassland", "Wetland", "Forest"))

mean_variation_cc<- meancover %>% mutate(
  landcover_names = fct_relevel(
    landcover_names,
    "Urban",
    "Farmland",
    "Grassland",
    "Wetland",
    "Forest"
  )
) %>% ggplot(aes(landcover_names, moremean, fill = landcover_names)) + geom_bar(stat = "identity", show.legend = F) + geom_pointrange(aes(x=landcover_names, ymin=moremean-meanse, ymax=moremean+meanse), size = 0.8, colour = "grey27", show.legend = F) + scale_fill_manual(values = landuseCols, labels = c(
  "Urban cover",
  "Farmland cover",
  "Grassland cover",
  "Wetland cover",
  "Forest cover"
)) + labs(x = "\nLand cover", y= "Mean community composition variation\n") + theme_minimal_grid() #+ theme(legend.text = element_text(size = 8), legend.title = element_blank(), legend.position = "bottom") + guides(colour=guide_legend(nrow=1)) # here we haven't removed the grassland outlier at 40% cover, which will be removed in models

figX_all_cc <- allInsects.long %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_500",
    "Open.uncultivated.land_1000",
    "Wetland_50",
    "Forest_500"
  )
) %>% ggplot(aes(cover, NMDS1, colour = landcover)) + geom_point() + scale_colour_manual(values = landuseCols, labels = c(
  "Urban cover (1000 m)",
  "Farmland cover (500 m)",
  "Grassland cover (100 m)",
  "Wetland cover (50 m)",
  "Forest cover (500 m)"
)) + geom_smooth(method=lm, aes(fill = landcover), alpha = 0.1, size =1.5, show.legend = F) + scale_fill_manual(values = landuseCols)  + scale_x_continuous(limits = c(0,100) , labels = function(x)
  paste0(x * 1, "%")) + theme_minimal_grid() + labs(x = "\nProportional land cover", y= "Community composition variation", colour = "Land cover") + theme(plot.subtitle = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 8), legend.position = "bottom", axis.title = element_text(size = 12, face = "bold")) + guides(colour=guide_legend(nrow=1))

#figX_cr <- plot_grid(figX_all_cc, figX_all_r,ncol=1, align = "h", labels = c("A", "B"))
#save_plot("plots/community_comp_richness_variation_all_landcovers.jpg", figX_cr, base_width = 12, base_height = 8, dpi = 400)

#Fig1 <- plot_grid(mean_variation_richness,ncol=1, align = "hv", labels = c("A", "B"))
cowplot::save_plot("plots/Figure1_barplot.jpg", mean_variation_cc, base_width = 12, base_height = 8, dpi = 400)

### Figure XX: the data - richness ##########################

# richness
allInsects.long <- allInsects %>% 
  select(richness.t, Urban_1000, Agriculture_500, Open.uncultivated.land_1000, Wetland_50, Forest_500) %>% pivot_longer(-c(richness.t), names_to = "landcover", values_to = "cover")

head(allInsects.long)

allInsects.long$cover_int <- floor((allInsects.long$cover - min(allInsects.long$cover)) / 10) + 1

# calculating the mean community variation across 5 land cover intervals
meancover <- allInsects.long %>%
  dplyr::group_by(cover_int, landcover) %>%
  dplyr::summarise(mean = mean(richness.t)) %>% ungroup(landcover)

withsd <- meancover %>% dplyr::group_by(landcover) %>%
  dplyr::mutate(Sp_SD = mapply(function(x)sd(mean[-x]), 1:dplyr::n()), se = mapply(function(x)sd(mean[-x]) / sqrt(length(mean[-x])), 1:dplyr::n()))

# calculate the mean variation for each land cover, disregarding each proportional interval
meancover <- withsd %>% dplyr::group_by(landcover) %>% dplyr::summarise(moremean = mean(mean), meansd = mean(Sp_SD), meanse = mean(se)) %>% ungroup(landcover)

meancover$landcover_names <- plyr::mapvalues(meancover$landcover, from = c("Urban_1000", "Agriculture_500", "Open.uncultivated.land_1000", "Wetland_50", "Forest_500"), to = c("Urban", "Farmland", "Grassland", "Wetland", "Forest"))

mean_variation_richness<- meancover %>% mutate(
  landcover_names = fct_relevel(
    landcover_names,
    "Urban",
    "Farmland",
    "Grassland",
    "Wetland",
    "Forest"
  )
) %>% ggplot(aes(landcover_names, moremean, fill = landcover_names)) + geom_bar(stat = "identity", show.legend = F) + geom_pointrange(aes(x=landcover_names, ymin=moremean-meanse, ymax=moremean+meanse), size = 0.8, colour = "grey27", show.legend = F) + scale_fill_manual(values = landuseCols, labels = c(
  "Urban cover",
  "Farmland cover",
  "Grassland cover",
  "Wetland cover",
  "Forest cover"
)) + labs(x = "\nLand cover", y= "Mean insect ASV richness (sqrt-transformed)\n") + theme_minimal_grid() #+ theme(legend.text = element_text(size = 8), legend.title = element_blank(), legend.position = "bottom") + guides(colour=guide_legend(nrow=1)) # here we haven't removed the grassland outlier at 40% cover, which will be removed in models

figX_all_r <- allInsects.long %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_500",
    "Open.uncultivated.land_1000",
    "Wetland_50",
    "Forest_1000"
  )
) %>% ggplot(aes(cover, richness, colour = landcover)) + geom_point() + scale_colour_manual(values = landuseCols, labels = c(
  "Urban cover (1000 m)",
  "Farmland cover (500 m)",
  "Grassland cover (1000 m)",
  "Wetland cover (50 m)",
  "Forest cover (1000 m)"
)) + geom_smooth(method=lm, aes(fill = landcover), alpha = 0.1, size =1.5, show.legend = F) + scale_fill_manual(values = landuseCols)  + scale_x_continuous(limits = c(0,100) , labels = function(x)
  paste0(x * 1, "%")) + theme_minimal_grid() + labs(x = "\nProportional land cover", y= "Richness", colour = "Land cover") + theme(plot.subtitle = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 8), legend.position = "bottom", axis.title = element_text(size = 12, face = "bold")) + guides(colour=guide_legend(nrow=1))

#figX_cr <- plot_grid(figX_all_cc, figX_all_r,ncol=1, align = "h", labels = c("A", "B"))
#save_plot("plots/community_comp_richness_variation_all_landcovers.jpg", figX_cr, base_width = 12, base_height = 8, dpi = 400)

Fig1 <- cowplot::plot_grid(mean_variation_cc, mean_variation_richness,ncol=1, align = "hv", labels = c("A", "B"))
cowplot::save_plot("plots/Figure1_barplot.jpg", Fig1, base_width = 16, base_height = 12)

### Linear Mixed Effects Model: Land covers (Table 1) #################

# jitter x and y slightly - fix later
allInsects$x2 <- allInsects$utm_x + rnorm(length(allInsects$utm_x),0,10)
allInsects$y2 <- allInsects$utm_y + rnorm(length(allInsects$utm_y),0,10)

#full and final model - comcomp
model1 <- lme(NMDS1 ~ Agriculture_500 + 
                Forest_1000 +
                Wetland_50 +
                Urban_1000 +
                Open.uncultivated.land_1000 +
                Time_band + 
                Time_band:cnumberTime +
                Temperature +
                Wind +
                Time_driven +
                cStops + 
                cyDay,
              random=~1|PilotID,
              correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
              data=subset(allInsects, Open.uncultivated.land_1000 < 20))
summary(model1)
model1$apVar
r.squaredGLMM(model1)
#R2m       R2c
#0.3492255 0.458991

#full and final model - richness
allInsects$richness.t <- sqrt(allInsects$richness)
model2 <- lme(richness.t ~ Agriculture_500 + 
                Forest_500 +
                Wetland_50 +
                Urban_1000 +
                Open.uncultivated.land_1000 +
                Time_band + 
                Time_band:cnumberTime + 
                cStops + 
                cyDay,
              random=~1|PilotID,
              correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
              data=subset(allInsects, Open.uncultivated.land_1000 < 20))
summary(model2)
model2$apVar
r.squaredGLMM(model2)
# R2m       R2c
#0.2860156 0.6400077

r2(model1)
r2(model2)

intervals(model2, which = "fixed")

tab_model(model1, p.style = "numeric_stars")
tab_model(model1, model2, dv.labels = c("Community composition", "Richness"), pred.labels = c("Intercept", "Farmland (500 m)", "Forest (1000 m)", "Wetland (50 m)", "Urban (1000 m)", "Grassland (50 m)", "Time band: midday vs evening", "Temperature (20-25 degree celsius)", "Temperature (25-30 degree celsius)", "Wind (light breeze)", "Wind (moderate breeze)", "Sampling duration", "Potential stops", "Day of year (Sampling within the month of June)", "Time within evening
(change in response per minute within time band)
", "Time within midday
(change in response per minute within time band)
"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value", show.icc = FALSE, p.style = "scientific")

#check variance inflation factor
vif(model1)
vif(model2)

get_variance(model1)
get_variance(model2)

### Figure X: effect plots ##########################
# comcomp
gls1.alleffects <- allEffects(model1)
plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model1)
plot(eall.lm1, lines=list(multiline=TRUE))
#plot(predictorEffects(model1, ~ Agriculture_500 + Urban_1000 + Open.uncultivated.land_50 + Forest_1000 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

### comcomp: ggplot effect plot ####
temp <- effectdata$Agriculture_500
temp$landcover <- "Agriculture_500"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_500
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Grassland_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_50
temp$landcover <- "Wetland_50"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_50
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_1000
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
model_time_cc <- lme(NMDS1 ~ Agriculture_500 + 
                       Forest_1000 +
                       Wetland_50 +
                       Urban_1000 +
                       Open.uncultivated.land_1000 +
                       Time_band + 
                       #Time_band:cnumberTime +
                       Temperature +
                       Wind +
                       Time_driven +
                       numberTime +
                       cStops + 
                       cyDay,
                     random=~1|PilotID,
                     correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
                     data=subset(allInsects, Open.uncultivated.land_50 < 20))

gls1.alleffects <- allEffects(model_time_cc)
#plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

temp <- effectdata$numberTime

#temp$landcover <- "Time_band:cnumberTime"
time <- temp

test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot_cc <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_500",
    "Grassland_1000",
    "Wetland_50",
    "Forest_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover (1000 m)",
      "Farmland cover (500 m)",
      "Grassland cover (1000 m)",
      "Wetland cover (50 m)",
      "Forest cover (1000 m)"
    )
  ) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x * 1, "%"), expand = c(0.01,0.12)) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "\nLand cover",
        y = "Community composition\n",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

### comcomp effect plot: time band ###############################

effectplot_time_cc <- time %>% ggplot(aes(x = numberTime, y = fit)) +
  geom_line(size = 2) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  ) + geom_ribbon(
    aes(
      ymin = fit-se,
      ymax = fit+se
    ),
    linetype = 2,
    alpha = 0.2,
    show.legend = F
  ) + labs(
    x = "Sampling time",
    y = "Community composition"
  ) + scale_x_continuous(breaks = c(700, 800, 900, 1000, 1100, 1200), labels = c("12:00", "13:30", "15:00", "17:00", "18:30", "20:00")) 

# richness
gls1.alleffects <- allEffects(model2)
plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model2)
plot(eall.lm1, lines=list(multiline=TRUE))
#plot(predictorEffects(model2, ~ Agriculture_500 + Urban_1000 + Open.uncultivated.land_50 + Forest_1000 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

### richness: ggplot effect plot ####
temp <- effectdata$Agriculture_500
temp$landcover <- "Agriculture_500"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_500
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Grassland_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_50
temp$landcover <- "Wetland_50"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_50
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_500
temp$landcover <- "Forest_500"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_500
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
model_time_r <- lme(richness.t ~ Agriculture_500 + 
                      Forest_500 +
                      Wetland_50 +
                      Urban_1000 +
                      Open.uncultivated.land_1000 +
                      Time_band + 
                      #Time_band:cnumberTime + 
                      cStops + 
                      numberTime +
                      cyDay,
                    random=~1|PilotID,
                    correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
                    data=subset(allInsects, Open.uncultivated.land_50 < 20))

gls1.alleffects <- allEffects(model_time_r)
#plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

temp <- effectdata$numberTime

#temp$landcover <- "Time_band:cnumberTime"
time <- temp

test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot_r <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_500",
    "Grassland_1000",
    "Wetland_50",
    "Forest_500"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover (1000 m)",
      "Farmland cover (500 m)",
      "Grassland cover (1000 m)",
      "Wetland cover (50 m)",
      "Forest cover (500 m)"
    )
  ) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x * 1, "%"), expand = c(0.01,0.12)) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "\nLand cover",
        y = "Insect ASV richness \n(sqrt-transformed)\n",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

effect_plots <- cowplot::plot_grid(effectplot_cc, effectplot_r, ncol = 1, align = "hv", labels = "AUTO")
cowplot::save_plot("plots/effect_plot_richness_t.png", effect_plots, base_width = 12, base_height = 8)

### richness effect plot: time band ###############################

effectplot_time_r <- time %>% ggplot(aes(x = numberTime, y = fit)) +
  geom_line(size = 2) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    axis.text = element_text(size = 8)
  ) + geom_ribbon(
    aes(
      ymin = fit-se,
      ymax = fit+se
    ),
    linetype = 2,
    alpha = 0.2,
    show.legend = F
  ) + labs(
    x = "Sampling time",
    y = "Insect ASV richness (sqrt-transformed)"
  ) + scale_x_continuous(breaks = c(700, 800, 900, 1000, 1100, 1200), labels = c("12:00", "13:30", "15:00", "17:00", "18:30", "20:00")) 

effectplot_time <- cowplot::plot_grid(effectplot_time_cc, effectplot_time_r, align = "hv", labels = "AUTO")
cowplot::save_plot("plots/effect_plot_time_notimeband.png", effectplot_time, base_width = 14, base_height = 7)

### Test timeband interactions  ################################

midday_cc <- allInsects %>% filter(Time_band == "midday")
evening_cc <- allInsects %>% filter(Time_band == "evening")

### comcomp time band ########
model_midday_cc <- lmer(NMDS1 ~ Agriculture_500 +
                         Forest_1000 +
                         Wetland_50 +
                         Urban_1000 +
                         Open.uncultivated.land_1000 +
                         cnumberTime +
                         Temperature +
                         Wind +
                         Time_driven +
                         cStops +
                         cyDay + 
                         (1|PilotID), data=subset(midday_cc, Open.uncultivated.land_1000 < 20))
summary(model_midday_cc)
r.squaredGLMM(model_midday_cc)
#            R2m       R2c
#[1,] 0.1136346 0.8820097

model_evening_cc <- lmer(NMDS1 ~ Agriculture_500 +
                          Forest_1000 +
                          Wetland_50 +
                          Urban_1000 +
                          Open.uncultivated.land_1000 +
                          cnumberTime +
                          Temperature +
                          Wind +
                          Time_driven +
                          cStops +
                          cyDay + 
                          (1|PilotID), data=subset(evening_cc, Open.uncultivated.land_1000 < 20))
summary(model_evening_cc)
r.squaredGLMM(model_evening_cc)
#           R2m       R2c
#[1,] 0.3824111 0.5388213

### richness time band ####### 

model_midday_r <- lmer(richness.t ~ Agriculture_500 +
                         Forest_500 +
                         Wetland_50 +
                         Urban_1000 +
                         Open.uncultivated.land_1000 +
                         cnumberTime +
                         cStops +
                         cyDay + 
                         (1|PilotID), data=subset(midday_cc, Open.uncultivated.land_1000 < 20))
summary(model_midday_r)
r.squaredGLMM(model_midday_r)
#R2m       R2c
#[1,]0.2195963 0.6446386

model_evening_r <- lmer(richness.t ~ Agriculture_500 +
                          Forest_500 +
                          Wetland_50 +
                          Urban_1000 +
                          Open.uncultivated.land_1000 +
                          cnumberTime +
                          cStops +
                          cyDay + 
                          (1|PilotID), data=subset(evening_cc, Open.uncultivated.land_1000 < 20))
summary(model_evening_r)
r.squaredGLMM(model_evening_r)
#R2m       R2c
#[1,0.2442609 0.6340512

tab_model(model_midday_cc, model_evening_cc, model_midday_r, model_evening_r, dv.labels = c("Community composition (midday)", "Community composition (evening)", "Richness (midday)", "Richness (evening)"), pred.labels = c("Intercept", "Farmland", "Forest (1000 m)", "Wetland", "Urban", "Grassland", "Sampling time", "Temperature (20-25)", "Temperature (25-30)", "Wind (light breeze)", "Wind (moderate breeze)", "Sampling duration", "Potential stops", "Day of year", "Forest (500 m)"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")

### pariwise test between slopes (general linear hypothesis test) ##########
library(multcomp)
coef(model1)
length(coef(model1))
K <- diag(16)
rownames(K) <- names(coef(model1))
model1.ht <- glht(model1, linfct = K)
summary(model1.ht)

# pairwise comparison to farmland
pair.ht <- glht(model1, linfct = c("Forest_1000 - Agriculture_500 = 0", "Wetland_50 - Agriculture_500 = 0", "Urban_1000 - Agriculture_500 = 0", "Open.uncultivated.land_1000 - Agriculture_500 = 0"))
confint(pair.ht)
summary(pair.ht)

# richness
coef(model2)
length(coef(model2))
K <- diag(11)
rownames(K) <- names(coef(model2))
model2.ht <- glht(model2, linfct = K)
summary(model2.ht)

# pairwise comparison to farmland
pair.ht <- glht(model2, linfct = c("Forest_500 - Agriculture_500 = 0", "Wetland_50 - Agriculture_500 = 0", "Urban_1000 - Agriculture_500 = 0", "Open.uncultivated.land_1000 - Agriculture_500 = 0"))
confint(pair.ht)
summary(pair.ht)

### other glht tests between land covers #########

# grassland comcomp and species richness
pair.ht <- glht(model1, linfct = c("Forest_1000 - Open.uncultivated.land_1000 = 0", "Wetland_50 - Open.uncultivated.land_1000 = 0", "Urban_1000 - Open.uncultivated.land_1000 = 0"))
confint(pair.ht)
summary(pair.ht) 

pair.ht <- glht(model2, linfct = c("Forest_500 - Open.uncultivated.land_1000 = 0", "Wetland_50 - Open.uncultivated.land_1000 = 0", "Urban_1000 - Open.uncultivated.land_1000 = 0"))
confint(pair.ht)
summary(pair.ht) # grassland has higher richness (non-significant for wetland comparison)

# urban comcomp and richness
pair.ht <- glht(model1, linfct = c("Forest_1000 - Urban_1000 = 0", "Wetland_50 - Urban_1000 = 0", "Open.uncultivated.land_1000 - Urban_1000 = 0"))
confint(pair.ht)
summary(pair.ht) # different species composition in urban areas compared to forest and farmland, and trend for grassland

pair.ht <- glht(model2, linfct = c("Forest_500 - Urban_1000 = 0", "Wetland_50 - Urban_1000 = 0"))
confint(pair.ht)
summary(pair.ht) # significantly higher richness in other land covers compared to urban, although non significant for wetland

### evenness ########
#Evenness is a measure of how homogeneous or even a community or ecosystem is in terms of the abundances of its species. A community in which all species are equally common is considered even and has a high degree of evenness.

# Pilou evenness (J)	compares the actual diversity value (such as the Shannon-Wiener Index, H') to the maximum possible diversity value (when all species are equally common, Hmax=ln s where S is the total number of species).

# Pilou evenness (J) is constrained between 0 and 1.0 and the more variation in abundances between different taxa within the community, the lower J. Unfortunately, Pilou's J is highly dependent on sample size (since S - the estimated number of species is dependent on sampling effort) and is also highly sensitive to rare taxa. 
tasvs <- t(asvs) # samples as rows, taxa as columns
tasvs <- decostand(tasvs, method = "pa") # transform to presence absence
tpa <- as.data.frame(tasvs)
tpa$richness <- rowSums(tpa) # get richness for each sample

richnessdata <- tpa %>% rownames_to_column(var = "PCRID") %>% select(PCRID, richness)

shannon <- diversity(tpa[-1], "shannon")
pilous_evenness <- shannon/log(specnumber(tpa))

shannon <- as.data.frame(shannon) %>% rownames_to_column(var = "PCRID")
pileous <- as.data.frame(pilous_evenness) %>% rownames_to_column(var = "PCRID")
data <- merge(allInsects, shannon, by.y = "PCRID", by.x = "SampleID")
data <- merge(data, pileous, by.y = "PCRID", by.x = "SampleID")

hist(data$shannon)
hist(data$pilous_evenness)
qqnorm(data$shannon)
qqnorm(data$pilous_evenness)

#full and final model - richness
modelpil <- lme(pilous_evenness ~ Agriculture_500 + 
                Forest_500 +
                Wetland_50 +
                Urban_1000 +
                Open.uncultivated.land_1000 +
                Time_band + 
                Time_band:cnumberTime + 
                cStops + 
                cyDay,
              random=~1|PilotID,
              correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
              data=subset(data, Open.uncultivated.land_1000 < 20))
summary(modelpil)
modelpil$apVar
r.squaredGLMM(modelpil)
# R2m       R2c
#0.1733493  0.4528909

r2(modelpil)
qqnorm(resid(modelpil))

intervals(modelpil, which = "fixed")

tab_model(modelpil, p.style = "numeric_stars", dv.labels = "Evenness", pred.labels = c("Intercept", "Farmland (500 m)", "Forest (1000 m)", "Wetland (50 m)", "Urban (1000 m)", "Grassland (50 m)", "Time band: midday vs evening", "Potential stops", "Day of year (Sampling within the month of June)", "Time within evening
(change in response per minute within time band)
", "Time within midday
(change in response per minute within time band)
"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value", show.icc = FALSE, digits = 4)

#check variance inflation factor
vif(modelpil)
get_variance(modelpil)

gls1.alleffects <- allEffects(modelpil)
plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(modelpil)
plot(eall.lm1, lines=list(multiline=TRUE))
#plot(predictorEffects(model2, ~ Agriculture_500 + Urban_1000 + Open.uncultivated.land_50 + Forest_1000 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

### richness: ggplot effect plot ####
temp <- effectdata$Agriculture_500
temp$landcover <- "Agriculture_500"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_500
  )%>% select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Grassland_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_1000
  ) %>% select(landcover, propcover, fit, se, lower, upper)

# Wetland
temp <- effectdata$Wetland_50
temp$landcover <- "Wetland_50"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_50
  )%>% select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_500
temp$landcover <- "Forest_500"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_500
  ) %>% select(landcover, propcover, fit, se, lower, upper)

test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot_e <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_500",
    "Grassland_1000",
    "Wetland_50",
    "Forest_500"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover (1000 m)",
      "Farmland cover (500 m)",
      "Grassland cover (1000 m)",
      "Wetland cover (50 m)",
      "Forest cover (500 m)"
    )
  ) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom"
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x * 1, "%"), expand = c(0.01,0.12)) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "\nLand cover",
        y = "Evenness",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)


effect_plots <- cowplot::plot_grid(effectplot_cc, effectplot_r, effectplot_e, ncol = 1, align = "hv", labels = "AUTO")
cowplot::save_plot("plots/effect_plots.png", effect_plots, base_width = 12, base_height = 14)

### pairwise comparison evenness #######
library(multcomp)
coef(modelpil)
length(coef(modelpil))
K <- diag(11)
rownames(K) <- names(coef(modelpil))
model2.ht <- glht(modelpil, linfct = K)
summary(model2.ht)

# pairwise comparison to farmland
pair.ht <- glht(modelpil, linfct = c("Forest_500 - Agriculture_500 = 0", "Wetland_50 - Agriculture_500 = 0", "Urban_1000 - Agriculture_500 = 0", "Open.uncultivated.land_1000 - Agriculture_500 = 0"))
confint(pair.ht)
summary(pair.ht)

pair.ht <- glht(modelpil, linfct = c("Forest_500 - Open.uncultivated.land_1000 = 0", "Wetland_50 - Open.uncultivated.land_1000 = 0", "Urban_1000 - Open.uncultivated.land_1000 = 0"))
confint(pair.ht)
summary(pair.ht) 

# urban comcomp and richness
pair.ht <- glht(modelpil, linfct = c("Urban_1000 - Forest_500 = 0", "Wetland_50 - Urban_1000 = 0", "Open.uncultivated.land_1000 - Urban_1000 = 0"))
confint(pair.ht)
summary(pair.ht)
