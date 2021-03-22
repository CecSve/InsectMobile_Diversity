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
allInsects_totasvs <- read.csv("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/cleaned-data/allInsects_totasvs.txt", row.names=1, sep="")

allInsects_totsample <- read.csv("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/cleaned-data/allInsects_totsample.txt", sep="")

taxonomy_Insecta <- read.delim("H:/Documents/Insektmobilen/Analysis/InsectMobile_Diversity/cleaned-data/DK_taxonomy_Insecta.txt")

### Covariation check #########################################

allInsects <- allInsects_totsample
asvs <- allInsects_totasvs

#check whether explanatory variables are strongly correlated
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_1000",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_500",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_250",names(allInsects))])])
cor(allInsects[,c("cStops","cTL",names(allInsects)[grepl("_50",names(allInsects))])])
#correlations between stops and urban cover...

### Figure XX: the data - richness ##########################

# richness
# alternative to figure X
allInsects.long <- allInsects %>% 
  select(richness, Urban_1000, Agriculture_500, Open.uncultivated.land_1000, Wetland_50, Forest_1000) %>% pivot_longer(-c(richness), names_to = "landcover", values_to = "cover")

head(allInsects.long)

allInsects.long$cover_int <- floor((allInsects.long$cover - min(allInsects.long$cover)) / 0.10) + 1

# calculating the mean community variation across 5 land cover intervals
meancover <- allInsects.long %>%
  dplyr::group_by(cover_int, landcover) %>%
  dplyr::summarise(mean = mean(richness)) %>% ungroup(landcover)

withsd <- meancover %>% dplyr::group_by(landcover) %>%
  dplyr::mutate(Sp_SD = mapply(function(x)sd(mean[-x]), 1:dplyr::n()), se = mapply(function(x)sd(mean[-x]) / sqrt(length(mean[-x])), 1:dplyr::n()))
# calculate the mean variation for each land cover, disregarding each proportional interval
meancover <- withsd %>% dplyr::group_by(landcover) %>% dplyr::summarise(moremean = mean(mean), meansd = mean(Sp_SD), meanse = mean(se)) %>% ungroup(landcover)

meancover$landcover_names <- plyr::mapvalues(meancover$landcover, from = c("Urban_1000", "Agriculture_500", "Open.uncultivated.land_1000", "Wetland_50", "Forest_1000"), to = c("Urban", "Farmland", "Grassland", "Wetland", "Forest"))

mean_variation_richness<- meancover %>% mutate(
  landcover_names = fct_relevel(
    landcover_names,
    "Urban",
    "Farmland",
    "Grassland",
    "Wetland",
    "Forest"
  )
) %>% ggplot(aes(landcover_names, moremean, fill = landcover_names)) + geom_bar(stat = "identity", show.legend = F) + geom_pointrange(aes(x=landcover_names, ymin=moremean-meanse, ymax=moremean+meanse), width = 0.5, size = 0.8, colour = "grey27", show.legend = F) + scale_fill_manual(values = landuseCols, labels = c(
  "Urban cover",
  "Farmland cover",
  "Grassland cover",
  "Wetland cover",
  "Forest cover"
)) + labs(x = "\nLand cover", y= "Mean richness\n") + theme_minimal_grid() #+ theme(legend.text = element_text(size = 8), legend.title = element_blank(), legend.position = "bottom") + guides(colour=guide_legend(nrow=1)) # here we haven't removed the grassland outlier at 40% cover, which will be removed in models

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
)) + geom_smooth(method=lm, aes(fill = landcover), alpha = 0.1, size =1.5, show.legend = F) + scale_fill_manual(values = landuseCols)  + scale_x_continuous(limits = c(0,1) , labels = function(x)
  paste0(x * 100, "%")) + theme_minimal_grid() + labs(x = "\nProportional land cover", y= "Richness", colour = "Land cover") + theme(plot.subtitle = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 8), legend.position = "bottom", axis.title = element_text(size = 12, face = "bold")) + guides(colour=guide_legend(nrow=1))

#figX_cr <- plot_grid(figX_all_cc, figX_all_r,ncol=1, align = "h", labels = c("A", "B"))
#save_plot("plots/community_comp_richness_variation_all_landcovers.jpg", figX_cr, base_width = 12, base_height = 8, dpi = 400)

#Fig1 <- plot_grid(mean_variation_richness,ncol=1, align = "hv", labels = c("A", "B"))
cowplot::save_plot("plots/Figure1_barplot.jpg", mean_variation_richness, base_width = 12, base_height = 8, dpi = 400)

### Linear Mixed Effects Model: Land covers (Table 1) #################
# used cStops instead of cTL for DK data

# for DK jitter x and y slightly - fix later
allInsects$x2 <- allInsects$utm_x + rnorm(length(allInsects$utm_x),0,10)
allInsects$y2 <- allInsects$utm_y + rnorm(length(allInsects$utm_y),0,10)

#full and final model - richness
model2 <- lme((richness*100) ~ Agriculture_500 + 
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
              random=~1|PilotID/RouteID,
              data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(model2)
model2$apVar
r.squaredGLMM(model2)
#           R2m       R2c
#[1,] 0.2899174   0.641468

r2(model2)

tab_model(model2, collapse.ci = TRUE, dv.labels = c("Richness"), pred.labels = c("Intercept", "Farmland (500 m)", "Forest (1000 m)", "Wetland (50 m)", "Urban (1000 m)", "Grassland (1000 m)", "Time band: midday vs evening", "Temperature (20-25 degree celsius)", "Temperature (25-30 degree celsius)", "Wind (light breeze)", "Wind (moderate breeze)", "Sampling time", "Potential stops", "Day of year", "Time within evening
(change in response per minute within time band)
", "Time within midday
(change in response per minute within time band)
"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")

#check variance inflation factor
vif(model2)

get_variance(model2)

### Figure X: effect plots ##########################
gls1.alleffects <- allEffects(model2)
plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model2)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(model2, ~ Agriculture_500 + Urban_1000 + Open.uncultivated.land_1000 + Forest_1000 + Wetland_50 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

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
temp <- effectdata$Forest_1000
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_1000
  ) %>% select(landcover, propcover, fit, se, lower, upper)

# Timeband
model_time_r <- lme((richness*100) ~ Agriculture_500 + 
                      Forest_1000 +
                      Wetland_50 +
                      Urban_1000 +
                      Open.uncultivated.land_1000 +
                      #Heathland_500 +
                      Time_band + 
                      #Time_band:cnumberTime + 
                      Temperature +
                      Wind +
                      Time_driven +
                      numberTime +
                      cStops + 
                      cyDay,
                    random=~1|PilotID/RouteID,
                    data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))

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
    limits = c(0, 1), expand=c(0,0.015),
    labels = function(x)
      paste0(x * 100, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "Land cover",
        y = "Richness",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

#effect_plots <- cowplot::plot_grid(effectplot_cc, effectplot_r, ncol = 1, align = "hv")
cowplot::save_plot("plots/effect_plot_richness.png", effectplot_r, base_width = 12, base_height = 8)

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
    y = "Richness"
  ) + scale_x_continuous(breaks = c(700, 800, 900, 1000, 1100, 1200), labels = c("12:00", "13:30", "15:00", "17:00", "18:30", "20:00")) 

#effectplot_time <- cowplot::plot_grid(effectplot_time_cc, effectplot_time_r, align = "hv", labels = "AUTO")
cowplot::save_plot("plots/effect_plot_time_notimeband_richness.png", effectplot_time_r, base_width = 14, base_height = 7)

### Test timeband interactions  ################################

midday_cc <- allInsects %>% filter(Time_band == "midday")
evening_cc <- allInsects %>% filter(Time_band == "evening")

### richness ####### examine whether the test

model_midday_r <- lmer((richness*100) ~ Agriculture_500 +
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
                         (1|PilotID), data=subset(midday_cc, Open.uncultivated.land_1000 < 0.2))
summary(model_midday_r)
r.squaredGLMM(model_midday_r)

model_evening_r <- lmer((richness*100) ~ Agriculture_500 +
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
                          (1|PilotID), data=subset(evening_cc, Open.uncultivated.land_1000 < 0.2))
summary(model_evening_r)
r.squaredGLMM(model_evening_r)

tab_model(model_midday_r, model_evening_r, collapse.ci = TRUE, dv.labels = c("Richness (midday)", "Richness (evening)"), pred.labels = c("Intercept", "Farmland", "Forest", "Wetland", "Urban", "Grassland", "Sampling time", "Temperature (20-25)", "Temperature (25-30)", "Wind (light breeze)", "Wind (moderate breeze)", "Sampling duration", "Potential stops", "Day of year"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")
