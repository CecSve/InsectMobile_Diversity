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

### Figure XX: the data - community composition ##########################

# community composition
allInsects.long <- allInsects %>% 
  select(distdif, Urban_1000, Agriculture_50, Open.uncultivated.land_1000, Wetland_250, Forest_1000) %>% pivot_longer(-c(distdif), names_to = "landcover", values_to = "cover")

head(allInsects.long)

allInsects.long$cover_int <- floor((allInsects.long$cover - min(allInsects.long$cover)) / 0.10) + 1

# calculating the mean community variation across 5 land cover intervals
meancover <- allInsects.long %>%
  dplyr::group_by(cover_int, landcover) %>%
  dplyr::summarise(mean = mean(distdif)) %>% ungroup(landcover)

withsd <- meancover %>% dplyr::group_by(landcover) %>%
  dplyr::mutate(Sp_SD = mapply(function(x)sd(mean[-x]), 1:dplyr::n()), se = mapply(function(x)sd(mean[-x]) / sqrt(length(mean[-x])), 1:dplyr::n()))

# calculate the mean variation for each land cover, disregarding each proportional interval
meancover <- withsd %>% dplyr::group_by(landcover) %>% dplyr::summarise(moremean = mean(mean), meansd = mean(Sp_SD), meanse = mean(se)) %>% ungroup(landcover)

meancover$landcover_names <- plyr::mapvalues(meancover$landcover, from = c("Urban_1000", "Agriculture_50", "Open.uncultivated.land_1000", "Wetland_250", "Forest_1000"), to = c("Urban", "Farmland", "Grassland", "Wetland", "Forest"))

ggplot(data = meancover,
       aes(x = landcover_names, y = moremean))+
  geom_bar(stat = "identity")

mean_variation_comcomp<- meancover %>% mutate(
  landcover_names = fct_relevel(
    landcover_names,
    "Urban",
    "Farmland",
    "Grassland",
    "Wetland",
    "Forest"
  )
) %>% ggplot(aes(landcover_names, moremean, fill = landcover_names)) + geom_bar(stat = "identity", show.legend = F) + geom_pointrange(aes(x=landcover_names, ymin=moremean-meanse, ymax=moremean+meanse), width = 0.5, size = 0.8, colour = "grey27", show.legend = F) + scale_fill_manual(values = landuseCols, labels = c(
  "Urban cover (1000 m)",
  "Farmland cover (50 m)",
  "Grassland cover (1000 m)",
  "Wetland cover (250 m)",
  "Forest cover (1000 m)"
)) + labs(x = "\nLand cover", y= "Mean variation in \ncommunity composition\n") + theme_minimal_grid() # here we haven't removed the grassland outlier at 40% cover, which will be removed in models

figX_all_cc <- allInsects.long %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_50",
    "Open.uncultivated.land_1000",
    "Wetland_250",
    "Forest_1000"
  )
) %>% ggplot(aes(cover, distdif, colour = landcover)) + geom_point() + scale_colour_manual(values = landuseCols, labels = c(
  "Urban cover (1000 m)",
  "Farmland cover (50 m)",
  "Grassland cover (1000 m)",
  "Wetland cover (250 m)",
  "Forest cover (1000 m)"
)) + geom_smooth(method=lm, aes(fill = landcover), alpha = 0.1, size =1.5, show.legend = F) + scale_fill_manual(values = landuseCols)  + scale_x_continuous(limits = c(0,1) , labels = function(x)
  paste0(x * 100, "%")) + theme_minimal_grid() + labs(x = "\nProportional land cover", y= "Community composition", colour = "Land cover") + theme(plot.subtitle = element_text(size = 16, face = "bold"),legend.title = element_blank(), legend.text = element_text(size = 8), legend.position = "bottom", axis.title = element_text(size = 12, face = "bold")) + guides(colour=guide_legend(nrow=1))

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

figX_cr <- plot_grid(figX_all_cc, figX_all_r,ncol=1, align = "h", labels = c("A", "B"))
save_plot("plots/community_comp_richness_variation_all_landcovers.jpg", figX_cr, base_width = 12, base_height = 8, dpi = 400)

Fig1 <- plot_grid(mean_variation_comcomp, mean_variation_richness,ncol=1, align = "hv", labels = c("A", "B"))
save_plot("plots/Figure1_barplots.jpg", Fig1, base_width = 12, base_height = 8, dpi = 400)

### Pie chart: community composition #####################################
routeMeans <- allInsects %>% 
  group_by(RouteID_JB) %>%
  dplyr::summarise(meanComposition = mean(distdif))

allInsects <- inner_join(allInsects,routeMeans,by="RouteID_JB")

#remove duplicates
allInsects_pie <- allInsects %>%
  select(RouteID_JB,meanComposition,Agriculture_1000,
         Forest_1000,Open.uncultivated.land_1000,
         Urban_1000,Wetland_1000) %>%
  distinct()

#fill in missing column
allInsects_pie$totalLand <- apply(allInsects_pie[,3:7],1,sum)
allInsects_pie$Other_1000 <- 1-allInsects_pie$totalLand

#divide up biomass into quantiles
allInsects_pie$CompositionCats <- cut_number(allInsects_pie$meanComposition,n=5)

#mean land cover per biomass cats
allInsects_cat <- allInsects_pie %>%
  group_by(CompositionCats) %>%
  dplyr::summarise(Agriculture_1000 = mean(Agriculture_1000),
                   Forest_1000 = mean(Forest_1000),
                   Open.uncultivated.land_1000 = mean(Open.uncultivated.land_1000),
                   Urban_1000 = mean(Urban_1000),
                   Wetland_1000 = mean(Wetland_1000),
                   Other_1000 = mean(Other_1000))

#plot data by biomass categories
allInsects_melt <- tidyr::gather(allInsects_cat, Land_cover, value, -CompositionCats)

#plot
allInsects_melt <- allInsects_melt %>% mutate(
  Land_cover = fct_relevel(
    Land_cover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated.land_1000",
    "Wetland_1000",
    "Forest_1000",
    "Other_1000"))  

composition.labs <- c("[-0.291,-0.0876]"="High", "(-0.0876,-0.0332]"="Medium", "(-0.0332,0.00841]"="Low", "(0.00841,0.0639]"="Medium", "(0.0639,0.47]"="High")
#names(biomass.labs) <- c("3-48.8 mg", "48.8-84.2 mg", "84.2-135 mg", "135-262 mg", "262-4360 mg")

fig_pie_cc <- ggplot(allInsects_melt,aes(x="",y=value,fill=Land_cover, order = Land_cover))+
  geom_bar(stat="identity")+
  facet_grid(~CompositionCats, labeller = labeller(CompositionCats=composition.labs))+
  coord_polar("y")+
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position="none", legend.title = element_blank(), plot.subtitle = element_text(size = 16, face = "bold"), legend.text = element_text(size = 8), strip.text.x = element_text(size = 12)) + scale_fill_manual(values = landuseCols, name = "Land cover", labels = c("Urban", "Farmland", "Grassland", "Wetland", "Forest", "Other")) +guides(fill=guide_legend(nrow=1,byrow=TRUE)) 

save_plot("plots/Landcover_communitycomp_proportions.png", fig_pie_cc, base_width = 12, base_height = 6)

### Pie chart: richness #####################################
routeMeans <- allInsects %>% 
  group_by(RouteID_JB) %>%
  dplyr::summarise(meanRichness = mean(richness))

allInsects <- inner_join(allInsects,routeMeans,by="RouteID_JB")

#remove duplicates
allInsects_pie <- allInsects %>%
  select(RouteID_JB,meanRichness,Agriculture_1000,
         Forest_1000,Open.uncultivated.land_1000,
         Urban_1000,Wetland_1000) %>%
  distinct()

#fill in missing column
allInsects_pie$totalLand <- apply(allInsects_pie[,3:7],1,sum)
allInsects_pie$Other_1000 <- 1-allInsects_pie$totalLand

#divide up biomass into quantiles
allInsects_pie$RichnessCats <- cut_number(allInsects_pie$meanRichness,n=5)

#mean land cover per biomass cats
allInsects_cat <- allInsects_pie %>%
  group_by(RichnessCats) %>%
  dplyr::summarise(Agriculture_1000 = mean(Agriculture_1000),
                   Forest_1000 = mean(Forest_1000),
                   Open.uncultivated.land_1000 = mean(Open.uncultivated.land_1000),
                   Urban_1000 = mean(Urban_1000),
                   Wetland_1000 = mean(Wetland_1000),
                   Other_1000 = mean(Other_1000))

#plot data by biomass categories
allInsects_melt <- tidyr::gather(allInsects_cat, Land_cover, value, -RichnessCats)

#plot
allInsects_melt <- allInsects_melt %>% mutate(
  Land_cover = fct_relevel(
    Land_cover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated.land_1000",
    "Wetland_1000",
    "Forest_1000",
    "Other_1000"))  

richness.labs <- c("[13.5,72.5]"="13.5-72.5", "(72.5,110]"="72.5-110", "(110,143]"="110-143", "(143,182]"="143-182", "(182,445]"="182-445")
#names(biomass.labs) <- c("3-48.8 mg", "48.8-84.2 mg", "84.2-135 mg", "135-262 mg", "262-4360 mg")

fig_pie_r <- ggplot(allInsects_melt,aes(x="",y=value,fill=Land_cover, order = Land_cover))+
  geom_bar(stat="identity")+
  facet_grid(~RichnessCats, labeller = labeller(RichnessCats=richness.labs))+
  coord_polar("y")+
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position="none", legend.title = element_blank(), plot.subtitle = element_text(size = 16, face = "bold"), legend.text = element_text(size = 8), strip.text.x = element_text(size = 12)) + scale_fill_manual(values = landuseCols, name = "Land cover", labels = c("Urban", "Farmland", "Grassland", "Wetland", "Forest", "Other")) +guides(fill=guide_legend(nrow=1,byrow=TRUE)) 

save_plot("plots/Landcover_richness_proportions.png", fig_pie_r, base_width = 12, base_height = 6)

pie_all <- plot_grid(fig_pie_cc, fig_pie_r, nrow = 1, labels = c("C", "D"))

Fig1 <- plot_grid(figX_cr, pie_all, nrow = 2, rel_heights = c(5,1))
save_plot("plots/Fig1.png", Fig1, base_width = 18, base_height = 12)

### Linear Mixed Effects Model: Land covers (Table 1) #################
# used cStops instead of cTL for DK data

# for DK jitter x and y slightly - fix later
allInsects$x2 <- allInsects$utm_x + rnorm(length(allInsects$utm_x),0,10)
allInsects$y2 <- allInsects$utm_y + rnorm(length(allInsects$utm_y),0,10)

#full and final model - community composition *100 to get 1% change
model1 <- lme((distdif*100) ~ Agriculture_50 +
      Forest_1000 +
      Wetland_250 +
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
    correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
    data=subset(allInsects, Open.uncultivated.land_1000 < 0.2), na.action = na.omit)

# data=subset(allInsects, Open.uncultivated.land_1000 < 0.2)
summary(model1)
model1$apVar
r.squaredGLMM(model1)
#           R2m       R2c
#[1,] 0.5824744  0.738737

#full and final model - richness
model2 <- lme((richness*100) ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              #Heathland_500 +
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

tab_model(model1, model2, collapse.ci = TRUE, dv.labels = c("Community composition", "Richness"), pred.labels = c("Intercept", "Farmland (50 m)", "Forest (1000 m)", "Wetland (250 m)", "Urban (1000 m)", "Grassland (1000 m)", "Time band: midday vs evening", "Temperature (20-25 degree celsius)", "Temperature (25-30 degree celsius)", "Wind (light breeze)", "Wind (moderate breeze)", "Sampling time", "Potential stops", "Day of year", "Time within evening
(change in response per minute within time band)
", "Time within midday
(change in response per minute within time band)
", "Farmland (500 m)", "Wetland (50 m)"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")

#check variance inflation factor
vif(model1)
vif(model2)

get_variance(model1)
get_variance(model2)

### AIC check ##############################################

# can't be performed on model 1 due to the spatial autocorrelation term, but it is  bit redundant anyways, because we have already made model selection. This is a bit more simple than the way we did it before, though

options(na.action = "na.fail")
dd <- dredge(model2)
subset(dd, delta < 2)

# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)

# Model average models with delta AICc < 4
#model.avg(dd, subset = delta < 2)

# best model with AIC 3858.909
lme1000 <- lme(richness ~ Agriculture_500 + 
                 Forest_1000 +
                 Wetland_50 +
                 Urban_1000 +
                 Open.uncultivated.land_1000 +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
# data=subset(allInsects, Open.uncultivated.land_1000 < 0.2)
summary(lme1000)
r.squaredGLMM(lme1000)

### Figure X: effect plots ##########################
gls1.alleffects <- allEffects(model1)
plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model1)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(model1, ~ Agriculture_50 + Urban_1000 + Open.uncultivated.land_1000 + Forest_1000 + Wetland_250 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

### community comp: ggplot effect plot ####
temp <- effectdata$Agriculture_50
temp$landcover <- "Agriculture_50"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_50
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
temp <- effectdata$Wetland_250
temp$landcover <- "Wetland_250"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_250
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Forest
temp <- effectdata$Forest_1000
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_1000
  ) %>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
model_time_cc <- lme((distdif*100) ~ Agriculture_50 +
                Forest_1000 +
                Wetland_250 +
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
              random=~1|PilotID/RouteID,
              correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
              data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))

gls1.alleffects <- allEffects(model_time_cc)
#plot(gls1.alleffects, 'Urban_1000', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

temp <- effectdata$numberTime
#temp$Time_band <- effectdata$Time_band
time <- temp

test <- rbind(urb, farm, grass, wet, forest)

# Visualization
effectplot_cc <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_50",
    "Grassland_1000",
    "Wetland_250",
    "Forest_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover (1000 m)",
      "Farmland cover (50 m)",
      "Grassland cover (1000 m)",
      "Wetland cover (250 m)",
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
        y = "Variation in community composition \n(+/- different from zero)",
        subtitle = "A: Community composition",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

### distdif effect plot: time band ###############################

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
        subtitle = "B: Richness",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

effect_plots <- cowplot::plot_grid(effectplot_cc, effectplot_r, ncol = 1, align = "hv")
cowplot::save_plot("plots/effect_plots.png", effect_plots, base_width = 11, base_height = 9)

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

effectplot_time <- cowplot::plot_grid(effectplot_time_cc, effectplot_time_r, align = "hv", labels = "AUTO")
cowplot::save_plot("plots/effect_plot_time_notimeband.png", effectplot_time, base_width = 14, base_height = 7)

### Test timeband interactions  ################################

midday_cc <- allInsects %>% filter(Time_band == "midday")
evening_cc <- allInsects %>% filter(Time_band == "evening")

### community comp ##########

model_midday_cc <- lmer((distdif*100) ~ Agriculture_50 +
                  Forest_1000 +
                  Wetland_250 +
                  Urban_1000 +
                  Open.uncultivated.land_1000 +
                  cnumberTime +
                  cStops +
                  cyDay + 
                  (1|PilotID), data=subset(midday_cc, Open.uncultivated.land_1000 < 0.2))
summary(model_midday_cc)
tab_model(model_midday_cc)
r.squaredGLMM(model_midday_cc)

model_evening_cc <- lmer((distdif*100) ~ Agriculture_50 +
                  Forest_1000 +
                  Wetland_250 +
                  Urban_1000 +
                  Open.uncultivated.land_1000 +
                  cnumberTime +
                  cStops +
                  cyDay + 
                  (1|PilotID), data=subset(evening_cc, Open.uncultivated.land_1000 < 0.2))
summary(model_evening_cc)
tab_model(model_midday_cc, model_evening_cc, collapse.ci = TRUE, dv.labels = c("Community composition (midday)", "Community composition (evening)"), pred.labels = c("Intercept", "Farmland", "Forest", "Wetland", "Urban", "Grassland", "Sampling time", "Potential stops", "Day of year"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")
r.squaredGLMM(model_evening_cc)

### richness #######

model_midday_r <- lmer((richness*100) ~ Agriculture_500 +
                          Forest_1000 +
                          Wetland_50 +
                          Urban_1000 +
                          Open.uncultivated.land_1000 +
                          cnumberTime +
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
                           cStops +
                           cyDay + 
                           (1|PilotID), data=subset(evening_cc, Open.uncultivated.land_1000 < 0.2))
summary(model_evening_r)
r.squaredGLMM(model_evening_r)

tab_model(model_midday_r, model_evening_r, collapse.ci = TRUE, dv.labels = c("Richness (midday)", "Richness (evening)"), pred.labels = c("Intercept", "Farmland", "Forest", "Wetland", "Urban", "Grassland", "Sampling time", "Potential stops", "Day of year"),
          string.pred = "Coefficient",
          string.ci = "Conf. Int (95%)",
          string.p = "P-Value")

### as separate models ###########
# community composition
lme1000 <- lme(distdif ~ Agriculture_50*Time_band +
                 Forest_1000 +
                 Wetland_250 +
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
               correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(distdif ~ Agriculture_50 +
                 Forest_1000*Time_band +
                 Wetland_250 +
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
               correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(distdif ~ Agriculture_50 +
                 Forest_1000 +
                 Wetland_250*Time_band +
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
               correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(distdif ~ Agriculture_50 +
                 Forest_1000 +
                 Wetland_250+
                 Urban_1000*Time_band  +
                 Open.uncultivated.land_1000 +
                 Time_band +
                 Time_band:cnumberTime +
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops +
                 cyDay,
               random=~1|PilotID/RouteID,
               correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(distdif ~ Agriculture_50 +
                 Forest_1000 +
                 Wetland_250 +
                 Urban_1000 +
                 Open.uncultivated.land_1000*Time_band +
                 Time_band +
                 Time_band:cnumberTime +
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops +
                 cyDay,
               random=~1|PilotID/RouteID,
               correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

# richness
lme1000 <- lme(richness ~ Agriculture_500*Time_band + 
                 Forest_1000 +
                 Wetland_50 +
                 Urban_1000 +
                 Open.uncultivated.land_1000 +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(richness ~ Agriculture_500 + 
                 Forest_1000*Time_band +
                 Wetland_50 +
                 Urban_1000 +
                 Open.uncultivated.land_1000 +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(richness ~ Agriculture_500 + 
                 Forest_1000 +
                 Wetland_50*Time_band +
                 Urban_1000 +
                 Open.uncultivated.land_1000 +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(richness ~ Agriculture_500 + 
                 Forest_1000 +
                 Wetland_50 +
                 Urban_1000*Time_band +
                 Open.uncultivated.land_1000 +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)

lme1000 <- lme(richness ~ Agriculture_500 + 
                 Forest_1000 +
                 Wetland_50 +
                 Urban_1000 +
                 Open.uncultivated.land_1000*Time_band +
                 #Heathland_500 +
                 Time_band + 
                 Time_band:cnumberTime + 
                 Temperature +
                 Wind +
                 Time_driven +
                 cStops + 
                 cyDay,
               random=~1|PilotID/RouteID,
               data=subset(allInsects, Open.uncultivated.land_1000 < 0.2))
summary(lme1000)
