# Simple regression models on richness

# Supplementary material analysis for InsectMobile diversity paper

#### Load required libraries ################################################################
library(cowplot)
library(ggplot2)
#library(wesanderson)
library(ggpubr)
library(scales)
library(ggpmisc)
library(grid)
library(gridExtra)
library(corrplot)
library(car)
library(lme4)
library(lmerTest)
library(nlme)
library(effects)

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified")

### get effect function #############################################################

#function to extract summary components of simpel buffer effect models
getEffect <- function(model){
  coefs <- summary(model)$coef[2,]
  temp <- confint(model)
  cis <- temp[5,]
  data.frame(t(coefs),t(cis))
}

# fitModels <- function(variable){
#   myformula()
#   
# }

### richness: simple regression plot land cover##############################################

### examine whether ASV richness is normally distibuted ##########
hist(allInsects_totsample$richness)
hist(sqrt(allInsects_totsample$richness))

qqnorm(allInsects_totsample$richness)
qqnorm(sqrt(allInsects_totsample$richness))

### buffer effect plots ######################
# NB! changed cTL to cStops since more data for DK 
str(allInsects_totsample)
#agriculture
hist(allInsects_totsample$richness)
hist(sqrt(allInsects_totsample$richness)) 
lme50 <- lmer(richness ~ (Agriculture_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness ~ (Agriculture_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness ~ (Agriculture_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness ~ (Agriculture_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban

lme50 <- lmer(richness ~ (Urban_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness ~ (Urban_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness ~ (Urban_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness ~ (Urban_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated

lme50 <- lmer(richness ~ (Open.uncultivated.land_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme250 <- lmer(richness ~ (Open.uncultivated.land_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme500 <- lmer(richness ~ (Open.uncultivated.land_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme1000 <- lmer(richness ~ (Open.uncultivated.land_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest

lme50 <- lmer(richness ~ (Forest_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness ~ (Forest_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness ~ (Forest_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness ~ (Forest_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland

lme50 <- lmer(richness ~ (Wetland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness ~ (Wetland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness ~ (Wetland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness ~ (Wetland_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated)
outAll$Land_use <- factor(outAll$Land_use,levels=landuseOrder)

ggplot(outAll)+
  geom_crossbar(aes(x=factor(Buffer),y=Estimate,
                    ymin=X2.5..,ymax=X97.5..,
                    fill=Land_use),
                width=0.5)+
  facet_wrap(~Land_use,scales="free",ncol=1)+
  scale_fill_manual(values=landuseCols)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,colour="black",linetype="dashed")+
  xlab("Buffer size (m)") + ylab("Effect of land cover on ASV richness")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black"))

ggsave("plots/DK_Landcover_buffer_richness.png",width=5,height=12)

### max land cover (1000 m) ###########
mini_data <- allInsects_totsample %>% select(PCRID, maxLand_use, maxareaProp, richness)
table(mini_data$maxLand_use)
variables <- c("Agriculture_1000", "Forest_1000", "Urban_1000")
mini_data <- mini_data %>% dplyr::filter(maxLand_use %in% variables) 
mini_data <- mini_data %>% dplyr::mutate(maxLand_use=dplyr::recode(maxLand_use, 
                                                                   "Agriculture_1000"="Farmland",
                                                                   "Forest_1000"="Forest",
                                                                   "Urban_1000"="Urban"))

hist(mini_data$maxareaProp)
qqnorm(mini_data$maxareaProp)
hist(mini_data$richness)
qqnorm(mini_data$richness)

# plot without regression line
u <- mini_data %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Farmland",
    "Forest",
    "Urban",
  )
) %>% ggplot() + geom_point(aes(x=maxareaProp,y=richness, colour = maxLand_use), show.legend = F, size = 3) + facet_wrap(. ~ maxLand_use, nrow = 3) +
  theme_pubclean() + xlab("Land cover at 1000 m") + ylab("Variation in community composition") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(
    limits = c(0, 1),
    labels = function(x)
      paste0(x * 100, "%")) + scale_colour_manual(values = landuseCols[c(2, 6, 1)])

# plot with regression line
upubr <- mini_data %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Farmland",
    "Forest",
    "Urban",
  )
) %>% ggscatter(x = "maxareaProp", y = "richness",
                color = "maxLand_use",
                palette = landuseCols[c(2,6,1)],
                shape = 19, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coeff.args = list(method = "spearman")
) + stat_cor(aes(label = paste(..rr.label..,
                               if_else(readr::parse_number(..p.label..) < 0.001, 
                                       "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0, label.y = 20) + xlab("Land cover at 1000 m") + ylab("Variation in richness (sqrt transformed)") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%"))  + facet_wrap(. ~ maxLand_use, nrow = 3) + theme_pubclean() + theme(legend.position = "none")

ggsave("plots/DK_LandcoverMax1000_richness_correlation.png",width=5,height=12)

