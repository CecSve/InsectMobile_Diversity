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

### load data ########
allInsects_totsample <- read.delim("cleaned-data/allInsects_totsample.txt", sep = " ")
totsample_asvs <- read.delim("cleaned-data/totsample_asvs.txt", sep = " ")

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

### comcomp: simple regression plot land cover##############################################
# change land covers to be 0-100 instead of 0-1
# allInsects_totsample[, c(30:141,143)] <- allInsects_totsample[, c(30:141,143)]*100

# farmland
lme50 <- lmer(NMDS1 ~ (Agriculture_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                 (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Agriculture_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Agriculture_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                  + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Agriculture_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban
lme50 <- lmer(NMDS1 ~ (Urban_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Urban_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Urban_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Urban_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated

lme50 <- lmer(NMDS1 ~ (Open.uncultivated.land_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme250 <- lmer(NMDS1 ~ (Open.uncultivated.land_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme500 <- lmer(NMDS1 ~ (Open.uncultivated.land_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme1000 <- lmer(NMDS1 ~ (Open.uncultivated.land_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest

lme50 <- lmer(NMDS1 ~ (Forest_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Forest_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Forest_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Forest_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland

lme50 <- lmer(NMDS1 ~ (Wetland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Wetland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Wetland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Wetland_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated)
outAll$Land_use <- factor(outAll$Land_use,levels=landuseOrder)

str(outAll)

effect_cc <- ggplot(outAll)+
  geom_crossbar(aes(x=factor(Buffer),y=Estimate,
                    ymin=Estimate-(Std..Error*1.96),ymax=Estimate+(Std..Error*1.96),
                    fill=Land_use),
                width=0.5)+
  facet_wrap(~Land_use,scales="free",ncol=1)+
  scale_fill_manual(values=landuseCols)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,colour="black",linetype="dashed")+
  xlab("Buffer size (m)") + ylab("Effect of land cover on community composition")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black"))


### richness: simple regression plot land cover##############################################

### examine whether ASV richness is normally distibuted ##########
hist(allInsects_totsample$richness)
hist(sqrt(allInsects_totsample$richness))

qqnorm(allInsects_totsample$richness)
qqnorm(sqrt(allInsects_totsample$richness))

### buffer effect plots ######################
# NB! changed cTL to cStops since more data for DK 
allInsects_totsample$richness.t <- sqrt(allInsects_totsample$richness) # make a response variable with sqrt transformed richness, for normal distribution

# farmland
lme50 <- lmer(richness.t ~ (Agriculture_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness.t ~ (Agriculture_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness.t ~ (Agriculture_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness.t ~ (Agriculture_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban

lme50 <- lmer(richness.t ~ (Urban_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness.t ~ (Urban_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness.t ~ (Urban_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness.t ~ (Urban_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated

lme50 <- lmer(richness.t ~ (Open.uncultivated.land_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme250 <- lmer(richness.t ~ (Open.uncultivated.land_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme500 <- lmer(richness.t ~ (Open.uncultivated.land_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme1000 <- lmer(richness.t ~ (Open.uncultivated.land_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest

lme50 <- lmer(richness.t ~ (Forest_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness.t ~ (Forest_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness.t ~ (Forest_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness.t ~ (Forest_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland

lme50 <- lmer(richness.t ~ (Wetland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(richness.t ~ (Wetland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(richness.t ~ (Wetland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(richness.t ~ (Wetland_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|PilotID), data=allInsects_totsample)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated)
outAll$Land_use <- factor(outAll$Land_use,levels=landuseOrder)

# somthing is of with the CIs so we will manually add them
effect_r <- ggplot(outAll)+
  geom_crossbar(aes(x=factor(Buffer),y=Estimate,
                    ymin=Estimate-(Std..Error*1.96),ymax=Estimate+(Std..Error*1.96),
                    fill=Land_use),
                width=0.5)+
  facet_wrap(~Land_use,scales="free",ncol=1)+
  scale_fill_manual(values=landuseCols)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none")+
  geom_hline(yintercept=0,colour="black",linetype="dashed")+
  xlab("Buffer size (m)") + ylab("Effect of land cover on insect ASV richness \n(sqrt-transformed)")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black"))

effectplots <- cowplot::plot_grid(effect_cc, effect_r, labels = "AUTO", align = "hv", ncol = 2)
ggsave("plots/DK_Landcover_buffer_effect.png",width=8,height=12)

