
# scatter plot with the car package
scatterplot(sqrt(richness)~log(Biomass+1) | hab50, data = allInsects, smooth = FALSE, grid = FALSE, frame = FALSE)

# scatterplot with the ggpubr package
library(ggpubr)

allInsects$Biomass.t <- log(allInsects$Biomass+1)

allInsects %>% ggscatter("Biomass.t", "richness.t", add = "reg.line", conf.int = T) + stat_cor(method = "pearson") + xlab("Insect biomass (+1 & log-transformed)") + ylab("Insect ASV richness (sqrt-transformed)")

allInsects %>% ggscatter("Biomass.t", "richness.t", add = "reg.line", conf.int = T, color = "hab50", palette = "jco",) + stat_cor(method = "pearson") 

allInsects_totsample$Biomass.t <- log(allInsects_totsample$Biomass+1)
allInsects_totsample %>% ggscatter("Biomass.t", "NMDS1", add = "reg.line", conf.int = T) + stat_cor(method = "pearson") + xlab("Insect biomass (+1 & log-transformed)") + ylab("Variation in community composition (NMDS1)")

# with only one nmds

# get proportions in 0-100 instead of 0-1
allInsects_totsample[, c(30:141, 143)] <- allInsects_totsample[, c(30:141, 143)] * 100

test <- allInsects_totsample

test$x2 <- test$utm_x + rnorm(length(test$utm_x),0,10)
test$y2 <- test$utm_y + rnorm(length(test$utm_y),0,10)

### simple regression plot ########
# farmland
lme50 <- lmer(NMDS1 ~ (Agriculture_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Agriculture_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Agriculture_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Agriculture_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban

lme50 <- lmer(NMDS1 ~ (Urban_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Urban_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Urban_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Urban_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated

lme50 <- lmer(NMDS1 ~ (Open.uncultivated.land_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme250 <- lmer(NMDS1 ~ (Open.uncultivated.land_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme500 <- lmer(NMDS1 ~ (Open.uncultivated.land_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
lme1000 <- lmer(NMDS1 ~ (Open.uncultivated.land_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.4))
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#forest

lme50 <- lmer(NMDS1 ~ (Forest_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Forest_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Forest_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Forest_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland

lme50 <- lmer(NMDS1 ~ (Wetland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(NMDS1 ~ (Wetland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(NMDS1 ~ (Wetland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(NMDS1 ~ (Wetland_1000) + Time_band + 
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
  xlab("Buffer size (m)") + ylab("Effect of land cover on insect ASV richness (sqrt-transformed)")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black"))



#test <- test %>% filter(Biomass > 0)

#full and final model - richness
#hist(test$Biomass)
#hist(log(test$Biomass))
hist(test$NMDS1)
hist(log(allInsects_totsample$NMDS1))

model2 <- lme(NMDS1 ~ Agriculture_500 + Urban_1000 + Open.uncultivated.land_50 + Forest_1000 + Wetland_50 +
                Time_band + 
                Time_band:cnumberTime + 
                Temperature +
                Wind +
                Time_driven +
                cStops + 
                cyDay,
              random=~1|PilotID/RouteID_JB,
              data=test)
summary(model2)
model2$apVar
r.squaredGLMM(model2)

tab_model(model2)

#check variance inflation factor
vif(model2)

get_variance(model2)

### Figure X: effect plots ##########################
gls1.alleffects <- allEffects(model2)
#plot(gls1.alleffects, 'Urban_1000:cover', ylab="Variation")
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model2)
plot(eall.lm1, lines=list(multiline=TRUE))
#plot(predictorEffects(model2, ~ landcover:cover + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))
