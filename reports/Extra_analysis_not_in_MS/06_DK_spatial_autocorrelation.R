# spatial autocorrelation model selection


### Load required libraries ########################################
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(lme4)# Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(ncf) # visualise spatial autocorrelation
library(DHARMa) # test spatial autocorrelation
library(MuMIn)

# for DK jitter x and y slightly - fix later
allInsects$x2 <- allInsects$utm_x + rnorm(length(allInsects$utm_x),0,10)
allInsects$y2 <- allInsects$utm_y + rnorm(length(allInsects$utm_y),0,10)

hist((allInsects$distdif)) # this is best - there's an outlier that we could consider to remove
hist(log(allInsects$distdif))

qqnorm(allInsects$distdif)
qqnorm(log(allInsects$distdif)) # this looks best, wil use log

### spatial autocorrelation ###########
#plot residuals
#full model
unique_coords <- distinct(allInsects, RouteID_JB, .keep_all = TRUE) # make dataframe with only one coordinate per row to test spatial autocorrelation for unique coordinates
unique_coords <- droplevels(unique_coords)

### community composition ##########
# since we have unique coordinates (which is required for the spatial autocorrelation test), route ID cannot be used as random effect, since they aren't < # of observations
lme1000 <- lme4::lmer(distdif ~ 
                        (Agriculture_1000) + 
                        (Urban_1000) +
                        (Open.uncultivated.land_1000)+
                        (Wetland_1000) +
                        (Forest_1000) +
                        Time_band + 
                        Time_band:cnumberTime + cTL + cyDay + (1|PilotID), data=unique_coords)
summary(lme1000)

unique_coords$resids <- as.numeric(residuals(lme1000))
unique_coords$resids_binary <- as.factor(ifelse(unique_coords$resids>0,1,-1))
#qplot(x,y,data=unique_coords,colour=resids)+ scale_colour_viridis_c()
#qplot(x,y,data=allInsects,colour=resids_binary)+ scale_colour_viridis_d()

# visualise it 
spline.autocorrelation_glm = spline.correlog(unique_coords$x2, unique_coords$y2, residuals(lme1000), latlon=T, resamp=100)
plot(spline.autocorrelation_glm)
summary(spline.autocorrelation_glm)

autocorrelation_glm = correlog(unique_coords$x2, unique_coords$y2, residuals(lme1000), increment = 1000, latlon=T, resamp=100)
plot(autocorrelation_glm)
autocorrelation_glm$correlation

#test it
simulationOutput <- simulateResiduals(fittedModel = lme1000, plot = T)
#residuals(simulationOutput)
#residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
testDispersion(lme1000)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)

res = simulateResiduals(lme1000)
testSpatialAutocorrelation(res, x =  unique_coords$x2, y = unique_coords$y2)

### richness ##################
lme1000 <- lme4::lmer(richness ~ 
                        (Agriculture_1000) + 
                        (Urban_1000) +
                        (Open.uncultivated.land_1000)+
                        (Wetland_1000) +
                        (Forest_1000) +
                        Time_band + 
                        Time_band:cnumberTime + cTL + cyDay + (1|PilotID), data=unique_coords)
summary(lme1000)

#unique_coords$resids <- as.numeric(residuals(lme1000))
#unique_coords$resids_binary <- as.factor(ifelse(unique_coords$resids>0,1,-1))
#qplot(x,y,data=unique_coords,colour=resids)+ scale_colour_viridis_c()
#qplot(x,y,data=allInsects,colour=resids_binary)+ scale_colour_viridis_d()

# visualise it 
spline.autocorrelation_glm = spline.correlog(unique_coords$x2, unique_coords$y2, residuals(lme1000), latlon=T, resamp=100)
plot(spline.autocorrelation_glm)
summary(spline.autocorrelation_glm)

autocorrelation_glm = correlog(unique_coords$x2, unique_coords$y2, residuals(lme1000), increment = 1000, latlon=T, resamp=100)
plot(autocorrelation_glm)
autocorrelation_glm$correlation

#test it
simulationOutput <- simulateResiduals(fittedModel = lme1000, plot = T)
#residuals(simulationOutput)
#residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
testDispersion(lme1000)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)

res = simulateResiduals(lme1000)
testSpatialAutocorrelation(res, x =  unique_coords$x2, y = unique_coords$y2)

### hertil ########
### model selection ###########
gls1 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
              Forest_250 +
              Wetland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID,
            correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=TRUE),
            data=allInsects)
summary(gls1)

#range     nugget 
#0.1668    0.1198  
AICc(gls1)#925.34

gls1 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
              Forest_250 +
              Wetland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
            data=allInsects)

summary(gls1)
#Parameter estimate(s):
#  range     nugget 
#0.1109   0.2492 
AICc(gls1)#923

gls1 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
              Forest_250 +
              Wetland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID,
            correlation=corExp(form=~x2+y2|PilotID/RouteID,nugget=FALSE),
            data=allInsects)
summary(gls1)

#range     nugget 
#0.165       
AICc(gls1)#923.2

gls1 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
              Forest_250 +
              Wetland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
            data=allInsects)
summary(gls1)
#range 
#0.1311   
AICc(gls1) # 923.3

gls1 <- lme(log(Biomass+1) ~ Agriculture_1000 + 
              Forest_250 +
              Wetland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID,
            data=allInsects)
summary(gls1)
AICc(gls1) # 921

#final model DK - use cStops instead of cTL - same story as without spatial correlation
gls1 <- lme(
  log(Biomass + 1) ~ #Urban_1000 + 
    Agriculture_1000 +
    Open.uncultivated.land_1000 +
    Wetland_50 +
    Forest_250 +
    Time_band +
    Time_band:cnumberTime +
    cStops +
    cyDay,
  random =  ~ 1 | PilotID / RouteID,
  data = allInsects
)
summary(gls1)
AICc(gls1) # 915 (best fit is 913, but with urban that should be removed due to corelation with stops
#keep in TL even if not significant