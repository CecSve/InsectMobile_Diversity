# spatial autocorrelation model selection


### Load required libraries ########################################
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(lme4)# Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(ncf) # visualise spatial autocorrelation
library(DHARMa) # test spatial autocorrelation
library(MuMIn)

# for DK jitter x and y slightly - fix later
allInsects_totsample$x2 <- allInsects_totsample$utm_x + rnorm(length(allInsects_totsample$utm_x),0,10)
allInsects_totsample$y2 <- allInsects_totsample$utm_y + rnorm(length(allInsects_totsample$utm_y),0,10)

hist((allInsects_totsample$NMDS1)) # this is best - there's an outlier that we could consider to remove
#hist(log(allInsects_totsample$distdif))

qqnorm(allInsects_totsample$NMDS1)
#qqnorm(log(allInsects_totsample$distdif)) 

# since the NMDS values are very low (and does not really hold any meaning) we multiply with 10000

allInsects_totsample$NMDS1 <- allInsects_totsample$NMDS1*10000

### spatial autocorrelation ###########
#plot residuals
#full model
unique_coords <- distinct(allInsects_totsample, RouteID_JB, .keep_all = TRUE) # make dataframe with only one coordinate per row to test spatial autocorrelation for unique coordinates
unique_coords <- droplevels(unique_coords)

### community composition ##########
# since we have unique coordinates (which is required for the spatial autocorrelation test), route ID cannot be used as random effect, since they aren't < # of observations
# buffer size with largest effect for each land cover
lme1000 <- lme4::lmer(NMDS1 ~ 
                        (Agriculture_500) + 
                        (Urban_1000) +
                        (Open.uncultivated.land_50)+
                        (Wetland_50) +
                        (Forest_1000) +
                        Time_band + 
                        Time_band:cnumberTime + cTL + cyDay + (1|PilotID), data=unique_coords)
summary(lme1000)

unique_coords$resids <- as.numeric(residuals(lme1000))
unique_coords$resids_binary <- as.factor(ifelse(unique_coords$resids>0,1,-1))
qplot(x2,y2,data=unique_coords,colour=resids)+ scale_colour_viridis_c()
#qplot(x2,y2,data=allInsects_totsample,colour=resids_binary)+ scale_colour_viridis_d()

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
hist(allInsects_totsample$richness)
hist(sqrt(allInsects_totsample$richness)) # sqrt transformation = normal distribution
unique_coords$richness.t <- sqrt(unique_coords$richness)
lme1000 <- lme4::lmer(richness.t ~ 
                        (Agriculture_500) + 
                        (Urban_1000) +
                        (Open.uncultivated.land_50)+
                        (Wetland_50) +
                        (Forest_1000) +
                        Time_band + 
                        Time_band:cnumberTime + cTL + cyDay + (1|PilotID), data=unique_coords)
summary(lme1000)

unique_coords$resids <- as.numeric(residuals(lme1000))
unique_coords$resids_binary <- as.factor(ifelse(unique_coords$resids>0,1,-1))
qplot(x2,y2,data=unique_coords,colour=resids)+ scale_colour_viridis_c()
#qplot(x,y,data=allInsects_totsample,colour=resids_binary)+ scale_colour_viridis_d()

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

### model selection community composition ###########
gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              Time_band + 
              Time_band:cnumberTime + 
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            correlation=corExp(form=~x2+y2|PilotID/RouteID_JB,nugget=TRUE),
            data=allInsects_totsample)
summary(gls1)

#range     nugget 
#6.126435e+00    4.426024e-07 
AICc(gls1)#  | with weather: temp, wind and time driven = 5624.944

gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              Time_band + 
              Time_band:cnumberTime + 
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
            data=allInsects_totsample)

summary(gls1)
#Parameter estimate(s):
#  range     nugget 
#6.127240e+00    5.011707e-07 
AICc(gls1)#  |  5622.692

gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              #Heathland_50 +
              Time_band + 
              Time_band:cnumberTime +
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            correlation=corExp(form=~x2+y2|PilotID/RouteID_JB,nugget=FALSE),
            data=allInsects_totsample)
summary(gls1)

#range     nugget 
#6.127248       
AICc(gls1)# |  5622.692

gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              #Heathland_50 +
              Time_band + 
              Time_band:cnumberTime +
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
            data=allInsects_totsample)
summary(gls1)
#range 
#6.127266   
AICc(gls1) #  | 5620.452

gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              #Heathland_50 +
              Time_band + 
              Time_band:cnumberTime + 
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            data=allInsects_totsample)
summary(gls1)
AICc(gls1) # | 5623.211

### community composition: final model ############
#use cStops instead of cTL - model AICs with or without spatial autocorrelation are very similar. Since we find spatial autocorrelation, we choose the model that accounts for this with the lowest AIC
gls1 <- lme(NMDS1 ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_50 +
              Time_band + 
              Time_band:cnumberTime +
              Temperature +
              Wind +
              Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
            data=allInsects_totsample)

summary(gls1)
AICc(gls1) 

### model selection richness ###########
allInsects_totsample$richness.t <- sqrt(allInsects_totsample$richness)

gls1 <- lme(richness.t ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Time_band + 
              Time_band:cnumberTime + 
              #Temperature +
              #Wind +
              #Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            correlation=corExp(form=~x2+y2|PilotID/RouteID_JB,nugget=TRUE),
            data=allInsects_totsample)
summary(gls1)

#range     nugget 
#0.1173881       0.1214073      
AICc(gls1)# 1747.406 | 1756.371 # better fit without weather, wind duration
qqnorm(resid(gls1))

gls1 <- lme(richness.t ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Time_band + 
              Time_band:cnumberTime + 
              #Temperature +
              #Wind +
              #Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
            data=allInsects_totsample)

summary(gls1)
#Parameter estimate(s):
#  range     nugget 
#78.2860456      0.7261387    
AICc(gls1)# 1745.19 | 1753.992
qqnorm(resid(gls1))

gls1 <- lme(richness.t ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Time_band + 
              Time_band:cnumberTime + 
              #Temperature +
              #Wind +
              #Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            correlation=corExp(form=~x2+y2|PilotID/RouteID_JB,nugget=FALSE),
            data=allInsects_totsample)
summary(gls1)

#range     nugget 
#0.1182003         
AICc(gls1)# 1745.218 | 
qqnorm(resid(gls1))

gls1 <- lme(richness.t ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Time_band + 
              Time_band:cnumberTime +
              #Temperature +
              #Wind +
              #Time_driven +
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=FALSE),
            data=allInsects_totsample)
summary(gls1)
#range 
#1.702399      
AICc(gls1) # 1747.29| 
qqnorm(resid(gls1))

gls1 <- lme(richness.t ~ Agriculture_500 + 
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
            random=~1|PilotID/RouteID_JB,
            data=allInsects_totsample)
summary(gls1)
AICc(gls1) # 1751.879 | 
#qqnorm(resid(gls1))

### richness: final model ############
#use cStops instead of cTL - model AICs with or without spatial autocorrelation are very similar. Since we find spatial autocorrelation, we choose the model that accounts for this with the lowest AIC
gls1 <- lme(richness.t ~ Agriculture_500 + 
              Forest_1000 +
              Wetland_50 +
              Urban_1000 +
              Open.uncultivated.land_1000 +
              Time_band + 
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID,
            correlation=corExp(form=~x2+y2|PilotID,nugget=TRUE),
            data=allInsects_totsample)
summary(gls1)
AICc(gls1) 
#qqnorm(resid(gls1))