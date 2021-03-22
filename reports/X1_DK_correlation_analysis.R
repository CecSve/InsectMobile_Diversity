# run script 02 first

### load libraries ###################
library(ggpubr)
library(tidyverse)
library(corrplot)
library(Hmisc)
library(ggfortify) # to plot PCA
library(ggbiplot)

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "chartreuse3", "#D55E00", "#56B4E9", "#009E73") # colour friendly, ordered by land cover with heathland

landuseOrder <- c("Urban","Farmland","Grassland", "Heathland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified")

### Denmark ###########################

### Correlation plot community composition ##################

# prepare plotting
par(mfrow = c(2, 2))

names(allInsects_totsample)

# select variables for PCA - we will only use land cover at 1000 m and cStops
biomass.pca <- prcomp(allInsects_totsample[,c(153,26,48:53)], center = TRUE,scale. = TRUE) #choose biomass, stops and land covers
summary(biomass.pca)
str(biomass.pca)

ggbiplot(biomass.pca)

# correlation plot for 1000 m buffer
someInsects <- allInsects_totsample[,c(153,26,48:53)]
colnames(someInsects)
colnames(someInsects) <- c("CommunityVariation", "Stops", "Farmland", "Forest", "Heathland", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

## add p-values on no significant coefficient
corrplot(p, p.mat = res1$p, method = "color", type = "upper",
         sig.level = c(.001, .01, .05), pch.cex = .9,
         insig = "label_sig", pch.col = "white", order = "AOE", tl.cex = 0.6, tl.col = "black", col = landuseCols)

# with correlation coefficient as well
corrplot.mixed(p, p.mat = res1$p, lower.col = "black", upper = "color",
               sig.level = c(.001, .01, .05), pch.cex = .9,
               insig = "label_sig", pch.col = "white", order = "AOE", tl.cex = 0.7, tl.col = "black", upper.col = landuseCols, number.cex = 0.7)

corrplot(p, p.mat = res1$p, method = "color", type = "upper",
         insig = "label_sig", pch.col = "white",
         pch = "p<.05", pch.cex = .8, order = "AOE", tl.cex = 0.6, tl.col = "black", col = landuseCols) # "AOE" is for the angular order of the eigenvectors

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = landuseCols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 1000 m", mar=c(0,0,1,0))

# correlation plot for 500 m buffer
names(allInsects_totsample)
someInsects <- allInsects_totsample[,c(153,26,42:47)]
colnames(someInsects)
colnames(someInsects) <- c("CommunityVariation", "Stops", "Farmland", "Forest", "Heathland", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = landuseCols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 500 m", mar=c(0,0,1,0))

# correlation plot for 250 m buffer
someInsects <- allInsects_totsample[,c(153,26,36:41)]
colnames(someInsects)
colnames(someInsects) <- c("CommunityVariation", "Stops", "Farmland", "Forest", "Heathland", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = landuseCols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 250 m", mar=c(0,0,1,0))

# correlation plot for 50 m buffer
someInsects <- allInsects_totsample[,c(153,26,30:35)]
colnames(someInsects)
colnames(someInsects) <- c("CommunityVariation", "Stops", "Farmland", "Forest", "Heathland", "Grassland", "Urban", "Wetland")

p <- cor(someInsects)

# add significance
res1 <- cor.mtest(someInsects, conf.level = .95)
res2 <- cor.mtest(someInsects, conf.level = .99)

# with correlation coefficient instead of p-values, coloured boxes = significant at a 0.05 level
corrplot(p, method = "color", col = landuseCols,
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal 
         diag = FALSE, title = "Correlation at 50 m", mar=c(0,0,1,0))


### PCA #########

# subset data for each land cover buffer prior to analysis
names(allInsects_totsample)
cor50 <- allInsects_totsample[,c(2,30:35)]
cor250 <- allInsects_totsample[,c(2,36:41)]
cor500 <- allInsects_totsample[,c(2,42:47)]
cor1000 <- allInsects_totsample[,c(2,48:53)]

# 50 m buffer correlation of land cover variables
fit <- princomp(cor50, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor250, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor500, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

# 250 m buffer correlation of land cover variables
fit <- princomp(cor1000, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

mydata <- allInsects_totsample[,c("cStops",names(allInsects_totsample)[grepl("_1000",names(allInsects_totsample))])]
names(mydata)
mydata <- mydata[,c(2:7)]
names(mydata)[names(mydata)=="Open.uncultivated.land_1000"] <- "Grassland_1000"
names(mydata)[names(mydata)=="Agriculture_1000"] <- "Farmland_1000"
names(mydata) <- gsub("_1000","",names(mydata))
allInsects_totsample$Land_use <- as.character(allInsects_totsample$Land_use)
allInsects_totsample$Land_use[allInsects_totsample$Land_use=="Dryland"] <- "Grassland"
allInsects_totsample$Land_use[allInsects_totsample$Land_use=="Open uncultivated land"] <- "Grassland"
landuseOrder
allInsects_totsample$Land_use <- factor(allInsects_totsample$Land_use, levels=landuseOrder)

fit <- princomp(mydata, cor=TRUE)

#pca with rotation
library(psych)
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
pca_rotated <- psych::principal(mydata, rotate="varimax", nfactors=2, scores=TRUE)
biplot(pca_rotated, main = "A: Denmark")
print(pca_rotated)

ggsave("plots/pca_with_rotation_1000_DK.png")

#add PCA axes scores to the dataset
allInsects_totsample$Urbanization_gradient <- pca_rotated$scores[,1]
allInsects_totsample$Forest_gradient <- pca_rotated$scores[,2]

#with ggplot
autoplot(fit)
dk_autoplot <- autoplot(fit, data = allInsects_totsample, colour = 'maxLand_use', 
                        loadings = TRUE, 
                        loadings.colour = 'black',
                        loadings.label = TRUE, 
                        loadings.label.size = 5) + 
  scale_colour_manual(values = landuseCols[1:6])+
  theme_bw() + labs(colour = "Land cover")

ggsave("plots/pca_1000_DK.png")
dk_autoplot <- plot_grid(dk_autoplot, labels = "AUTO")
save_plot("plots/pca_1000_DK_numbered.png", dk_autoplot, base_height = 8, base_width = 12)

lme1000 <- lmer(log(Biomass+1) ~ 
                  Urbanization_gradient + 
                  Forest_gradient +
                  Time_band + 
                  Time_band:cnumberTime + cStops + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
summary(lme1000)

library(MuMIn)
r.squaredGLMM(lme1000)


### other stuff - not used ##################
# examine interaction betweeen DNA conc, volume and biomass of insects
labdata$Qubit <- as.numeric(labdata$Qubit)
str(labdata)
plot(log(labdata$Biomass), log(labdata$Qubit))
plot(log(labdata$Qubit), log(labdata$Biomass))
plot(log(labdata$InsectVolume), log(labdata$Qubit))
plot(log(labdata$Biomass), log(labdata$InsectVolume))

# basic regression plots with pearson correlation
ggscatter(
  labdata,
  x = "Qubit",
  y = "Biomass",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray")
) +  stat_cor(method = "pearson") + scale_y_log10() +scale_x_log10() # DNA conc increases with increasing biomass

ggscatter(
  labdata,
  x = "Qubit",
  y = "InsectVolume",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray")
) +  stat_cor(method = "pearson") + scale_y_log10() +scale_x_log10() # DNA conc increases with more insects OR larger insects (volume is based on how much space the insects take up in a 15 ml centrifuge tube)

ggscatter(
  labdata,
  x = "Qubit",
  y = "DNABufferVolumeAdded",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray")
) +  stat_cor(method = "pearson") +scale_x_log10() # More buffer (which is increasing with insectvolume) increases DNA conc

ggscatter(
  labdata,
  x = "DNABufferVolumeAdded",
  y = "InsectVolume",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray")
) +  stat_cor(method = "pearson") + scale_y_log10() # insect volume is (and should be) associated with how much buffer is added, although it is not a linear relationship, but rather asymptotic

test <- merge(data, labdata)

ggscatter(
  test,
  x = "yDay",
  y = "Biomass",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "darkgrey", fill = "lightgray"), color = "roughLand_use"
) +  stat_cor(method = "pearson") + scale_y_log10()
