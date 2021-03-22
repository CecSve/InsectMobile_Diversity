
### load libraries ##########
library(tidyverse)
library(vegan)
library(ggrepel)
library(ape)
library(reshape2)

### load data ##########
asvs <- read.delim("cleaned-data/DK_asvtable_2018_data.txt")
taxonomy_insecta <- read.delim("cleaned-data//DK_taxonomy_Insecta.txt",sep="\t")
taxonomy_99 <- read.delim("cleaned-data//DK_taxonomy_Insecta_99.txt",sep="\t")
data <- read.delim("cleaned-data/allInsects.txt", sep = " ")

### aligning data for analysis ###################
names(asvs) # only 2018 samples, with no negatives and blanks, german samples, or 2017 samples

keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset asvtable to only contain samples with metadata

keep <- names(asvs)
data <- data %>% filter(PCRID %in% keep)

# do all samples have seuqnces?
colSums(asvs) # yes
#asvs <- asvs[,colSums(asvs) > 0]

# remove empty asvs
rowSums(asvs)
min(rowSums(asvs)) # no empty asvs or samples

length(unique(data[["RouteID_JB"]])) # count how many routes were sampled - but notice some have received new routes
length(unique(data[["PilotID"]])) # count how many pilots that carried out the sampling
length(unique(data[["SampleID"]])) # count how many samples
table(data$maxLand_use)

### merging total samples from size sorted samples ###########
# 
data_unique <- data %>% distinct(SampleID, .keep_all = TRUE)
keep <- data %>% select(PCRID, SampleID)

t.asvs <- t(asvs)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- full_join(keep, t.asvs, by = "PCRID")

test2 <- test %>% select(-PCRID) %>% group_by(SampleID) %>% summarise_all(list(sum))

totsample_asvs <- test2 %>% column_to_rownames(var = "SampleID")
totsample_asvs <- as.data.frame(t(totsample_asvs))

# remove outlier
totsample_asvs <- totsample_asvs %>% rownames_to_column(var = "PCRID") %>% select(-P153.1B) %>% column_to_rownames(var = "PCRID")

### Community composition ##################
tasvs <- t(totsample_asvs) # samples as rows, taxa as columns
#bray.data<-vegdist(tasvs, Type = "bray") # 
# the Bray-Curtis index is based on abundance data, while the Sorensen index is based on presence/absence data. Both indices have similarity and dissimilarity (or distance) versions. 
tasvs <- decostand(tasvs, method = "pa") # transform to presence absence
sorensen.data <- designdist(tasvs, "(A+B-2*J)/(A+B)") # for presence absence

# due to issues in the ape package (negative eigen values comes before positive eigen values and causes an error), we have to use NMDS for ordiation vizualisation 
set.seed(41)
meta_mds <- metaMDS(sorensen.data, k = 1, trymax = 100)
meta_mds$points
#MDS <- cmdscale(vegdist(bray.data, method = "bray"), k = 2, eig = T, add = T )
#round(MDS$eig*100/sum(MDS$eig),1) # % of variance explained by the MDS axes
stressplot(meta_mds)
plot(meta_mds)
#plot(MDS$points)

class(meta_mds)
#class(MDS)
methods(class="metaMDS")
meta_mds$stress
Y <- scores(meta_mds, display="sites")
#G <- scores(MDS, display = "sites")
#plot(Y, type="n")
#text(Y[,1], Y[,2], rownames(Y), col="red")
#plot(G, type = "n")
#text(G[,1], G[,2], rownames(G), col="red")
all(rownames(tasvs) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data

Y <- data.frame(Y)
Y <- Y %>% rownames_to_column(var = "SampleID") 
Y$richness <- rowSums(tasvs) # change to relative?

#G <- data.frame(G)
#G <- G %>% rownames_to_column(var = "SampleID") 
#G$readcount <- rowSums(tasvs) # change to relative?

# add differences in distance column
#Y$distdif <- Y$NMDS1 - Y$NMDS2 # not working because it does not yield a result that show similarity between samples

allInsects <- merge(data, Y, by = "SampleID")
allInsects_totsample <- merge(data_unique, Y, by = "SampleID")

#allinsects_mds <- merge(data_unique, G, by = "SampleID")

# plot all samples if both NMDS are included only
allInsects %>% ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size = 5, aes(colour = NMDS1), show.legend = T) + 
  theme_minimal() + labs(colour = "dissimilarity") + scale_colour_viridis_c(option = "D")

allInsects %>% ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(size = 5, aes(colour = maxLand_use), show.legend = T) + 
  theme_minimal() + labs(colour = "dissimilarity") + stat_ellipse( aes(colour = maxLand_use))

allinsects_mds %>% ggplot(aes(x=Dim1, y=Dim2)) +
  geom_point(size = 5, aes(colour = maxLand_use), show.legend = T) + 
  theme_minimal() + labs(colour = "dissimilarity") + stat_ellipse( aes(colour = maxLand_use))

adonis2(sorensen.data~ Agriculture_1000 + Urban_1000 + Open.uncultivated.land_1000 + Forest_1000 + Wetland_1000 + Temperature + Time_driven + Wind + cStops + Time_band + cyDay, strata = c("RouteID_JB", "PilotID"), data = allInsects)

library(ecodist)

MRM(dist(sorensen.data)~dist(Agriculture_1000) + dist(Urban_1000) + dist(Open.uncultivated.land_1000) + dist(Forest_1000) + dist(Wetland_1000) + dist(Time_driven) + dist(cStops) + dist(cyDay), data = allInsects)

### Diversity  not used#############
pa <- decostand(totsample_asvs, method = "pa")
tpa <- t(pa)
tpa <- as.data.frame(tpa)
tpa$richness <- rowSums(tpa)

richnessdata <- tpa %>% rownames_to_column(var = "SampleID") %>% select(SampleID, richness)

allInsects <- merge(allInsects, richnessdata, by = "SampleID") # data frame based on size sorted samples
allInsects_totsample <- merge(data_unique, richnessdata, by = "SampleID")
allInsects_totsample <- merge(allInsects_totsample, Y, by = "SampleID") # data frame based on total samples

#write.table(allInsects_totsample, file = "cleaned-data/allInsects_totsample.txt")
#write.table(totsample_asvs, file = "cleaned-data/allInsects_totasvs.txt")

### envfit community comp  - not used ###########
env_var <- allInsects_totsample %>% select(Agriculture_1000, Urban_1000, Open.uncultivated.land_1000, Wetland_1000, Forest_1000, maxLand_use, distdif)
fit <- envfit(meta_mds, env_var, perm = 999, na.rm = T)
scores(fit, "vectors")
plot(meta_mds)
plot(fit)
plot(fit, p.max = 0.05, col = "deeppink")
