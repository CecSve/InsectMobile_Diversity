
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
data_unique <- data %>% dplyr::distinct(SampleID, .keep_all = TRUE)
keep <- data %>% dplyr::select(PCRID, SampleID)

t.asvs <- t(asvs)
t.asvs <- as.data.frame(t.asvs) %>% rownames_to_column(var = "PCRID") 
test <- full_join(keep, t.asvs, by = "PCRID")

test2 <- test %>% select(-PCRID) %>% group_by(SampleID) %>% summarise_all(list(sum))

totsample_asvs <- test2 %>% column_to_rownames(var = "SampleID")
totsample_asvs <- as.data.frame(t(totsample_asvs))

# remove outlier P153.1B if running community composition analysis

### Diversity #############
pa <- decostand(totsample_asvs, method = "pa") # transform matrix into presence absence
tpa <- t(pa) # transpose
tpa <- as.data.frame(tpa)
tpa$richness <- rowSums(tpa) # get richness for eack sample

richnessdata <- tpa %>% rownames_to_column(var = "SampleID") %>% select(SampleID, richness)

allInsects <- merge(data, richnessdata, by = "SampleID") # data frame based on size sorted samples
allInsects_totsample <- merge(data_unique, richnessdata, by = "SampleID")

### community comp ########
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
all(rownames(tasvs) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data

Y <- data.frame(Y)
Y <- Y %>% rownames_to_column(var = "SampleID") 

allInsects <- merge(allInsects, Y, by = "SampleID")
allInsects_totsample <- merge(allInsects_totsample, Y, by = "SampleID")

write.table(allInsects_totsample, file = "cleaned-data/allInsects_totsample.txt")
write.table(totsample_asvs, file = "cleaned-data/totsample_asvs.txt")

### indicator species analysis ###########

# change land covers to be 0-100 instead of 0-1
allInsects_totsample[, c(30:141,143)] <- allInsects_totsample[, c(30:141,143)]*100

# Create new column that picks out samples that are at the extreme in representing a single land cover type (>60% or >80% of a single land cover type) at the buffer zone with the largest effect size (script 05)

allInsects_totsample$hab50 = 'Mix' # samples that do not have more than 60% of one specific land type
allInsects_totsample$hab50[allInsects_totsample$Agriculture_500>=50]<-'Agriculture50'
allInsects_totsample$hab50[allInsects_totsample$Forest_1000>=50]<-'Forest50'
allInsects_totsample$hab50[allInsects_totsample$Urban_1000>=50]<-'Urban50'
table(allInsects_totsample$hab50) # notice the variation in sample size

library(indicspecies)
abund <- t(pa) # based on presence absence
land <- allInsects_totsample$hab50

# analysis with indicator value - explained in De Cáceres et al. (2010) - the accounts for unequal sample sizes
indval <-  multipatt(abund, land, func = "IndVal.g", duleg = TRUE, control = how(nperm=999)) # number of permutations affect the precision of the p-value
summary(indval, indvalcomp=TRUE)

# Component A: Component `A' is sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found. This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group.

#Component B: Component `B' is sample estimate of the probability of finding the species in sites belonging to the site group. This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group.

# Pearson's phi coefficient of association [Chytrý et al., 2002]. This coefficient is a measure of the correlation between two binary vectors. It is a good practice to correct the phi coefficient for the fact that some groups have more sites than others [Tichý and Chytrý, 2006]. To do that, we need to use func = "r.g" instead of func = "r":
phi <- multipatt(abund, land, func = "r.g", duleg = TRUE, control = how(nperm=999))
summary(phi)

# Correlation indices are used for determining the ecological preferences of species among a set of alternative site groups or site group combinations. Indicator value indices are used for assessing the predictive values of species as indicators of the conditions prevailing in site groups, e.g. for field determination of community types or ecological monitoring.

# finding indicators for forest, allowing for trio species combinations. We can discard species combinations with low indicator values by setting thresholds for components A and B.
sc <- indicators(X=abund, func = "IndVal.g", cluster=land, group="Forest50", max.order = 3, verbose=TRUE, At=0.8, Bt=0.5) # at least present in 80% of all forest sites and found in 50% of all samples for the land cover 
print(sc, sqrtIVt = 0.75) # only show the most useful indicators
plotcoverage(sc)

# urban
sc <- indicators(X=abund, func = "IndVal.g", cluster=land, group="Urban50", max.order = 3, verbose=TRUE, At=0.8, Bt=0.2) # at least present in 80% of all forest sites and found in 50% of all samples for the land cover 
print(sc) # only show the most useful indicators
plotcoverage(sc)

# Farmland
sc <- indicators(X=abund, func = "IndVal.g", cluster=land, group="Agriculture50", max.order = 3, verbose=TRUE, At=0.8, Bt=0.2) # at least present in 80% of all forest sites and found in 50% of all samples for the land cover 
print(sc) # only show the most useful indicators

# Farmland
sc <- indicators(X=abund, func = "IndVal.g", cluster=land, group="Mix", max.order = 3, verbose=TRUE, At=0.8, Bt=0.2) # at least present in 80% of all forest sites and found in 50% of all samples for the land cover 
print(sc) # only show the most useful indicators


