# run 03 script

### load libraries ##########
library(tidyverse)
library(vegan)
library(ggrepel)
library(ape)

### load data ##########
asvs <- read.delim("cleaned-data/DK_asvtable.txt")
taxonomy <- read.delim("cleaned-data/DK_taxonomy.txt",sep="\t")
data <- read.delim("cleaned-data/allInsects.txt", sep = " ")

### aligning data for analysis ###################
names(asvs) = gsub(pattern = "X*", replacement = "", x = names(asvs)) # remove the Xs that have emerged in the column headers that are numerical
asvs <- asvs[, grepl("IM18_*", names(asvs))] # only keep Danish samples from 2018 - be aware this removes blanks and negatives as well
names(asvs)

keep <- data$PCRID
asvs <- asvs[, (names(asvs) %in% keep)] # subset asvtable to only contain samples with metadata

keep <- names(asvs)
data <- data %>% filter(PCRID %in% keep)

# do all samples have seuqnces?
colSums(asvs) # yes
#asvs <- asvs[,colSums(asvs) > 0]

# remove empty asvs
rowSums(asvs)
asvs <- asvs[rowSums(asvs) > 0, ]

# remove the asvs from the taxonomy
keep <- rownames(asvs)
taxonomy <- taxonomy[rownames(taxonomy) %in% keep, ]

# rmeove asvs that does not fit with taxonomy of class = insecta, 99% identity score threshold etc
keep <- rownames(taxonomy)
asvs <- asvs[rownames(asvs) %in% keep, ]

# do all samples have seuqnces?
colSums(asvs) # no
asvs <- asvs[,colSums(asvs) > 0]

# subset data to the samples with reads
keep <- colnames(asvs)
data <- data %>% filter(PCRID %in% keep)

length(unique(data[["RouteID_JB"]])) # count how many routes were sampled - but notice some have received new routes
length(unique(data[["PilotID"]])) # count how many pilots that carried out the sampling
length(unique(data[["SampleID"]])) # count how many samples

### NMDS ##################
tasvs <- t(asvs) # samples as rows, taxa as columns
bray.data<-vegdist(tasvs, Type = "bray", binary = T) # presence absence data with dissimilarity index bray - use jaccard if presence absence

# due to issues in the ape package (negative eigen values comes before positive eigen values and causes an error), we have to use NMDS for ordiation vizualisation 
set.seed(41)
meta_mds <- metaMDS(bray.data, k = 2, trymax = 50)
stressplot(meta_mds)
plot(meta_mds)

class(meta_mds)
methods(class="metaMDS")

Y <- scores(meta_mds, display="sites")
plot(Y, type="n")
text(Y[,1], Y[,2], rownames(Y), col="red")
all(rownames(tasvs) == rownames(Y)) # make sure that rows of Y are in the same order as rows of norm.data

Y <- data.frame(Y)
Y$abundance <- rowSums(tasvs)
Y$labels <- data$maxLand_use
Y$colour <- data$Time_band

# plot all samples
ggplot(Y, aes(x=NMDS1, y=NMDS2, size = abundance, colour = labels, group = colour)) +
  geom_point() + 
  theme_minimal() + scale_size_continuous(range = c(3,8)) + labs(colour = "Highest proportion land cover on route")

# # The function envfit will add the environmental variables as vectors to the ordination plot
env_var <- data %>% select(Biomass, Time_band, yDay, Time_driven, Agriculture_1000, Forest_1000, Urban_1000, Heathland_1000, Wetland_1000, Open.uncultivated.land_1000, cnumberTime)

ef <- envfit(meta_mds, env_var, permu = 10000)
ef

# The two last columns are of interest: the squared correlation coefficient and the associated p-value
# Plot the vectors of the significant correlations and interpret the plot
plot(meta_mds, type = "t", display = "sites")
plot(ef, p.max = 0.05) # no explanatory variables are significant and have really low R2

# # Homogeneity of dispersion test - betadiversity
beta <- betadisper(bray.data, data$maxLand_use) # ## Calculate multivariate dispersions
plot(beta)
boxplot(beta) ## Draw a boxplot of the distances to centroid for each group
beta.HSD <- TukeyHSD(beta) ## Tukey's Honest Significant Differences. A Tukey's test can be done to see if and which groups differ in relation to their variances 
plot(beta.HSD, las = 1)
beta.HSD[["group"]]
anova(beta) # use group dispersions to perform an ANOVA test
permutest(beta, permutations = 999) ## Permutation test for F
adonis(bray.data ~ maxLand_use, data = data)
