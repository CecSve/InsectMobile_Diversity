# run script 04 DK ordination community composition richness per sample script first

# phyloseq and ampvis

### load libraries ###############

library(phyloseq)
library(ampvis2)
library(metagMisc) # for extracting phyloseq objects into dataframes
library(ggpubr)
library(effects)
library(lme4)
library(lmerTest)
library(nlme)
library(cowplot)
library(tidyverse)

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

### Most abundant ASVs ###############

# run script 04 prior to analysis

OTU <-
  otu_table(pa, taxa_are_rows = T) # changed taxa_are_rows to true for phyloseq to work

# prepare metadata
rownames(allInsects) <- NULL
data_meta <-
  allInsects %>% column_to_rownames(var = "SampleID")
sampledata <- sample_data(data_meta)

# prepare taxonomic data
names(taxonomy_Insecta)
taxon <-
  taxonomy_Insecta %>% rownames_to_column(var = "otuid") %>% dplyr::select(otuid, kingdom, phylum, class, order, family, genus, species) %>% column_to_rownames(var = "otuid")
taxon <- as.matrix(taxon) #Makes a dataframe
TAX <- tax_table(taxon)

setdiff(taxa_names(OTU), taxa_names(TAX)) # looks right

#Create create phyloseq object named physeq
physeq <- phyloseq(OTU, sampledata,TAX)
physeq <- merge_phyloseq(physeq, sampledata)

OTUsub = names(sort(taxa_sums(physeq), TRUE)[1:20]) # select the most abundant taxa
bsub=prune_taxa(OTUsub, physeq) # only select the 50 most abundant ASVs in the physeq object

#plot_bar(bsub, "maxLand_use", fill = "order") # plot on all land covers - note that sample size has not been accounted for at this stage

physeqdata <- phyloseq_to_df(bsub, addtax = T, addtot = F, addmaxrank = T,
                             sorting = "abundance") # convert the phyloseq object into a dataframe for plotting
table(physeqdata$LowestTaxRank)

names(physeqdata)
#mydata <- mutate_each(test, funs(ifelse(. > 1, 1, 0)), -c(OTU, kingdom, phylum, class, order, family, genus, species, Total)) # convert into presence absence
data_long <- tidyr::gather(physeqdata, samples, abundance, P1.2A:P99.1B) # reorder data so samples are in one column instead of multiple columns
#data_long

allphyseqdata <- merge(data_long, allInsects, by.x = "samples", by.y = "SampleID") # merge the data set so we can visualise based on max land use

#allphyseqdata$genusF <- as.factor(allphyseqdata$genus)
subphy <- allphyseqdata %>% tidyr::replace_na(list(species = "sp.")) %>% mutate(order_species = paste(order, family, genus, species, sep = ":"))

modelphy <-
  lmerTest::lmer(
    abundance ~ Agriculture_50:order_species +  Forest_50:order_species + Wetland_50:order_species + Urban_50:order_species + Open.uncultivated.land_50:order_species + (1 | PilotID / RouteID),
    na.action = na.omit,
    data = subphy
  )

summary(modelphy)
r.squaredGLMM(modelphy)
#R2m       R2c
#0.07020821 0.1541349
confint_modelphy<- confint.merMod(modelphy)
test <- as.data.frame(confint_modelphy)
write.table(test, file = "confidence_modelphy.txt", sep = "\t")
