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

#landuseCols <- c("#CC79A7", "#E69F00", "chartreuse3", "#D55E00", "#56B4E9", "#009E73") # colour friendly, ordered by land cover with heathland

#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified")

### Most abundant ASVs ###############

#Load data, taxnomy file and metadata into phyloseq formats

# prepare otu table
data_asvs <-
  allInsects_totasvs %>% rownames_to_column(var = "otuid") %>% mutate_each(funs(ifelse(. > 1, 1, 0)),-otuid) %>% column_to_rownames(var = "otuid") # convert into presence absence

OTU <-
  otu_table(data_asvs, taxa_are_rows = T) # changed taxa_are_rows to true for phyloseq to work

# prepare metadata
rownames(allInsects_totsample) <- NULL
data_meta <-
  allInsects_totsample %>% column_to_rownames(var = "SampleID")
sampledata <- sample_data(data_meta)

# prepare taxonomic data
names(taxonomy_Insecta)
taxon <-
  taxonomy_Insecta %>% rownames_to_column(var = "otuid") %>% select(otuid, kingdom, phylum, class, order, family, genus, species) %>% column_to_rownames(var = "otuid")
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

# for DK jitter x and y slightly - fix later
allInsects_totsample$x2 <- allInsects_totsample$utm_x + rnorm(length(allInsects_totsample$utm_x),0,10)
allInsects_totsample$y2 <- allInsects_totsample$utm_y + rnorm(length(allInsects_totsample$utm_y),0,10)

allphyseqdata <- merge(data_long, allInsects_totsample, by.x = "samples", by.y = "SampleID") # merge the data set so we can visualise based on max land use

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

# load package to visualise table in html format
#library(sjPlot)
#library(sjmisc)
#library(sjlabelled)

#tab_model(modelphy)

#subphy <- allphyseqdata %>% select(maxLand_use, order, species, abundance.x) #%>% drop_na(species) # perhaps choose to drop species with no assigned name (if we only want named species, for now it is interesting to see if the most abundant ASVs are poorly assigned)
#subphy <- subphy %>% replace_na(list(species = "sp.")) %>% mutate(order_species = paste(order, species, sep = ":"))

# account for sample size and only include the three max land cover with most samples
#table(subphy$maxLand_use)

#agr <- subphy %>% filter(maxLand_use == "Agriculture_1000")
#agr$relabun <- agr$abundance.x/13950 

#forest <- subphy %>% filter(maxLand_use == "Forest_1000")
#forest$relabun <- forest$abundance.x/2250

#urb <- subphy %>% filter(maxLand_use == "Urban_1000")
#urb$relabun <- urb$abundance.x/1800

#dat <- rbind(agr, forest, urb)

# summarise by land cover and species
#dat <- dat %>% select(maxLand_use, order_species, relabun) # only select the relevant columns to plot

#dat <- dat %>% group_by(maxLand_use, order_species) %>% summarise_all(sum)

#dat %>%  group_by(order_species) %>% arrange(desc(relabun)) %>% slice(1:10) %>% ggplot(aes(x=maxLand_use, y=order_species, fill= relabun)) + geom_tile() + scale_fill_viridis_c()  

#dat <- as_tibble(dat)
# rename land cover variables
plotdata <- dat %>% mutate(
  maxLand_use = dplyr::recode(
    maxLand_use,
    "Agriculture_1000" = "Farmland",
    "Urban_1000" = "Urban",
    "Forest_1000" = "Forest"
  )
)

plotdata %>% 
  ggplot(aes(x=maxLand_use, y=order_species, fill= relabun)) + 
  geom_tile() + scale_fill_viridis_c(guide = guide_legend(
    direction = "horizontal",
    title.position = "top",
    label.position = "bottom"
  )) + theme_pubclean() + theme(legend.position = "bottom",legend.title = element_text(size = 10), axis.title = element_text(size = 12, face = "bold")) + labs(fill = "Abundance in samples \n(sample size corrected)") + xlab("Highest propotional land cover at 1000 m") + ylab("Insect order and species") + 
  scale_x_discrete(expand=c(0,0.1)) + 
  scale_y_discrete(expand=c(0,0.1))

mod <- glm(relabun~maxLand_use, data = dat)
summary(mod)

qqnorm(resid(mod))

gls1.alleffects <- allEffects(mod)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(mod)
plot(eall.lm1, lines=list(multiline=F))

### GLS on all data with random effect ######
physeqdata <- phyloseq_to_df(physeq, addtax = T, addtot = T, addmaxrank = F,
                             sorting = "abundance") # convert the phyloseq object into a dataframe for plotting
names(physeqdata)
#mydata <- mutate_each(test, funs(ifelse(. > 1, 1, 0)), -c(OTU, kingdom, phylum, class, order, family, genus, species, Total)) # convert into presence absence
data_long <- tidyr::gather(physeqdata, samples, abundance, P1.2A:P99.1B)

test <- merge(data_long, allInsects_totsample, by.x = "samples", by.y = "SampleID")

gls1 <- lme(abundance.x ~ Agriculture_1000 + Urban_1000 + Open.uncultivated.land_1000 + Heathland_1000 + Forest_1000 + Wetland_1000 + 
              Time_band +
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            data=test)

summary(gls1)
qqnorm(resid(gls1))

# plot relative effects of each variable - very slow
gls1.alleffects <- allEffects(gls1)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(gls1)
#plot(eall.lm1, lines=list(multiline=TRUE))

### ggplot effect plot ####
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

temp <- effectdata$Open.uncultivated.land_1000
temp$landcover <- "Open.uncultivated.land_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Open.uncultivated.land_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

temp <- effectdata$Heathland_1000 
temp$landcover <- "Heathland_1000"
heath <- temp %>% 
  dplyr::rename(
    propcover = Heathland_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

temp <- effectdata$Forest_1000 
temp$landcover <- "Forest_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Forest_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

temp <- effectdata$Wetland_1000 
temp$landcover <- "Wetland_1000"
wet <- temp %>% 
  dplyr::rename(
    propcover = Wetland_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

effect_data <- rbind(urb, farm, grass, heath, wet, forest)

# Visualization
effectplot <- effect_data %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000",
    "Open.uncultivated.land_1000",
    "Heathland_1000",
    "Wetland_1000",
    "Forest_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban", "Farmland",
      "Grassland",
      "Heathland",
      "Wetland",
      "Forest" 
    )
  ) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.position = "bottom"
  ) + scale_x_continuous(
    limits = c(0, 1),
    labels = function(x)
      paste0(x * 100, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        x = "Land cover at 1000 m",
        y = "Predicted effect change in ASV abundance",
        colour = "Land cover type"
      ) + scale_fill_manual(values =  landuseCols) + guides(colour = guide_legend(nrow = 1))

save_plot("plots/DK_predictedeffect_physeqabun_allcovers.png", effectplot, base_width = 10, base_height = 6)

### ampvis ################

taxa <- taxonomy_Insecta %>% rownames_to_column(var = "SampleID") %>% select(SampleID, kingdom, phylum, class, order, family, genus, species)

pa.otus <- rownames_to_column(pa, var = "SampleID")

otutable <- left_join(pa.otus, taxa, by = "SampleID")

otutable  <- column_to_rownames(otutable, var = "SampleID") # set sequenceids as rownames
# metadata must have sample data in the first column and column classes matter 
otutable <- otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

str(allInsects_totsample)
envvar <- c("Forest_1000", "Agriculture_1000", "Urban_1000")
env_data <- allInsects_totsample %>% filter(maxLand_use %in% envvar)
table(env_data$maxLand_use)

amp <- amp_load(otutable = otutable, metadata = env_data)

# heatmap (normalise gives the proportion of the sample)
show_col(viridis_pal()(20))

amp_heatmap(
  data = amp,
  group_by = "maxLand_use",
  normalise = F,
  tax_aggregate = "Species",
  tax_add = "Order",
  plot_values = T,
  #plot_colorscale = "log10",
  color_vector = c("#FDE725FF","#29AF7FFF","#32648EFF","#440154FF"),
  tax_empty = "remove",
  tax_show = 50,
)+ guides(fill = guide_legend(title = "Relative ASV abundance", reverse = T)) 

# extracting the heatmap object data and normalising by sample size prior to plotting - probably not the best way
dat <- test$data

ggplot(data = dat, aes(x=Group, y=Display, fill= Abundance)) + 
  geom_tile() + scale_fill_viridis_c()

agr <- dat %>% filter(Group == "Agriculture_1000")
agr$relabun <- agr$Sum/279

forest <- dat %>% filter(Group == "Forest_1000")
forest$relabun <- forest$Sum/45

urb <- dat %>% filter(Group == "Urban_1000")
urb$relabun <- urb$Sum/36

dat <- rbind(agr, forest, urb)

ggplot(data = dat, aes(x=Group, y=Display, fill= relabun)) + 
  geom_tile() + scale_fill_viridis_c()

# extracting the data from the amp object

d_long <- amp_export_long(amp, metadata_vars = "maxLand_use", tax_levels = c("OTU", "Order", "Species"))
d_long$Species <- sub("^$", NA, d_long$Species)
sub_d_long <- d_long %>% drop_na(Species) %>% select(maxLand_use, Order,Species, count)

dat <- sub_d_long %>%
  group_by(maxLand_use, Order, Species) %>%   
  summarise_all(sum)

table(dat$maxLand_use)
agr <- dat %>% filter(maxLand_use == "Agriculture_1000")
agr$relabun <- agr$count/588969

forest <- dat %>% filter(maxLand_use == "Forest_1000")
forest$relabun <- forest$count/94995

urb <- dat %>% filter(maxLand_use == "Urban_1000")
urb$relabun <- urb$count/75996

dat <- rbind(agr, forest, urb)

tjek <- dat %>% filter_all(all_vars(!is.na(.)))
tjek2 <- tjek %>% arrange(desc(relabun)) %>%
  slice(1:30)

ggplot(data = dat, aes(x=maxLand_use, y=Order, fill= count)) + 
  geom_tile() + scale_fill_viridis_c()

dat %>% arrange(desc(count))%>% slice(1:10) %>% 
  ggplot(., aes(x=maxLand_use, y=Species, fill= count)) + geom_tile() 

### venn diagram: https://madsalbertsen.github.io/ampvis2/reference/amp_venn.html #####
max_venn <- amp_venn(amp,
                     group_by = "maxLand_use",
                     cut_a = 0.1,
                     cut_f = 1)

# boxplot/point plot
amp_boxplot(
  amp,
  group_by = "maxLand_use",
  sort_by = "mean",
  plot_type = "boxplot",
  tax_aggregate = "Family",
  tax_show = 10,
  tax_empty = "remove",
  plot_log = TRUE
) + ylab('Relative abundance') + guides(colour = guide_legend(title = "?", reverse = T)) + scale_color_manual(values = landuseCols)

# richness indices
amp_alphadiv <-
  amp_alphadiv(
    amp,
    measure = c("observed", "shannon", "simpson", "invsimpson"),
    richness = TRUE,
    #rarefy = 200
  )

amp_alphadiv_select <- amp_alphadiv %>% select(maxLand_use, ObservedOTUs, Shannon, Simpson, Chao1)

# plot richness indices on log10 transformed y-axis
amp_alphadiv_select %>% tidyr::gather("id", "Richness", 2:5) %>% ggplot(., aes(maxLand_use, Richness, fill = maxLand_use)) + geom_boxplot(show.legend = F) + facet_wrap(~ id, ncol = 5) + scale_y_log10() + scale_fill_manual(values = landuseCols[c(2,5,1)]) + theme_pubclean() + labs(fill = "Dominant land cover") + xlab("Dominant land cover") + ylab("Richness (log-transformed)")

amp_ordinate(amp,
             type = "pcoa",
             filter_species = 0,
             transform = "none",
             distmeasure = "bray",
             sample_color_by = "maxLand_use",
             sample_colorframe = "maxLand_use",
             sample_plotly = "maxareaProp")

### rarefaction curve on read abundance ##############
# join otus and taxonomy
otus <- rownames_to_column(totsample_asvs, var = "otuid")
taxa <- rownames_to_column(taxonomy_insecta, var = "otuid")
taxa <- taxa %>% select(otuid, kingdom, phylum, class, order, family, genus, species)
otutable <- left_join(otus, taxa, by = "otuid")

otutable <- otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)
# revert otuids back to rownames if needed
otus <- column_to_rownames(otus, var = "otuid")
otutable  <- column_to_rownames(otutable, var = "otuid") # set sequenceids as rownames

# metadata must have sample data in the first column and column classes matter 
str(allInsects_totsample)
#metadata <- rownames_to_column(metadata, var = "otuid")

# load amp data
#otutable <- otutable %>% dplyr::rename(Kingdom = kingdom, Phylum = phylum, Class = class, Order = order, Family = family, Genus = genus, Species = species)

amp <- amp_load(otutable = otutable, metadata = env_data)

# ordination

amp_ordinate(amp,
             type = "PCOA",
             filter_species = 0,
             transform = "none",
             distmeasure = "jsd",
             sample_color_by = "maxLand_use",
             sample_colorframe = "maxLand_use",
             sample_plotly = "maxareaProp")

# rarefaction curve
rarecurve <-
  amp_rarecurve(
    amp,
    facet_by = "maxLand_use",
    stepsize = 100,
    color_by = "maxLand_use",
    #facet_scales = "free"
  )

min(colSums(otus[, 2:364]))
max(colSums(otus[, 2:364]))

rarecurve + geom_line(size = 1, na.rm = TRUE) + scale_x_continuous(
  labels = scales::number,
  limits = c(0, 100000),
  breaks = seq(from = 0, to = 100000, by = 5000)
) + theme(legend.position = "bottom") + scale_colour_manual(values = landuseCols[c(2,6,1)], labels = c("Farmland", "Forest", "Urban")) + theme_pubclean() + 
  theme(
    strip.text.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) + ylab("Number of observed ASVs") + xlab("Land cover at 1000 m") + guides(x =  guide_axis(angle = 90))

# coverage
mean(allInsects$Forest_1000)
mean(allInsects$Forest_500)
mean(allInsects$Forest_250)
mean(allInsects$Forest_50)
(0.1585302+0.1556052+0.1492961+0.1429663)/4

mean(allInsects$Urban_1000)
mean(allInsects$Urban_500)
mean(allInsects$Urban_250)
mean(allInsects$Urban_50)
(0.1274353+0.1431824+0.1691167+0.2395598)/4

mean(allInsects$Agriculture_1000)
mean(allInsects$Agriculture_500)
mean(allInsects$Agriculture_250)
mean(allInsects$Agriculture_50)
(0.536347+0.5293759+0.5190831+0.4282332)/4

