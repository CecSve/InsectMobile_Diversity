# make sure data is already loaded (script 07)

### vegan analysis #######

library(vegan)
### how many species in the samples? #######
envvar <- c("Forest_1000", "Agriculture_1000", "Urban_1000")
site_type <- allInsects %>% filter(maxLand_use %in% envvar) %>% select(SampleID, maxLand_use)
keep <- site_type$SampleID
otus <- subset(asvs, select=keep)

pa <- decostand(otus, method = "pa")
tpa <- t(pa)
sppr <- specnumber(tpa)

### analysis of variance ###########
# takes the same form as the usual models you'd see in R: response ~ dependent, data = environmental grouping
sppr_aov <- aov(sppr ~ maxLand_use, data = site_type)
summary(sppr_aov)

sppr_df <- sppr %>% 
  enframe() %>% 
  full_join(site_type, by = c("name" = "SampleID"))

sppr_df <- sppr_df %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Urban_1000",
    "Agriculture_1000",
    "Forest_1000")) 

table(sppr_df$maxLand_use)

plot_sppr <- ggplot(sppr_df, aes(x = maxLand_use, y = value, fill = maxLand_use)) +
  geom_boxplot() +
  scale_fill_manual(values = landuseCols[c(1,2,5)]) +
  scale_x_discrete(labels = c("Urban \n (n = 36)", "Farmland \n (n = 279)", "Forest \n (n = 45)")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("white"),
        panel.grid = element_line("grey90"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Land cover",
       y = "Number of species per site",
       title = "Richness")
plot_sppr

#### rarefy ##########
#pa <- decostand(asvs, method = "pa")
tpa <- t(otus)
sppr <- specnumber(tpa)
S <- specnumber(tpa) # observed number of species
raremax <- min(rowSums(tpa))
Srare <- rarefy(tpa, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
#rarecurve(tpa, step = 20, sample = raremax, col = "blue", cex = 0.6)

### how diverse are the communities ##########
shannondiv <- diversity(tpa)
head(shannondiv)

sppdiv_aov <- aov(shannondiv ~ maxLand_use, data = site_type)
summary(sppdiv_aov)

shandiv_df <- shannondiv %>% 
  # put all those calculations into a data frame
  enframe() %>% 
  # rename columns for ease of joining
  rename(SampleID = name,
         shan_div = value)

div_plot_df <- shandiv_df %>% 
  # join with site_type
  full_join(site_type, ., by = "SampleID") %>% 
  # group by landtype
  group_by(maxLand_use) %>% 
  # calculate mean and standard error of diversity
  summarize(mean = round(mean(shan_div), 2),
            err = sd(shan_div)/sqrt(length(shan_div))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))

plot_shandiv <- ggplot(div_plot_df, aes(x = maxLand_use, y = mean, fill = maxLand_use)) +
  geom_col(color = "black") +
  scale_fill_manual(values = landuseCols[c(1,2,5)]) +
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = maxLand_use, y = mean + err + 0.07, label = mean_label)) +
  scale_x_discrete(labels = c("Urban \n (n = 36)", "Farmland \n (n = 279)", "Forest \n (n = 45)")) +
  scale_y_continuous(limits = c(0, 3.5), expand = c(0,0)) +
  clean_background + 
  theme(legend.position = "none") +
  labs(x = "Land cover",
       y = "Mean Shannon diversity",
       title = "Shannon diversity")
plot_shandiv

### how different are the communities in species composition? ########
site_type <- allInsects %>% filter(maxLand_use %in% envvar) %>% select(SampleID, maxLand_use)
keep <- site_type$SampleID
otus <- subset(asvs, select=keep)

#pa <- decostand(otus, method = "pa")
tpa <- t(otus)

ins_perm <- adonis(tpa ~ maxLand_use, data = site_type)
ins_perm

### PCA ########
insPCA <- rda(tpa)
insPCA

PCAscores <- scores(insPCA, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  full_join(site_type, by = "SampleID")

PCAvect <- scores(insPCA, display = "species") %>% 
  as.data.frame()

plot_PCA <- ggplot() +
  geom_point(data = PCAscores, aes(x = PC1, y = PC2, color = maxLand_use)) +
  scale_color_manual(values = landuseCols[c(2,5,1)]) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = PCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  #geom_text(data = PCAvect, aes(x = PC1, y = PC2, label = rownames(PCAvect))) +
  clean_background +
  labs(x = "PC1 (29.71%)",
       y = "PC2 (23.11%)",
       title = "Principal Components Analysis") + xlim(-5, 25) + ylim(-100, 100)
plot_PCA

# trying the same with Laurens method
met <- allInsects
met$hab60 = 'Mix'
met$hab60[met$Agriculture_1000>=.6]<-'Agriculture60'
met$hab60[met$Forest_1000>=.6]<-'Forest60'
met$hab60[met$Urban_1000>=.6]<-'Urban60'

pa <- decostand(asvs, method = "pa")
tpa <- t(pa)
sppr <- specnumber(tpa)

site_type <- met %>% select(SampleID, hab60)

sppr_aov <- aov(sppr ~ hab60, data = site_type)
summary(sppr_aov)

sppr_df <- sppr %>% 
  enframe() %>% 
  full_join(site_type, by = c("name" = "SampleID"))

sppr_df <- sppr_df %>% mutate(
  hab60 = fct_relevel(
    hab60,
    "Urban60",
    "Agriculture60",
    "Forest60",
    "Mix")) 

table(sppr_df$hab60)

plot_sppr_lauren <- ggplot(sppr_df, aes(x = hab60, y = value, fill = hab60)) +
  geom_boxplot() +
  scale_fill_manual(values = landuseCols[c(1,2,5,6)]) +
  scale_x_discrete(labels = c("Urban \n (n = 10)", "Farmland \n (n = 182)", "Forest \n (n = 6)", "Mix \n (n = 166)")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("white"),
        panel.grid = element_line("grey90"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Land cover",
       y = "Number of species per site",
       title = "Richness")
plot_sppr_lauren

pa <- decostand(asvs, method = "pa")
tpa <- t(pa)
sppr <- specnumber(tpa)

site_type <- met %>% select(SampleID, hab60)

shannondiv <- diversity(tpa)
head(shannondiv)

sppdiv_aov <- aov(shannondiv ~ hab60, data = site_type)
summary(sppdiv_aov)

shandiv_df <- shannondiv %>% 
  # put all those calculations into a data frame
  enframe() %>% 
  # rename columns for ease of joining
  rename(SampleID = name,
         shan_div = value)

div_plot_df <- shandiv_df %>% 
  # join with site_type
  full_join(site_type, ., by = "SampleID") %>% 
  # group by landtype
  group_by(hab60) %>% 
  # calculate mean and standard error of diversity
  summarize(mean = round(mean(shan_div), 2),
            err = sd(shan_div)/sqrt(length(shan_div))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))

plot_shandiv <- ggplot(div_plot_df, aes(x = hab60, y = mean, fill = hab60)) +
  geom_col(color = "black") +
  scale_fill_manual(values = landuseCols[c(1,2,5,6)]) +
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = hab60, y = mean + err + 0.07, label = mean_label)) +
  scale_x_discrete(labels = c("Urban \n (n = 36)", "Farmland \n (n = 279)", "Forest \n (n = 45)", "Mix")) +
  scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  clean_background + 
  theme(legend.position = "none") +
  labs(x = "Land cover",
       y = "Mean Shannon diversity",
       title = "Shannon diversity")
plot_shandiv
