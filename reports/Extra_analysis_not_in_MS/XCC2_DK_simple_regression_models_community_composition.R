# Simple regression models on community composition

# Supplementary material analysis for InsectMobile diversity paper

#### Load required libraries ################################################################
library(cowplot)
library(ggplot2)
#library(wesanderson)
library(ggpubr)
library(scales)
library(ggpmisc)
library(grid)
library(gridExtra)
library(corrplot)
library(car)
library(lme4)
library(lmerTest)
library(nlme)
library(effects)

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseCols <- c("#CC79A7", "#E69F00", "chartreuse3", "#D55E00", "#56B4E9", "#009E73") # colour friendly, ordered by land cover with heathland

landuseOrder <- c("Urban","Farmland","Grassland", "Heathland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified")

### get effect function #############################################################

#function to extract summary components of simpel buffer effect models
getEffect <- function(model){
  coefs <- summary(model)$coef[2,]
  temp <- confint(model)
  cis <- temp[5,]
  data.frame(t(coefs),t(cis))
}

# fitModels <- function(variable){
#   myformula()
#   
# }

### community composition: simple regression plot land cover##############################################

# urban
fit = lm((allInsects_totsample$distdif) ~ allInsects_totsample$Urban_1000)
coef(summary(fit))[2,4] 

qU <- ggplot(allInsects_totsample,aes(x=(Urban_1000),y=(distdif)))+
  geom_point(col=landuseCols[1])+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
  xlab("Urban cover") +ylab("") + labs(subtitle = "A") + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%")) 

fit = lm((allInsects_totsample$distdif) ~ (allInsects_totsample$Agriculture_1000))
coef(summary(fit))[2,4] 

qF <- ggplot(allInsects_totsample,aes(x=Agriculture_1000,y=(distdif)))+
  geom_point(col=landuseCols[2])+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
  xlab("Farmland cover") +ylab("") + labs(subtitle = "B") + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%")) 

fit = lm((allInsects_totsample$distdif) ~ allInsects_totsample$Open.uncultivated.land_1000)
coef(summary(fit))[2,4] # non-significant

qD <- allInsects_totsample %>% filter(Open.uncultivated.land_1000 < 0.4)  %>% ggplot(aes(x=Open.uncultivated.land_1000, y=distdif))+
  geom_point(col=landuseCols[3])+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
  xlab("Grassland cover") +ylab("") + labs(subtitle = "C") + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%")) 

fit = lm((allInsects_totsample$distdif) ~ allInsects_totsample$Heathland_1000)
coef(summary(fit))[2,4] # non-significant

qH <- allInsects_totsample  %>% ggplot(aes(x=Heathland_1000, y=distdif))+
  geom_point(col=landuseCols[4])+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
  xlab("Heathland cover") +ylab("") + labs(subtitle = "C") + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%")) 

fit = lm((allInsects_totsample$distdif) ~ allInsects_totsample$Wetland_1000)
coef(summary(fit))[2,4] # non-significant

qW <- ggplot(allInsects_totsample,aes(x=Wetland_1000,y=distdif))+
               geom_point(col=landuseCols[5])+
               scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
               theme_bw() +
               geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
               xlab("Wetland cover") +ylab("") + labs(subtitle = "D") + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%")) 

fit = lm(allInsects_totsample$distdif ~ allInsects_totsample$Forest_1000)
coef(summary(fit))[2,4] # non-significant

qFo <- ggplot(allInsects_totsample,aes(x=Forest_1000,y=distdif))+
  geom_point(col=landuseCols[6])+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_smooth(method="lm", color="darkgrey", linetype = "solid")+
  xlab("Forest cover") +ylab("") + labs(subtitle = "E")+ scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x*100, "%"))

# will not include unspecified land cover in the plot
#qUns <- ggplot(allInsects_totsample,aes(x=sqrt(Unspecified.land.cover_1000),y=(Biomass+1)))+geom_point(col="darkgrey")+scale_y_log10() +theme_bw() +geom_smooth(method="lm",color="grey70")+xlab("Unspecified land cover") +ylab("") + labs(subtitle = "F")+ scale_x_continuous(labels = function(x) paste0(x*100, "%"))

y.grob <- textGrob("Variation in community composition (distance from zero)", 
                   gp=gpar(fontface="bold", col="black", fontsize=10), rot=90)

plot <- plot_grid(qU,qF,qD,qH,qW,qFo,ncol=1)
plot <- grid.arrange(arrangeGrob(plot, left = y.grob))
#ggsave("plots/DK_Landcover_percent.png",width=3,height=8)
save_plot("plots/DK_Landcover_percent_community_comp.png", plot, base_width = 3, base_height = 8)

### DK plot buffers#################################################

### urban ####
mini_data <- allInsects_totsample %>% select(PCRID,  Urban_50, Urban_250, Urban_500, Urban_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
dt_long

# renmae
landcover_names <- list(
  'Urban_50'="50 m",
  'Urban_250'="250 m",
  'Urban_500'="500 m",
  'Urban_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
u <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Urban_50",
    "Urban_250",
    "Urban_500",
    "Urban_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[1]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Urban land cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(
    limits = c(0, 1),
    labels = function(x)
      paste0(x * 100, "%")) 

# plot with regression line
upubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Urban_50",
    "Urban_250",
    "Urban_500",
    "Urban_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[1], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.5, label.y = 0.5) + xlab("Urban land cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### farmland ############
mini_data <- allInsects_totsample %>% select(PCRID,  Agriculture_50, Agriculture_250, Agriculture_500, Agriculture_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
dt_long

# renmae
landcover_names <- list(
  'Agriculture_50'="50 m",
  'Agriculture_250'="250 m",
  'Agriculture_500'="500 m",
  'Agriculture_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
a <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Agriculture_50",
    "Agriculture_250",
    "Agriculture_500",
    "Agriculture_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[2]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Farmland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%")) 

# plot with regression line
apubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Agriculture_50",
    "Agriculture_250",
    "Agriculture_500",
    "Agriculture_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[2], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.5, label.y = 0.5) + xlab("Farmland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### grassland ###############
mini_data <- allInsects_totsample %>% select(PCRID,  Open.uncultivated.land_50, Open.uncultivated.land_250, Open.uncultivated.land_500, Open.uncultivated.land_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
#dt_long

# rename
landcover_names <- list(
  'Open.uncultivated.land_50'="50 m",
  'Open.uncultivated.land_250'="250 m",
  'Open.uncultivated.land_500'="500 m",
  'Open.uncultivated.land_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
g <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Open.uncultivated.land_50",
    "Open.uncultivated.land_250",
    "Open.uncultivated.land_500",
    "Open.uncultivated.land_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[3]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Grassland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%")) 

# plot with regression line
gpubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Open.uncultivated.land_50",
    "Open.uncultivated.land_250",
    "Open.uncultivated.land_500",
    "Open.uncultivated.land_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[3], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.2, label.y = 0.5) + xlab("Grassland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### heathland ################
mini_data <- allInsects_totsample %>% select(PCRID,  Heathland_50, Heathland_250, Heathland_500, Heathland_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
#dt_long

# renmae
landcover_names <- list(
  'Heathland_50'="50 m",
  'Heathland_250'="250 m",
  'Heathland_500'="500 m",
  'Heathland_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
h <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Heathland_50",
    "Heathland_250",
    "Heathland_500",
    "Heathland_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[4]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Heathland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%")) 

# plot with regression line
hpubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Heathland_50",
    "Heathland_250",
    "Heathland_500",
    "Heathland_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[4], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.25, label.y = 0.5) + xlab("Heathland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### wetland #############
mini_data <- allInsects_totsample %>% select(PCRID,  Wetland_50, Wetland_250, Wetland_500, Wetland_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
#dt_long

# renmae
landcover_names <- list(
  'Wetland_50'="50 m",
  'Wetland_250'="250 m",
  'Wetland_500'="500 m",
  'Wetland_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
w <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Wetland_50",
    "Wetland_250",
    "Wetland_500",
    "Wetland_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[5]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Wetland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%")) 

# plot with regression line
wpubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Wetland_50",
    "Wetland_250",
    "Wetland_500",
    "Wetland_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[5], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.2, label.y = 0.5) + xlab("Wetland cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### forest ###########################
mini_data <- allInsects_totsample %>% select(PCRID,  Forest_50, Forest_250, Forest_500, Forest_1000, distdif)

# wide to long format
dt_long <- mini_data %>% 
  tidyr::gather(landcover_type, propcover, -PCRID, -distdif)
#dt_long

# renmae
landcover_names <- list(
  'Forest_50'="50 m",
  'Forest_250'="250 m",
  'Forest_500'="500 m",
  'Forest_1000' = "1000 m"
)

# create labels for plotting
landcover_labeller <- function(variable,value){
  return(landcover_names[value])
}

# plot without regression line
f <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Forest_50",
    "Forest_250",
    "Forest_500",
    "Forest_1000"
  )
) %>% ggplot() +
  geom_point(aes(x=propcover,y=distdif), show.legend = F, size = 3, colour = landuseCols[6]) +
  facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) +
  theme_pubclean() + xlab("Forest cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+ scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%")) 

# plot with regression line
fpubr <- dt_long %>% mutate(
  landcover_type = fct_relevel(
    landcover_type,
    "Forest_50",
    "Forest_250",
    "Forest_500",
    "Forest_1000"
  )) %>% ggscatter(x = "propcover", y = "distdif",
                   color = landuseCols[6], shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.45, label.y = 0.5) + xlab("Forest cover") + ylab("") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(     limits = c(0, 1),     labels = function(x)       paste0(x * 100, "%"))  + facet_wrap(. ~ landcover_type, labeller = landcover_labeller, ncol = 2) + theme_pubclean()

### merge plots ####################
# plot with no regression line
#plot <- plot_grid(u,a,g,h,w,f,ncol=2, labels = c("A", "B", "C", "D", "E", "F"), align = "hv")

plot <- plot_grid(
  # column 1
  plot_grid(u,a,g, ncol = 1, labels = c('A', 'B', 'C')) #+ theme(plot.background = element_rect(color = "black"))
  
  # column 2
  , plot_grid(h,w,f, ncol = 1, labels = c('D', 'E', 'F')) #+ theme(plot.background = element_rect(color = "black"))
  , ncol = 2, align = "hv")

y.grob <- textGrob("Variation in community composition", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

test <- grid.arrange(arrangeGrob(plot, left = y.grob))
grid.newpage()
grid.draw(test)
#ggdraw(add_sub(plot, "Label", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))
ggsave('plots/Landcover_percent_buffer_noreg.png', plot = test, width=50,height=35, units = "cm", dpi=300)
#ggsave("plots/Landcover_percent_buffer_noreg.png",width=25,height=15)

# plot with regression line
plot <- plot_grid(
  # column 1
  plot_grid(upubr,apubr,gpubr, ncol = 1, labels = c('A', 'B', 'C')) #+ theme(plot.background = element_rect(color = "black"))
  
  # column 2
  , plot_grid(hpubr, wpubr,fpubr, ncol = 1, labels = c('D', 'E', 'F')) #+ theme(plot.background = element_rect(color = "black"))
  , ncol = 2, align = "hv")

#plot_grid(upubr,apubr,gpubr, hpubr, wpubr,fpubr,ncol=2, labels = c("A", "B", "C", "D", "E", "F"))
test <- grid.arrange(arrangeGrob(plot, left = y.grob))
grid.newpage()
grid.draw(test)
#ggdraw(add_sub(plot, "Label", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5))
ggsave('plots/Landcover_percent_with_stats.png', plot = test, width=50,height=35, units = "cm", dpi=300)
#ggsave("plots/Landcover_percent_with_stats.png",width=6,height=20)

### intensive vs organic farming predicted community composition ###############################
# Intensive
hist(allInsects_totsample$distdif) # no transformation needed
hist(allInsects_totsample$richness)
hist(log(allInsects_totsample$richness))
hist(sqrt(allInsects_totsample$richness)) # normal distribution with square root transformation

lme1000 <- lmer(distdif ~ 
                  Intensiv_1000 +
                  Ekstensiv_1000 +
                  Semi.intensiv_1000 +
                  #Markblok_1000 +
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_totsample)
summary(lme1000)
newData = data.frame(Intensiv_1000=0.5, Ekstensiv_1000 = 0.5, Semi.intensiv_1000 = 0.5, Markblok_1000 = 0.5, cStops=0, cyDay = 0, Time_band = "midday",cnumberTime = 0)

#make predictions
Intensiv1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

predFun <- function(fit) {
  predict(fit,newData,re.form=NA)
}

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
conventional <- bb[["data"]]
Intensiv2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

Intensiv <- cbind(Intensiv1, Intensiv2)
Intensiv <- as.data.frame(Intensiv)
colnames(Intensiv)
names(Intensiv)[1] <- "predCommunityVariation"
names(Intensiv)[2] <- "lowCI"
names(Intensiv)[3] <- "highCI"
row.names(Intensiv) <- "Intensive"

# propOeko
lme1000 <- lmer(distdif ~ 
                  Intensiv_organic_1000 +
                  Ekstensiv_organic_1000 +
                  Semi.intensiv_organic_1000 +
                  #Markblok_organic_1000 +  
                  Time_band + 
                  Time_band:cnumberTime +
                  log(cStops+1) + cyDay + 
                  (1|RouteID_JB) + (1|PilotID), 
                data=allInsects_totsample)
summary(lme1000)
newData = data.frame(Intensiv_organic_1000=0.5, Ekstensiv_organic_1000 = 0.5, Semi.intensiv_organic_1000 = 0.5, Markblok_organic_1000 = 0.5,
                     cStops=0,
                     cyDay = 0,
                     Time_band = "midday",
                     cnumberTime = 0)

#make predictions
propOeko1 <- t(as_tibble(exp(predict(lme1000,newdata=newData,re.form=NA))))

predFun <- function(fit) {
  predict(fit,newData,re.form=NA)
}

bb <- bootMer(lme1000,nsim=1000,FUN=predFun,seed=101)
organic <- bb[["data"]]
propOeko2 <- t(as_tibble(exp(quantile(bb$t,c(0.025,0.975)))))

propOeko <- cbind(propOeko1, propOeko2)
propOeko <- as.data.frame(propOeko)
colnames(propOeko)
names(propOeko)[1] <- "predCommunityVariation"
names(propOeko)[2] <- "lowCI"
names(propOeko)[3] <- "highCI"
row.names(propOeko) <- "Organic"

predConfData <- rbind(Intensiv, propOeko)
predConfData <- rownames_to_column(predConfData, var = "Farmland_type")

# plot
p <- predConfData %>% ggplot(aes(Farmland_type, predCommunityVariation, colour = Farmland_type))
finalplot <- p + geom_pointrange(aes(ymin = lowCI, ymax = highCI), size =1.5) + scale_colour_manual(values = c("#F09018", "#E3B622" )) + theme_minimal_grid() + theme(legend.title = element_blank(), legend.key = element_rect(size = 0.1), legend.key.size = unit(1, 'cm')) + labs(x = "\nFarming system", y = "Predicted community variation and 95% CIs\n", subtitle = "A") + theme(plot.subtitle = element_text(size = 20, face = "bold")) #+ scale_y_log10()

save_plot("plots/DK_predicted_community_variation_farmtype.png", finalplot, base_width = 8, base_height = 5)

### hertil###
### combining the predicted data to re-run model and calculate effects ####
predeffect <- merge(conventional, organic)

#predeffect <- predeffect %>% rename(Biomass = `log(Biomass + 1)`) # be mindful that biomass is +1 and logtransformed here, the sme for stops
predeffect <- predeffect %>% dplyr::rename(cStops = `log(cStops + 1)`) 

# run model (marklok data not useable)
gls1 <- lme(distdif ~ Intensiv_1000 +
              Ekstensiv_1000 +
              Semi.intensiv_1000 +
              #Markblok_1000 + 
              Intensiv_organic_1000 +
              Ekstensiv_organic_1000 +
              Semi.intensiv_organic_1000 +
              #Markblok_organic_1000 + 
              Time_band +
              Time_band:cnumberTime + 
              cStops + 
              cyDay,
            random=~1|PilotID/RouteID_JB,
            data=predeffect)

summary(gls1)

gls1.alleffects <- allEffects(gls1)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(gls1)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(gls1, ~ Intensiv_1000 + Ekstensiv_1000 + Semi.intensiv_1000 + Intensiv_organic_1000 + Ekstensiv_organic_1000 + Semi.intensiv_organic_1000 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

### ggplot effect plot ####
temp <- effectdata$Intensiv_1000
temp$landcover <- "Intensiv_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Intensiv_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Ekstensiv_1000
temp$landcover <- "Ekstensiv_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Ekstensiv_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

# Open.uncultivated.land
temp <- effectdata$Semi.intensiv_1000
temp$landcover <- "Semi.intensiv_1000"
grass <- temp %>% 
  dplyr::rename(
    propcover = Semi.intensiv_1000
  ) %>% select(landcover, propcover, fit, se, lower, upper)

# 
temp <- effectdata$Intensiv_organic_1000
temp$landcover <- "Intensiv_organic_1000"
wet <- temp %>% 
  dplyr::rename(
    propcover = Intensiv_organic_1000
  )%>% select(landcover, propcover, fit, se, lower, upper)

# 
temp <- effectdata$Ekstensiv_organic_1000
temp$landcover <- "Ekstensiv_organic_1000"
forest <- temp %>% 
  dplyr::rename(
    propcover = Ekstensiv_organic_1000
  ) %>% select(landcover, propcover, fit, se, lower, upper)

# 
temp <- effectdata$Semi.intensiv_organic_1000
temp$landcover <- "Semi.intensiv_organic_1000"
semiint <- temp %>% 
  dplyr::rename(
    propcover = Semi.intensiv_organic_1000
  ) %>% select(landcover, propcover, fit, se, lower, upper)


test <- rbind(urb, farm, grass, wet, forest, semiint)

# Visualization
effectplot <- test %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Intensiv_1000",
    "Intensiv_organic_1000",
    "Semi.intensiv_1000",
    "Semi.intensiv_organic_1000",
    "Ekstensiv_1000",
    "Ekstensiv_organic_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_grey(
    labels = c(
      "Intensive cover", "Intensive organic cover",
      "Semi-intensive",
      "Semi-intensive organic cover",
      "Exstensive cover",
      "Exstensive organic cover" 
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
        x = "Land cover",
        y = "Predicted effect change in community composition",
        subtitle = "A",
        colour = "Land cover type"
      ) + scale_fill_grey()

#save_plot("plots/DK_predictedeffect_communitycomp_farming_landcover.png", effectplot, base_width = 10, base_height = 6)

### buffer effect plots ######################
# NB! changed cTL to cStops since more data for DK 
str(allInsects_totsample)
#agriculture
hist(allInsects_totsample$Agriculture_500)
hist(sqrt(allInsects_totsample$Agriculture_500)) # does not help
lme50 <- lmer(distdif ~ (Agriculture_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(distdif ~ (Agriculture_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(distdif ~ (Agriculture_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(distdif ~ (Agriculture_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outAgri <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outAgri <- as.data.frame(outAgri)
outAgri$Buffer <- c(50,250,500,1000)
outAgri$Land_use <- "Farmland"

#urban
hist(allInsects_totsample$Urban_1000)#should we log it?
hist(sqrt(allInsects_totsample$Urban_1000)) #sqrt is bettwe
lme50 <- lmer(distdif ~ (Urban_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(distdif ~ (Urban_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(distdif ~ (Urban_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(distdif ~ (Urban_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outUrban <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUrban <- as.data.frame(outUrban)
outUrban$Buffer <- c(50,250,500,1000)
outUrban$Land_use <- "Urban"

#Open.uncultivated
hist(allInsects_totsample$Open.uncultivated.land_250)#log??
lme50 <- lmer(distdif ~ (Open.uncultivated.land_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.16))
lme250 <- lmer(distdif ~ (Open.uncultivated.land_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.16))
lme500 <- lmer(distdif ~ (Open.uncultivated.land_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.16))
lme1000 <- lmer(distdif ~ (Open.uncultivated.land_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Open.uncultivated.land_50 < 0.16))
outOpen.uncultivated <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outOpen.uncultivated <- as.data.frame(outOpen.uncultivated)
outOpen.uncultivated$Buffer <- c(50,250,500,1000)
outOpen.uncultivated$Land_use <- "Grassland"

#Heathland
hist(allInsects_totsample$Heathland_250)#log??
lme50 <- lmer(distdif ~ (Heathland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Heathland_50 < 0.16))
lme250 <- lmer(distdif ~ (Heathland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Heathland_50 < 0.16))
lme500 <- lmer(distdif ~ (Heathland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Heathland_50 < 0.16))
lme1000 <- lmer(distdif ~ (Heathland_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=subset(allInsects_totsample, Heathland_50 < 0.16))
outHeathland <- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outHeathland <- as.data.frame(outHeathland)
outHeathland$Buffer <- c(50,250,500,1000)
outHeathland$Land_use <- "Heathland"

#forest
hist(allInsects_totsample$Forest_250)#log??
lme50 <- lmer(distdif ~ (Forest_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay + 
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(distdif ~ (Forest_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(distdif ~ (Forest_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(distdif ~ (Forest_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outForest<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outForest <- as.data.frame(outForest)
outForest$Buffer <- c(50,250,500,1000)
outForest$Land_use <- "Forest"

#wetland
hist(allInsects_totsample$Wetland_1000)#log??
lme50 <- lmer(distdif ~ (Wetland_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(distdif ~ (Wetland_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(distdif ~ (Wetland_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(distdif ~ (Wetland_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outWetland<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outWetland <- as.data.frame(outWetland)
outWetland$Buffer <- c(50,250,500,1000)
outWetland$Land_use <- "Wetland"

#unspecified
hist(allInsects_totsample$Unspecified.land.cover_1000)#log??
hist(log(allInsects_totsample$Unspecified.land.cover_1000))
lme50 <- lmer(log(Biomass+1) ~ log(Unspecified.land.cover_50) + Time_band + 
                Time_band:cnumberTime + cStops + cyDay +  
                (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme250 <- lmer(log(Biomass+1) ~ log(Unspecified.land.cover_250) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme500 <- lmer(log(Biomass+1) ~ log(Unspecified.land.cover_500) + Time_band + 
                 Time_band:cnumberTime + cStops + cyDay +  
                 (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
lme1000 <- lmer(log(Biomass+1) ~ log(Unspecified.land.cover_1000) + Time_band + 
                  Time_band:cnumberTime + cStops + cyDay +  
                  (1|RouteID_JB) + (1|PilotID), data=allInsects_totsample)
outUnspecified<- rbind(getEffect(lme50),getEffect(lme250),getEffect(lme500),getEffect(lme1000))
outUnspecified <- as.data.frame(outUnspecified)
outUnspecified$Buffer <- c(50,250,500,1000)
outUnspecified$Land_use <- "Unspecified"

#combine all effects
outAll <- rbind(outForest,outAgri,outUrban,outWetland,outOpen.uncultivated,outHeathland)
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
  xlab("Buffer size (m)") + ylab("Effect of land cover on community composition")+ theme(plot.subtitle=element_text(size=18, face="bold", color="black")) + labs(subtitle = "A")

ggsave("plots/DK_Landcover_buffer_community_composition.png",width=5,height=12)

### max land cover (1000 m) ###########
mini_data <- allInsects_totsample %>% select(PCRID, maxLand_use, maxareaProp, distdif)
table(mini_data$maxLand_use)
variables <- c("Agriculture_1000", "Forest_1000", "Urban_1000")
mini_data <- mini_data %>% dplyr::filter(maxLand_use %in% variables) 
mini_data <- mini_data %>% dplyr::mutate(maxLand_use=dplyr::recode(maxLand_use, 
                                                "Agriculture_1000"="Farmland",
                                                "Forest_1000"="Forest",
                                                "Urban_1000"="Urban"))

hist(mini_data$maxareaProp)
qqnorm(mini_data$maxareaProp)
hist(mini_data$distdif)
qqnorm(mini_data$distdif)

# plot without regression line
u <- mini_data %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Farmland",
    "Forest",
    "Urban",
  )
) %>% ggplot() + geom_point(aes(x=maxareaProp,y=distdif, colour = maxLand_use), show.legend = F, size = 3) + facet_wrap(. ~ maxLand_use, nrow = 3) +
  theme_pubclean() + xlab("Land cover at 1000 m") + ylab("Variation in community composition") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(
    limits = c(0, 1),
    labels = function(x)
      paste0(x * 100, "%")) + scale_colour_manual(values = landuseCols[c(2, 6, 1)])

# plot with regression line
upubr <- mini_data %>% mutate(
  maxLand_use = fct_relevel(
    maxLand_use,
    "Farmland",
    "Forest",
    "Urban",
  )
) %>% ggscatter(x = "maxareaProp", y = "distdif",
                  color = "maxLand_use",
                palette = landuseCols[c(2,6,1)],
                   shape = 19, size = 3, # Points color, shape and size
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "Darkgrey", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
                   cor.coeff.args = list(method = "spearman")
  ) + stat_cor(aes(label = paste(..rr.label..,
                                 if_else(readr::parse_number(..p.label..) < 0.001, 
                                         "p<0.001", ..p.label..), sep = "~`,   `~")), label.x = 0.5, label.y = 0.5) + xlab("Land cover at 1000 m") + ylab("Variation in community composition") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%"))  + facet_wrap(. ~ maxLand_use, nrow = 3) + theme_pubclean() + theme(legend.position = "none")

ggsave("plots/DK_LandcoverMax1000_community_composition_correlation.png",width=5,height=12)
