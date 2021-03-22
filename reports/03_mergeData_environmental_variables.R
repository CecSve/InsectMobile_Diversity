### load libraries ##########
library(reshape2)
library(tidyverse)

### load data ###############
insectsDK <- read.table("cleaned-data/DK_rough_landuse_biomass.txt")

# load pilotripids and routeids (prepared by Jesper)
tripids <- read.delim("cleaned-data/DK_pilotTripIdToRouteID.txt")

# load centroid coordinates for each route
coords <- read.delim("covariate-data/DK_ruter2018_pkt_koordinater.txt")

# load traffic light counts per route
tfstops <- read.delim("covariate-data/DK_TrafficLightsCount.txt")

# add all stops
allstops <- read.delim("covariate-data/DK_ruter2018_countStops.txt")

# environmental data
environData <- read.delim("covariate-data/environData_DK.txt",as.is=T)

# add land use intensity for urban and agriculture
landuseUrban <- read.delim("covariate-data/urban_landuse_intensity_DK.txt",as.is=T)
landuseFarmland <- read.delim("covariate-data/farmland_landuse_intensity_DK.txt",as.is=T)
landuseWetland <- read.delim("covariate-data/wetland_landuse_intensity_DK.txt",as.is=T)

### merge data #################
# merging data by new routeIDs (RouteID_JB) so other data can be merged as well
allInsects <- insectsDK

mergedData <-
  merge(allInsects, tripids, by.x = "SampleID", by.y = "PilotTripID")
setdiff(allInsects$SampleID, tripids$PilotTripID) # all metadatasamples are included - yay!

# adding stopping effect (proxy = count of traffic lights on route)
with_tfstops <-
  merge(mergedData,
        tfstops,
        by.x = "RouteID_JB",
        by.y = "routeID",
        all = T)

# removing samples from stops that don't have biomass
mergedData <- with_tfstops %>% drop_na(SampleID)
mergedData <-
  mergedData %>% mutate(Num_trafficLights = replace(Num_trafficLights, is.na(Num_trafficLights), 0)) # recode NAs to zeros for number of traffic ligths on the routes

# add all stops
allstops <- read.delim("covariate-data/DK_ruter2018_countStops.txt")
with_allstops <-
  merge(mergedData,
        allstops,
        by.x = "RouteID_JB",
        by.y = "routeId2018",
        all = T)

# removing samples from stops that don't have biomass
mergedData <- with_allstops %>% drop_na(SampleID)

# adding utm coordinates for route centroids
mergedData <-
  merge(mergedData, coords, by.x = "RouteID_JB", by.y = "routeID")
mergedData <-
  select(mergedData,-OBJECTID) # remove objectid column since it is not needed  

# first replace commas with points for decimals
mergedData$utm_x <- sapply(mergedData$utm_x, gsub, pattern = ",", replacement= ".")
mergedData$utm_y <- sapply(mergedData$utm_y, gsub, pattern = ",", replacement= ".")
str(mergedData)

# change from character to numeric
mergedData$utm_x <- as.numeric(mergedData$utm_x)
mergedData$utm_y <- as.numeric(mergedData$utm_y)

mergedData <-
  mergedData %>% mutate(
    Land_use = recode(
      roughLand_use,
      "Forest" = "Forest",
      "Urban" = "Urban",
      "Dryland" = "Open uncultivated land",
      "Wetland" = "Wetland",
      "Farmland" = "Agriculture"
    )
  )

# rename stop columns
mergedData <- plyr::rename(mergedData, c("Num_trafficLights" = "tr_signals"))
mergedData <- plyr::rename(mergedData, c("COUNT_STOPS" = "stops"))

#add on environmental data
# merge land cover data with merged data
allInsects <- merge(mergedData,environData,by.x="RouteID_JB",by.y="routeID",all.x=T)

# add land use intensity to allInsects data
allInsects <- merge(allInsects,landuseUrban,by.x="RouteID_JB",by.y="routeID",all.x=T)
allInsects <- merge(allInsects,landuseFarmland,by.x="RouteID_JB",by.y="routeID",all.x=T)
allInsects <- merge(allInsects,landuseWetland,by.x="RouteID_JB",by.y="routeID",all.x=T)

# add two columns for the 1000 buffer for land use intensity analysis, where the land use with highest proportion is added and the corresponding areaProp is listed
allInsects <- allInsects %>% 
  rownames_to_column('id') %>%
  left_join(
    allInsects %>% 
      rownames_to_column('id') %>%
      tidyr::gather(maxLand_use, maxareaProp, Agriculture_1000:Wetland_1000) %>% 
      group_by(id) %>% 
      slice(which.max(maxareaProp)), 
    by = 'id'
  )

# joining introduced .y and .x to headers and should be removed
str(allInsects)
allInsects <- allInsects[, -grep(".y$", colnames(allInsects))]
names(allInsects) <- gsub(".x","",names(allInsects),fixed = TRUE)
allInsects <- column_to_rownames(allInsects, var = "id")

#sort time data to standard each around the time band
allInsects$numberTime <- as.numeric(hms(allInsects$StartTime))#Denmark
tail(allInsects)
table(allInsects$maxLand_use)

#centering
allInsects$cyDay <- allInsects$yDay - median(allInsects$yDay)
allInsects$cStops <- log(allInsects$stops+1) - median(log(allInsects$stops+1))
allInsects$cTL <- log(allInsects$tr_signals+1) - median(log(allInsects$tr_signals+1))

#transform to minutes
allInsects$numberTime <- allInsects$numberTime/60 

#centre time around each time band
middayMean <- median(allInsects$numberTime[allInsects$Time_band=="midday"],na.rm=T)# 795
eveningMean <- median(allInsects$numberTime[allInsects$Time_band=="evening"],na.rm=T)# 1087

allInsects$cnumberTime <- NA
allInsects$cnumberTime[allInsects$Time_band=="midday"] <- allInsects$numberTime[allInsects$Time_band=="midday"] -middayMean
allInsects$cnumberTime[allInsects$Time_band=="evening"] <- allInsects$numberTime[allInsects$Time_band=="evening"] -eveningMean

#set missing values to mean of time band
allInsects$cnumberTime[is.na(allInsects$cnumberTime)]<- 0

write.table(allInsects, "cleaned-data/allInsects.txt")
