# metadata preparation

### load libraries ######################
library(tidyverse)
library(data.table)
library(plyr)
library(lubridate)
library(reshape2)

### lab data ############################
#sampleids <- read.table("data/labdata/SampleID_PCRID_sizes.txt", sep = "\t", header = T) # need to reload data - this is the correct PCRIDs and they didn't correpsond to PCRIDSamplesize before, check this and re-run
pcrids <- read.table("data/labdata/SampleID_size_SampleID_PCRID.txt", sep = "\t", header = T) # this is now the correct PCRIDS 
biomass <- read.table("data/labdata/SampleID_Biomass_sizes.txt", sep = "\t", header = T)
setdiff(pcrids$SampleID, biomass$SampleID) # the difference is extraction blanks and test samples
setdiff(biomass$SampleID, pcrids$SampleID) # unprocessed samples (may need updating)
labdata <- merge(pcrids, biomass, by = c("SampleID", "PID")) 

buffer <- read.table("data/labdata/PID_SampleID_size_insectvolume_buffer.txt", sep = "\t", header = T)
image <- read.delim("data/labdata/imagerecognition.txt")

labdata <- merge(labdata, buffer, by = c("SampleID_size", "PID"))
str(image)
str(labdata)

test <- merge(labdata, image, by.x = c("PCRID3", "SampleID", "SampleID_size", "InsectVolume", "DNABufferVolumeAdded"), by.y = c("PCRID", "SampleID", "SampleID_size", "InsectVolume", "DNABufferVolumeAdded"), all.x = T)

labdata <- test %>% select(PCRID3, PID, SampleID, SampleID_size, InsectVolume, DNABufferVolumeAdded, Qubit.Final.QIASymphony, DryMassSmallNumber..mg., DryMassLargeNumber..mg., SampleBiomassNumber..mg.)

# trying to merge data with duplicate sampleids due to size sorting
str(labdata)

# change column headers 
labdata <-
  labdata %>% dplyr::rename(
    Biomass_small = DryMassSmallNumber..mg.,
    Biomass_large = DryMassLargeNumber..mg.,
    Biomass = SampleBiomassNumber..mg.,
    Qubit = Qubit.Final.QIASymphony,
    PCRID = PCRID3
  )

# change too low to zero and too high readings to 100
test <- labdata %>% 
  mutate(Qubit = ifelse(as.character(Qubit) == "too low", "0", as.character(Qubit))) 

labdata <- test %>% 
  mutate(Qubit = ifelse(as.character(Qubit) == "Too low", "0", as.character(Qubit)))

labdata <- labdata %>% 
  mutate(Qubit = ifelse(as.character(Qubit) == "too high", "100", as.character(Qubit)))
  
### sampling event data ############################
SamplingEvent <- read.csv("data/samplingevent/SamplingEvent.csv", sep=";")
routeID <- read.delim("data/samplingevent/DK_routeID.txt")
metadata <- merge(SamplingEvent, routeID, by = "SampleID")
setdiff(labdata$SampleID, metadata$SampleID) # as far as I remember, those samples were not collected, but they at least don't have any data in the master sheet either, so we can ignore them
metadata <- merge(metadata, labdata, by = c("SampleID", "PID"))
str(metadata) 
metadata <-
  metadata %>% select(
    PCRID,
    SampleID,
    SampleID_size,
    PID,
    TripID,
    RouteID,
    Biomass_small,
    Biomass_large,
    Biomass,
    LandUSeType,
    subLandUseType,
    Date,
    StartTime,
    EndTime,
    Wind,
    Temperature,
    FullySampled
  ) # only select relevant columns for analysis

table(metadata$LandUSeType)

# Rename the values from Danish to English
metadata$LandUSeType <- mapvalues(
  metadata$LandUSeType,
  from = c(
    'mark',
    'skov',
    'skov, tør',
    'tør',
    'tør, mark',
    'tør, våd',
    'urban',
    'våd',
    'våd, tør'
  ),
  to = c(
    'Farmland',
    'Forest',
    'Forest_dry',
    'Dryland',
    'Dry_agriculture',
    'Dry_wet',
    'Urban',
    'Wetland',
    'Dry_wet'
  )
)

table(metadata$FullySampled) # all are fully sampled

# get summaries of how many samples there is for each variable and their levels
length(unique(metadata[["RouteID"]])) # count how many routes were sampled - but notice some have received new routes
length(unique(metadata[["PID"]])) # count how many pilots that carried out the sampling
length(unique(metadata[["SampleID"]])) # count how many samples
data.frame(table(metadata$Wind)) # how often were the different wind categories registered
data.frame(table(metadata$Temperature)) # how many samples were collected at different temperature intervals
data.frame(table(metadata$Date)) # how many samples per day

# change column header to macth DE variable names and drop google maps route links
data <-
  metadata %>% dplyr::rename(
    roughLand_use = LandUSeType,
    PilotID = PID
  )

head(data)
table(data$roughLand_use)

# make a column for whether sampling was midday or evening
time <- as.POSIXct(strptime(data$StartTime, "%H:%M"), "UTC")
x = as.POSIXct(strptime(c("120000", "150000", "170000", "200000"), "%H%M%S"), "UTC")
data$Time_band <-
  case_when(between(time, x[1], x[2]) ~ "midday", between(time, x[3], x[4]) ~
              "evening")

# adding a column for route length
data$Route_length <- '5000'

# Add Distance_driven which is 10000 for DK data
data$Distance_driven <- '10000'

# Change the wind types to Light, Gentle, Moderate, so without breeze and numbers
data <- data %>% mutate(Wind=recode(Wind, 
                                    "Light Breeze 1.6-3.3"="Light breeze",
                                    "Gentle breeze 3.4-5.5"="Gentle breeze",
                                    "Moderate breeze 5.5-7.9"="Moderate breeze"))

#format Date
data$Date <- as.Date(data$Date, "%d-%m-%Y")
str(data)
data$yDay <- yday(data$Date)

# Add time driven column Time_driven
# convert the time columns to datetimes
test <- data
str(test)
test$StartTime <- as.POSIXct(data$StartTime,
                             format='%H:%M:%S')
test$EndTime   <- as.POSIXct(data$EndTime,
                             format='%H:%M:%S')

data$Time_driven <- difftime(test$EndTime, test$StartTime, units = "mins") 

# setting route_length and distance_driven as numeric values
data$Route_length <- as.double(data$Route_length)
data$Distance_driven <- as.double(data$Distance_driven)
data$Time_driven <- as.double(data$Time_driven)

# Add Velocity (Route_length*2)/Time_driven - we think it could account for some of the variation between samples, especially urban (many stops)
str(data)
data$Velocity <- (data$Route_length*2)/data$Time_driven

write.table(data, file = "cleaned-data/DK_rough_landuse_biomass.txt")
