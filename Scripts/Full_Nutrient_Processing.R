#### Process August nutrient data and save as csv
#### Created by Danielle Barnas
#### Created on 8/25/2022

## Load Libraries
library(tidyverse)
library(lubridate)
library(here)
library(curl) # pull data from url


## Read in data
AugChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv')) %>% mutate(Season = "Dry")
turb1<-read_csv(here("Data","Biogeochem","August2021","Turb_NC.csv"))
MarchChemData <- read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data//March2022/CarbonateChemistry/pHProbe_Data_calculated_POcorrect.csv')) %>% mutate(Season = "Wet")
turb2<-read_csv(here("Data","Biogeochem","March2022","Turb_NC.csv"))
depth <- read_csv(here("Data","Adj_Sandwich_Depth.csv"))


## Clean Data
## remove unnecessary/redundant turbinaria data to not include in cluster analysis
turb1 <- turb1 %>%
  select(CowTagID, del15N, C_N, N_percent) %>%
  mutate(Season = "Dry")
turb2 <- turb2 %>%
  select(CowTagID, del15N, C_N, N_percent) %>%
  mutate(Season = "Wet")
turb <- full_join(turb1, turb2) %>% distinct()

## Join august and march data


## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.
# Remove outliers / irrelevant data points
removeSite1 <- AugChemData %>%
  filter(Season == "Dry",
         CowTagID == "V2",
         Tide == " Low",
         Day_Night == "Day",
         Date == ymd("2021-08-08"))

removeSite2 <- AugChemData %>%
  filter(Season == "Dry",
         CowTagID == "C4",
         Tide =="Low",
         Day_Night == "Night",
         Date == ymd("2021-08-09"))

## Filter out redundant low tide day sample, first 'low tide' was super high
removeSite3 <- AugChemData %>%
  filter(Season == "Dry",
         Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

removeSite5 <- AugChemData %>%
  filter(CowTagID %in% c("VSPRING", "Varari_Well", "CSPRING"))

## Remove outlier (samples run on different day from rest) from March
removeSite4 <- MarchChemData %>%
  filter(Season == "Wet",
         CowTagID %in% c("VSPRING", # sample locations unnecessary for my analysis
                         "V17",
                         "VRC",
                         "CRC",
                         "CPIT",
                         "CPIT_Bottom",
                         "CSPRING_ROAD",
                         "CSPRING_BEACH",
                         "CSPRING_BEACH2"))

## Pull out Location, lat, lon of CowTagIDs
gps <- AugChemData %>% select(Location, CowTagID, lat, lon) %>% distinct()

## Remove unnecessary/redundant data
AugChem <- AugChemData %>%
  select(-c(Date,
            Time,
            DateTime,
            Plate_Seep,
            Top_Plate_ID,
            Bottom_Plate_ID,
            Jamie_Plate_ID,
            Temperature,
            # Craig recommends to remove "because we don't hypothesize them to be
            # orthogonal to any of the other fDOM we're using"
            MarineHumic_Like, Lignin_Like)) %>%
  anti_join(removeSite1) %>% # remove outlier/irrelevant data
  anti_join(removeSite2) %>%
  anti_join(removeSite3) %>%
  anti_join(removeSite5) %>%
  left_join(turb, by = c('CowTagID','Season')) %>% # join with T. ornata nutrient loading data
  left_join(gps, by = c('Location','CowTagID','lat','lon'))

MarchChem <- MarchChemData %>%
  select(-c(Date,
            SeepCode,
            SamplingTime,
            TempInSitu,
            Notes)) %>%
  anti_join(removeSite4) %>%
  left_join(turb, by = c('CowTagID','Season')) %>%
  left_join(gps, by = c('CowTagID'))

ReducedChemData <- full_join(AugChem, MarchChem) %>%
  relocate(Season, .after = Day_Night)


##### Summarise: Max and Min of parameters across seasons ####
maxmin_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # get parameter max and min across sampling periods by site and plate
  summarise_at(vars(Salinity:N_percent), .funs = range, na.rm=T) %>% # returns two rows containing max and min value per CowTagID
  #summarise_at(vars(Salinity:Tyrosine_Like), .funs = diff, na.rm=T) %>%  # returns difference between rows per CowTagID (yields range at each location)
  ungroup()

max_data <- maxmin_data %>%
  group_by(Location, CowTagID, Season) %>%
  summarise_at(vars(Salinity:N_percent), .funs = max, na.rm = T) %>%  # select for max values
  mutate(MaxMin = "Maximum")

min_data <- maxmin_data %>%
  group_by(Location, CowTagID, Season) %>%
  summarise_at(vars(Salinity:N_percent), .funs = min, na.rm = T) %>%  # select for max values
  mutate(MaxMin = "Minimum")
# join max and min values and other data sets
maxmin_data <- full_join(max_data, min_data)

## Join with depth
maxmin_data <- depth %>%
  select(Location, CowTagID, adj_CT_depth_cm) %>%
  right_join(maxmin_data)
## Join with GPS
maxmin_data <- gps %>%
  right_join(maxmin_data)
## move Max Min notation to front of df
maxmin_data <- maxmin_data %>%
  relocate(MaxMin, .before = adj_CT_depth_cm)

## Write csv ####
write_csv(maxmin_data, here("Data","Biogeochem","Nutrient_Processed_MaxMin.csv"))

##### Summarise: Coefficient of Variation of parameters - group by Seasons ####
mean_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "MeanVal")

sd_data <- ReducedChemData %>%
  group_by(Location, CowTagID, Season) %>%
  # Calculate SD values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "SDVal")

Full_data <- full_join(mean_data,sd_data) %>%
  mutate(CVVal = SDVal / MeanVal) %>% # calculate CV
  select(-c(MeanVal, SDVal)) %>% # remove intermediate columns
  pivot_wider(names_from = Parameters, values_from = CVVal) %>%  # pivot back to wide
  left_join(turb) # add t ornata data (only variation across seasons, not within)

## Join depth to Full_data
Full_data <- depth %>%
  select(Location, CowTagID, adj_CT_depth_cm) %>%
  right_join(Full_data)

## Join GPS to Full_data
Full_data <- gps %>%
  right_join(Full_data) %>%
  relocate(Season, .after = CowTagID)


## Write csv ####
write_csv(Full_data, here("Data","Biogeochem","Nutrient_Processed_CV_seasons.csv"))


##### Summarise: Coefficient of Variation of parameters ####

# The coefficient of variation (CV) is the ratio of the standard deviation to the mean and shows the extent of variability
#in relation to the mean of the population. The higher the CV, the greater the dispersion.
mean_data_all <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "MeanVal")

sd_data_all <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate SD values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "SDVal")

rangeTurb <- turb %>%
  group_by(CowTagID) %>%
  summarise_at(vars(del15N:N_percent), .funs = diff, na.rm = T) # get difference between two turb samples

Full_data_all <- full_join(mean_data_all,sd_data_all) %>%
  mutate(CVVal = SDVal / MeanVal) %>% # calculate CV
  select(-c(MeanVal, SDVal)) %>% # remove intermediate columns
  pivot_wider(names_from = Parameters, values_from = CVVal) %>%  # pivot back to wide
  left_join(rangeTurb) # add t ornata data (only variation across seasons, not within)

## Join depth to Full_data_all
Full_data_all <- depth %>%
  select(Location, CowTagID, adj_CT_depth_cm) %>%
  right_join(Full_data_all)

## Join GPS to Full_data_all
Full_data_all <- gps %>%
  right_join(Full_data_all)

## Write csv ####
write_csv(Full_data_all, here("Data","Biogeochem","Nutrient_Processed_CV.csv"))


