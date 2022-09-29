#### Process August nutrient data and save as csv
#### Created by Danielle Barnas
#### Created on 8/25/2022

## Load Libraries
library(tidyverse)
library(lubridate)
library(here)
library(curl) # pull data from url


## Read in data
AllChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))
turb1<-read_csv(here("Data","Biogeochem","August2021","Turb_NC.csv"))
gps <- read_csv(here("Data","Sandwich_Locations_Final.csv"))
depth <- read_csv(here("Data","Adj_Sandwich_Depth.csv"))


## Clean Data
## remove unnecessary/redundant turbinaria data to not include in cluster analysis
turb <- turb1 %>%
  select(CowTagID, del15N, C_N, N_percent)


## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.  Remove this point
removeSite1 <- AllChemData %>%
  filter(CowTagID == "V2",
         Tide == " Low",
         Day_Night == "Day",
         Date == ymd("2021-08-08"))

removeSite2 <- AllChemData %>%
  filter(CowTagID == "C4",
         Tide =="Low",
         Day_Night == "Night",
         Date == ymd("2021-08-09"))

## Filter out redundant low tide day sample, first 'low tide' was super high
removeSite3 <- AllChemData %>%
  filter(Tide == "Low",
         Day_Night == "Day",
         Date == ymd("2021-08-06"))

## Remove unnecessary/redundant data
ReducedChemData <- AllChemData %>%
  # remove outlier data and bad sampling date
  anti_join(removeSite1) %>%
  anti_join(removeSite2) %>%
  anti_join(removeSite3) %>%
  # only keep Plate and Seep data
  filter(Plate_Seep %in% c("Plate", "Seep")) %>%  # ignores springs and well
  # remove unnecessary columns
  select(-c(Date,
            Time,
            DateTime,
            Plate_Seep,
            Top_Plate_ID,
            Bottom_Plate_ID,
            Jamie_Plate_ID,
            Temperature,
            MarineHumic_Like, Lignin_Like)) # Craig recommends to remove "because we don't hypothesize them to be orthogonal to any of the other fDOM we're using"


## Summarise: Range of parameters
# Full_data <- ReducedChemData %>%
#   group_by(Location, CowTagID) %>%
#   # get parameter ranges across sampling periods by site and plate
#   summarise_at(vars(Salinity:Tyrosine_Like), .funs = range, na.rm=T) %>% # returns two rows containing max and min value per CowTagID
#   summarise_at(vars(Salinity:Tyrosine_Like), .funs = diff, na.rm=T) %>%  # returns difference between rows per CowTagID (yields range at each location)
#   ungroup() %>%
#   left_join(turb, by = "CowTagID")

## Summarise: Coefficient of Variation of parameters
mean_data <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate mean values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = mean, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "MeanVal")

sd_data <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # Calculate SD values and pivot longer to join for CV calculation
  summarise_at(vars(Salinity:Tyrosine_Like), .funs = sd, na.rm = T) %>%
  ungroup() %>%
  pivot_longer(cols = c(Salinity:Tyrosine_Like), names_to = "Parameters", values_to = "SDVal")

Full_data <- full_join(mean_data,sd_data) %>%
  mutate(CVVal = SDVal / MeanVal * 100) %>% # calculate CV
  select(-c(MeanVal, SDVal)) %>% # remove intermediate columns
  pivot_wider(names_from = Parameters, values_from = CVVal) %>% # pivot back to wide
  left_join(turb, by = 'CowTagID') # join with turb

## Join depth to Full_data
Full_data <- depth %>%
  select(Location = Site, CowTagID, adj_CT_depth_cm) %>%
  right_join(Full_data)

## Join GPS to Full_data
Full_data <- gps %>%
  select(Location,CowTagID,lat,lon) %>%
  right_join(Full_data)

## Write csv
write_csv(Full_data, here("Data","Biogeochem","AugNutrient_Processed_CV.csv"))



