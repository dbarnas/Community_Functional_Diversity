#### Process August nutrient data and save as csv
#### Created by Danielle Barnas
#### Created on 8/25/2022

## Load Libraries
library(tidyverse)
library(lubridate)
library(here)
library(curl) # pull data from url


## Read in data
AllChemData <- read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data//March2022/CarbonateChemistry/pHProbe_Data_calculated_POcorrect.csv'))
turb1<-read_csv(here("Data","Biogeochem","March2022","Turb_NC.csv"))
gps <- read_csv(here("Data","Sandwich_Locations_Final.csv"))


## Clean Data
## remove unnecessary/redundant turbinaria data to not include in cluster analysis
turb <- turb1 %>%
  select(CowTagID, del15N, N_percent)


## Filter out unnecessary/redundant data
removeSite <- AllChemData %>%
  select(CowTagID) %>%
  filter(CowTagID %in% c("VSPRING", # sample locations unnecessary for my analysis
                         "VRC",
                         "CRC",
                         "CPIT",
                         "CPIT_Bottom",
                         "CSPRING_ROAD",
                         "CSPRING_BEACH",
                         "CSPRING_BEACH2"))


ReducedChemData <- AllChemData %>%
  # remove unnecessary columns
  select(-c(TempInSitu
            #MarineHumic_Like,
            #Lignin_Like # Craig recommends to remove "because we don't hypothesize them to be orthogonal to any of the other fDOM we're using"
  )) %>%
  mutate(Location = if_else(str_detect(CowTagID, "^V") == T, "Varari", "Cabral")) # if CowTagID starts with V, label location as varari, else (if starts with C) label as Cabral


## Summarise: Range of parameters

Full_data <- ReducedChemData %>%
  group_by(Location, CowTagID) %>%
  # get parameter ranges across sampling periods by site and plate
  summarise_at(vars(Salinity:Phosphate_umolL), .funs = range, na.rm=T) %>% # returns two rows containing max and min value per CowTagID
  summarise_at(vars(Salinity:Phosphate_umolL), .funs = diff, na.rm=T) %>%  # returns difference between rows per CowTagID (yields range at each location)
  ungroup() %>%
  left_join(turb, by = 'CowTagID') %>%  # join with turb
  # only keep Plate and Seep data
  anti_join(removeSite) # ignores springs, pits, and offshore samples

## Join GPS to Full_data
Full_data <- gps %>%
  select(Location,CowTagID,lat,lon) %>%
  right_join(Full_data)

write_csv(Full_data, here("Data","Biogeochem","MarchNutrient_Processed.csv"))



