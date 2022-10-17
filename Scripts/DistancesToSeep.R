

#### LIBRARIES ####
library(tidyverse)
library(here)
library(curl)
library(geosphere)


#### READ IN DATA ####
survey<-read_csv(here("Data","Surveys","Varari_CC_Survey_reduced_taxon.csv"))
traits<-read_csv(here("Data","Functional_Traits.csv"))
turb<-read_csv(here("Data","Biogeochem","Turb_NC.csv"))
AllChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))


### CALCULATE DISTANCES TO SGD TREE

# isolate seep lat and lon at Varari
seepData <- AllChemData %>%
  filter(Plate_Seep == 'Seep',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon, Plate_Seep)

# isolate single numeric value for lat and lon
seepLat <- as.numeric(seepData$lat[1])
seepLon <- as.numeric(seepData$lon[1])

# select distinct points for each plate location to calculate distances
distData <- AllChemData %>%
  filter(Plate_Seep == 'Plate',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon, Top_Plate_ID, Plate_Seep) %>%
  distinct() %>%
  mutate(lat_seep = seepLat,
         lon_seep = seepLon) %>%
  # find Haversine distance
  mutate(dist_to_seep_m = distHaversine(cbind(lon_seep, lat_seep), cbind(lon, lat))) %>%
  # group by sample Site Number
  group_by(CowTagID) %>%
  # choose only minimum distances
  slice(which.min(dist_to_seep_m)) %>%
  select(-c(lat_seep, lon_seep))

# isolate V13, which is in ambient upcurrent of SGD
V13dist <- distData %>%
  filter(CowTagID == 'V13') %>%
  mutate(dist_to_seep_m = -dist_to_seep_m) # get negative value because of opposite direction from other locations
# remove V13 from distData then rejoin with new value from above
distData <- distData %>%
  filter(CowTagID != 'V13') %>%
  rbind(V13dist)

write_csv(distData, here("Data", "Plate_Distance_to_Seep.csv"))

# associate distance order to Top Plate ID order
orderPlates <- distData %>%
  ungroup() %>%
  select(dist_to_seep_m, CowTagID) %>%
  distinct() %>%
  arrange(dist_to_seep_m) %>%
  # as_factor creates levels based on current position
  mutate(CowTagID = as_factor(as.character(CowTagID)))

# sort factor levels by distance
# distLevels <- paste(sort(as.numeric(levels(full_data$dist_to_seep_m))))
distLevels <- paste(levels(orderPlates$CowTagID))

# assign order to factor levels by distance
distData$CowTagID <- factor(distData$CowTagID, levels = distLevels)
levels(distData$CowTagID) # check
