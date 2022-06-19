### Summer 2022 Benthic Community Composition Surveys

#### LIBRARIES ####
library(tidyverse)
library(here)
library(tidytext)
library(geosphere)
library(curl)


#### READ IN DATA ####
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Surveys", "Survey_Metadata.csv"))
locations<-read_csv(here("Data","Sandwich_Locations_Final.csv"))
AllChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))


### CALCULATE DISTANCES TO SGD TREE

# isolate seep lat and lon at Varari
seepData <- locations %>%
  filter(Plate_Seep == 'Seep',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon)

# isolate single numeric value for lat and lon
seepLat <- as.numeric(seepData$lat[1])
seepLon <- as.numeric(seepData$lon[1])

# select distinct points for each plate location to calculate distances
distData <- locations %>%
  filter(Plate_Seep == 'Plate',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon) %>%
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


# associate distance order to Top Plate ID order
orderDistance <- distData %>%
  ungroup() %>%
  select(dist_to_seep_m, CowTagID) %>%
  distinct() %>%
  arrange(dist_to_seep_m) %>%
  # as_factor creates levels based on current position
  mutate(CowTagID = as_factor(as.character(CowTagID)))


### FACET ORDER CowTagID BY SILICATE RANGE

chem <- AllChemData %>%
  filter(Plate_Seep == 'Plate') %>%
  select(CowTagID,
         Date, Time, DateTime,
         Tide, Day_Night,
         Salinity, Silicate_umolL) %>%
  #unite(Tide, Day_Night, col = Tide_DayNight, sep = "_") %>% # create unique ID for tide and time; remove columns
  group_by(CowTagID) %>%
  mutate(rangeSi = max(Silicate_umolL, na.rm = T) - min(Silicate_umolL, na.rm = T)) %>% # range of Silicate experienced across sample times/tides
  mutate(rangeSal = max(Salinity, na.rm = T) - min(Salinity, na.rm = T)) %>% # range of salinities experienced across sample times/tides
  distinct(CowTagID, rangeSi, rangeSal)

# associate salinity range order to Top Plate ID order
orderSalinity <- chem %>%
  ungroup() %>%
  select(rangeSal, CowTagID) %>%
  distinct() %>%
  arrange(rangeSal) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID))) # as_factor creates levels based on current position

# associate salinity range order to Top Plate ID order
orderSilicate <- chem %>%
  ungroup() %>%
  select(rangeSi, CowTagID) %>%
  distinct() %>%
  arrange(rangeSi) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID))) # as_factor creates levels based on current position



### Processing

percent <- survey %>%
  select(CowTagID, Taxa, SpeciesCounts) %>%
  group_by(CowTagID) %>%
  mutate(TotalCounts = sum(SpeciesCounts)) %>%
  ungroup() %>%
  mutate(PercentTaxa = SpeciesCounts / TotalCounts * 100) %>%
  left_join(distData) %>%
  left_join(chem)


# sort factor levels by distance
distLevels <- paste(levels(orderDistance$CowTagID))
# sort factor levels by salinity
salLevels <- paste(levels(orderSalinity$CowTagID))
# sort factor levels by silicate
siLevels <- paste(levels(orderSilicate$CowTagID))

# assign order to factor levels by distance
percent$CowTagID <- factor(percent$CowTagID, levels = distLevels)
# assign order to factor levels by salinity
percent$CowTagID <- factor(percent$CowTagID, levels = salLevels)
# assign order to factor levels by silicate
percent$CowTagID <- factor(percent$CowTagID, levels = siLevels)

levels(percent$CowTagID) # check


### Plotting

NotIn <- percent %>%
  filter(Taxa %in% c('coral1', 'coral2', 'coral3', 'coral4', 'Halimeda2', # remove unidentified taxa
                     'Bare Rock', 'Rubble', 'Sand')) # remove abiotic substrates

percent %>%
  anti_join(NotIn) %>%
  ggplot(aes(x = CowTagID, y = PercentTaxa, fill = Taxa)) +
  geom_col()


