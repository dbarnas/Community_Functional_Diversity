### Summer 2022 Benthic Substrate Surveys

#### LIBRARIES ####
library(tidyverse)
library(here)
library(tidytext)
library(geosphere)
library(curl)


#### READ IN DATA ####
survey <- read_csv(here("Data","Surveys","Substrate_2022.csv"))
meta <- read_csv(here("Data", "Surveys", "Survey_Metadata.csv"))
locations<-read_csv(here("Data","Sandwich_Locations_Final.csv"))
AugustChem<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))
MarchChem <-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Nutrients/Nutrients_watersampling_Mar22.csv'))



### CALCULATE MEAN VALUES (SALINITY AND SILICATE)

AugChem <- AugustChem %>%
  #filter(Plate_Seep == 'Plate') %>%
  filter(Location == 'Varari') %>%
  select(CowTagID,
         Date, Time, DateTime,
         Tide, Day_Night,
         Salinity, Silicate_umolL) %>%
  #unite(Tide, Day_Night, col = Tide_DayNight, sep = "_") %>% # create unique ID for tide and time; remove columns
  group_by(CowTagID) %>%
  mutate(rangeSi = max(Silicate_umolL, na.rm = T) - min(Silicate_umolL, na.rm = T)) %>% # range of Silicate experienced across sample times/tides
  mutate(rangeSal = max(Salinity, na.rm = T) - min(Salinity, na.rm = T)) %>% # range of salinities experienced across sample times/tides
  distinct(CowTagID, rangeSi, rangeSal) %>%
  filter(CowTagID != 'Varari_Well')

MarchChem <- MarchChem %>%
  #filter(Plate_Seep == 'Plate') %>%
  #filter(Location == 'Varari') %>%
  select(CowTagID,
         Date,
         Tide, Day_Night,
         Silicate_umolL) %>%
  #unite(Tide, Day_Night, col = Tide_DayNight, sep = "_") %>% # create unique ID for tide and time; remove columns
  group_by(CowTagID) %>%
  mutate(rangeSi = max(Silicate_umolL, na.rm = T) - min(Silicate_umolL, na.rm = T)) %>% # range of Silicate experienced across sample times/tides
  distinct(CowTagID, rangeSi)

remove <- MarchChem %>% # remove sample locations not in my surveys
  filter(CowTagID %in% c('VSPRING', 'VRC', 'CRC', 'CPIT', 'CPIT_Bottom',
                         'CSPRING_ROAD', 'CSPRING_BEACH', 'CSPRING_BEACH2'))
MarchChem <- MarchChem %>%
  anti_join(remove)

### FACET ORDER CowTagID BY SILICATE or SALINITY RANGE

# associate salinity range order to Top Plate ID order
# orderSalinity <- AugChem %>%
#   ungroup() %>%
#   select(rangeSal, CowTagID) %>%
#   distinct() %>%
#   arrange(rangeSal) %>%
#   mutate(CowTagID = as_factor(as.character(CowTagID))) # as_factor creates levels based on current position

# associate silicate range order to Top Plate ID order
orderSilicate <- MarchChem %>%
  ungroup() %>%
  select(rangeSi, CowTagID) %>%
  distinct() %>%
  arrange(rangeSi) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID))) # as_factor creates levels based on current position



### Processing: percent cover and ordering by biogeochemistry

# isolate latitude and longitude to join
locations <- locations %>%
  select(CowTagID, lat, lon)

# Calculate percent cover per quad
percent <- survey %>%
  #left_join(locations) %>% # join with lat and lon data
  group_by(CowTagID) %>%
  mutate(TotalCounts = sum(LiveCoral, Rubble, DeadCoral, Sand)) %>%
  ungroup() %>%
  pivot_longer(cols = c(LiveCoral:Sand), names_to = 'Substrate', values_to = 'SubstrateCounts') %>%
  mutate(PercentSub = SubstrateCounts / TotalCounts * 100)


# sort factor levels by Salinity and Silicate
#SalLevels <- paste(levels(orderSalinity$CowTagID))
SiLevels <- paste(levels(orderSilicate$CowTagID))

# assign order to factor levels by Salinity or Silicate
#percent$CowTagID <- factor(percent$CowTagID, levels = SalLevels)
percent$CowTagID <- factor(percent$CowTagID, levels = SiLevels)
levels(percent$CowTagID) # check


### Merge Silicate range with
MarchChem <- MarchChem %>%
  left_join(survey) %>%
  select(Site, CowTagID, rangeSi) %>%
  drop_na()


### Plotting

percent %>%
  drop_na() %>%
  ggplot(aes(x = CowTagID, y = PercentSub, fill = Substrate)) +
  geom_col() +
  facet_wrap(~Substrate, scales = "fixed") # same scale on x and y axes across plots

MarchChem %>%
  filter(CowTagID != 'CSEEP' & CowTagID != 'VSEEP') %>%
  ggplot(aes(y = rangeSi, x = fct_reorder(.f = CowTagID, .x = rangeSi))) +
  geom_point() +
  theme_bw() +
  labs(x = 'CowTagID')

