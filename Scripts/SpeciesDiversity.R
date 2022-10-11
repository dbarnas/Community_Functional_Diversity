### CALCULATE DIVERSITY OF SPECIES AND FUNCTIONAL TRAITS
### Created by Danielle Barnas
### Created on 9/27/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(vegan) # diversity()
library(ggmap) # mapping
library(ggrepel) # repel labels on map
#library(BiodiversityR)
library(stats) # lm() and anova()
library(emmeans)
library(agricolae) # HSD.test()


#### READ IN DATA ####
survey <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
gps <- read_csv(here("Data", "Sandwich_Locations_Final.csv"))


#### PROCESS SPECIES DATA FOR DIVERSITY INDEX ####
div <- survey %>%
  filter(Location != "Cabral") %>%
  select(CowTagID, Taxa, SpeciesCounts) %>% # only relevant columns
  pivot_wider(names_from = Taxa, values_from = SpeciesCounts) %>% # go wide to associate all taxa with all survey sites
  pivot_longer(cols = 2:66, names_to = 'Taxa', values_to = 'SpeciesCounts') %>% # go back to long form for is.na test
  mutate(SpeciesCounts = if_else(is.na(SpeciesCounts), 0, SpeciesCounts)) %>% # all NA = 0
  pivot_wider(names_from = Taxa, values_from = SpeciesCounts)  # wide form for diversity index
siteID <- div[,1] # select CowTagID for later cbind
div <- div %>%
  select(-CowTagID) # only numerical data

# calculate Shannon diversity H'
index <- diversity(div, index = "shannon")
as_tibble(index) # make index a df

# add diversity values to df with CowTagIDs
index <- cbind(siteID, index)
index <- index %>%
  rename(ShannonDivSpecies = index)

# save as csv
write_csv(index, here("Data", "Surveys", "Species_Diversity.csv"))

# read in same csv
index <- read_csv(here("Data", "Surveys", "Species_Diversity.csv")) %>%
  filter(Location == "Varari") # until I get gps points for Maya's locations

#### MAP SPECIES DIVERSITY ####

# mean lat and long for the maps
MeanGPS <- gps %>%
  filter(Location == "Varari") %>%
  summarise(lon = median(lon, na.rm = TRUE) + 0.0002, # add 0.0002 to fit all on map at zoom 19
            lat = median(lat, na.rm = TRUE))
SiteGPS <- gps %>%
  filter(Location == "Varari") %>%
  group_by(CowTagID) %>%
  summarise(lon = mean(lon, na.rm = TRUE),
            lat = mean(lat, na.rm = TRUE))

# add lat and lon to index
index_loc <- left_join(index, SiteGPS) %>%
  mutate(ShannonDivSpecies = round(ShannonDivSpecies, 2))

# base map
VarariBaseMap<-get_map(MeanGPS %>%
                         #filter(Location == "Varari") %>%
                         select(lon,lat),
                       maptype = 'satellite', zoom = 19)

plotmap <- ggmap(VarariBaseMap) +
  geom_label_repel(data = index_loc,
             aes(x = lon,
                 y = lat,
                 label = paste0(CowTagID,"\n",ShannonDivSpecies)),
             size = 3,
             max.overlaps = 21) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major = element_line(color = 'white', linetype = "dashed",size = 0.5),
        plot.background=element_rect(fill='white')) +
  labs(title = "Varari Species Diversity",
       subtitle = "Shannon Diverity Index (H')")

ggsave(here("Output","V_Diversity_Map.png"), plotmap, height = 10, width = 10, device = "png")



