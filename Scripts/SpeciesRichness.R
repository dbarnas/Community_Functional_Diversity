### CALCULATE RICHNESS OF SPECIES AND FUNCTIONAL TRAITS
### Created by Danielle Barnas
### Created on 8/29/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(ggmap) # mapping
library(ggrepel) # repel labels on map
library(stats) # lm() and anova()
library(emmeans)
library(agricolae) # HSD.test()

spcomp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))

### Clean data ###
spcomp <- spcomp  %>%
  select(-c(Notes, PhotoNum))


### calculate species richness ###
richness <- spcomp %>%
  select(Location, CowTagID, Taxa) %>%
  distinct() %>%
  group_by(Location, CowTagID) %>%
  count(name = 'spRichness', Taxa) %>% # value of 1 for each distinct taxa per CowTagID
  mutate(spRichness = sum(spRichness)) %>%
  ungroup() %>%
  select(-Taxa) %>%
  distinct()

### write csv
write_csv(richness, here("Data", "Surveys", "Species_Richness.csv"))

