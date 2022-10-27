### CALCULATE RICHNESS OF SPECIES AND FUNCTIONAL TRAITS
### Created by Danielle Barnas
### Created on 8/29/2022

#### INSTALL NECESSARY LIBRARIES ####
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")}
if("here" %in% rownames(installed.packages()) == FALSE){install.packages("here")}
if("ggmap" %in% rownames(installed.packages()) == FALSE){install.packages("ggmap")}
if("ggrepel" %in% rownames(installed.packages()) == FALSE){install.packages("ggrepel")}
if("stats" %in% rownames(installed.packages()) == FALSE){install.packages("stats")}

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(ggmap) # mapping
library(ggrepel) # repel labels on map
library(stats) # lm() and anova()

#### READ IN DATA ####
spcomp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))

### Clean data ###

# remove non taxonomic / substrate classifications
noTaxa <- spcomp %>%
  filter(Taxa %in% c("Bare Rock", "Sand", "Rubble"))
spcomp <- spcomp  %>%
  anti_join(noTaxa) %>%
  select(-c(Notes, PhotoNum))

# Cabral 0 species richness
cRich <- tibble('Location' = as.character("Cabral"), 'CowTagID' = as.character("CSEEP"), 'spRichness' = as.numeric(0))


### calculate species richness ###
richness <- spcomp %>%
  select(Location, CowTagID, Taxa) %>%
  distinct() %>%
  group_by(Location, CowTagID) %>%
  count(name = 'spRichness', Taxa) %>% # value of 1 for each distinct taxa per CowTagID
  mutate(spRichness = sum(spRichness)) %>%
  ungroup() %>%
  select(-Taxa) %>%
  distinct() %>%
  bind_rows(cRich)


### write csv
write_csv(richness, here("Data", "Surveys", "Species_Richness.csv"))

