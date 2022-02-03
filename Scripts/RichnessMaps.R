### Map of community richness ###
### Created by: Nyssa Silbiger ####
### Edited by: Danielle Barnas ####
### Created on: 10/22/2021 ####
### Edited on: 10/26/2021 ####

## load library#####
library(here)
library(tidyverse)
library(patchwork)
library(ggmap)
library(viridis)
library(maptools)
library(kriging)
library(ggnewscale)
library(wql)
library(glue)
library(gridExtra)

### READ IN DATA
survey<-read_csv(here("Data","Varari_CC_Survey_reduced_taxon.csv"))
biogeo<-read_csv(here("Data","Allbiogeochemdata_QC.csv"))

### PERCENT COVER OF DOMINANT SUBSTRATE
sp_count<-survey %>% 
  pivot_longer(cols = DS.1:DS.25, names_to = 'quadID', values_to = 'Taxa') %>% 
  group_by(Top_Plate_ID,plate.num, Taxa) %>% 
  count(Taxa) %>% 
  rename(spCount = n) %>% 
  mutate(percentCover = spCount) # percent cover 

halfsurvey<-sp_count %>%  # two sandwiches only surveyed 0.5m2
  filter(Top_Plate_ID == 'V10' | Top_Plate_ID == 'V8') %>% 
  mutate(percentCover = percentCover * 2)

sp_count<-sp_count %>%
  filter(Top_Plate_ID != 'V10' & Top_Plate_ID !='V8') %>% 
  rbind(halfsurvey) 

## Join species percent cover data with all biogeochem data
biogeo_sp<-sp_count %>% 
  left_join(biogeo, by="Top_Plate_ID")

## Join turbinaria, nutrient, trait, carbonate chemistry, and gps data
turb<-turb %>% 
  rename(Top_Plate_ID = 'CowTagID')
carbonate<-carbonate %>% 
  rename(Top_Plate_ID = 'CowTagID')

full_count <-sp_count %>% 
  left_join(nutrient, by = 'Top_Plate_ID') %>% 
  left_join(traits, by = 'Taxa') %>% 
  left_join(gps, by = 'Top_Plate_ID') %>% 
  left_join(turb, by = 'Top_Plate_ID') %>% 
  left_join(carbonate, by = 'Top_Plate_ID')

### PLOT PERCENT COVER AS FUNCTION OF BIOGEOCHEM DATA
biogeo_sp %>% 
  filter(Tide=="Low") %>% 
  filter(Day_Night=="Day") %>% 
  ggplot(aes(x = Ammonia_umolL,
           y = percentCover))+
  geom_point()
# silicate could be a good proxy Silicate_umolL

full_count %>% 
  ggplot(aes(x=spCount,
             y))+
  geom_point()
  
  filter(Tide=="Low") %>% 
  filter(Day_Night=="Night") %>% 
  ggplot(aes(x = Phosphate,
             y = percentCover))+
  geom_point()

### PLOT RICHNESS AS A CONTINUOUS VARIABLE ALONG POLYGON (Use AugustMaps.R)


