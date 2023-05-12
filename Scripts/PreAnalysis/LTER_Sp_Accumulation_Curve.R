### Species Accumulation Curve for MCR LTER Suveys ###
### Created by Danielle Barnas
### Created on November 13, 2021

library(tidyverse)
library(here)
rm(list = ls())

### READ IN DATA
data<-read_csv(here("Data","McrLter","MCR_LTER_Annual_Survey_Benthic_Cover_20201221.csv"))

### ANALYSIS

## Total quads by species/taxa counts
## species accumulation curve shows how many new species you get by doing more surveys

# Return list of distinct species by quad survey
quad.data<-data %>% 
  filter(Taxonomy_Substrate_Functional_Group != "Sand" &
           #Year == "2019" &
           Habitat == "Fringing" &
           Site == "LTER 6" &
           Percent_Cover > 0) %>% 
  unite(Transect,Quadrat, col = "Location", sep = "_") %>% 
  group_by(Location) %>% 
  distinct(Taxonomy_Substrate_Functional_Group) %>%
  ungroup()

sp.data<-quad.data %>% 
  select(Taxonomy_Substrate_Functional_Group) %>% 
  distinct() %>% 
  arrange(Taxonomy_Substrate_Functional_Group)

# remove any duplicate species from dataframe to only display one of each across full survey
quad.data.order <- quad.data[order(quad.data[,'Location'],quad.data[,'Taxonomy_Substrate_Functional_Group']),]
quad.data.order <- quad.data.order[!duplicated(quad.data.order$Taxonomy_Substrate_Functional_Group),]

# rejoin with df containing all quads
quad.data<-quad.data %>% 
  select(Location) %>% 
  left_join(quad.data.order) %>% 
  distinct()

# create column "num" containing a value of 1 next to each species
quad.data<-quad.data %>% 
  group_by(Location) %>% 
  count(Taxonomy_Substrate_Functional_Group) %>% 
  rename(num = n) %>% 
  ungroup()

# replace values of 1 with values of 0 for all NA species (indicating no new species for a later survey)
quad.data<-quad.data %>% 
  mutate(num = if_else(is.na(Taxonomy_Substrate_Functional_Group) == TRUE, true = (num=0), false = (num=1)))

# add total species per quadrat survey
quad.data<-quad.data %>% 
  group_by(Location) %>% 
  summarise(sp.sum = sum(num))

# maximum richness in a quadrat survey
quad.data %>% 
  summarise(max(sp.sum))

# add column to enumerate each quadrat survey
quad.data$row.num=seq.int(nrow(quad.data)) # add column of row numbers aka quad number

# sequentially add species richness as number of surveys increases
quad.data<-quad.data %>% 
  mutate(species.sum = cumsum(sp.sum))

quad.data %>% 
  ggplot(aes(y = species.sum, x = row.num)) +
  geom_point() +
  labs(x = "Total Quadrat Surveys",
       y = "Species Richness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_smooth()



