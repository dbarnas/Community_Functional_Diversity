### Species Accumulation Curve ###
### Created by Danielle Barnas
### Created on November 13, 2021

library(tidyverse)
library(here)
rm(list = ls())

### READ IN DATA
quad.data<-read_csv(here("Data","Varari_CC_Survey.csv"))

### ANALYSIS

## Total quads by species/taxa counts
## species accumulation curve shows how many new species you get by doing more surveys

# Return list of distinct species by quad survey
quad.data<-quad.data %>% 
  pivot_longer(cols = DS.1:DS.25, names_to = "quad.unit", values_to = "dom.species") %>% 
  filter(dom.species != "Sand/Bare Rock") %>% 
  unite(Top_Plate_ID, quad, sep = "_", col = "quadID") %>% 
  group_by(quadID) %>% 
  distinct(dom.species) %>% 
  ungroup()

sp.data<-quad.data %>% 
  select(dom.species) %>% 
  distinct()

# remove any duplicate species from dataframe to only display one of each across full survey
quad.data.order <- quad.data[order(quad.data[,'quadID'],quad.data[,'dom.species']),]
quad.data.order <- quad.data.order[!duplicated(quad.data.order$dom.species),]

# rejoin with df containing all quads
quad.data<-quad.data %>% 
  select(quadID) %>% 
  left_join(quad.data.order) %>% 
  distinct()

# create column "num" containing a value of 1 next to each species
quad.data<-quad.data %>% 
  group_by(quadID) %>% 
  count(dom.species) %>% 
  rename(num = n) %>% 
  ungroup()

# replace values of 1 with values of 0 for all NA species (indicating no new species for a later survey)
quad.data<-quad.data %>% 
  mutate(num = if_else(is.na(dom.species) == TRUE, true = (num=0), false = (num=1)))

# add total species per quadrat survey
quad.data<-quad.data %>% 
  group_by(quadID) %>% 
  summarise(sp.sum = sum(num))

# maximum richness in a quadrat survey
quad.data %>% 
  summarise(max(sp.sum))

# add column to enumerate each quadrat survey
quad.data$row.num=seq.int(nrow(quad.data)) # add column of row numbers aka quad number

# sequentially add species richness as number of surveys increases
quad.data<-quad.data %>% 
  mutate(dom.species.sum = cumsum(sp.sum))

quad.data %>% 
  ggplot(aes(y = dom.species.sum, x = row.num)) +
  geom_point() +
  labs(x = "Total Quadrat Surveys",
       y = "Species Richness") +
  theme_bw() +
  geom_smooth()

  

