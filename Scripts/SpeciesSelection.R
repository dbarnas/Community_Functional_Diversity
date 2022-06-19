### Summer 2022 Benthic Community Composition Surveys
### To determine species selection for community physiology

#### LIBRARIES ####
library(tidyverse)
library(here)
library(tidytext)
library(geosphere)


#### READ IN DATA ####
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Surveys", "Survey_Metadata.csv"))

# total richness of each survey location
full <- survey %>%
  select(CowTagID, Taxa, SpeciesCounts) %>%
  full_join(meta) %>%
  select(CowTagID, RelativeSGD, Taxa, SpeciesCounts) %>%
  distinct() %>%
  filter(CowTagID %in% c('V19', 'V9', 'V16', 'V13', 'V8', 'V18', 'V15', 'V11')) # remove sites that are not on extreme ends of H/L

full

richness <- full %>%
  group_by(RelativeSGD, Taxa) %>%
  summarise(sum = sum(SpeciesCounts)) %>%
  filter(Taxa != 'Sand' & Taxa != 'Bare Rock', Taxa != 'Rubble')

richness%>%
  filter(sum > 1) %>%
  ggplot(aes(y = sum,
             x = reorder_within(Taxa,sum,RelativeSGD), # reorder within facets
             fill = RelativeSGD)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ RelativeSGD, scales = 'free') + # remove NA values
  scale_x_reordered() # reorder within facets

Hal <- full %>%
  filter(Taxa %in% c('Halimeda opuntia', 'Halimeda incrassata', 'Halimeda minima', 'Halimeda distorta', 'Halimeda macroloba'))
full <- full %>%
  anti_join(Hal)
Hal <- Hal %>%
  group_by(CowTagID, RelativeSGD) %>%
  summarise(SpeciesCounts = sum(SpeciesCounts)) %>%
  mutate(Taxa = 'Halimeda')
full <- full %>%
  full_join(Hal)


AVG <- full %>%
  pivot_wider(names_from = CowTagID, values_from = SpeciesCounts) %>%
  pivot_longer(cols = c(V16:V11) ,names_to = 'CowTagID', values_to = 'SpeciesCounts') %>%
  mutate(SpeciesCounts = replace_na(SpeciesCounts, replace = 0)) %>%
  group_by(RelativeSGD, Taxa) %>%
  summarise(avg = mean(SpeciesCounts)) %>%
  filter(Taxa != 'Sand' &
           Taxa != 'Bare Rock' &
           Taxa != 'Rubble' &
           Taxa != 'Anemone')

AVG %>%
  # filter(avg > 1) %>%
  ggplot(aes(y = avg,
             x = reorder_within(Taxa,avg,RelativeSGD), # reorder within facets
             fill = RelativeSGD)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ RelativeSGD, scales = 'free') + # remove NA values
  scale_x_reordered() # reorder within facets

View(AVG)

