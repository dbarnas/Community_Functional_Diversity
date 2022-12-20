### Top Species in High and Low SGD from Summer 2022 Surveys###
### Created by Danielle Barnas
### Created on December 6, 2022

### LOAD LIBRARIES
library(tidyverse)
library(here)
library(ggrepel)


### READ IN DATA
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrient_Processed_CV.csv")) %>% select(-c(del15N:N_percent))

### PROCESS SPECIES
# high groundwater influence
survey <- survey %>%
  filter(Location == "Varari") %>%
  filter(Taxa != "Sand" & Taxa != "Bare Rock" & Taxa != "Rubble") %>%
  left_join(meta) %>%
  select(CowTagID, Taxa, SpeciesCounts, del15N:N_percent,dist_to_seep_m) %>%
  mutate(dist_to_seep_m = if_else(CowTagID == "V13", dist_to_seep_m*-1, dist_to_seep_m))
# view total richness across sites
plot <- survey %>%
  filter(CowTagID != "V13") %>%
  group_by(CowTagID) %>%
  count(Taxa) %>%
  summarise(total = sum(n)) %>%
  left_join(chem) %>%
  left_join(meta) %>%
  arrange(dist_to_seep_m) %>%
  ggplot(aes(x = Silicate_umolL, y = total)) +
  geom_point() +
  geom_text_repel(aes(label = CowTagID)) +
  theme_bw()
plot

# High SGD: 14, 17, 15, 11, 18
# Low SGD: 1, 2, 3, 4, 16
highsgd <- survey %>%
  filter(CowTagID == "V14" | CowTagID == "V17" | CowTagID == "V15" |
           CowTagID == "V11" | CowTagID == "V18") %>%
  group_by(Taxa) %>%
  summarise(total = sum(SpeciesCounts)) %>%
  arrange(desc(total)) %>%
  filter(Taxa != "Turf", Taxa != "Crustose Corallines", Taxa != "Heteractis magnifica") %>%
  rename(HighTaxa = Taxa)
highsgd[1:10,]


# low groundwater influence
# check other outer sites
lowsgd <- survey %>%
  filter(CowTagID == "V2" | CowTagID == "V3" |
           CowTagID == "V4" | CowTagID == "V16") %>%
  group_by(Taxa) %>%
  summarise(total = sum(SpeciesCounts)) %>%
  arrange(desc(total)) %>%
  filter(Taxa != "Turf", Taxa != "Crustose Corallines", Taxa != "Heteractis magnifica") %>%
  rename(LowTaxa = Taxa)
lowsgd[1:10,]

compareHighLow <- cbind(highsgd[1:10,], lowsgd[1:10,])





highsgd2 <- survey %>%
  filter(CowTagID == "V14" | CowTagID == "V17" | CowTagID == "V15" |
           CowTagID == "V11" | CowTagID == "V18" | CowTagID == "V8") %>%
  group_by(Taxa) %>%
  summarise(total = sum(SpeciesCounts)) %>%
  arrange(desc(total)) %>%
  filter(Taxa != "Turf", Taxa != "Crustose Corallines", Taxa != "Heteractis magnifica") %>%
  rename(HighTaxa = Taxa)
highsgd2[1:10,]


# low groundwater influence
# check other outer sites
lowsgd2 <- survey %>%
  filter(CowTagID == "V1" | CowTagID == "V2" | CowTagID == "V3" |
           CowTagID == "V4" | CowTagID == "V16" | CowTagID == "V13") %>%
  group_by(Taxa) %>%
  summarise(total = sum(SpeciesCounts)) %>%
  arrange(desc(total)) %>%
  filter(Taxa != "Turf", Taxa != "Crustose Corallines", Taxa != "Heteractis magnifica") %>%
  rename(LowTaxa = Taxa)
lowsgd2[1:10,]

compareHighLow2 <- cbind(compareHighLow, highsgd2[1:10,], lowsgd2[1:10,])





