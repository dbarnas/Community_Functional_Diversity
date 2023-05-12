### Plotting Taxonomic groups and substrate as bar plots, relative abundance

### Created by Danielle Barnas
### Created on March 1, 2023


#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(ggrepel)
library(PNWColors)
library(vegan)
library(patchwork)


#### READ IN DATA ####
meta <- read_csv(here("Data", "Full_Metadata.csv"))
comp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
taxa <- read_csv(here("Data","Surveys","Distinct_Taxa.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrient_Processed_CV.csv"))


mypalette <- pnw_palette(name = "Bay", n = 10)

dist <- meta %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
   arrange(dist_to_seep_m) %>%
   select(CowTagID)
dist <- dist$CowTagID[1:20]

#### Calculate relative abundance of taxon_groups ####

relCover <- comp %>%
  left_join(taxa) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Taxa, Taxon_Group, SpeciesCounts) %>%
  group_by(CowTagID) %>%
  mutate(Total = sum(SpeciesCounts)) %>%
  ungroup() %>%
  group_by(CowTagID,Taxon_Group) %>%
  mutate(pCover = SpeciesCounts / Total * 100,
         pCover = sum(pCover))

nspecies <- relCover %>%
  select(CowTagID,Taxa,Taxon_Group) %>%
  distinct() %>%
  ungroup()
nspecies$CowTagID <- factor(nspecies$CowTagID, levels = dist)

totalTaxa <- nspecies %>% # count total taxa in each taxon group across community
  group_by(Taxon_Group) %>%
  distinct(Taxa) %>%
  dplyr::count(Taxon_Group) %>%
  rename(total = n)

# display percentages of species observed in taxonomic groups
nspecies %>% # count total taxa in each taxon group per cowtag
  group_by(CowTagID, Taxon_Group) %>%
  dplyr::count(Taxa) %>%
  summarise(nTaxa = sum(n)) %>%
  left_join(totalTaxa) %>%
  mutate(pTaxa = nTaxa / total * 100) %>%
  ggplot(aes(x = Taxon_Group,
             y = pTaxa,
             fill = Taxon_Group)) +
  geom_col() +
  facet_wrap(~CowTagID) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))



listspecies <- nspecies %>% # list all species present in taxon groups
  select(Taxon_Group,Taxa) %>%
  arrange(Taxon_Group) %>%
  distinct()

listcowtag <- nspecies %>% # list all species present in taxon groups by cowtag
  relocate(Taxon_Group, .after = CowTagID) %>%
  arrange(CowTagID, Taxon_Group) %>%
  distinct()


relative_taxon_sub_abundance <- relCover %>% # add any similar taxon groups
  distinct(CowTagID, Taxon_Group, pCover) %>%  # remove duplicate taxon groups
  ungroup()

# order taxon-group
mylevels <- unique(relative_taxon_sub_abundance$Taxon_Group)[3:12] # remove hard substrate and sand
mylevels <- c(mylevels, "Hard Substrate", "Sand")
relative_taxon_sub_abundance$Taxon_Group <- factor(relative_taxon_sub_abundance$Taxon_Group, levels = mylevels)

# order cowtags
ctlevels <- dist
relative_taxon_sub_abundance$CowTagID <- factor(relative_taxon_sub_abundance$CowTagID, levels = ctlevels)

mypalette[11] <- "grey" # hard substrate
mypalette[12] <- "azure4" # sand

relative_taxon_sub_abundance %>%
  ggplot(aes(x = CowTagID, y = pCover, fill = Taxon_Group)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = mypalette)

# calculate biotic vs hard substrate vs sand
biotic_sub <- relative_taxon_sub_abundance %>%
  mutate(Taxon_Group = if_else(Taxon_Group == "Sand", "Sand",
                       if_else(Taxon_Group == "Hard Substrate", "Hard Substrate", "Biotic"))) %>%
  group_by(CowTagID, Taxon_Group) %>%
  summarise(pCover = sum(pCover)) %>%
  pivot_wider(names_from = Taxon_Group, values_from = pCover) %>%
  mutate_all(.funs = ~if_else(is.na(.), 0, .))

sand_sub <- relative_taxon_sub_abundance %>%
  mutate(Taxon_Group = if_else(Taxon_Group == "Sand", "Sand", "Hard or Biotic Substrate")) %>%
  group_by(CowTagID, Taxon_Group) %>%
  summarise(pCover = sum(pCover)) %>%
  pivot_wider(names_from = Taxon_Group, values_from = pCover) %>%
  mutate_all(.funs = ~if_else(is.na(.), 0, .))

 #### Plot relative abundance of taxon_groups ####


ferich <- comp %>%
  left_join(taxa) %>%
  filter(Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Location, CowTagID, Taxa,
         Morph2,
         Life_Span,
         #Max_size,
         #Growth_rate,
         Symb,
         Calc,
         ER,
         FM) %>%
  unite(Morph,
        Life_Span,
        #Max_size,
        #Growth_rate,
        Symb,
        Calc,
        ER,
        FM,
        sep = ",", remove = T, col = "FE") %>%
  distinct(Location, CowTagID, FE) %>%
  group_by(Location, CowTagID) %>%
  dplyr::count(FE) %>%
  mutate(n = 1,
         fer = sum(n)) %>%
  distinct(Location, CowTagID, fer)


## select only distance from seep from site data
meta <- meta %>%
  select(Location, CowTagID, dist_to_seep_m, meanRugosity)

# Full_data <- full_join(sprich, ferich) %>%
#   left_join(meta) %>%
#   filter(Location == "Varari" & CowTagID != "V13")
