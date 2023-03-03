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

#### Calculate relative abundance of taxon_groups ####

relative_taxon_sub_abundance <- comp %>%
  left_join(taxa) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Taxa, Taxon_Group, SpeciesCounts) %>%
  group_by(CowTagID) %>%
  mutate(Total = sum(SpeciesCounts)) %>%
  ungroup() %>%
  group_by(CowTagID,Taxon_Group) %>%
  mutate(pCover = SpeciesCounts / Total * 100,
         pCover = sum(pCover)) %>% # add any similar taxon groups
  distinct(CowTagID, Taxon_Group, pCover) %>%  # remove duplicate taxon groups
  ungroup()

# order taxon-group
mylevels <- unique(relative_taxon_sub_abundance$Taxon_Group)[3:12] # remove hard substrate and sand
mylevels <- c(mylevels, "Hard Substrate", "Sand")
relative_taxon_sub_abundance$Taxon_Group <- factor(relative_taxon_sub_abundance$Taxon_Group, levels = mylevels)

# order cowtags
ctlevels <- dist$CowTagID
relative_taxon_sub_abundance$CowTagID <- factor(relative_taxon_sub_abundance$CowTagID, levels = ctlevels)

mypalette[11] <- "grey" # hard substrate
mypalette[12] <- "azure4" # sand

relative_taxon_sub_abundance %>%
  ggplot(aes(x = CowTagID, y = pCover, fill = Taxon_Group)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = mypalette)


#### Plot relative abundance of taxon_groups ####



relative_taxon_sub_abundance %>%
  ggplot(aes(x = CowTagID))






ferich <- comp %>%
  left_join(taxa) %>%
  filter(Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Location, CowTagID, Taxa,
         Morphology,
         #Life_Span,
         #Max_size,
         #Growth_rate,
         Zooxanthellate,
         Calcification,
         Energetic_Resource,
         Feeding_Mode) %>%
  unite(Morphology,
        #Life_Span,
        #Max_size,
        #Growth_rate,
        Zooxanthellate,
        Calcification,
        Energetic_Resource,
        Feeding_Mode,
        sep = ", ", remove = T, col = "FE") %>%
  distinct(Location, CowTagID, FE) %>%
  group_by(Location, CowTagID) %>%
  count(FE) %>%
  mutate(n = 1,
         fer = sum(n)) %>%
  distinct(Location, CowTagID, fer)


## select only distance from seep from site data
meta <- meta %>%
  select(Location, CowTagID, dist_to_seep_m, meanRugosity)

Full_data <- full_join(sprich, ferich) %>%
  left_join(meta) %>%
  filter(Location == "Varari" & CowTagID != "V13")
