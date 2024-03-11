### Comparing functional presence and volume across each quadrat

### Created by Danielle Barnas
### Created on December 14, 2022
### Modified March 5, 2023

##### LOAD LIBRARIES #####

library(tidyverse)
library(here)
library(FD)
library(tripack) # Triangulation of Irregularly Spaced Data
library(geometry) # Mesh Generation and Surface Tessellation
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(patchwork)
library(PNWColors)
library(ggrepel)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


##### READ IN DATA #####

alphatag <- read_csv(here("Data","CowTag_to_AlphaTag.csv")) # rename CowTags with alphabetical lettering

traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
comp <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv")) # richness, % richness of community pool, and % volume of community pool

meta <- read_csv(here("Data", "Full_Metadata.csv"))
shore <- read_csv(here("Data", "Shore_distance.csv")) %>% select(-Location)
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv"))

# PCoA axes of traits
fd.coord.sgd <- read.csv(here("Data","FE_4D_coord_dmb.csv"), row.names = 1) # class data.frame
# FE with trait groups
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
# species abundances (%) wide format
myspecies <- read_csv(here("Data", "Species_Abundances_wide.csv"))
# species and functional entities
species_entities <- read_csv(here("Data", "Species_FE.csv"))



##### CLEAN AND ANALYSIZE #####

chem <- chem %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         #CowTagID != "VSEEP" &
         CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
redchem <- chem %>% select(CowTagID, Phosphate_umolL, NN_umolL)
# have one df with the CV, mean, mediam, august, march, both data, look back at min and max values for pH and Silicate


comp <- comp %>%
  filter(Location == "Varari") %>%  # only analyze varari for now
  filter(CowTagID != "V13")


# only cowtag ID's
quad.label <- chem %>%
  #filter(Location == "Varari",
  #CowTagID != "VSEEP",
  #CowTagID != "V13") %>%
  distinct(CowTagID)


mylongspecies <- myspecies %>%
  pivot_longer(cols = 2:52, names_to = "Taxa", values_to = "pCover") %>%
  left_join(meta) %>%
  select(Location, CowTagID, Taxa, pCover)


# df with unique functional entities for each row
entity <- fes_traits.sgd

#  CowTagIDs as factors, to be used as relative SGD identifiers along gradient
relative.sgd <- quad.label$CowTagID



## Percent Cover

# CowTagIDs as rownames
sgd.sp <- column_to_rownames(.data = myspecies, var = "CowTagID")
sgd.sp <- as.data.frame(sgd.sp)





## Plot convex hull (modified Teixido script)
### Bar plot

# order CowTagID's by distance from the seepage point
tagOrder <- meta %>%
  filter(Location == "Varari",
         #CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  arrange(dist_to_seep_m) %>% # set arrange factor
  select(CowTagID) %>%
  left_join(alphatag)

# set cowtag order as arrange factor order
tagOrder <- tagOrder$AlphaTag[1:20] # exclude maya's sites



############################################# plot convex hull

# color palette assignment
cols <- pnw_palette("Bay",20,type="continuous")
cols <- rev(cols) # reverse color pattern so high sgd gets red
names(cols) <- tagOrder

# add raw richness values to Fric to put value above plot bars
Fric_rich <- Fric %>%
  select(CowTagID, NbSp, NbFEs) %>%
  rename(Sp = NbSp, FE = NbFEs) %>%
  pivot_longer(cols = Sp:FE, names_to = 'Parameters', values_to = 'richness')
# change legend to show distances from seep values rather than survey location (or just remove?)


# All volumes in distinct plots
p <- Fric %>%
  left_join(alphatag) %>%
  left_join(meta) %>%
  select(Sp = NbSpP, FE = NbFEsP, Vol4D = Vol8D, AlphaTag) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  mutate(AlphaTag = factor(AlphaTag)) %>%
  mutate(Parameters = factor(Parameters, levels = c("Sp", "FE", "Vol4D"))) %>%
  # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = AlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative % of Community") +
  geom_text(aes(x = Parameters, label = round(Values,0)),
            size = 3, vjust = -0.4)
p


### Raw data figure for species and functional richness at each plot - before relative plot above (1a, 1b)
Richness <- Fric_rich %>%
  pivot_wider(names_from = 'Parameters', values_from = 'richness')

# plot same format as above but raw richness values
richness_plot <- Richness %>%
  left_join(alphatag) %>%
  mutate(CowTagID = factor(AlphaTag)) %>%
  pivot_longer(cols = 2:3, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Parameters, y = Values, fill = AlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  geom_text(aes(x = Parameters, label = Values),
            size = 3, vjust = -0.4) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Richness")
richness_plot


