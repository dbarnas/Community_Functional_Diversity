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


### Functional space using PCoA (modified Teixido script)

q <- list()
All.ch.tib <- tibble(x = as.numeric(),
                     y = as.numeric(),
                     CowTagID = as.character())
All.m.sgd <- tibble(PC1 = as.numeric(),
                    PC2 = as.numeric(),
                    PC3 = as.numeric(),
                    PC4 = as.numeric(),
                    CowTagID = as.character(),
                    FE = as.character())

# needs to be df before running for loop
species_entities <- as.data.frame(column_to_rownames(species_entities, var = 'Taxa'))

cowtagOrder <- alphatag %>%
  arrange(AlphaTag) %>%
  select(CowTagID)
cowtagOrder <- cowtagOrder$CowTagID[2:20] # removes seep


for(i in cowtagOrder) {

  tag = i # use for rbinding data below

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[i,] > 0)]
  # only species present in each treatment

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]
  m.sgd <- data.matrix(m.sgd) # parse from data frame to matrix array

  mid.m.sgd <- as_tibble(m.sgd) %>%
    mutate(CowTagID = tag, FE = rownames(m.sgd))
  All.m.sgd <- All.m.sgd %>% rbind(mid.m.sgd)

  tr.sgd <- tri.mesh(m.sgd[,1],m.sgd[,2], duplicate = "remove") # duplicate: default = "error", "strip" = removes all duplicate points, "remove" = leaves one point of duplicate points

  ch.sgd <- convex.hull(tr.sgd)

  ch.tib <- cbind(ch.sgd$x, ch.sgd$y, ch.sgd$i) # parse as tibble df
  colnames(ch.tib) <- c("x", "y", "i")
  ch.tib <- as_tibble(ch.tib) %>%
    select(x,y) %>%
    mutate(CowTagID = tag)

  All.ch.tib <- All.ch.tib %>% rbind(ch.tib)
}

# add VSEEP coordinates
seep.species.sgd <- colnames(sgd.sp)[which(sgd.sp['VSEEP',] > 0)] # select present species from seep
seep.fes_cond.sgd <- species_entities[rownames(species_entities) %in% seep.species.sgd, ]
seep.m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% seep.fes_cond.sgd, ]
seep.mid.m.sgd <- as_tibble(seep.m.sgd) %>%
  mutate(CowTagID = "VSEEP", FE = rownames(seep.m.sgd))
sub.ch.tib <- seep.mid.m.sgd %>%
  select(CowTagID, x = PC1, y = PC2)

All.ch.tib <- All.ch.tib %>%
  rbind(sub.ch.tib) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
  mutate(AlphaTag = factor(AlphaTag))
All.m.sgd <- All.m.sgd %>%
  rbind(seep.mid.m.sgd) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
  mutate(AlphaTag = factor(AlphaTag))

# graph faceted polygons showing functional volume

qAll <- ggplot(data = All.ch.tib, aes(x = x, y = y)) +
  geom_polygon(aes(fill = AlphaTag, color = AlphaTag), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = AlphaTag)) +
  theme_bw() +
  theme(#legend.position = "none",
    strip.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~AlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))


### COMBINE VOLUME FIGURE FROM Presence SCRIPT WITH ABOVE TO INCLUDE POLYGON

mypalette <- rev(pnw_palette(name = "Bay", n = 20))
alphapalette <- c(paste0(mypalette[1],"70"), paste0(mypalette[2],"70"), paste0(mypalette[3],"70"),
                  paste0(mypalette[4],"70"), paste0(mypalette[4],"70"), paste0(mypalette[6],"70"),
                  paste0(mypalette[7],"70"), paste0(mypalette[8],"70"), paste0(mypalette[9],"70"),
                  paste0(mypalette[10],"70"), paste0(mypalette[11],"70"), paste0(mypalette[12],"70"),
                  paste0(mypalette[13],"70"), paste0(mypalette[14],"70"), paste0(mypalette[15],"70"),
                  paste0(mypalette[16],"70"), paste0(mypalette[17],"70"), paste0(mypalette[18],"70"),
                  paste0(mypalette[19],"70"), paste0(mypalette[20], "70"))

#load data
ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID')) # move tag names to rownames and make data.frame class


spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))

alphatag <- read_csv(here("Data","CowTag_to_AlphaTag.csv"))

################################## Data manipulation and arrangements

ab.conditions.sgd <- ab.sgd

################################# compute abundance of FEs for the three conditions

fes.sgd <- levels(as_factor(spe_fes.sgd$FE))

ab.conditions.sgd <- rownames_to_column(ab.conditions.sgd, var = "CowTagID")
ab.conditions.sgd2 <- ab.conditions.sgd %>%
  pivot_longer(names_to = "Taxa", values_to = "pCover", cols = 2:ncol(ab.conditions.sgd)) %>%
  left_join(spe_fes.sgd) %>%
  group_by(CowTagID, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup()
ab.conditions.sgd <- ab.conditions.sgd2 %>%
  pivot_wider(names_from = FE, values_from = pCover)




######################


# Figure 2. Overall distribution of FE abundance across the functional space

names(alphapalette) <- alphatag$AlphaTag[1:20]

## relative abundance in ggplot
fig5.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  arrange(AlphaTag)

fig5a <- fig5.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = AlphaTag,
                 fill = AlphaTag),
             shape = 21) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib,
               aes(x = x, y = y,
                   color = AlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig5a


ggsave(here("Output", "PaperFigures", "Fig5a_FEVol_Abund_PCoA_distance.png"), fig4b, height = 5, width = 7)

