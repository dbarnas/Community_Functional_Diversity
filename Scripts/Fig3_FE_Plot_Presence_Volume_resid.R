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
library(pairwiseAdonis)

##### READ IN DATA #####

# rename CowTags with alphabetical lettering
alphatag <- read_csv(here("Data","CowTag_to_AlphaTag.csv"))

traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
comp <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv")) %>%
  mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         #CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
redchem <- chem %>% select(CowTagID, Phosphate_umolL, NN_umolL)
# have one df with the CV, mean, mediam, august, march, both data, look back at min and max values for pH and Silicate

# richness, % richness of community pool, and % volume of community pool
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))

# PCoA axes of traits
fd.coord.sgd <- read.csv(here("Data","FE_4D_coord_dmb.csv"), row.names = 1) # class data.frame
# FE with trait groups
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
# species abundances (%) wide format
myspecies <- read_csv(here("Data", "Species_Abundances_wide.csv"))
# species and functional entities
species_entities <- read_csv(here("Data", "Species_FE.csv"))





##### CLEAN AND ANALYSIS #####
resFric <- Fric %>%
  left_join(meta)
resSp <- residuals(lm(data = resFric, NbSp ~ meanRugosity))
resSpp <- residuals(lm(data = resFric, NbSpP ~ meanRugosity))
resFE <- residuals(lm(data = resFric, NbFEs ~ meanRugosity))
resFEp <- residuals(lm(data = resFric, NbFEsP ~ meanRugosity))
resVol <- residuals(lm(data = resFric, Vol8D ~ meanRugosity))
resFric <- resFric %>%
  cbind(resSp) %>%
  cbind(resSpp) %>%
  cbind(resFE) %>%
  cbind(resFEp) %>%
  cbind(resVol)


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



#### RESIDUALS AND DISTANCE FROM SEEP
summary(lm(data = resFric %>% left_join(meta), # %>% filter(CowTagID != "VSEEP"),
           resSpp ~ dist_to_seep_m))
summary(lm(data = resFric %>% left_join(meta), # %>% filter(CowTagID != "VSEEP"),
           resFEp ~ dist_to_seep_m))
# summary(lm(data = resFric %>% left_join(meta), # %>% filter(CowTagID != "VSEEP"),
#            resVol ~ dist_to_seep_m))

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
Fric_rich <- resFric %>%
  select(CowTagID, resSp, resFE) %>%
  rename(Sp = resSp, FE = resFE) %>%
  pivot_longer(cols = Sp:FE, names_to = 'Parameters', values_to = 'richness')
# change legend to show distances from seep values rather than survey location (or just remove?)

# Figure 1. Species and functional diversity changes along SGD gradient
# All volumes in distinct plots
pdist <- resFric %>%
  left_join(alphatag) %>%
  left_join(meta) %>%
  select(Sp = resSpp, FE = resFEp, Vol4D = Vol8D, AlphaTag) %>%
  pivot_longer(cols = 1:2, names_to = "Parameters", values_to = "Values") %>%
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
  geom_hline(yintercept = 0) +
  #ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative % of Community (residuals)")
# geom_text(aes(x = Parameters, label = round(Values,0)),
#           size = 3, vjust = -0.4)
pdist


ggsave(here("Output", "PaperFigures", "Figure1barplot_res_distance.png"), pdist, width = 6, height = 5)


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

qAll <- All.ch.tib %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = AlphaTag, color = AlphaTag), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = AlphaTag)) +
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~AlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
qAll


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

# ab.conditions.sgd2 <- ab.conditions.sgd2 %>%
#   left_join(shore) %>%
#   mutate(rcover = pCover / shore_dist_m)



######################


# Figure 2A. Overall distribution of FE abundance across the functional space colored by NN_umolL

disc.redchem <- redchem %>%
  mutate(Phosphate_umolL = as.character(round(Phosphate_umolL,3)),
         NN_umolL = as.character(round(NN_umolL,3)))

#names(alphapalette) <- alphatag$AlphaTag[1:20]
names(alphapalette) <- (disc.redchem %>% left_join(alphatag) %>% arrange(NN_umolL))$NN_umolL

## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  arrange(AlphaTag) %>%
  left_join(disc.redchem)

fig2bdist <- fig2.fd.sgd %>%
  left_join(disc.redchem) %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 # color = AlphaTag,
                 # fill = AlphaTag),
                 color = NN_umolL,
                 fill = NN_umolL),
             shape = 21,
             show.legend = FALSE) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib %>% left_join(alphatag) %>% left_join(disc.redchem),
               aes(x = x, y = y,
                   #color = AlphaTag),
                   color = NN_umolL),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2", color = "N+N (umol L-1)") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig2bdist


ggsave(here("Output", "PaperFigures", "Plot_Vol_Abund_PCoA_res_dist_NN.png"), fig2bdist, height = 5, width = 7)



# Figure 2B. Overall distribution of FE abundance across the functional space colored by Phosphate_umolL

disc.redchem <- redchem %>%
  mutate(Phosphate_umolL = as.character(round(Phosphate_umolL,3)),
         NN_umolL = as.character(round(NN_umolL,3)))

#names(alphapalette) <- alphatag$AlphaTag[1:20]
names(alphapalette) <- (disc.redchem %>% left_join(alphatag) %>% arrange(Phosphate_umolL))$Phosphate_umolL

## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  arrange(AlphaTag) %>%
  left_join(disc.redchem)

fig2bdist <- fig2.fd.sgd %>%
  left_join(disc.redchem) %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 # color = AlphaTag,
                 # fill = AlphaTag),
                 color = Phosphate_umolL,
                 fill = Phosphate_umolL),
             shape = 21,
             show.legend = FALSE) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib %>% left_join(alphatag) %>% left_join(disc.redchem),
               aes(x = x, y = y,
                   #color = AlphaTag),
                   color = Phosphate_umolL),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2", color = "Phosphate (umol L-1)") +
  facet_wrap(~AlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig2bdist


ggsave(here("Output", "PaperFigures", "Plot_Vol_Abund_PCoA_res_dist_phos.png"), fig2bdist, height = 5, width = 7)




###########################################################################
## RESIDUALS AND PHOSPHATE
###########################################################################


# prepare chem data for joining other df for graphing
redchem <- redchem %>%
  left_join(alphatag) %>%
  arrange(Phosphate_umolL) %>%
  mutate(PAlphaTag = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T")) %>%
  unite(PAlphaTag, CowTagID, col = "PAlphaTag", sep = "-", remove=F) %>%
  arrange(NN_umolL) %>%
  mutate(NNAlphaTag = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T")) %>%
  unite(NNAlphaTag, CowTagID, col = "NNAlphaTag", sep = "-", remove=F)


summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resSp ~ poly(Phosphate_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resFEp ~ poly(Phosphate_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resVol ~ poly(Phosphate_umolL, 2)))

# color palette assignment
cols <- pnw_palette("Bay",20,type="continuous")
cols <- rev(cols) # reverse color pattern so high sgd gets red
tagorder <- redchem %>% arrange(PAlphaTag)
names(cols) <- tagorder$PAlphaTag

# Figure 1. Species and functional diversity changes along SGD gradient
# All volumes in distinct plots
pphos <- resFric %>%
  left_join(redchem) %>%
  select(Sp = resSpp, FE = resFEp, Vol4D = Vol8D, PAlphaTag) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  mutate(PAlphaTag = factor(PAlphaTag)) %>%
  mutate(Parameters = factor(Parameters, levels = c("Sp", "FE", "Vol4D"))) %>%
    # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = PAlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~PAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  #ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative % of Community (residuals)")
  # geom_text(aes(x = Parameters, label = round(Values,0)),
  #           size = 3, vjust = -0.4)
pphos

ggsave(here("Output", "PaperFigures", "Figure1barplot_res_phosphate.png"), pphos, width = 6, height = 5)



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
  left_join(redchem) %>%
  mutate(PAlphaTag = factor(PAlphaTag))
All.m.sgd <- All.m.sgd %>%
  rbind(seep.mid.m.sgd) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
  left_join(redchem) %>%
  mutate(PAlphaTag = factor(PAlphaTag))

# graph faceted polygons showing functional volume

qAll <- All.ch.tib %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = PAlphaTag, color = PAlphaTag), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = PAlphaTag)) +
  theme_bw() +
  theme(#legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~PAlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
qAll
# I want to view this on a map. Can I put polygons as map points?

#ggsave(here("Output", "PaperFigures", "Teixido_Figure1volume_dmb_CowTags.png"), qAll, width = 8, height = 5)



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
redchem <- redchem %>% arrange(PAlphaTag)
names(alphapalette) <- redchem$PAlphaTag[1:20]

## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  left_join(redchem) %>%
  arrange(PAlphaTag)

fig2bphos <- fig2.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = PAlphaTag,
                 fill = PAlphaTag),
             shape = 21) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib,
               aes(x = x, y = y,
                   color = PAlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  facet_wrap(~PAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig2bphos


ggsave(here("Output", "PaperFigures", "Plot_Vol_Abund_PCoA_res_Phosphate.png"), fig2bphos, height = 5, width = 7)




###########################################################################
## RESIDUALS AND NITRATES + NITRITES
###########################################################################


# prepare chem data for joining other df for graphing
redchem <- redchem %>%
  left_join(alphatag) %>%
  arrange(Phosphate_umolL) %>%
  mutate(PAlphaTag = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T")) %>%
  unite(PAlphaTag, CowTagID, col = "PAlphaTag", sep = "-", remove=F) %>%
  arrange(NN_umolL) %>%
  mutate(NNAlphaTag = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T")) %>%
  unite(NNAlphaTag, CowTagID, col = "NNAlphaTag", sep = "-", remove=F)


summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resSp ~ poly(NN_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resFEp ~ poly(NN_umolL, 2)))
summary(lm(data = resFric %>% left_join(redchem) %>% filter(CowTagID != "VSEEP"),
           resVol ~ poly(NN_umolL, 2)))

# color palette assignment
cols <- pnw_palette("Bay",20,type="continuous")
cols <- rev(cols) # reverse color pattern so high sgd gets red
tagorder <- redchem %>% arrange(NNAlphaTag)
names(cols) <- tagorder$NNAlphaTag

# Figure 1. Species and functional diversity changes along SGD gradient
# All volumes in distinct plots
pnn <- resFric %>%
  left_join(redchem) %>%
  select(Sp = resSpp, FE = resFEp, Vol4D = Vol8D, NNAlphaTag) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  mutate(NNAlphaTag = factor(NNAlphaTag)) %>%
  mutate(Parameters = factor(Parameters, levels = c("Sp", "FE", "Vol4D"))) %>%
  # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = NNAlphaTag)) +
  geom_col(color = "black") +
  facet_wrap(~NNAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  #ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative % of Community (residuals)")
# geom_text(aes(x = Parameters, label = round(Values,0)),
#           size = 3, vjust = -0.4)
pnn

ggsave(here("Output", "PaperFigures", "Figure1barplot_res_NN.png"), pnn, width = 6, height = 5)



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
  left_join(redchem) %>%
  mutate(NNAlphaTag = factor(NNAlphaTag))
All.m.sgd <- All.m.sgd %>%
  rbind(seep.mid.m.sgd) %>%
  left_join(alphatag) %>%
  select(-CowTagID) %>%
  left_join(redchem) %>%
  mutate(NNAlphaTag = factor(NNAlphaTag))

# graph faceted polygons showing functional volume

qAll <- All.ch.tib %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = NNAlphaTag, color = NNAlphaTag), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = NNAlphaTag)) +
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~NNAlphaTag) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
qAll
# I want to view this on a map. Can I put polygons as map points?

#ggsave(here("Output", "PaperFigures", "Teixido_Figure1volume_dmb_CowTags.png"), qAll, width = 8, height = 5)



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
redchem <- redchem %>% arrange(NNAlphaTag)
names(alphapalette) <- redchem$NNAlphaTag[1:20]

## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(alphatag) %>%
  left_join(redchem) %>%
  arrange(NNAlphaTag)

fig2bnn <- fig2.fd.sgd %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = NNAlphaTag,
                 fill = NNAlphaTag),
             shape = 21) + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib,
               aes(x = x, y = y,
                   color = NNAlphaTag),
               alpha = 0.5,
               fill = NA) + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  facet_wrap(~NNAlphaTag) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values = alphapalette) +
  scale_color_manual(values = mypalette)

fig2bnn


ggsave(here("Output", "PaperFigures", "Plot_Vol_Abund_PCoA_res_NN.png"), fig2bnn, height = 5, width = 7)

##########################################################################################
##########################################################################################

## checking things with nyssa

# bring taxonomic group up to Phyla and rename group to Phyla
# normalize to rugosity and use residuals and check patterns above again
# model everything with and without normalizing to rugosity and send nyssa a Rmd of paired relationships with response variables
# - normalized and not

# create a mega plot with all the variables for means and CV and min and max, etc. for reference





### View location of each functional entity:
fd.coord.sgd.tibble <- as_tibble(rownames_to_column(as.data.frame(fd.coord.sgd))) %>%
  rename(FE = "rowname")

FE_pca_plot <- fd.coord.sgd.tibble %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = FE),
                  size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank())
FE_pca_plot
ggsave(here("Output", "PaperFigures", "FE_pca_labeled.png"), FE_pca_plot, width = 5, height = 5)


## View functional trait faceted figures
plot_fe_group_pcoa <- fd.coord.sgd.tibble %>%
  separate(FE, into = c('Phyla','Morphology','Calcification','Energetic Resource'),
           sep = ",", remove = F) %>%
  pivot_longer(cols = 'Phyla':'Energetic Resource', names_to = "Group", values_to = "Trait")
plot_fe_group_pcoa$Group <- factor(plot_fe_group_pcoa$Group,
                                      levels = c("Phyla", "Morphology", "Calcification", "Energetic Resource"))
plot_fe_group_pcoa <- plot_fe_group_pcoa %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = Trait),
                  size = 3,
                  max.overlaps = 18) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Group)
plot_fe_group_pcoa
ggsave(here("Output", "PaperFigures", "FE_grouped_pcoa.png"), plot_fe_group_pcoa, width = 6, height = 6)


### View representative species for each functional entity
FE_representatives <- as_tibble(rownames_to_column(species_entities)) %>%
  rename(Species = "rowname") %>%
  left_join(fd.coord.sgd.tibble) %>%
  group_by(FE) %>%
  filter(row_number()==1)

FE_reps_pca_plot <- FE_representatives %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = Species),
                  size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank())
FE_reps_pca_plot
#ggsave(here("Output", "PaperFigures", "FE_pca_labeled_representatives.png"), FE_reps_pca_plot, width = 10, height = 10)


### Show all trait points possible as a blank diagram for visualization of full volume

FE_pca_plot_allPoints <- fig2.fd.sgd %>%
  filter(pCover > 0) %>%  # to show outline diagram of polygon and V2 hits every outer point
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() + # shape of a fillable circle. lets us fill with alpha values
  geom_polygon(data = All.ch.tib %>% filter(AlphaTag == "S"),
               aes(x = x, y = y),
               alpha = 0.5,
               fill = NA,
               color = "black") + # no fill on the polygon
  labs(x = "PCoA1", y = "PCoA2") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white"))

FE_pca_plot_allPoints
ggsave(here("Output", "PaperFigures","Example_Polygon.png"), FE_pca_plot_allPoints, width = 6, height = 6)

mylm <- fd.coord.sgd.tibble %>%
  separate(FE, into = c('Taxon_Group', 'Morph', 'Calc', 'ER'),
           sep = ",", remove = F)
  # unite(Taxon_Group, Morph, col = "Taxon_Morph", remove = F) %>%
  # unite(Taxon_Group, Calc, col = "Taxon_Calc", remove = F) %>%
  # unite(Morph, Calc, col = "Morph_Calc", remove = F) %>%
  #pivot_longer(cols = c('Taxon_Morph', 'Taxon_Calc', 'Morph_Calc'), names_to = "Group", values_to = "Trait")
summary(lm(data = mylm,PC2 ~ Taxon_Group))
summary(lm(data = mylm,PC1 ~ ER))


## Can use the three values above (SpR, FER, Vol4D), and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness

#
#
# ### relative abundance
# FE_nmds_data <- myspecies %>%
#   filter(CowTagID != "VSEEP") %>%  # remove seep for nMDS for now
#   pivot_longer(cols = Turf:'Caulerpa racemosa', names_to = "Taxa", values_to = "pCover") %>%
#   filter(pCover > 0) %>%
#   left_join(as_tibble(rownames_to_column(species_entities, var = "Taxa"))) %>%
#   group_by(CowTagID, FE) %>% # get relative abundance of FE (pCvoer is already percent, so just add percentages of FE)
#   mutate(pCoverFE = sum(pCover)) %>%
#   distinct(CowTagID, FE, pCoverFE) %>%
#   drop_na(FE) %>%
#   pivot_wider(names_from = FE, values_from = pCoverFE) %>% # longform for the nmds and will establish absence through NAs
#   mutate_at(vars(2:ncol(.)), .funs = ~if_else(is.na(.), 0, .)) # zero for NA's 1's for presence
# # will cbind cowtags later
# # set levels as numerical order of plates
# CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')
# FE_nmds_data$CowTagID <- factor(FE_nmds_data$CowTagID, levels = CTlevels)
# # arrange by cowtag and then remove for nmds
# FE_nmds_data <- FE_nmds_data %>%
#   arrange(CowTagID) %>%
#   ungroup() %>%
#   select(-CowTagID)
#
#
# ord1 <- metaMDS(FE_nmds_data, k=2, distance='bray')
#
# # stress with k=2 dimensions. Is it < 0.3?
# ord1$stress
#
# # stress plot - want to minimize scatter
# stressplot(ord1)
#
# #param_mds <- nMDS_species(ord1) # MDS1 and MDS2 for FEs
# # get points for species
# Group <- rownames(ord1$species) # get characteristic names
# MDS1 <- c(ord1$species[,1]) # MDS1 for characteristics
# MDS2 <- c(ord1$species[,2]) # MDS2 for characteristics
# Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
#   mutate(MDS1 = as.numeric(MDS1), # as numeric
#          MDS2 = as.numeric(MDS2)) %>%
#   #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
#   select(MDS1, MDS2, Group)
#
#
# #param_mds_cat <- nMDS_points(ord1, meta, c('CowTagID', 'dist_to_seep_m')) # MDS1 and MDS2 for CowTagID
# Groupb <- as.character(CTlevels) # assign CowTagID
# MDS1b <- ord1$points[,1] # MDS1 for CowTagID
# MDS2b <- ord1$points[,2] # MDS2 for CowTagID
# Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
#   mutate(MDS1b = as.numeric(MDS1b), # as numeric
#          MDS2b = as.numeric(MDS2b)) %>%
#   rename(CowTagID = Groupb)
#
# joinDF <- meta %>%
#   select(CowTagID, dist_to_seep_m)
#
# Datab <- Datab %>%
#   left_join(joinDF)
#
#
# ## plot
# nMDSplot <- ggplot(data = Data,
#                    aes(x = MDS1,
#                        y = MDS2)) +
#   geom_point(color = "black") +
#   geom_point(data = Datab,
#              aes(x = MDS1b,
#                  y = MDS2b,
#                  color = (dist_to_seep_m)),
#              size = 3) +
#   geom_text_repel(data = Data, # site characteristics
#                   aes(x = MDS1,
#                       y = MDS2,
#                       label = Group),
#                   size = 2) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         legend.position = "right") +
#   geom_label_repel(data = Datab, # site characteristics
#                    aes(x = MDS1b,
#                        y = MDS2b,
#                        label = CowTagID),
#                    size = 5,
#                    max.overlaps = 16) + # increase from 10 because too many FEs overlapping
#   scale_color_gradient(low = "red", high = "yellow")
# nMDSplot
#
# ggsave(here("Output", "PaperFigures", "FE_nmds_plot.png"), nMDSplot, width = 10, height = 10)
#
# ### PERMANOVA
# richPermFull <- cbind(Groupb, FE_nmds_data) %>%  # bind cowTagIDs
#   rename(CowTagID = Groupb) %>%
#   left_join(joinDF) %>%
#   mutate(relDist = if_else(dist_to_seep_m <= 26, "Near", if_else(dist_to_seep_m > 100, "Far", "Mid")))
#
# # dist < 26 > 100 ***0.001 V20, V17, V14 near vs mid 0.012, near vs far 0.045
# # dist < 47 > 100 **0.003 V20, V17, V14, V9 near vs mid 0.006
# # dist < 50 > 100 insignif. V20, V17, V14, V9, V10
#
# # dist < 26 > 120 **0.002 near vs mid 0.009
# # dist < 47 > 120 **0.001 near vs mid 0.021
#
# permanovamodel<-adonis2(richPermFull[,2:25]~relDist, richPermFull, permutations = 999,
#                         method="bray") # should change out cowtagid with some grouping name
# permanovamodel
#
# #If we are to trust the results of the permanova, then we have to assume that the dispersion among
# #data is the same in each group. We can test with assumption with a PermDisp test:
# disper<-vegdist(richPermFull[,2:25])
# betadisper(disper, richPermFull$relDist)
# #Look at the Average distance to median...these numbers should be reasonably similar
# #A rule of thumb is that one number should not be twice as high as any other
#
# pairwise.adonis(richPermFull[2:25], richPermFull$relDist, perm=999)
#
# #Get coefficients to see which species are most important in explaining site differences:
# #permanovamodel$coefficients
#
#
#

#
#
# #################################################################################
# #################################################################################
# #################################################################################
#
#
#
# ### Intersecting functional space using PCoA (modified Teixido script)
#
#
# # Supplementary Figure 1. Intersection of the three functional volumes among pH zones.
# # all volumes in 1 figure,
#
#
# r <- All.ch.tib %>%
#   mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%  # if not already run above
#   ggplot(aes(x = x, y = y)) +
#   geom_polygon(aes(color = CowTagID, fill = CowTagID), alpha = 0.3) + # create polygon using product of convex.hull(tri.mesh)
#   labs(x = "PCoA 1", y = "PCoA 2", fill = "Survey Location", color = "Survey Location") +
#   geom_point(data = All.m.sgd, aes(x = PC1, y = PC2, color = CowTagID)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2])) +
#   scale_color_manual(values = cols) +
#   scale_fill_manual(values = cols)
# r
#
# ggsave(here("Output", "Teixido", "Teixido_S1_dmb_CowTags.png"), r, width = 6, height = 5)
#
#
#
# ## Null Model of Functional Richness among SGD zones
#
#
# ############### null model of Functional richness among pH zones
#
# sgd.sp <- data.matrix(sgd.sp) # replace ab.conditions - must be matrix class
#
# n_perm = 100
# min.relative.sgd <- relative.sgd[c(1:3,5,8:19)]
# Fric_perm <- lapply(min.relative.sgd, function (x) { # condition
#
#   species.sgd <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions
#
#   perm <- sapply((1:n_perm), function (z) {
#
#     species_entities$FE <- sample(species_entities$FE) # spe_fes
#
#     fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]
#
#     m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]
#
#     ch.sgd <- convhulln(m.sgd, options = "FA")
#
#     chg.sgd <- convhulln(fd.coord.sgd, options = "FA")
#
#     c(length(species.sgd), length(species.sgd)/ncol(sgd.sp)*100, dim(m.sgd)[1], dim(m.sgd)[1]/dim(fd.coord.sgd)[1]*100, ch.sgd$vol/chg.sgd$vol*100)
#
#   })#eo sapply
#
#   rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")
#
#
#   perm
#
# })#eo lapply
#
# names(Fric_perm) = min.relative.sgd # condition
#
#
#
# Fric_perm_Q <- lapply(Fric_perm, function (x) {
#
#   rowQuantiles(x, probs=c(0.05, 0.95))
#
# })#eo lapply
#
#
#
# Fric = as.data.frame(Fric)
#
# Fric$lowerFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply # condition
# Fric$upperFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
# Fric$lowerVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
# Fric$upperVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
# Fric$cond <- relative.sgd
# relative.sgd <- factor(relative.sgd, levels = tagOrder)
#
# Fric$cond <- as.factor(relative.sgd)
# levels(Fric$cond)
# colnames(Fric) <- c("NbSp", "NbSpP", "NbFE","NbFEP", "Vol8D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")
#
#
# #Plot the null model
# #Supplementary Figure 2. Null model of functional richness (functional volume) among pH zones.
#
# tiff(filename="Output/Teixido/Figure_S2_dmb_CowTags.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)
#
#
# s <- Fric %>%
#   ggplot(aes(x = cond, y = Vol8D)) +
#   geom_point(aes(color = cond), size = 3) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_color_manual(values = cols) +
#   labs(x = "Survey Locations", y = "Relative Richness (%)") +
#   ylim(0,100) # set y axis scale
# s
#
# ggsave("Output/Teixido/Figure_S2_dmb_CowTags.png", s, width = 6, height = 5)
#
#
#
#
#
#
#
#
#
#
#
#
# ################################################ Convex Hull Intersect
#
# #load intersect function to compute convex hull (vertices + volume) of two set of points and their intersection
#
#
#
# # source(here("Scripts","Teixido","intersect.R"))
# #
# #
# # mat_int <- Fric <- lapply(relative.sgd, function (x) { # condition
# #
# #   species <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions
# #
# #   fes_cond <- species_entities[rownames(species_entities) %in% species, ] # spe_fes
# #
# #   m <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond,]
# #
# #   return(m)
# #
# # })#eo lapply
# #
# # names(mat_int) = relative.sgd # condition
# #
# # ###############intersect Low with Moderate
# #
# # Low_int_Mod <- CHVintersect(mat_int[["Low"]],mat_int[["Moderate"]])
# #
# # # percentage of the Moderate volume within Low
# # Low_int_Mod$vol[3]/Low_int_Mod$vol[2]
# #
# #
# #
# #
# # ###############intersect Low with High
# #
# # Low_int_High <- CHVintersect(mat_int[["Low"]],mat_int[["High"]])
# #
# # #pergcentage of the Extreme Low volume within Ambient
# # Low_int_High$vol[3]/Low_int_High$vol[2]
# #
# #
# #
# # ###############intersect Moderate with High
# #
# # Mod_int_High <- CHVintersect(mat_int[["Moderate"]],mat_int[["High"]])
# #
# # #percentage of the High volume within Low
# # Mod_int_High$vol[3]/Mod_int_High$vol[2]
# #
#
#
#
#
#
#
#
#
#
#
#
#
# ########################################### BETA DIVERSITY
#
# library('betapart')
#
#
# ###### taxonomic (Jaccard)
# sgd.sp[which(sgd.sp>0)] = 1 # ab.conditions
# bata.taxo <- beta.pair(sgd.sp, index.family="jaccard")
#
#
# ###### functional (Jaccard like)
#
# # Compute abundances of FEs for the three conditions
#
# # Load again the spe_fes matrix, 2 column variables
#
# #spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=",")
# species_entities <- rownames_to_column(species_entities)
# colnames(species_entities) <- c("species", "FE")
# species_entities$FE <- as_factor(species_entities$FE)
#
#
#
# # tidy version of computing abundances of FE's for three conditions
# ab.fe.conditions <- sgd.sp %>%
#   as_tibble(sgd.sp) %>%
#   mutate(CowTagID = relative.sgd) %>% # effectively makes rownames a column
#   relocate(CowTagID, .before = Turf) %>%
#   pivot_longer(cols = 2:ncol(.), names_to = "species", values_to = "presence") %>%
#   left_join(species_entities) %>%
#   drop_na() %>%  # for now removes cyanobacteria
#   group_by(CowTagID, FE) %>%
#   summarise(presence = sum(presence)) %>%
#   mutate(presence = if_else(presence > 0, 1, 0), # binary for presence-absence
#          FE = as.character(FE)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = FE, values_from = presence)
#
# ab.fe.conditions <- column_to_rownames(ab.fe.conditions, var = "CowTagID")
#
#
# # true functional beta; essential: colnames(ab.fe.conditions) == rownames(fd.coord.sgd)
# beta.fun <- functional.beta.pair(data.matrix(ab.fe.conditions), fd.coord.sgd, index.family="jaccard")
#
# #######Plot categories of the 15 functional traits across the functional space
#
# #get data to plot the traits
# #]fes <- read.csv2("Data/Teixido/Data_FEs.csv", sep=",", dec=".", row.names=1)
# fes <- column_to_rownames(entity, var = "FE")
#
# #spe_fes <- read.csv2("Data/Teixido/Data_Species_FEs.csv", sep=";", dec=",", row.names=1)
# spe_fes <- column_to_rownames(species_entities, var = "species")
#
# ###### Supplementary Figure 6. Distribution of functional trait categories across the functional space
#
# tiff(filename="Output/Teixido/Figure_S6_dmb_CowTags.tif", height=20, width=30, units="cm", compression = c("lzw"), res=300, pointsize=10)
#
# ftr <- colnames(fes)
#
# par(mfrow=c(3,5))
#
# for (i in ftr) {
#
#   lab <- as.factor(sort(unique(fes[,i])))
#
#   plot(fd.coord.sgd[,1], fd.coord.sgd[,2], pch=16, cex=1.2, col = as.numeric(fes[,i]), xlim = c(-0.3, 0.7),
#        main=gsub("."," ", i, fixed = T), xlab="PCoA 1", ylab="PCoA 2")
#   legend(x=-0.25, y=0.2, legend=lab, pch=rep(16, length(lab)), col=as.numeric(as.factor(lab)), bty = "n")
#
# }
#
# dev.off()
#
#
#
#
#
#
#
# ########################################### VOLUME MAPPED ALONG REEF
#
# library(ggmap)
# library(maptools)
#
# # mean lat and long for the maps
# LocationGPS <- meta %>%
#   filter(Location == "Varari",
#          CowTagID != "V13") %>%
#   group_by(Location) %>% # varari vs cabral
#   summarise(lon = median(lon, na.rm = TRUE),
#             lat = median(lat, na.rm = TRUE))
#
#
# VarariBaseMap <- get_map(LocationGPS %>% filter(Location == "Varari") %>% select(lon,lat), maptype = 'satellite', zoom = 19)
#
# VmapSites <- ggmap(VarariBaseMap) +
#   geom_point(data = Full_data %>% filter(CowTagID != "V13"),
#              aes(x = lon, y = lat),
#              color = "white",
#              size = 2) +
#   labs(x = "Longitude", y = "Latitude",  #label x and y axes
#        title = "Varari Sample Locations") +
#   geom_label(data = Full_data %>% filter(CowTagID != "V13"),
#              aes(x = lon, y = lat,
#                  label = CowTagID),
#              size = 1.5)
#
# VmapSites
#
#
# volSites <- VmapSites + ggplot(data = All.ch.tib, aes(x = x, y = y)) +
#   geom_polygon(aes(fill = CowTagID, color = CowTagID), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
#   labs(x = "PCoA 1", y = "PCoA 2") +
#   geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = CowTagID)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_fill_manual(values = cols) +
#   scale_color_manual(values = cols) +
#   facet_wrap(~CowTagID) +
#   xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
#
# #ggsave(here("Output","Teixido", "SiteMap_Volume.png"), volSites, width = 10, height = 9)
#
#
#
#
#
#
#
#
#
#
#
