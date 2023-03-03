### Comparing functional presence and volume across each quadrat


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



traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
comp <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv"))
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
Fric <- data.matrix(column_to_rownames(Fric, var = "CowTagID")) # convert Fric from tibble to matrix array


comp <- comp %>%
  filter(Location == "Varari", # only analyze varari for now
         CowTagID != "VSEEP",
         CowTagID != "V13")

# rename for separate analyses
myspecies <- comp


# only cowtag ID's
quad.label <- myspecies %>%
  select(CowTagID) %>%
  distinct()


# create attribute for Functional Entity of each species
traits <- traits %>%
  full_join(myspecies) %>% # apply traits to species composition df
  filter(Location == "Varari",
         CowTagID != "VSEEP", CowTagID != "V13",
         Identified == "yes",
         Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Taxa, Taxon_Group, Morph,
         #Life_Span,MS_cat,GR_cat, # still have NAs
         Zoox:FM) %>%
  unite(col = "fun_entity", Morph:FM, sep = ",", remove = F) %>%
  relocate(fun_entity, .after = FM) %>%
  distinct()

# write csv for all distinct FE combinations
fes_traits.sgd <- traits %>%
  select(FE = fun_entity, Morph:FM) %>%
  distinct()

write_csv(fes_traits.sgd, here("Data", "Distinct_FE.csv"))


# species abundance as wide format species
myspecies <- myspecies %>%
  filter(Taxa != "Bare Rock" & # only include biological data
           Taxa != "Sand" &
           Taxa != "Rubble") %>%
  select(CowTagID, Taxa, SpeciesCounts) %>%
  group_by(CowTagID) %>%
  mutate(totalCount = sum(SpeciesCounts), # calculate percent cover of taxa
         cover = SpeciesCounts / totalCount * 100) %>%
  ungroup() %>%
  select(-c(SpeciesCounts,totalCount)) %>%
  group_by(CowTagID, Taxa) %>% # make sure all distinct taxa are summed together
  mutate(cover = sum(cover)) %>%
  distinct() %>%  # include only one of each taxa category. may have some repeats after ID'ing previously unknown organisms
  pivot_wider(names_from = Taxa, values_from = cover) %>%
  mutate_all(.funs = ~if_else(is.na(.), 0, .)) %>%  # zero percent for any NA values
  ungroup()

# write csv for species abundances
write_csv(myspecies, here("Data", "Species_Abundances_wide.csv"))

# select only for species and functional entities
species_entities <- traits %>%
  select(Taxa, fun_entity)

# write csv for species entities tied to FE
write_csv(species_entities, here("Data", "Species_FE.csv"))

# df with unique functional entities for each row
entity <- traits %>%
  select(-c(Taxa, Taxon_Group)) %>%
  distinct() %>%
  relocate(fun_entity, .before = Morph) %>%
  mutate_all(.funs = as_factor) # groups need to be factors to run quality_funct_space()

#  CowTagIDs as factors, to be used as relative SGD identifiers along gradient
relative.sgd <- quad.label$CowTagID




## Percent Cover
# Set categories of relative SGD (based on clustering or PERMANOVA from biogeochemistry)
# OR Set categories of CowTagIDs (will sort later based on biogeochemistr or nutrients)
sgd.sp <- myspecies %>%
  left_join(quad.label) %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  group_by(CowTagID, Taxa) %>% # make sure cover is totaled by taxa
  mutate(cover = sum(cover)) %>%
  ungroup() %>%
  distinct() %>% # remove duplicate values after summing up percent cover
  group_by(CowTagID) %>%
  pivot_wider(names_from = Taxa, values_from = cover)  # pivot species wide by relative SGD

# relative sgd as rownames
# OR CowTagIDs as rownames
sgd.sp <- column_to_rownames(.data = sgd.sp, var = "CowTagID")
sgd.sp <- as.data.frame(sgd.sp)

# species names as rownames to parse to dataframe (if not already)
# species_entities <- column_to_rownames(.data = species_entities, var = "Taxa")
# species_entities <- as.data.frame(species_entities)



## Plot convex hull (modified Teixido script)
### Bar plot

# order CowTagID's by distance from the seepage point
tagOrder <- meta %>%
  filter(Location == "Varari",
         CowTagID != "VSeep",
         CowTagID != "V13") %>%
  arrange(desc(dist_to_seep_m)) %>% # set arrange factor
  select(CowTagID)
# set cowtag order as arrange factor order
tagOrder <- tagOrder$CowTagID[1:19] # exclude maya's sites



############################################# plot convex hull
# cols <- c("#3A5FCD", "#FFA500", "#CD2626")
# colstr <- c("#3A5FCD70", "#FFA50070", "#CD262670") # sets same color but at 70% opacity
#
# names(cols) <- c("Low", "Moderate", "High") # name the colors for the ggplot for loop
# names(colstr) <- c("Low", "Moderate", "High")

cols <- pnw_palette("Bay",19,type="continuous")
names(cols) <- tagOrder

# add raw richness values to Fric to put value above plot bars
Fric_rich <- Fric %>%
  as_tibble() %>%
  mutate(CowTagID = relative.sgd) %>%
  select(CowTagID, NbSp, NbFEs) %>%
  rename(Sp = NbSp, FE = NbFEs) %>%
  pivot_longer(cols = Sp:FE, names_to = 'Parameters', values_to = 'richness')

# change legend to show distances from seep values rather than survey location (or just remove?)


# Figure 1. Species and functional diversity changes along SGD gradient
# All volumes in distinct plots
p <- Fric %>%
  as_tibble() %>%
  mutate(CowTagID = relative.sgd) %>%
  select(Sp = NbSpP, FE = NbFEsP, Vol4D = Vol8D, CowTagID) %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  #left_join(Fric_rich) %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%
    # plot facet_wrapped
  ggplot(aes(x = Parameters, y = Values, fill = CowTagID)) +
  geom_col(color = "black") +
  facet_wrap(~CowTagID) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative richness (%)") +
  geom_text(aes(x = Parameters, label = round(Values,0)),
            size = 3, vjust = -0.4)
p

ggsave(here("Output", "PaperFigures", "Teixido_Figure1barplot_dmb_CowTags.png"), p, width = 6, height = 5)


### Raw data figure for species and functional richness at each plot - before relative plot above (1a, 1b)
Richness <- Fric_rich %>%
  pivot_wider(names_from = 'Parameters', values_from = 'richness')

# plot same format as above but raw richness values
richness_plot <- Richness %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%
  pivot_longer(cols = 2:3, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Parameters, y = Values, fill = CowTagID)) +
  geom_col(color = "black") +
  facet_wrap(~CowTagID) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  geom_text(aes(x = Parameters, label = Values),
            size = 3, vjust = -0.4) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Richness")
richness_plot

ggsave(here("Output", "PaperFigures", "Raw_richness_CowTags.png"), richness_plot, width = 6, height = 5)


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

for(i in tagOrder) {

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

All.ch.tib <- All.ch.tib %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder))
All.m.sgd <- All.m.sgd %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder))

# graph faceted polygons showing functional volume
qAll <- ggplot(data = All.ch.tib, aes(x = x, y = y)) +
  geom_polygon(aes(fill = CowTagID, color = CowTagID), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.spacing = unit(1, "lines")) + # increase facet wrap panel spacing
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~CowTagID) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
qAll
# I want to view this on a map. Can I put polygons as map points?

ggsave(here("Output", "PaperFigures", "Teixido_Figure1volume_dmb_CowTags.png"), qAll, width = 8, height = 5)


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
ggsave(here("Output", "PaperFigures", "FE_pca_labeled.png"), FE_pca_plot, width = 10, height = 10)


## View subset figures
fd.coord.sgd.tibble %>%
  separate(FE, into = c('Morphology', 'Zooxanthellate', 'Calcification', 'Energetic_Resource', 'Feeding_Mode'),
           sep = ",", remove = F) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = Morphology),
                  size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank())

### View representative species for each functional entity
FE_representatives <- as_tibble(rownames_to_column(species_entities)) %>%
  rename(Species = "rowname",
         FE = "fun_entity") %>%
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
ggsave(here("Output", "PaperFigures", "FE_pca_labeled_representatives.png"), FE_reps_pca_plot, width = 10, height = 10)


## Can use the three values above (SpR, FER, Vol4D), and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness



### relative abundance
FE_nmds_data <- myspecies %>%
  pivot_longer(cols = Turf:'Caulerpa racemosa', names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(as_tibble(rownames_to_column(species_entities, var = "Taxa"))) %>%
  rename(FE = fun_entity) %>%
  group_by(CowTagID, FE) %>% # get relative abundance of FE (pCvoer is already percent, so just add percentages of FE)
  mutate(pCoverFE = sum(pCover)) %>%
  distinct(CowTagID, FE, pCoverFE) %>%
  drop_na(FE) %>%
  pivot_wider(names_from = FE, values_from = pCoverFE) %>% # longform for the nmds and will establish absence through NAs
  mutate_at(vars(2:24), .funs = ~if_else(is.na(.), 0, .)) # zero for NA's 1's for presence
# will cbind cowtags later
# set levels as numerical order of plates
CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')
FE_nmds_data$CowTagID <- factor(FE_nmds_data$CowTagID, levels = CTlevels)
# arrange by cowtag and then remove for nmds
FE_nmds_data <- FE_nmds_data %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  select(-CowTagID)


ord1 <- metaMDS(FE_nmds_data, k=2, distance='bray')

# stress with k=2 dimensions. Is it < 0.3?
ord1$stress

# stress plot - want to minimize scatter
stressplot(ord1)

#param_mds <- nMDS_species(ord1) # MDS1 and MDS2 for FEs
# get points for species
Group <- rownames(ord1$species) # get characteristic names
MDS1 <- c(ord1$species[,1]) # MDS1 for characteristics
MDS2 <- c(ord1$species[,2]) # MDS2 for characteristics
Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
  mutate(MDS1 = as.numeric(MDS1), # as numeric
         MDS2 = as.numeric(MDS2)) %>%
  #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
  select(MDS1, MDS2, Group)


#param_mds_cat <- nMDS_points(ord1, meta, c('CowTagID', 'dist_to_seep_m')) # MDS1 and MDS2 for CowTagID
Groupb <- as.character(CTlevels) # assign CowTagID
MDS1b <- ord1$points[,1] # MDS1 for CowTagID
MDS2b <- ord1$points[,2] # MDS2 for CowTagID
Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
  mutate(MDS1b = as.numeric(MDS1b), # as numeric
         MDS2b = as.numeric(MDS2b)) %>%
  rename(CowTagID = Groupb)

joinDF <- meta %>%
  select(CowTagID, dist_to_seep_m)

Datab <- Datab %>%
  left_join(joinDF)


## plot
nMDSplot <- ggplot(data = Data,
                   aes(x = MDS1,
                       y = MDS2)) +
  geom_point(color = "black") +
  geom_point(data = Datab,
             aes(x = MDS1b,
                 y = MDS2b,
                 color = (dist_to_seep_m)),
             size = 3) +
  geom_text_repel(data = Data, # site characteristics
                  aes(x = MDS1,
                      y = MDS2,
                      label = Group),
                  size = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  geom_label_repel(data = Datab, # site characteristics
                   aes(x = MDS1b,
                       y = MDS2b,
                       label = CowTagID),
                   size = 5,
                   max.overlaps = 16) + # increase from 10 because too many FEs overlapping
  scale_color_gradient(low = "red", high = "yellow")
nMDSplot

ggsave(here("Output", "PaperFigures", "FE_nmds_plot.png"), nMDSplot, width = 10, height = 10)

### PERMANOVA
richPermFull <- cbind(Groupb, FE_nmds_data) %>%  # bind cowTagIDs
  rename(CowTagID = Groupb) %>%
  left_join(joinDF) %>%
  mutate(relDist = if_else(dist_to_seep_m <= 50, "Near", if_else(dist_to_seep_m > 100, "Far", "Mid")))

# dist < 26 > 100 ***0.001 V20, V17, V14 near vs mid 0.012, near vs far 0.045
# dist < 47 > 100 **0.003 V20, V17, V14, V9 near vs mid 0.006
# dist < 50 > 100 insignif. V20, V17, V14, V9, V10

# dist < 26 > 120 **0.002 near vs mid 0.009
# dist < 47 > 120 **0.001 near vs mid 0.021

permanovamodel<-adonis2(richPermFull[,2:24]~relDist, richPermFull, permutations = 999,
                        method="bray") # should change out cowtagid with some grouping name
permanovamodel

#If we are to trust the results of the permanova, then we have to assume that the dispersion among
#data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(richPermFull[,2:24])
betadisper(disper, richPermFull$relDist)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

pairwise.adonis(richPermFull[2:24], richPermFull$relDist, perm=999)

#Get coefficients to see which species are most important in explaining site differences:
permanovamodel$coefficients










### Intersecting functional space using PCoA (modified Teixido script)


# Supplementary Figure 1. Intersection of the three functional volumes among pH zones.
# all volumes in 1 figure,


r <- All.ch.tib %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%  # if not already run above
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(color = CowTagID, fill = CowTagID), alpha = 0.3) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2", fill = "Survey Location", color = "Survey Location") +
  geom_point(data = All.m.sgd, aes(x = PC1, y = PC2, color = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2])) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)
r

ggsave(here("Output", "Teixido", "Teixido_S1_dmb_CowTags.png"), r, width = 6, height = 5)



## Null Model of Functional Richness among SGD zones


############### null model of Functional richness among pH zones

sgd.sp <- data.matrix(sgd.sp) # replace ab.conditions - must be matrix class

n_perm = 100
min.relative.sgd <- relative.sgd[c(1:3,5,8:19)]
Fric_perm <- lapply(min.relative.sgd, function (x) { # condition

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions

  perm <- sapply((1:n_perm), function (z) {

    species_entities$fun_entity <- sample(species_entities$fun_entity) # spe_fes

    fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

    m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]

    ch.sgd <- convhulln(m.sgd, options = "FA")

    chg.sgd <- convhulln(fd.coord.sgd, options = "FA")

    c(length(species.sgd), length(species.sgd)/ncol(sgd.sp)*100, dim(m.sgd)[1], dim(m.sgd)[1]/dim(fd.coord.sgd)[1]*100, ch.sgd$vol/chg.sgd$vol*100)

  })#eo sapply

  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")


  perm

})#eo lapply

names(Fric_perm) = min.relative.sgd # condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {

  rowQuantiles(x, probs=c(0.05, 0.95))

})#eo lapply



Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply # condition
Fric$upperFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- relative.sgd
relative.sgd <- factor(relative.sgd, levels = tagOrder)

Fric$cond <- as.factor(relative.sgd)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFE","NbFEP", "Vol8D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model
#Supplementary Figure 2. Null model of functional richness (functional volume) among pH zones.

tiff(filename="Output/Teixido/Figure_S2_dmb_CowTags.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)


s <- Fric %>%
  ggplot(aes(x = cond, y = Vol8D)) +
  geom_point(aes(color = cond), size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = cols) +
  labs(x = "Survey Locations", y = "Relative Richness (%)") +
  ylim(0,100) # set y axis scale
s

ggsave("Output/Teixido/Figure_S2_dmb_CowTags.png", s, width = 6, height = 5)












################################################ Convex Hull Intersect

#load intersect function to compute convex hull (vertices + volume) of two set of points and their intersection



# source(here("Scripts","Teixido","intersect.R"))
#
#
# mat_int <- Fric <- lapply(relative.sgd, function (x) { # condition
#
#   species <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions
#
#   fes_cond <- species_entities[rownames(species_entities) %in% species, ] # spe_fes
#
#   m <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond,]
#
#   return(m)
#
# })#eo lapply
#
# names(mat_int) = relative.sgd # condition
#
# ###############intersect Low with Moderate
#
# Low_int_Mod <- CHVintersect(mat_int[["Low"]],mat_int[["Moderate"]])
#
# # percentage of the Moderate volume within Low
# Low_int_Mod$vol[3]/Low_int_Mod$vol[2]
#
#
#
#
# ###############intersect Low with High
#
# Low_int_High <- CHVintersect(mat_int[["Low"]],mat_int[["High"]])
#
# #pergcentage of the Extreme Low volume within Ambient
# Low_int_High$vol[3]/Low_int_High$vol[2]
#
#
#
# ###############intersect Moderate with High
#
# Mod_int_High <- CHVintersect(mat_int[["Moderate"]],mat_int[["High"]])
#
# #percentage of the High volume within Low
# Mod_int_High$vol[3]/Mod_int_High$vol[2]
#












########################################### BETA DIVERSITY

library('betapart')


###### taxonomic (Jaccard)
sgd.sp[which(sgd.sp>0)] = 1 # ab.conditions
bata.taxo <- beta.pair(sgd.sp, index.family="jaccard")


###### functional (Jaccard like)

# Compute abundances of FEs for the three conditions

# Load again the spe_fes matrix, 2 column variables

#spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=",")
species_entities <- rownames_to_column(species_entities)
colnames(species_entities) <- c("species", "fun_entity")
species_entities$fun_entity <- as_factor(species_entities$fun_entity)

# fes <- levels(species_entities$fun_entity) # spe_fes
# ab.fe.conditions <- lapply(relative.sgd, function (z) { # condition
#   abund.fes <-  sapply(fes, function (x) {
#     spec <- as.character(species_entities[which(species_entities$fun_entity == x),]$species) # spe_fes
#     sum(sgd.sp[z,spec]) # ab.conditions
#   })#eo sapply
#   abund.fes
# })#eo lapply
#
# names(ab.fe.conditions) = relative.sgd # condition
#
# ab.fe.conditions <- do.call(rbind, ab.fe.conditions)
# ab.fe.conditions[which(ab.fe.conditions>0)] = 1

# tidy version of computing abundances of FE's for three conditions
ab.fe.conditions <- sgd.sp %>%
  as_tibble(sgd.sp) %>%
  mutate(CowTagID = relative.sgd) %>% # effectively makes rownames a column
  relocate(CowTagID, .before = Turf) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "species", values_to = "presence") %>%
  left_join(species_entities) %>%
  drop_na() %>%  # for now removes cyanobacteria
  group_by(CowTagID, fun_entity) %>%
  summarise(presence = sum(presence)) %>%
  mutate(presence = if_else(presence > 0, 1, 0), # binary for presence-absence
         fun_entity = as.character(fun_entity)) %>%
  ungroup() %>%
  pivot_wider(names_from = fun_entity, values_from = presence)

ab.fe.conditions <- column_to_rownames(ab.fe.conditions, var = "CowTagID")


# true functional beta; essential: colnames(ab.fe.conditions) == rownames(fd.coord.sgd)
beta.fun <- functional.beta.pair(data.matrix(ab.fe.conditions), fd.coord.sgd, index.family="jaccard")

#######Plot categories of the 15 functional traits across the functional space

#get data to plot the traits
#]fes <- read.csv2("Data/Teixido/Data_FEs.csv", sep=",", dec=".", row.names=1)
fes <- column_to_rownames(entity, var = "fun_entity")

#spe_fes <- read.csv2("Data/Teixido/Data_Species_FEs.csv", sep=";", dec=",", row.names=1)
spe_fes <- column_to_rownames(species_entities, var = "species")

###### Supplementary Figure 6. Distribution of functional trait categories across the functional space

tiff(filename="Output/Teixido/Figure_S6_dmb_CowTags.tif", height=20, width=30, units="cm", compression = c("lzw"), res=300, pointsize=10)

ftr <- colnames(fes)

par(mfrow=c(3,5))

for (i in ftr) {

  lab <- as.factor(sort(unique(fes[,i])))

  plot(fd.coord.sgd[,1], fd.coord.sgd[,2], pch=16, cex=1.2, col = as.numeric(fes[,i]), xlim = c(-0.3, 0.7),
       main=gsub("."," ", i, fixed = T), xlab="PCoA 1", ylab="PCoA 2")
  legend(x=-0.25, y=0.2, legend=lab, pch=rep(16, length(lab)), col=as.numeric(as.factor(lab)), bty = "n")

}

dev.off()







########################################### VOLUME MAPPED ALONG REEF

library(ggmap)
library(maptools)

# mean lat and long for the maps
LocationGPS <- chem %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  group_by(Location) %>% # varari vs cabral
  summarise(lon = median(lon, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))

# Full_data <- Fric %>%
#   as_tibble() %>%
#   mutate(CowTagID = relative.sgd) %>%
#   relocate(CowTagID, .before = NbSp) %>%
#   left_join(chem) %>%
#   select(CowTagID:Vol8D, lat, lon)

VarariBaseMap <- get_map(LocationGPS %>% filter(Location == "Varari") %>% select(lon,lat), maptype = 'satellite', zoom = 19)

VmapSites <- ggmap(VarariBaseMap) +
  geom_point(data = Full_data %>% filter(CowTagID != "V13"),
             aes(x = lon, y = lat),
             color = "white",
             size = 2) +
  labs(x = "Longitude", y = "Latitude",  #label x and y axes
       title = "Varari Sample Locations") +
  geom_label(data = Full_data %>% filter(CowTagID != "V13"),
             aes(x = lon, y = lat,
                 label = CowTagID),
             size = 1.5)

VmapSites


volSites <- VmapSites + ggplot(data = All.ch.tib, aes(x = x, y = y)) +
  geom_polygon(aes(fill = CowTagID, color = CowTagID), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~CowTagID) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))

#ggsave(here("Output","Teixido", "SiteMap_Volume.png"), volSites, width = 10, height = 9)











