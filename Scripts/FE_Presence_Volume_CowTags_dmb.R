### Comparing functional presence and volume across each quadrat


library(tidyverse)
library(here)
library(FD)
library(tripack) # Triangulation of Irregularly Spaced Data
library(geometry) # Mesh Generation and Surface Tessellation
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(patchwork)
library(PNWColors)



traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
myspecies <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv"))

myspecies <- myspecies %>%
  filter(Location == "Varari", # only analyze varari for now
         CowTagID != "VSEEP",
         CowTagID != "V13")

quad.label <- myspecies %>%
  select(CowTagID) %>%
  distinct()

meta %>%
  filter(Location == "Varari",
         CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  select(CowTagID, N_percent, dist_to_seep_m)


# create attribute for Functional Entity of each species
traits <- traits %>%
  full_join(myspecies) %>%
  filter(Location == "Varari",
         CowTagID != "VSEEP", CowTagID != "V13",
         Identified == "yes",
         Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Taxa, Taxon_Group, Morphology,
         #Life_Span,MS_cat,GR_cat, # still have NAs
         Zooxanthellate:Feeding_Mode) %>%
  # unite(col = "taxon_entity", Taxon_Group:Feeding_Mode, sep = ",", remove = F) %>%
  # relocate(taxon_entity, .after = Feeding_Mode) %>%
  unite(col = "fun_entity", Morphology:Feeding_Mode, sep = ",", remove = F) %>%
  relocate(fun_entity, .after = Feeding_Mode) %>%
  distinct()


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



# select only for species and functional entities
species_entities <- traits %>%
  select(Taxa, fun_entity)

# df with unique functional entities for each row
entity <- traits %>%
  select(-c(Taxa, Taxon_Group)) %>%
  distinct() %>%
  relocate(fun_entity, .before = Morphology) %>%
  mutate_all(.funs = as_factor) # groups need to be factors to run quality_funct_space()

# relative SGD levels
#relative.sgd <- c("Low", "Moderate", "High")
# OR CowTagIDs as factors
relative.sgd <- quad.label$CowTagID

## Computing multidimensional functional space

# source function from Teixido: function for computing the quality of functional dendrogramm and multidimensional functional spaces
source(here("Scripts", "Teixido", "quality_funct_space.R"))

qfs <- quality_funct_space(mat_funct = entity, # distinct functional trait entities (NAs not allowed. Traits can be different types: numeric, ordinal, nominal)
                           traits_weights = NULL, # default = same weight for all traits
                           nbdim = 14, # default = 7 max number of dimensions
                           metric = "Gower", # other option is "Euclidean" for cluster::daisy dissimilarity matrix calculation; when using Gower's distance, number of dimensions tested should be lower than number of species
                           dendro = FALSE, # whether the best functional dendrogram should be looked for. default = FALSE
                           plot = "DMB_quality_funct_space") # set the name of the jpeg file for plots illustrating the quality of functional space. NA means no plot

# low meanSD = high quality space
round(qfs$meanSD, 4)

# keep coordinates on 4 dimensions, where meanSD < 0.004
# WHY < 0.004? but also 11-14 have the lowest meanSD, so why are 1:4 chosen?
fd.coord.sgd <- qfs$details_funct_space$mat_coord[,1:4] %>%
  cbind(entity) %>% # rejoin functional entity values for each PC row
  select(PC1:fun_entity)
rownames(fd.coord.sgd) <- fd.coord.sgd$fun_entity # provide rownames
fd.coord.sgd <- fd.coord.sgd[,1:4]
fd.coord.sgd <- data.matrix(fd.coord.sgd) # paste from df to matrix array

# write_csv(fd.coord, here("Data","FE_4D_coord.csv))

# see variance explained by the PCoA axes
gower <- qfs$details_funct_space$mat_dissim

fit <- cmdscale(gower, eig = TRUE, k = 4) # PCoA

# variance explained by the axes
cumsum(fit$eig[fit$eig >= 0]) / sum(fit$eig[fit$eig > 0])



# Functional Richness

## Calculate species and function entity richness

# Species Richness within quadrats
sprich <- myspecies %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  distinct() %>%
  filter(cover > 0) %>%
  mutate(CowTagID = as_factor(CowTagID)) %>%
  dplyr::count(CowTagID, Taxa) %>%
  group_by(CowTagID) %>%
  summarise(counts = sum(n)) %>%
  left_join(quad.label)


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

# species names as rownames
species_entities <- column_to_rownames(.data = species_entities, var = "Taxa")
species_entities <- as.data.frame(species_entities)


## Calculate convex hull (modified Teixido script)

#### Calculate convex hull (Teixido script after I gave up at convhulln)

Fric <- lapply(relative.sgd, function (x) {

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)]

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]

  ch.sgd <- convhulln(m.sgd, options = "FA") # FA: generalized areas and volumes

  chg.sgd <- convhulln(fd.coord.sgd, options = "FA")

  c(length(species.sgd), length(species.sgd)/ncol(sgd.sp)*100, dim(m.sgd)[1], dim(m.sgd)[1]/dim(fd.coord.sgd)[1]*100, ch.sgd$vol/chg.sgd$vol*100)
  #  72 is Teixido's total number of species, so I am dividing by my total number of species (52 at Varari)

})#eo lapply


names(Fric) = relative.sgd

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ) , and the volume among the 3 pH zones
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol8D")
# Number of species per condition
# percentage of species present per condition
# number of functional entities present per condition
# percentage of functional entities present per condition
# percentage of functional space volume occupied




## Plot convex hull (modified Teixido script)
### Bar plot

# order CowTagID's by distance from the seepage point
tagOrder <- meta %>%
  filter(Location == "Varari",
         CowTagID != "Seep",
         CowTagID != "V13") %>%
  arrange(desc(dist_to_seep_m)) %>% # set arrange factor
  select(CowTagID)
tagOrder <- tagOrder$CowTagID[1:19] # set cowtag order as arrange factor order


############################################# plot convex hull
# cols <- c("#3A5FCD", "#FFA500", "#CD2626")
# colstr <- c("#3A5FCD70", "#FFA50070", "#CD262670") # sets same color but at 70% opacity
#
# names(cols) <- c("Low", "Moderate", "High") # name the colors for the ggplot for loop
# names(colstr) <- c("Low", "Moderate", "High")

cols <- pnw_palette("Bay",19,type="continuous")
names(cols) <- tagOrder


# Figure 1. Species and functional diversity changes among pH zones.
# All volumes in distinct plots
p <- Fric %>%
  as_tibble() %>%
  mutate(CowTagID = relative.sgd) %>%
  select(Sp = NbSpP, FE = NbFEsP, Vol4D = Vol8D, CowTagID) %>%
  mutate(CowTagID = factor(CowTagID, levels = tagOrder)) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Parameters, y = Values, fill = CowTagID)) +
  geom_col(color = "black") +
  facet_wrap(~CowTagID) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylim(0,100) +
  scale_fill_manual(values = cols) +
  labs(fill = "Survey Location", x = "", y = "Relative richness (%)")
p

ggsave(here("Output", "Teixido", "Teixido_Figure1barplot_dmb_CowTags.png"), p, width = 6, height = 5)

### Functional space using PCoA (modified Teixido script)

## so far unable to populate plots that show differences across the three sgd levels. need to reevaluate the for loop

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

qAll <- ggplot(data = All.ch.tib, aes(x = x, y = y)) +
  geom_polygon(aes(fill = CowTagID, color = CowTagID), alpha = 0.5) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2") +
  geom_point(data = as_tibble(All.m.sgd), aes(x = PC1, y = PC2, color = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~CowTagID) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
qAll
# I want to view this on a map. Can I put polygons as map points?

ggsave(here("Output", "Teixido", "Teixido_Figure1volume_dmb_CowTags.png"), qAll, width = 6, height = 5)




### Linear Regression Sp and FE and Vol4D ~ SGD
reg.Fric <- Fric %>%
  as_tibble() %>%
  mutate(CowTagID = relative.sgd) %>%
  relocate(CowTagID, .before = NbSp) %>%
  left_join(meta)

summary(lm(data = reg.Fric, NbSp ~ dist_to_seep_m)) # **
summary(lm(data = reg.Fric, NbFEs ~ dist_to_seep_m)) # **
summary(lm(data = reg.Fric, Vol8D ~ dist_to_seep_m)) # not significant

summary(lm(data = reg.Fric, NbSp ~ N_percent)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ N_percent)) # not significant
summary(lm(data = reg.Fric, Vol8D ~ N_percent)) # not significant

summary(lm(data = reg.Fric, NbSp ~ del15N)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ del15N)) # . 0.099 not significant
summary(lm(data = reg.Fric, Vol8D ~ del15N)) # not significant

summary(lm(data = reg.Fric, NbSp ~ meanRugosity)) # **
summary(lm(data = reg.Fric, NbFEs ~ meanRugosity)) # *
summary(lm(data = reg.Fric, Vol8D ~ meanRugosity)) # **

summary(lm(data = reg.Fric, NbSp ~ adj_CT_depth_cm)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ adj_CT_depth_cm)) # not significant
summary(lm(data = reg.Fric, Vol8D ~ adj_CT_depth_cm)) # not significant




### Calculate regressions again with Residuals (remove structure/substrate)

#fit model: linear relationship
resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
resModVol <- lm(Vol8D ~ meanRugosity, data=reg.Fric) # entity volume

#view model summary
summary(resModSpR) # strong significance
summary(resModFER) # significance
summary(resModVol) # strong significance

#calculate the standardized residuals
resSpR <- residuals(resModSpR)
resFER <- residuals(resModFER)
resVol <- residuals(resModVol)

#column bind standardized residuals back to original data frame
res_data <- reg.Fric %>%
  cbind(resSpR) %>%
  cbind(resFER) %>%
  cbind(resVol)



#plot predictor variable vs. standardized residuals
rug_res_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = res)) +
  geom_point() +
  geom_smooth(method="lm", formula = y~poly(x,2)) +
  geom_text_repel(aes(label = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Residuals (% Cover Filamentous Algae ~ Rugosity)",
       x = "Distance to seep (m)")
rug_res_plot
summary(lm(data = res_data, res ~ poly(dist_to_seep_m,2)))
summary(lm(data = res_data, resFER ~ dist_to_seep_m))

reg.Fric %>%
  ggplot(aes(x = dist_to_seep_m, y = NbFEs/NbSp)) +
  geom_point(aes(color = dist_to_seep_m)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  geom_label(aes(label = CowTagID))

summary(lm(data = reg.Fric, (NbFEs/NbSp) ~ dist_to_seep_m))

# ratio of FE:Sp richness ~ distance to seep
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity


## Can use the three values above, and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness







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

ggsave(here("Output","Teixido", "SiteMap_Volume.png"), volSites, width = 10, height = 9)











