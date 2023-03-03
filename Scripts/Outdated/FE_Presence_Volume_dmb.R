### Clustering as  Low, Moderate, or High levels of SGD


library(tidyverse)
library(here)
library(FD)
library(tripack) # Triangulation of Irregularly Spaced Data
library(geometry) # Mesh Generation and Surface Tessellation
library(matrixStats) # Functions that Apply to Rows and Columns of Matrices (and to Vectors)
library(patchwork)



traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
myspecies <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))



sgd.label <- myspecies %>%
  select(CowTagID) %>%
  mutate(relSGD =
           if_else(CowTagID == "VSEEP" |
                     CowTagID == "V14",
                   # CowTagID == "V17" |
                   # CowTagID == "V15" |
                   # CowTagID == "V11",
                   "High",
                   if_else(CowTagID == "V1" |
                             CowTagID == "V2" |
                             CowTagID == "V3" |
                             CowTagID == "V4" |
                             CowTagID == "V16", "Low", "Moderate"))) # any unlisted cow tag id's will be labeled as Moderate



# create attribute for Functional Entity of each species
traits <- traits %>%
  full_join(myspecies) %>%
  filter(Location == "Varari",
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
  filter(Location == "Varari", # only analyze Varari sites for now
         Taxa != "Bare Rock" & # only include biological data
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



# select only for species and funcitonal entities
species_entities <- traits %>%
  select(Taxa, fun_entity)

# df with unique functional entities for each row
entity <- traits %>%
  select(-c(Taxa, Taxon_Group)) %>%
  distinct() %>%
  relocate(fun_entity, .before = Morphology) %>%
  mutate_all(.funs = as_factor) # groups need to be factors to run quality_funct_space()

# relative SGD levels
relative.sgd <- c("Low", "Moderate", "High")


## Computing multidimensional functional space

# source function from Teixido: function for computing the quality of functional dendrogramm and multidimensional functional spaces
source(here("Scripts", "Teixido", "quality_funct_space.R"))

qfs <- quality_funct_space(mat_funct = entity, # distinct functional trait entities (NAs not allowed. Traits can be different types: numeric, ordinal, nominal)
                           traits_weights = NULL, # defualt = same weight for all traits
                           nbdim = 14, # default = 7 max number of dimensions
                           metric = "Gower", # other option is "Euclidean" for cluster::daisy dissimilarity matrix calculation
                           dendro = FALSE, # whenter the best functional dendrogram should be looked for. default = FALSE
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

# Species Richness across sgd levels
sprich <- myspecies %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  distinct() %>%
  filter(cover > 0) %>%
  mutate(CowTagID = as_factor(CowTagID)) %>%
  dplyr::count(CowTagID, Taxa) %>%
  group_by(CowTagID) %>%
  summarise(counts = sum(n)) %>%
  left_join(sgd.label)

# anova of species richness
mod <- aov(data = sprich, counts ~ relSGD)
anova(mod)
TukeyHSD(mod)

# Set categories of relative SGD (based on clustering or PERMANOVA from biogeochemistry)
sgd.sp <- myspecies %>%
  left_join(sgd.label) %>%
  relocate(relSGD, .after = CowTagID) %>%
  pivot_longer(cols = Turf:ncol(.), names_to = "Taxa", values_to = "cover") %>% # reset to pivot by SGD
  select(-CowTagID) %>%
  group_by(relSGD, Taxa) %>% # need to sum up percent cover by relative sgd and taxa
  mutate(cover = sum(cover)) %>%
  ungroup() %>%
  distinct() %>% # remove duplicate values after summing up percent cover
  group_by(relSGD) %>%
  pivot_wider(names_from = Taxa, values_from = cover)  # pivot species wide by relative SGD

# relative sgd as rownames
sgd.sp <- column_to_rownames(.data = sgd.sp, var = "relSGD")
sgd.sp <- as.data.frame(sgd.sp)

# species names as rownames
species_entities <- column_to_rownames(.data = species_entities, var = "Taxa")
species_entities <- as.data.frame(species_entities)



## Calculate convex hull (modified Teixido script)

#### Calculate convex hull (Teixido script after I gave up at convhulln)

Fric <- lapply(relative.sgd, function (x) {

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)]

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd,]

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

############################################# plot convex hull
cols <- c("#3A5FCD", "#FFA500", "#CD2626")
colstr <- c("#3A5FCD70", "#FFA50070", "#CD262670")

names(cols) <- c("Low", "Moderate", "High") # name the colors for the ggplot for loop
names(colstr) <- c("Low", "Moderate", "High")


# Figure 1. Species and functional diversity changes among pH zones.
# All volumes in distinct plots

p <- Fric %>%
  as_tibble() %>%
  mutate(relSGD = relative.sgd) %>%
  select(Sp = NbSpP, FE = NbFEsP, Vol4D = Vol8D, relSGD) %>%
  mutate(relSGD = factor(relSGD, levels = c("Low", "Moderate", "High"))) %>%
  pivot_longer(cols = 1:3, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Parameters, y = Values, fill = relSGD)) +
  geom_col(color = "black") +
  facet_wrap(~relSGD, nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = cols) +
  labs(fill = "Relative SGD", x = "", y = "Relative richness (%)")


### Functional space using PCoA (modified Teixido script)

## so far unable to populate plots that show differences across the three sgd levels. need to reevaluate the for loop

q <- list()

for(i in relative.sgd) {

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[i,] > 0)]
  # only species present in each treatment

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]
  m.sgd <- data.matrix(m.sgd) # parse from data frame to matrix array

  tr.sgd <- tri.mesh(m.sgd[,1],m.sgd[,2], duplicate = "strip") # duplicate: default = "error", "strip" = removes all duplicate points, "remove" = leaves one point of duplicate points

  ch.sgd <- convex.hull(tr.sgd)

  ch.tib <- cbind(ch.sgd$x, ch.sgd$y, ch.sgd$i) # parse as tibble df
  colnames(ch.tib) <- c("x", "y", "i")
  ch.tib <- as_tibble(ch.tib)

  if(i == "Low"){
  q[[i]] <- ggplot(data = ch.tib, aes(x = x, y = y)) +
    geom_polygon(fill = colstr[i], color = cols[i]) + # create polygon using product of convex.hull(tri.mesh)
    labs(x = "PCoA 1", y = "PCoA 2") +
    geom_point(data = as_tibble(m.sgd), aes(x = PC1, y = PC2),
               color = cols[i]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))
  } else {
    q[[i]] <- ggplot(data = ch.tib, aes(x = x, y = y)) +
      geom_polygon(fill = colstr[i], color = cols[i]) + # create polygon using product of convex.hull(tri.mesh)
      labs(x = "PCoA 1", y = "") +
      geom_point(data = as_tibble(m.sgd), aes(x = PC1, y = PC2),
                 color = cols[i]) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2]))

  }
}


pq <- (p) / (q[["Low"]] + q[["Moderate"]] + q[["High"]])
ggsave(here("Output", "Teixido", "Teixido_Figure1_dmb.png"), pq, width = 6, height = 5)


### Intersecting functional space using PCoA (modified Teixido script)


# Supplementary Figure 1. Intersection of the three functional volumes among pH zones.
# all volumes in 1 figure,

ch.full <- tibble(x = as.numeric(),
                  y = as.numeric(),
                  i = as.numeric(),
                  relSGD = as.character())
m.full <- tibble(PC1 = as.numeric(),
                 PC2 = as.numeric(),
                 PC3 = as.numeric(),
                 PC4 = as.numeric(),
                 relSGD = as.character())

for (i in relative.sgd) {

  sgd <- i

  species.sgd <- colnames(sgd.sp)[which(sgd.sp[i,] > 0)]

  fes_cond.sgd <- species_entities[rownames(species_entities) %in% species.sgd, ]

  m.sgd <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond.sgd, ]

  tr.sgd <- tri.mesh(m.sgd[,1],m.sgd[,2], duplicate = "strip")

  ch.sgd <- convex.hull(tr.sgd)

  ch.tib <- cbind(ch.sgd$x, ch.sgd$y, ch.sgd$i) # parse as tibble df
  colnames(ch.tib) <- c("x", "y", "i")

  ch.tib <- as_tibble(ch.tib) %>%
    mutate(relSGD = sgd)

  m.sgd <- as_tibble(m.sgd) %>%
    mutate(relSGD = sgd)
  m.full <- m.full %>%
    rbind(m.sgd)

  ch.full <- ch.full %>% rbind(ch.tib)
}

r <- ch.full %>%
  mutate(relSGD = factor(relSGD, levels = c("Low", "Moderate", "High"))) %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(color = relSGD, fill = relSGD), alpha = 0.3) + # create polygon using product of convex.hull(tri.mesh)
  labs(x = "PCoA 1", y = "PCoA 2", fill = "Relative SGD", color = "Relative SGD") +
  geom_point(data = m.full, aes(x = PC1, y = PC2, color = relSGD)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(min(fd.coord.sgd[,1]), max(fd.coord.sgd[,1])) + ylim(min(fd.coord.sgd[,2]), max(fd.coord.sgd[,2])) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = colstr)

r
ggsave(here("Output", "Teixido", "Teixido_S1_dmb.png"), r, width = 6, height = 5)


## Null Model of Functional Richness among SGD zones


############### null model of Functional richness among pH zones

sgd.sp <- data.matrix(sgd.sp) # replace ab.conditions - must be matrix class

n_perm = 100

spe_fes_r = species_entities # spe_fes

Fric_perm <- lapply(relative.sgd, function (x) { # condition

  species <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions


  perm <- sapply((1:n_perm), function (z) {


    spe_fes_r$fun_entity <- sample(species_entities$fun_entity) # spe_fes

    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]

    m <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond,]

    ch <- convhulln(m, options = "FA")

    chg <- convhulln(fd.coord.sgd, options = "FA")

    c(length(species), length(species)/ncol(sgd.sp)*100, dim(m)[1], dim(m)[1]/dim(fd.coord.sgd)[1]*100, ch$vol/chg$vol*100)



  })#eo sapply

  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")


  perm

})#eo lapply

names(Fric_perm) = relative.sgd # condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {

  rowQuantiles(x, probs=c(0.05, 0.95))

})#eo lapply



Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply # condition
Fric$upperFE <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(relative.sgd, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- relative.sgd
relative.sgd <- factor(relative.sgd, levels = c("Low", "Moderate", "High"))

Fric$cond <- as.factor(relative.sgd)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFE","NbFEP", "Vol8D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")




#Plot the null model
#Supplementary Figure 2. Null model of functional richness (functional volume) among pH zones.


tiff(filename="Output/Teixido/Figure_S2_dmb.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)


plot(Vol8D ~ cond, data=Fric, border="white", xlab="", ylab="Relative Richness (%)", ylim=c(0,100))
points(Vol8D ~ cond, data=Fric, pch=16, col=cols[c(1,2,3)], cex=2)

Low <- rbind(Fric[1,], Fric[1,])
Low$Vol8D <- c(Low$lowerVol[1], Low$upperVol[1])
lines(Vol8D ~ cond, data=Low, lwd=3, col=cols["Low"])

Moderate <- rbind(Fric[2,], Fric[2,])
Moderate$Vol8D <- c(Moderate$lowerVol[1], Moderate$upperVol[1])
lines(Vol8D ~ cond, data=Moderate, lwd=3, col=cols["Moderate"])

High<- rbind(Fric[3,], Fric[3,])
High$Vol8D <- c(High$lowerVol[1], High$upperVol[1])
lines(Vol8D ~ cond, data=High, lwd=3, col=cols["High"])

dev.off()










################################################ Convex Hull Intersect

#load intersect function to compute convex hull (vertices + volume) of two set of points and their intersection



source(here("Scripts","Teixido","intersect.R"))


mat_int <- Fric <- lapply(relative.sgd, function (x) { # condition

  species <- colnames(sgd.sp)[which(sgd.sp[x,] > 0)] # ab.conditions

  fes_cond <- species_entities[rownames(species_entities) %in% species, ] # spe_fes

  m <- fd.coord.sgd[rownames(fd.coord.sgd) %in% fes_cond,]

  return(m)

})#eo lapply

names(mat_int) = relative.sgd # condition

###############intersect Low with Moderate

Low_int_Mod <- CHVintersect(mat_int[["Low"]],mat_int[["Moderate"]])

# percentage of the Moderate volume within Low
Low_int_Mod$vol[3]/Low_int_Mod$vol[2]




###############intersect Low with High

Low_int_High <- CHVintersect(mat_int[["Low"]],mat_int[["High"]])

#pergcentage of the Extreme Low volume within Ambient
Low_int_High$vol[3]/Low_int_High$vol[2]



###############intersect Moderate with High

Mod_int_High <- CHVintersect(mat_int[["Moderate"]],mat_int[["High"]])

#percentage of the High volume within Low
Mod_int_High$vol[3]/Mod_int_High$vol[2]













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
  mutate(relSGD = relative.sgd) %>% # effectively makes rownames a column
  relocate(relSGD, .before = Turf) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "species", values_to = "presence") %>%
  left_join(species_entities) %>%
  drop_na() %>%  # for now removes cyanobacteria
  group_by(relSGD, fun_entity) %>%
  summarise(presence = sum(presence)) %>%
  mutate(presence = if_else(presence > 0, 1, 0), # binary for presence-absence
         fun_entity = as.character(fun_entity)) %>%
  ungroup() %>%
  pivot_wider(names_from = fun_entity, values_from = presence)

ab.fe.conditions <- column_to_rownames(ab.fe.conditions, var = "relSGD")


# true functional beta; essential: colnames(ab.fe.conditions) == rownames(fd.coord.sgd)
beta.fun <- functional.beta.pair(data.matrix(ab.fe.conditions), fd.coord.sgd, index.family="jaccard")

#######Plot categories of the 15 functional traits across the functional space

#get data to plot the traits
#]fes <- read.csv2("Data/Teixido/Data_FEs.csv", sep=",", dec=".", row.names=1)
fes <- column_to_rownames(entity, var = "fun_entity")

#spe_fes <- read.csv2("Data/Teixido/Data_Species_FEs.csv", sep=";", dec=",", row.names=1)
spe_fes <- column_to_rownames(species_entities, var = "species")

###### Supplementary Figure 6. Distribution of functional trait categories across the functional space

tiff(filename="Output/Teixido/Figure_S6_dmb.tif", height=20, width=30, units="cm", compression = c("lzw"), res=300, pointsize=10)


ftr <- colnames(fes)

par(mfrow=c(3,5))



for (i in ftr) {

  lab <- as.factor(sort(unique(fes[,i])))

  plot(fd.coord.sgd[,1], fd.coord.sgd[,2], pch=16, cex=1.2, col = as.numeric(fes[,i]), xlim = c(-0.3, 0.7),  main=gsub("."," ", i, fixed = T), xlab="PCoA 1", ylab="PCoA 2")
  legend(x=0.45, y=0.35, legend=lab, pch=rep(16, length(lab)), col=as.numeric(as.factor(lab)), bty = "n")


}

dev.off()



















