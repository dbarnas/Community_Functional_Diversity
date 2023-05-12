### IDENTIFY WHICH SPECIES DRIVE PATTERNS ALONG SGD
### In response to Pete's suggestions about reducing how many species need to be functionally categorized
### Created by Danielle Barnas
### Created on March 14, 2023

##### Load libraries #####
library(tidyverse)
library(here)
library(vegan)
library(ggrepel)
library(patchwork)
library(PNWColors)
library(pairwiseAdonis)


##### Read in data #####
meta <- read_csv(here("Data", "Full_Metadata.csv"))
spec <- read_csv(here("Data", "Species_FE.csv"))
fun <- read_csv(here("Data", "Distinct_FE.csv"))
comp <- read_csv(here("Data", "Species_Abundances_wide.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv"))
shore <- read_csv(here("Data", "Shore_distance.csv"))


##### Prepare data #####
spfun <- comp %>%
  pivot_longer(cols = 2:ncol(comp), names_to = "Taxa", values_to = "pcover") %>%
  left_join(spec) %>%
  left_join(fun) %>%
  drop_na(FE) %>%
  select(CowTagID, Taxa, pcover, FE:ER)
  # calculate percent cover

meta <- chem %>%
  left_join(meta) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  left_join(shore)


##### Create functions for analysis #####
nMDS_species <- function(mdsord){

  # get points for species
  Group <- rownames(mdsord$species) # get characteristic names
  MDS1 <- c(mdsord$species[,1]) # MDS1 for characteristics
  MDS2 <- c(mdsord$species[,2]) # MDS2 for characteristics
  Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
    mutate(MDS1 = as.numeric(MDS1), # as numeric
           MDS2 = as.numeric(MDS2)) %>%
    #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
    select(MDS1, MDS2, Group)

  return(Data)
}

nMDS_points <- function(mdsord, param_df, param_select){
  # get points for
  Groupb <- as.character(CTlevels) # assign CowTagID
  MDS1b <- mdsord$points[,1] # MDS1 for CowTagID
  MDS2b <- mdsord$points[,2] # MDS2 for CowTagID
  Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
    mutate(MDS1b = as.numeric(MDS1b), # as numeric
           MDS2b = as.numeric(MDS2b)) %>%
    rename(CowTagID = Groupb)

  joinDF <- param_df %>%
    select(param_select)

  Datab <- Datab %>%
    left_join(joinDF)

  return(Datab)
}


##### Run nMDS #####

# CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20','VSEEP')
CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')

# use species percent cover to look at differences across sites
# perm1 <- comp
perm1 <- comp %>% filter(CowTagID!= "VSEEP")

# arrange data in order of CowTagID for ordination plot later
perm1$CowTagID <- factor(perm1$CowTagID, levels = CTlevels)
perm1 <- perm1 %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  select(-CowTagID)

ord1 <- metaMDS(perm1, k=2, distance='bray')

# stress with k=2 dimensions. Is it < 0.3?
ord1$stress

# stress plot - want to minimize scatter
stressplot(ord1)

#param_mds <- nMDS_species(ord1)
#param_mds_cat <- nMDS_points(ord1, chem, allChem)

# get points for species
Group <- rownames(ord1$species) # get characteristic names
MDS1 <- c(ord1$species[,1]) # MDS1 for characteristics
MDS2 <- c(ord1$species[,2]) # MDS2 for characteristics
Data <- as_tibble(cbind(Group, MDS1, MDS2)) %>%  # bind all cols into tibble
  mutate(MDS1 = as.numeric(MDS1), # as numeric
         MDS2 = as.numeric(MDS2)) %>%
  #mutate(Taxon_Group = if_else(Taxa == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
  select(MDS1, MDS2, Group)

# get points cowtagid's
Groupb <- as.character(CTlevels) # assign CowTagID
MDS1b <- ord1$points[,1] # MDS1 for CowTagID
MDS2b <- ord1$points[,2] # MDS2 for CowTagID
Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
  mutate(MDS1b = as.numeric(MDS1b), # as numeric
         MDS2b = as.numeric(MDS2b)) %>%
  rename(CowTagID = Groupb)

joinDF <- meta %>%
  select(CowTagID, dist_to_seep_m, shore_dist_m)

Datab <- Datab %>%
  left_join(joinDF)

# based on permanova classifications, classify sites by relative distance from seep
Datab <- Datab %>%
  mutate(relDist = if_else(dist_to_seep_m <= 26, "Near", if_else(dist_to_seep_m > 80, "Far", "Mid")))

ggplot(data = Data,
                   aes(x = MDS1,
                       y = MDS2)) +
  geom_point(color = "black") +
  geom_point(data = Datab,
             aes(x = MDS1b,
                 y = MDS2b,
                 color = relDist),
             size = 3) +
  geom_text_repel(data = Data, # site characteristics
                  aes(x = MDS1,
                      y = MDS2,
                      label = Group),
                  size = 3) +
  geom_text_repel(data = Datab, # site characteristics
                  aes(x = MDS1b,
                      y = MDS2b,
                      label = CowTagID),
                  size = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  ggforce::geom_mark_ellipse(
    data = Datab,
    aes(x = MDS1b,
        y = MDS2b,
        fill = relDist, label = relDist),
    alpha = .15, show.legend = FALSE,
    #label.buffer = unit(1, "mm")
    )
  ##coord_cartesian(xlim = c(-4, 7), ylim = c(-6, 4)) +
  ##scale_shape_manual(values = c(22,16))+
  ##scale_colour_hue(l = 45)+ # l: luminance/lightness 0-100
  ##scale_fill_hue(l = 45)
  #scale_color_gradient(low = "red", high = "yellow") # set color scale


### PERMANOVA
PermFull <- cbind(Groupb, perm1) %>%  # bind cowTagIDs
  rename(CowTagID = Groupb) %>%
  left_join(joinDF) %>%
  mutate(relDist = if_else(dist_to_seep_m <= 26, "Near", if_else(dist_to_seep_m > 80, "Far", "Mid")))


permanovamodel<-adonis2(PermFull[,2:52]~relDist, PermFull, permutations = 999,
                        method="bray") # should change out cowtagid with some grouping name
permanovamodel

pairwise.adonis(PermFull[,2:52], PermFull$relDist, perm=999)
#significant diff between Near and all farther sites (V14, V17, and V20 cluster away from all else)


### Let's look at PCA

#species composition at each site
# arrange data in order of CowTagID for ordination plot later
comp1 <- comp
comp1$CowTagID <- factor(comp1$CowTagID, levels = CTlevels)
comp1 <- comp1 %>%
  arrange(CowTagID) %>%
  ungroup()
perm1 # numerical species comp only, ordered by site numerically


# Run the PCA
pca_full <- prcomp(perm1, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_full <- as_tibble(pca_full$x[,1:2]) # scores for each sample location
PC_loadings_full <- as_tibble(pca_full$rotation) %>% # loadings for each parameter
  bind_cols(labels = rownames(pca_full$rotation))

# Put it with all the original data
pca_fulldata_all <- comp %>%
  filter(CowTagID != "VSEEP") %>%
  bind_cols(PC_scores_full) %>%
  left_join(joinDF) %>%
  mutate(relDist = if_else(dist_to_seep_m <= 26, "Near", if_else(dist_to_seep_m > 80, "Far", "Mid")))

## for easier plot viewing, let's see only species with >=5% cover (10 points)
species5 <- pca_fulldata_all %>%
  pivot_longer(cols = 2:52 , names_to = "species", values_to = "pcover") %>%
  filter(pcover >= 5) %>%
  select(CowTagID, labels = species, pcover) %>%
  left_join(PC_loadings_full) %>%
  select(CowTagID, labels, pcover, PC1, PC2)


## scores plot

p1_full <- pca_fulldata_all %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_point(shape = 21,
             color = "black",
             aes(fill = dist_to_seep_m),
             size = 4) +
  scale_fill_gradient(low = "red", high = "yellow") +
  #scale_fill_manual(values = mypalette) +
  theme_bw()+
  theme(legend.position = "top",
        #panel.background = element_rect(fill = "azure4"),
        panel.grid = element_blank()) +
  labs(title = "Varari Scores Plot", fill = "Distance to seep (m)") +
  geom_text_repel(aes(label = CowTagID),
                  size = 4)
p1_full


## loadings plot

# loadings plot
p1_full_loading <- ggplot(data = species5,  # PC_loadings_full
         aes(x = PC1,
             y = PC2,
             #label = labels
             )) + # labels are the parameters (species in this case)
  geom_segment(data = species5,
               aes(x = 0,
                   y = 0,
                   xend = PC1*15,
                   yend = PC2*15,
                   color = labels),
               arrow = arrow(length = unit(0.1,"cm")),
               #color = "grey"
               ) +
  annotate("text",
           x = species5$PC1*15+0.1,
           y = species5$PC2*15+0.1,
           label = species5$labels,
           size = 3) +
  #coord_cartesian(xlim = c(-4, 5), ylim = c(-6, 5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Varari Loadings") +
  geom_point(data = pca_fulldata_all,
             shape = 21,
             color = "black",
             size = 4,
             aes(fill = relDist)) +
  geom_text_repel(data = pca_fulldata_all,
                  aes(x = PC1, y = PC2,
                      label = CowTagID),
                  size = 2)
  #scale_fill_gradient(low = "red", high = "yellow")

p1_full_loading





### Patch PCA plots
VarariPCA <- p1_full + p1_full_loading +
  patchwork::plot_annotation("Varari PCA and Loadings (w/ Residuals)",
                             theme = theme(plot.title = element_text(size = rel(1.5),
                                                                     face = "bold",
                                                                     hjust = 0.5,
                                                                     margin = margin(t = 10, b = 20,
                                                                                     unit = "pt"))))
VarariPCA




