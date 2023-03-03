### FE Redundancy: bar plot and nMDS comparing FE redundancy across distance from seep

### Created by Danielle Barnas
### Crated on March 1, 2023

library(tidyverse)
library(here)
library(patchwork)
library(PNWColors)
library(ggrepel)
library(stats)



traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
comp <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv"))
#Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
#Fric <- data.matrix(column_to_rownames(Fric, var = "CowTagID")) # convert Fric from tibble to matrix array

mypalette <- pnw_palette(name = "Bay", n = 20)
scatterpalette <- mypalette[c(13, 15, 29, 14, 28, 13, 27, 12, 26, 11, 25, 10, 24, 9, 23, 8, 22, 7, 21, 6, 20, 5, 19, 4, 18, 3, 17, 2, 16, 1)]


comp <- comp %>%
  filter(Location == "Varari", # only analyze varari for now
         #CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  select(Location:SpeciesCounts)

dist <- meta %>%
  filter(Location == "Varari") %>%
  select(CowTagID, dist_to_seep_m)

# only cowtag ID's
quad.label <- comp %>%
  left_join(meta) %>%
  select(Location, CowTagID, dist_to_seep_m) %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  distinct() %>%
  arrange(dist_to_seep_m) %>%
  select(CowTagID)
dist_order <- quad.label$CowTagID

# create attribute for Functional Entity of each species
traits <- comp %>%
  full_join(traits) %>% # apply traits to species composition df
  filter(Location == "Varari",
         #CowTagID != "VSEEP",
         CowTagID != "V13",
         #Identified == "yes",
         Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Taxa, Taxon_Group,
         #Morphology,
         #Life_Span,Max_size,Growth_Rate, # still have NAs
         #Zooxanthellate:Feeding_Mode
         Morph:FM) # use reduced names trait groups
fe_traits <- traits %>%
  unite(col = "FE", Morph2:FM, sep = ",", remove = T) %>%
  relocate(FE, .after = Taxon_Group) %>%
  distinct()
fe_traits


FER_pres_abs <- comp %>%
  left_join(fe_traits) %>%
  drop_na(FE) %>% # remove abiotic 'taxa' and also cyanobacteria
  group_by(Location, CowTagID) %>%
  dplyr::count(FE) %>%
  arrange(FE)
FER_pres_abs$CowTagID <- parse_factor(FER_pres_abs$CowTagID, levels = dist_order)

FER_pres_abs_plot <- FER_pres_abs %>%
  arrange(CowTagID) %>%
  ggplot(aes(x = CowTagID, y = n, fill = FE)) +
  geom_col(position = "stack", color = "white") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = mypalette) +
  labs(y = "FER", x = "Survey Location")
FER_pres_abs_plot

ggsave(here("Output", "PaperFigures", "FER_pres_abs.png"), FER_pres_abs_plot, width = 9, height = 6)


### analyze redundancy over distance to seep

FER_counts <- comp %>%
  left_join(fe_traits) %>%
  drop_na(FE) %>% # remove abiotic 'taxa' and also cyanobacteria
  # group_by(Location, CowTagID) %>%
  # dplyr::count(FE) %>%
  arrange(FE) %>%
  left_join(meta) %>%
  select(Location, CowTagID, dist_to_seep_m, FE) %>%
  mutate(PA = 1) %>%
  filter(Location == "Varari",
         CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  arrange(CowTagID, dist_to_seep_m) %>%
  mutate(FE = as_factor(FE))
FER_counts$CowTagID <- parse_factor(FER_counts$CowTagID, levels = dist_order)
summary(FER_counts)
## regression: FE counts ~ distance to seep (via CowTags)


FER_counts <- FER_pres_abs %>%
  mutate(FE = as_factor(FE)) %>%
  rename(Counts = n) %>%
  left_join(dist)
mod <- summary(lm(Counts ~ dist_to_seep_m*FE, data = FER_counts))
as_tibble(rownames_to_column(as.data.frame(mod$coefficients), var = "Variables")) %>%
  select(Variables, pval= 'Pr(>|t|)') %>%
  filter(pval < 0.06)

FER_counts %>%
  ggplot(aes(x = dist_to_seep_m, y = Counts, color = FE)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm") +
  facet_wrap(~FE)



facet_fe_redundancy_plot <- FER_counts %>%
  arrange(CowTagID, dist_to_seep_m) %>%
  ggplot(aes(x = Counts, fill = FE)) +
  geom_bar(color = "white") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect("white")) +
  facet_wrap(~CowTagID) +
  scale_fill_manual(values = mypalette) +
  labs(y = "FE Count Tally", x = "FE Counts (Redundancy)")
facet_fe_redundancy_plot
ggsave(here("Output", "PaperFigures", "fe_redundacy_bar_plot.png"), facet_fe_redundancy_plot, width = 8, heigh = 7)












###########################################################
### redundancy / presence.absence nMDS of FE
###########################################################


FE_nmds_data <- FER_pres_abs %>%
  left_join(meta) %>%
  select(Location, CowTagID, dist_to_seep_m, FE, Counts = n) %>%
  filter(Location == "Varari",
         CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  arrange(CowTagID, dist_to_seep_m) %>%
  mutate(Counts = as.double(Counts)) %>%
  pivot_wider(names_from = FE, values_from = Counts) %>%
  mutate_at(vars(4:23), .funs = ~if_else(is.na(.), 0, .))

# will cbind cowtags later
# set levels as numerical order of plates
CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')
FE_nmds_data$CowTagID <- factor(FE_nmds_data$CowTagID, levels = CTlevels)

# arrange by cowtag and then remove for nmds
FE_nmds_data <- FE_nmds_data %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  select(-c(Location,CowTagID, dist_to_seep_m))


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
  scale_color_gradient(low = "red", high = "yellow") +
  labs(color = "Distance to seep (m)")
nMDSplot

ggsave(here("Output", "PaperFigures", "FE_redundancy_nmds_plot.png"), nMDSplot, width = 10, height = 10)

### PERMANOVA
richPermFull <- cbind(Groupb, FE_nmds_data) %>%  # bind cowTagIDs
  rename(CowTagID = Groupb) %>%
  left_join(joinDF) %>%
  mutate(relDist = if_else(dist_to_seep_m <= 50, "Near", if_else(dist_to_seep_m > 120, "Far", "Mid")))

# dist < 26 > 100 **0.003 V20, V17, V14 near vs mid 0.006
# dist < 47 > 100 **0.01 V20, V17, V14, V9 near vs mid 0.024
# dist < 50 > 100 insignif. V20, V17, V14, V9, V10

# dist < 26 > 120 **0.004 near vs mid 0.003
# dist < 47 > 120 **0.006 near vs mid 0.003


permanovamodel<-adonis2(richPermFull[,2:21]~relDist, richPermFull, permutations = 999,
                        method="bray") # should change out cowtagid with some grouping name
permanovamodel


#If we are to trust the results of the permanova, then we have to assume that the dispersion among
#data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(richPermFull[,2:21])
betadisper(disper, richPermFull$relDist)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

pairwise.adonis(richPermFull[2:21], richPermFull$relDist, perm=999)

#Get coefficients to see which species are most important in explaining site differences:
#permanovamodel$coefficients


###########################################################
### redundancy of individual FE trait groups
###########################################################


red_traits <- comp %>%
  left_join(traits) %>%
  select(-Date) %>%
  distinct() %>% # remove accidental duplicate rows
  drop_na() %>%  # remove abiotic 'taxa'
  group_by(Location, CowTagID) # stay grouped for further processing
red_traits$CowTagID <- parse_factor(red_traits$CowTagID, levels = dist_order)


FE_Trait_Fxn <- function(mytrait){
  fer_trait <- red_traits %>%
    dplyr::count({{mytrait}}) %>%
    arrange(CowTagID)

  myfill <- enquo(mytrait)

  myplot <- fer_trait %>%
    arrange(CowTagID) %>%
    ggplot(aes(x = CowTagID, y = n, fill = !!myfill)) +
    geom_col(position = "stack", color = "white") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = mypalette)  +
    labs(y = "FER", x = "Survey Location")

  return(myplot)

}


FE_Trait_Fxn(Morph)

fer_morph <- red_traits %>%
  dplyr::count(Morph) %>%
  arrange(Morph)

fer_morph %>%
  arrange(CowTagID) %>%
  ggplot(aes(x = CowTagID, y = n, fill = Morph)) +
  geom_col(position = "stack", color = "white") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = mypalette) +
  labs(y = "FER", x = "Survey Location")




















###########################################################
### redundancy / presence.absence nMDS of FE
###########################################################

# will cbind cowtags later
# set levels as distance from seep order of plates
CTlevels <- dist_order
#CTlevels <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V14','V15','V16','V17','V18','V19','V20')
FE_nmds_dist_data <- FER_pres_abs %>%
  left_join(meta) %>%
  select(Location, CowTagID, dist_to_seep_m, FE, Counts = n) %>%
  filter(Location == "Varari",
         CowTagID != "VSEEP",
         CowTagID != "V13") %>%
  distinct() %>%
  arrange(CowTagID) %>%
  mutate(Counts = as.double(Counts)) %>%
  pivot_wider(names_from = FE, values_from = Counts) %>%
  mutate_at(vars(4:23), .funs = ~if_else(is.na(.), 0, .)) %>%
  ungroup()
# order cowtags by distance to seep as factor levels
FE_nmds_dist_data$CowTagID <- factor(FE_nmds_dist_data$CowTagID, levels = CTlevels)
FE_nmds_dist_data <- FE_nmds_dist_data %>%
  arrange(CowTagID) %>%
  pivot_longer(cols = 4:23, names_to = "FE", values_to = "Counts") %>%
  pivot_wider(names_from = CowTagID, values_from = Counts, id_cols = FE) %>%  # id_cols sets the unique identifier for data
  arrange(FE)

# get rownames
fe_order <- FE_nmds_dist_data$FE

# remove FE for nmds
FE_nmds_dist_data <- FE_nmds_dist_data %>%
  select(-c(FE)) # remove non-numerical column


ord1 <- metaMDS(FE_nmds_dist_data, k=2, distance='bray')

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
Groupb <- as.character(fe_order) # assign fe from rownames
MDS1b <- ord1$points[,1] # MDS1 for CowTagID
MDS2b <- ord1$points[,2] # MDS2 for CowTagID
Datab <- as_tibble(cbind(Groupb, MDS1b, MDS2b)) %>%  # bind all cols into tibble
  mutate(MDS1b = as.numeric(MDS1b), # as numeric
         MDS2b = as.numeric(MDS2b)) %>%
  rename(FE = Groupb)

joinDF <- meta %>%
  select(Group = CowTagID, dist_to_seep_m)

Data <- Data %>%
  left_join(joinDF)


## plot
nMDSplot <- ggplot(data = Data,
                   aes(x = MDS1,
                       y = MDS2)) +
  geom_point(aes(color = dist_to_seep_m)) +
  geom_point(data = Datab,
             aes(x = MDS1b,
                 y = MDS2b),
             size = 3) +
  geom_text_repel(data = Data, # site characteristics
                  aes(x = MDS1,
                      y = MDS2,
                      label = Group),
                  size = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  geom_label_repel(data = Datab, # site characteristics
                   aes(x = MDS1b,
                       y = MDS2b,
                       label = FE),
                   size = 2,
                   max.overlaps = 16) + # increase from 10 because too many FEs overlapping
  scale_color_gradient(low = "red", high = "yellow") +
  labs(color = "Distance to seep (m)")
nMDSplot

### PERMANOVA
richPermFull <- cbind(Groupb, FE_nmds_dist_data) %>%  # bind cowTagIDs
  rename(FE = Groupb) #%>%
  #left_join(joinDF) %>%
  #mutate(relDist = if_else(dist_to_seep_m <= 50, "Near", if_else(dist_to_seep_m > 120, "Far", "Mid")))


permanovamodel<-adonis2(richPermFull[,2:20]~FE, richPermFull, permutations = 999,
                        method="bray") # should change out cowtagid with some grouping name
permanovamodel


#If we are to trust the results of the permanova, then we have to assume that the dispersion among
#data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(richPermFull[,2:20])
betadisper(disper, richPermFull$FE)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

pairwise.adonis(richPermFull[2:20], richPermFull$FE, perm=999)

#Get coefficients to see which species are most important in explaining site differences:
#permanovamodel$coefficients

