#### Cluster significant difference analysis ####
### Created by Danielle Barnas
### Created on Aug 29, 2022

#### Load Libraries ####
library(tidyverse)
library(here)
library(stats) # lm() and anova()
library(emmeans)
library(agricolae) # HSD.test()
library(vegan) # permanova
library(pairwiseAdonis)

#### Read in Data ####
cluster <- read_csv(here("Data", "Biogeochem", "Cluster_metadata_FullSeasons_noseep.csv"))
#turbcluster <- read_csv(here("Data", "Surveys", "Cluster_turb_metadata_noSeep.csv"))

#### CLUSTERS FROM FULL BIOGEOCHEM ####

### parse V_clusters to factor
cluster <- cluster %>%
  mutate(cluster = as.factor(cluster))
# turbcluster <- turbcluster %>%
#   rename(V_cluster = V_cluster_turb) %>%
#   mutate(V_cluster = as.factor(V_cluster))

Vclust <- cluster %>% filter(Location == "Varari")
Cclust <- cluster %>% filter(Location == "Cabral")


# manova with all three clusters
modelV <- manova(data = Vclust,
                 cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                       Silicate_umolL, NN_umolL, Ammonia_umolL,
                       M_C, HIX, VisibleHumidic_Like, Tryptophan_Like, Tyrosine_Like,
                       del15N, C_N, N_percent) ~ cluster)
summary(modelV)

modelC <- manova(data = Cclust,
                 cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                       Silicate_umolL, NN_umolL, Ammonia_umolL,
                       M_C, HIX, VisibleHumidic_Like, Tryptophan_Like, Tyrosine_Like,
                       del15N, C_N, N_percent) ~ cluster)
summary(modelC)


#### VARARI
# manova by cluster duo:
VclusterHL <- Vclust %>% filter(cluster %in% c("High", "Low"))
VclusterHM <- Vclust %>% filter(cluster %in% c("High", "Mid"))
VclusterLM <- Vclust %>% filter(cluster %in% c("Low", "Mid"))

modelHLV <- manova(data = VclusterHL,
                   cbind(pH, NN_umolL, Salinity,
                         Silicate_umolL, del15N,
                         Tyrosine_Like, M_C) ~ cluster) # cannot run model when response variables outnumber sample size

modelHMV <- manova(data = VclusterHM,
                   cbind(pH, NN_umolL, Salinity,
                         Silicate_umolL, del15N,
                         Tyrosine_Like, M_C) ~ cluster)

modelLMV <- manova(data = VclusterLM,
                   cbind(pH, NN_umolL, Salinity,
                         Silicate_umolL, del15N,
                         Tyrosine_Like, M_C) ~ cluster)

summary(modelHLV) # p > 0.3
summary(modelHMV) # p > 0.05
summary(modelLMV) # ***


#### CABRAL
# manova by cluster duo:
CclusterHL <- Cclust %>% filter(cluster %in% c("High", "Low"))
CclusterHM <- Cclust %>% filter(cluster %in% c("High", "Mid"))
CclusterLM <- Cclust %>% filter(cluster %in% c("Low", "Mid"))

modelHLC <- manova(data = CclusterHL,
                   cbind(Tyrosine_Like, VisibleHumidic_Like,
                         Phosphate_umolL) ~ cluster) # , Tryptophan_Like, M_C, N_percent, NN_umolL

modelHMC <- manova(data = CclusterHM,
                   cbind(Tyrosine_Like, VisibleHumidic_Like,
                         Phosphate_umolL) ~ cluster) # cannot run model when response variables outnumber sample size

modelLMC <- manova(data = CclusterLM,
                   cbind(Tyrosine_Like, VisibleHumidic_Like,
                         Phosphate_umolL) ~ cluster) # , Tryptophan_Like, M_C, N_percent, NN_umolL

summary(modelHLC) # ***
summary(modelHMC) # p > 0.1
summary(modelLMC) # **





#### CLUSTERS FROM TURB ONLY ####

# # manova with all three clusters
# model2 <- manova(data = turbcluster,
#                  cbind(del15N, C_N, N_percent) ~ V_cluster)
# summary(model2) # ***
#
# # manova by cluster duo:
# clusterHLt <- turbcluster %>% filter(V_cluster %in% c("High", "Low"))
# clusterHMt <- turbcluster %>% filter(V_cluster %in% c("High", "Mid"))
# clusterLMt <- turbcluster %>% filter(V_cluster %in% c("Low", "Mid"))
#
# modelHLt <- manova(data = clusterHLt,
#                   cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables
#
# modelHMt <- manova(data = clusterHMt,
#                    cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables
#
# modelLMt <- manova(data = clusterLMt,
#                    cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables
#
#
# summary(modelHLt) # ***
# summary(modelHMt) # **
# summary(modelLMt) # ***


###################################
###### nMDS AND perMANOVA
###################################

# perMANOVA: is there a difference in response variables among groups?

# 3 steps in using classification functions:
# 1. Get the coefficients and constant of the classification equation for each group
# 2. Get a classification score for each observation for each group
#     a. Calculated by using the actual values for each variable to solve the classification for that group
# 3. Classify each observation into the group it matches most closely
#     a. This may or may not be the actual group from which the observation came
#     b. You have to know what group the observation goes into.
#        You have to know the classifications already, and then you can see how good the DFA is at fitting the data

#### PERMANOVA WITH ALL BIOGEOCHEMICAL DATA ####


# select numerical and single category data only
V_reduced <- cluster %>%
  filter(Location == "Varari") %>%
  select(pH, NN_umolL, Salinity, # only including discriminitive parameters from LCM
         Silicate_umolL, del15N,
         Tyrosine_Like, M_C, cluster)

C_reduced <- cluster %>%
  filter(Location == "Cabral") %>%
  select(Tyrosine_Like, VisibleHumidic_Like, # only including discriminitive parameters from LCM
         Phosphate_umolL, cluster)

# numerical data only
V_mds_data <- V_reduced %>%
  select(pH:M_C)

C_mds_data <- C_reduced %>%
  select(Tyrosine_Like:Phosphate_umolL)

#create the ordination output using bray curtis dissimilarity matrix
# numerical data only
Vord <- metaMDS(V_mds_data, k=2, distance='bray')
Cord <- metaMDS(C_mds_data, k=2, distance='bray')

#let's look at the stress with k=2 dimensions. Is it < 0.3?
Vord$stress
Cord$stress

# Let's look at the stress plot
# want to minimize scatter
stressplot(Vord)
stressplot(Cord)

# add back on the cluster IDs
V_mds_data <- V_mds_data %>% cbind(V_reduced$cluster) %>% rename(cluster = 'V_reduced$cluster')
C_mds_data <- C_mds_data %>% cbind(C_reduced$cluster) %>% rename(cluster = 'C_reduced$cluster')

# numeric data only
V_mds <- cluster %>%
  filter(Location == "Varari") %>%
  select(pH, NN_umolL, Salinity, # only including discriminitive parameters from LCM
         Silicate_umolL, del15N,
         Tyrosine_Like, M_C, Ammonia_umolL,
         VisibleHumidic_Like, cluster)
model2V <- adonis2(V_mds[,1:9]~cluster, V_mds, permutations = 999,
                       method="bray")
model2V

C_mds <- cluster %>%
  filter(Location == "Cabral") %>%
  select(Tyrosine_Like, VisibleHumidic_Like, # only including discriminitive parameters from LCM
         Phosphate_umolL, Tryptophan_Like,
         M_C, N_percent, NN_umolL, C_N, TA, cluster)
model2C <- adonis2(C_mds[,1:9]~cluster, C_mds, permutations = 999,
                 method="bray")
model2C

# If we are to trust the results of the permanova,
# then we have to assume that the dispersion among
# data is the same in each group. We can test with
# assumption with a PermDisp test:

# Look at the Average distance to median...these numbers should be
# reasonably similar
# A rule of thumb is that one number should *not* be twice as high as any other

disperV<-vegdist(V_mds[,1:9])
betadisper(disperV, V_mds$cluster)

disperC<-vegdist(C_mds[,1:9], na.rm=T)
betadisper(disperC, C_mds$cluster)

### Dispersion is all over the place, so cannot trust

pairwise.adonis(V_mds[,1:9], V_mds$cluster, perm=999)
pairwise.adonis(C_mds[,1:9], C_mds$cluster, perm=999)

#Get coefficients to see which species are most important in explaining site differences:
model2V$coefficients
model2C$coefficients





