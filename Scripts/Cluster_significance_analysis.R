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
cluster <- read_csv(here("Data", "Surveys", "Cluster_metadata.csv"))
turbcluster <- read_csv(here("Data", "Surveys", "Cluster_turb_metadata.csv"))

#### CLUSTERS FROM FULL BIOGEOCHEM ####

### parse V_clusters to factor
cluster <- cluster %>%
  mutate(V_cluster = as.factor(V_cluster))
turbcluster <- turbcluster %>%
  rename(V_cluster = V_cluster_turb) %>%
  mutate(V_cluster = as.factor(V_cluster))


# manova with all three clusters
model1 <- manova(data = cluster,
              cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                    Silicate_umolL, NN_umolL, del15N,
                    C_N, N_percent) ~ V_cluster)
summary(model1) # **


# manova by cluster duo:
clusterHL <- cluster %>% filter(V_cluster %in% c("High", "Low"))
clusterHM <- cluster %>% filter(V_cluster %in% c("High", "Mid"))
clusterLM <- cluster %>% filter(V_cluster %in% c("Low", "Mid"))

modelHL <- manova(data = clusterHL,
                         cbind(TA, Phosphate_umolL, # bind response variables
                               Silicate_umolL, NN_umolL, del15N,
                               C_N, N_percent) ~ V_cluster) # cannot run model when response variables outnumber sample size
# removed pH and Salinity (lowest discrimitive power)

modelHM <- manova(data = clusterHM,
                  cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                        Silicate_umolL, NN_umolL, del15N,
                        C_N, N_percent) ~ V_cluster)

modelLM <- manova(data = clusterLM,
                  cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                        Silicate_umolL, NN_umolL, del15N,
                        C_N, N_percent) ~ V_cluster)

summary(modelHL) # p > 0.06
summary(modelHM) # p > 0.08
summary(modelLM) # **



#### CLUSTERS FROM TURB ONLY ####

# manova with all three clusters
model2 <- manova(data = turbcluster,
                 cbind(del15N, C_N, N_percent) ~ V_cluster)
summary(model2) # ***

# manova by cluster duo:
clusterHLt <- turbcluster %>% filter(V_cluster %in% c("High", "Low"))
clusterHMt <- turbcluster %>% filter(V_cluster %in% c("High", "Mid"))
clusterLMt <- turbcluster %>% filter(V_cluster %in% c("Low", "Mid"))

modelHLt <- manova(data = clusterHLt,
                  cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables

modelHMt <- manova(data = clusterHMt,
                   cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables

modelLMt <- manova(data = clusterLMt,
                   cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables


summary(modelHLt) # ***
summary(modelHMt) # **
summary(modelLMt) # ***


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
cluster_reduced <- cluster %>%
  select(-c(Location, CowTagID, lat, lon)) %>%
  select(Salinity, TA, pH, Phosphate_umolL, # only including discriminitive parameters from LCM
         Silicate_umolL, NN_umolL, del15N,
         C_N, N_percent, V_cluster)

# numerical data only
mds_data <- cluster_reduced %>%
  select(Salinity:N_percent)

#create the ordination output using bray curtis dissimilarity matrix
# numerical data only
ord<-metaMDS(mds_data,k=2, distance='bray')

#let's look at the stress with k=2 dimensions. Is it < 0.3?
ord$stress

# Let's look at the stress plot
# want to minimize scatter
stressplot(ord)

# basic plot
# dots represent sites (sandwich locations) and + represents parameters
ordiplot(ord)

# add text
ordiplot(ord, type = 'text')

# add back on the cluster IDs
mds_data <- cluster_reduced

# numeric data only
model3 <- adonis(mds_data[,1:9]~V_cluster, mds_data, permutations = 999,
                       method="bray")
model3

# If we are to trust the results of the permanova,
# then we have to assume that the dispersion among
# data is the same in each group. We can test with
# assumption with a PermDisp test:

# Look at the Average distance to median...these numbers should be
# reasonably similar
# A rule of thumb is that one number should *not* be twice as high as any other

disper<-vegdist(mds_data[,1:9])
betadisper(disper, mds_data$V_cluster)

### Dispersion is all over the place, so cannot trust

pairwise.adonis(mds_data[,1:9], mds_data$V_cluster, perm=999)

#Get coefficients to see which species are most important in explaining site differences:
model3$coefficients





#### PERMANOVA WITH ALL TURBINARIA DATA ####

# select numerical and single category data only
turbcluster_reduced <- turbcluster %>%
  select(del15N, C_N, N_percent, V_cluster)

# numerical data only
mds_turbdata <- turbcluster_reduced %>%
  select(del15N:N_percent)

#create the ordination output using bray curtis dissimilarity matrix
# numerical data only
ord<-metaMDS(mds_turbdata,k=2, distance='bray')

#let's look at the stress with k=2 dimensions. Is it < 0.3?
ord$stress

# Let's look at the stress plot
# want to minimize scatter
stressplot(ord)

# basic plot
# dots represent sites (sandwich locations) and + represents parameters
ordiplot(ord)

# add text
ordiplot(ord, type = 'text')

# add back on the cluster IDs
mds_turbdata <- turbcluster_reduced

# numeric data only
model4 <- adonis(mds_turbdata[,1:3]~V_cluster, mds_turbdata, permutations = 999,
                 method="bray")
model4

# If we are to trust the results of the permanova,
# then we have to assume that the dispersion among
# data is the same in each group. We can test with
# assumption with a PermDisp test:

# Look at the Average distance to median...these numbers should be
# reasonably similar
# A rule of thumb is that one number should *not* be twice as high as any other

disper<-vegdist(mds_data[,1:3])
betadisper(disper, mds_turbdata$V_cluster)

### Dispersion is all over the place, so cannot trust

# post-hoc
pairwise.adonis(mds_turbdata[,1:3], mds_turbdata$V_cluster, perm=999)

#Get coefficients to see which parameters are most important in explaining site differences:
model4$coefficients



