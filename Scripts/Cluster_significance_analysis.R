#### Cluster significant difference analysis ####
### Created by Danielle Barnas
### Created on Aug 29, 2022

#### Load Libraries ####
library(tidyverse)
library(here)
library(stats) # lm() and anova()
library(emmeans)
library(agricolae) # HSD.test()

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
cluster12 <- cluster %>% filter(V_cluster %in% c(1, 2))
cluster13 <- cluster %>% filter(V_cluster %in% c(1, 3))
cluster23 <- cluster %>% filter(V_cluster %in% c(2, 3))

model12 <- manova(data = cluster12,
                         cbind(TA, Phosphate_umolL, # bind response variables
                               Silicate_umolL, NN_umolL, del15N,
                               C_N, N_percent) ~ V_cluster) # cannot run model when response variables outnumber sample size
# removed pH and Salinity (lowest discrimitive power)

model13 <- manova(data = cluster13,
                  cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                        Silicate_umolL, NN_umolL, del15N,
                        C_N, N_percent) ~ V_cluster)

model23 <- manova(data = cluster23,
                  cbind(Salinity, TA, pH, Phosphate_umolL, # bind response variables
                        Silicate_umolL, NN_umolL, del15N,
                        C_N, N_percent) ~ V_cluster)

summary(model12) # p > 0.06
summary(model13) # p > 0.08
summary(model23) # **



#### CLUSTERS FROM TURB ONLY ####

# manova with all three clusters
model2 <- manova(data = turbcluster,
                 cbind(del15N, C_N, N_percent) ~ V_cluster)
summary(model2) # ***

# manova by cluster duo:
cluster12t <- turbcluster %>% filter(V_cluster %in% c(1, 2))
cluster13t <- turbcluster %>% filter(V_cluster %in% c(1, 3))
cluster23t <- turbcluster %>% filter(V_cluster %in% c(2, 3))

model12t <- manova(data = cluster12t,
                  cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables

model13t <- manova(data = cluster13t,
                   cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables

model23t <- manova(data = cluster23t,
                   cbind(del15N, C_N, N_percent) ~ V_cluster) # bind response variables


summary(model12t) # ***
summary(model13t) # **
summary(model23t) # ***

