### Plotting Species and Functional entity richness by distance to seep

### Created by Danielle Barnas
### Created on February 27, 2023


#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(ggrepel)
library(PNWColors)
library(vegan)
library(pairwiseAdonis)
library(patchwork)


#### READ IN DATA ####
meta <- read_csv(here("Data", "Full_Metadata.csv")) %>% # input seep rugosity
  mutate_at(vars(meanRugosity), .funs = ~if_else(is.na(meanRugosity)==TRUE, 0.975, .))
comp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
taxa <- read_csv(here("Data","Surveys","Distinct_Taxa.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrient_Processed_CV.csv"))



#### Calculate richness ####

sprich <- comp %>%
  filter(Taxa != "Sand" & Taxa != "Bare Rock" & Taxa != "Rubble" & Taxa != "Bare Rock - exposed") %>%
  group_by(Location, CowTagID) %>%
  count(Taxa) %>%
  mutate(n = 1,
         spr = sum(n)) %>%
  distinct(Location, CowTagID, spr) %>%
  ungroup()

ferich <- comp %>%
  left_join(taxa) %>%
  filter(Taxon_Group != "Hard Substrate" &
           Taxon_Group != "Sand") %>%
  select(Location, CowTagID, Taxa,
         Morphology,
         #Life_Span,
         #Max_size,
         #Growth_rate,
         Zooxanthellate,
         Calcification,
         Energetic_Resource,
         Feeding_Mode) %>%
  unite(Morphology,
        #Life_Span,
        #Max_size,
        #Growth_rate,
        Zooxanthellate,
        Calcification,
        Energetic_Resource,
        Feeding_Mode,
        sep = ", ", remove = T, col = "FE") %>%
  distinct(Location, CowTagID, FE) %>%
  group_by(Location, CowTagID) %>%
  count(FE) %>%
  mutate(n = 1,
         fer = sum(n)) %>%
  distinct(Location, CowTagID, fer)


## select only distance from seep from site data
meta <- meta %>%
  select(Location, CowTagID, dist_to_seep_m, meanRugosity)

Full_data <- full_join(sprich, ferich) %>%
  left_join(meta) %>%
  filter(Location == "Varari" & CowTagID != "V13")



#### CALCULATE RESIDUALS (richness ~ distance to seep) ####

#fit models
spResMod <- lm(spr ~ meanRugosity, data=Full_data %>% filter(CowTagID != "V13"))
feResMod <- lm(fer ~ meanRugosity, data=Full_data %>% filter(CowTagID != "V13"))

#view model summary
summary(spResMod)
summary(feResMod)

#calculate the standardized residuals
spres <- residuals(spResMod)
feres <- residuals(feResMod)

#view the standardized residuals
spres
feres

#column bind standardized residuals back to original data frame
res_data <- Full_data %>%
  filter(CowTagID != "V13") %>%
  cbind(spres, feres) %>%
  relocate(spres, .after = CowTagID) %>%
  relocate(feres, .after = spres)

#### PLOT RESIDUALS ####

#plot predictor variable vs. standardized residuals
p1 <- res_data %>%
  ggplot(aes(x = meanRugosity,
             y = spres)) +
  geom_point(color = if_else(spres > 0, "black", "red"),
             size = 3) +
  geom_text_repel(aes(label = CowTagID)) +
  theme_bw()
p2 <- res_data %>%
  ggplot(aes(x = meanRugosity,
             y = feres)) +
  geom_point(color = if_else(feres > 0, "black", "red"),
             size = 3) +
  geom_text_repel(aes(label = CowTagID)) +
  theme_bw()
p1+p2


### plot residuals against distance parameter

sprespval <- anova(lm(data = res_data, spres ~ dist_to_seep_m))[5]$'Pr(>F)'[1]
spresr2val <- summary(lm(data = res_data, spres ~ dist_to_seep_m))$'r.squared'

ferespval <- anova(lm(data = res_data, feres ~ dist_to_seep_m))[5]$'Pr(>F)'[1]
feresr2val <- summary(lm(data = res_data, feres ~ dist_to_seep_m))$'r.squared'

p3 <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = spres)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", fill = "darkseagreen") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(y = "Residuals (Species Richness ~ Mean Rugosity)",
       x = "Distance to seepage point (m)",
       caption = paste("p-value:",round(sprespval,3),'\n',
                       "r-squared:", round(spresr2val,3))) +
  geom_text_repel(aes(label = CowTagID), size = 3)
p4 <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = feres)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", fill = "darkseagreen") +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(y = "Residuals (FE Richness ~ Mean Rugosity)",
       x = "Distance to seepage point (m)",
       caption = paste("p-value:",round(ferespval,3),'\n',
                       "r-squared:", round(feresr2val,3))) +
  geom_text_repel(aes(label = CowTagID), size = 3)
p34 <- p3+p4
p34

ggsave(here("Output", "PaperFigures","plot_rich_residuals_dist.png"), p34, width = 7, height = 6)




