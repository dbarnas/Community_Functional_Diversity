#### Linear Regressions of Fric and Fric residuals
### species richness, functional entity richness, and functional volume

### Created by Danielle Barnas
### Created on February 26, 2023


library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)


Fric <- read_csv(here("Data", "Sp_FE_Vol_subset.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv")) %>%
  select(Location, CowTagID, Salinity:Tyrosine_Like)
shore <- read_csv(here("Data", "Shore_distance.csv"))
meta <- meta %>%
  left_join(chem)


### Join Sp and FE and Vol4D with metadata
reg.Fric <- Fric %>%
  as_tibble() %>%
  left_join(meta) %>%
  mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity)) %>%
  left_join(shore)




### Calculate regressions again with Residuals (remove structure/substrate)

#fit model: linear relationship
resModSpR <- lm(NbSp ~ shore_dist_m, data=reg.Fric) # species richness
resModFER <- lm(NbFEs ~ shore_dist_m, data=reg.Fric) # entity richness
resModSpRp <- lm(NbSpP ~ shore_dist_m, data=reg.Fric) # relative species richness
resModFERp <- lm(NbFEsP ~ shore_dist_m, data=reg.Fric) # relative entity richness

#check other parameters if no relationship to rugosity
colnames(reg.Fric)
#resModSpR <- lm(NbSp ~  shore_dist_m, data=reg.Fric)
# > 10% cover * C_N, dist_to_seep_m, del15N, Salinity, HIX, N_percent
# no relationship: Silicate, NN_umolL, Ammonia_umolL, M_C, VisibleHumidic_Like, Tryptophan_Like
# > 5% cover * dist_to_seep_m, Salinity, TA, pH, Phosphate_umolL, shore_dist_m
# no relationship: del15N, C_N, N_percent, meanRugosity (0.053), Temp, Silicate (0.93), NN (0.3), Ammonia (0.98),
# M_C, HIX, VisibleHumidic_Like, Tryptophan_Like, Tyrosine_Like
#
# ggplot(data = reg.Fric,
#        aes(x = dist_to_seep_m, y = NbSp)) +
#   geom_point() +
#   geom_smooth(method = "lm")

#view model summary
summary(resModSpR) #
summary(resModFER) #

#calculate the standardized residuals
resSpR <- residuals(resModSpR)
resFER <- residuals(resModFER)
resSpRp <- residuals(resModSpRp)
resFERp <- residuals(resModFERp)

#column bind standardized residuals back to original data frame
res_data <- reg.Fric %>%
  cbind(resSpR) %>%
  cbind(resFER) %>%
  cbind(resSpRp) %>%
  cbind(resFERp)



# plot predictor variable vs. standardized residuals

# isolate seep poing
ptseep <- res_data %>%
  filter(CowTagID == "VSEEP")

## Raw richness residuals
rug_res_SpR_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpR)) +
  geom_point(size = 2.5) +
  geom_point(data = ptseep,
             aes(x = dist_to_seep_m,
                 NbSp),
                 #y = resSpR),
             shape = 23,
             size = 5,
             color = "black",
             fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Dist offshore-normalized)",
       x = "Distance to seep (m)")
rug_res_SpR_plot

rug_res_FER_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFER)) +
  geom_point(size = 2.5) +
  geom_point(data = ptseep,
             aes(x = dist_to_seep_m,
                 y = NbFEs),
                 #y = resSpR),
             shape = 23,
             size = 5,
             color = "black",
             fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FER (Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_FER_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpRp)) +
  geom_point(size = 2.5) +
  geom_point(data = ptseep,
             aes(x = dist_to_seep_m,
                 y = resSpR),
             shape = 23,
             size = 5,
             color = "black",
             fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_SpRp_plot
summary(lm(data = res_data, resSpRp ~ Silicate_umolL))


rug_res_FERp_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFERp)) +
  geom_point(size = 2.5) +
  geom_point(data = ptseep,
             aes(x = dist_to_seep_m,
                 y = resSpR),
             shape = 23,
             size = 5,
             color = "black",
             fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_FERp_plot


### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1
ggsave(here("Output", "PaperFigures", "Subset", "lm_rich_residuals_dist_subset.png"), plot1, height = 6, width = 10)


summary(lm(data = res_data, resSpR ~ dist_to_seep_m))
summary(lm(data = res_data, resSpR ~ Ammonia_umolL))


summary(lm(data = res_data %>% filter(CowTagID != "VSEEP"), resSpR ~ dist_to_seep_m))
summary(lm(data = res_data %>% filter(CowTagID != "VSEEP"), resFER ~ dist_to_seep_m))
summary(lm(data = res_data %>% filter(CowTagID != "VSEEP"), resSpRp ~ dist_to_seep_m))
summary(lm(data = res_data %>% filter(CowTagID != "VSEEP"), resFERp ~ dist_to_seep_m))





res_data %>%
  ggplot(aes(x = dist_to_seep_m, y = resFERp/resSpRp)) +
  geom_point(aes(color = dist_to_seep_m)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  geom_label(aes(label = CowTagID))
summary(lm(data = res_data, (NbFEs/NbSp) ~ dist_to_seep_m))
summary(lm(data = res_data, (resFERp/resSpRp) ~ dist_to_seep_m))


# ratio of FE:Sp richness ~ distance to seep
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity


## Can use the three values above, and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness


