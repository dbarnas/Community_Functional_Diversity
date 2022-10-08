### Broad functional categorizations across SGD

#### LIBRARIES ####
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)


#### READ IN DATA ####
taxa <- read_csv(here("Data", "Surveys","Distinct_Varari_Taxa.csv"))
survey <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv")) %>%
  filter(Location == "Varari")
nutrient <- read_csv(here("Data", "Biogeochem", "AugNutrient_Processed_CV.csv"))


#### CLEAN DATA ####

# group macroalgae
funTrait <- taxa %>%
  mutate(coral_algae = if_else(Taxon_Group == "Rhodophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Chlorophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Phaeophyta", true = "Macroalgae", false = Taxon_Group))))

# calculate percent cover of primary taxon groups
pCover <- funTrait %>%
  full_join(survey) %>%
  group_by(CowTagID) %>%
  count(coral_algae) %>%
  pivot_wider(names_from = coral_algae, values_from = n) %>%
  mutate(total = sum(across(Abiotic:Cyanobacteria), na.rm = T)) %>% # calculate total individuals per category
  mutate(across(Abiotic:Cyanobacteria, ~ . / total * 100), .keep = c("all")) %>%  # calculate % cover for benthic categories per cowtagid
  select(-total) %>%   # remove intermediate metric
  pivot_longer(Abiotic:Cyanobacteria, names_to = "Taxa", values_to = "values") %>% # pivot to replace na with zero
  mutate(values = if_else(is.na(values), 0, values)) %>%
  pivot_wider(names_from = Taxa, values_from = values) # pivot back to wide

# binning calcifiers and fleshy+turf algae
pCover <- pCover %>%
  mutate(TotalCalc = CCA + Coral,
         TotalAlgae = Macroalgae + Turf)

# join with biogeochemical data
pCover <- pCover %>%
  filter(CowTagID != "VSEEP") %>%
  left_join(nutrient)


#### PLOTTING ####

# Macroalgae vs Coral

#lm(data = pCover, TotalAlgae ~ TotalCalc) # observe intercept
pCover %>%
  ggplot(aes(x = TotalCalc, y = TotalAlgae)) +
  geom_point() +
  #geom_abline(slope = -1, intercept = 66.86) +
  geom_text(aes(label = ifelse(TotalCalc < 5, CowTagID, NA)),
            nudge_y = -3, nudge_x = 2, size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #scale_y_continuous(limits = c(0,70))


#### Comparing coral/algae ratios across sgd ####

# del15N
p1 <- pCover %>%
  ggplot(aes(y = log((TotalCalc+1)/(Macroalgae+1)), x = del15N)) +
  geom_point(aes(size = N_percent)) +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)),
  #           nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, lty = 2)+
  scale_size_binned("%N (Nutrient Loading)") +
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "del15N (Nutrient Source)")+
  annotate("text", x = 5, y = 0.5, label = "Calcifier-dominated")+
  annotate("text", x = 5, y = -0.7, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.direction = "horizontal",
        legend.position = c(.18, .1))
ggsave(here("Output","CoralAlgae_15N.pdf"), p1, width = 10, height = 10)

# percent N
p2 <- pCover %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = N_percent))+
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)),
  #           nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_hline(yintercept = 0, lty = 2) + #lty: line type, 2: dotted line
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "%N (Nutrient Loading)") +
  # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
  # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none")
ggsave(here("Output","CoralAlgae_Npercent.pdf"), p2, width = 10, height = 10)

# C:N ratio
p3 <- pCover %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = C_N))+
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)),
  #           nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_hline(yintercept = 0, lty = 2) + #lty: line type, 2: dotted line
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "C:N Ratio (Nutrient Loading)") +
  # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
  # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none")
ggsave(here("Output","CoralAlgae_CN.pdf"), p3, width = 10, height = 10)

modBenthic15N<-lm(log((TotalCalc+1)/(TotalAlgae+1))~del15N, data = pCover)
anova(modBenthic15N)

modBenthicNpercent<-lm(log((TotalCalc+1)/(TotalAlgae+1))~N_percent, data = pCover)
anova(modBenthicNpercent)

p1+p2 / p3


## without log transformation
pCover %>%
  ggplot(aes(y = (TotalCalc/TotalAlgae), x = C_N))+
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)),
  #           nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_hline(yintercept = 1, lty = 2) + #lty: line type, 2: dotted line
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "ratio of calcifiers to fleshy algae",
       color = "",
       x = "C:N Ratio (Nutrient Loading)") +
  # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
  # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none")






