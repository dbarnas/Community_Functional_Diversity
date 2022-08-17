###########################################################
### GROWTH
###########################################################

library(tidyverse)
library(here)
library(stats)


### READ IN DATA
bw <- read_csv(here("Data","Growth","buoyant_weights.csv"))
ww <- read_csv(here("Data", "Growth", "wet_weights.csv"))
disp <- read_csv(here("Data", "Growth", "water_displacement.csv"))

############################################################
### BUOYANT WEIGHTS
############################################################

bw_red <- bw %>%
  filter(Before_After == 'yes') %>% # only keep organisms with before and after pairings in high and low environment
  select(Individual_ID, Weight_g, pre_post, Species_ID) %>% # remove unecessary columns
  pivot_wider(names_from = 'pre_post', values_from = 'Weight_g') %>% # wide format for pre and post weights
  mutate(growth_g = post - pre) %>% # calculate growth
  arrange(growth_g) %>%
  filter(Species_ID != "Montipora efflorescans") %>%  # too low sample size; remove
  filter(Species_ID != 'Lithophyllum kotchyanum') # just for now, logistics of H species
# separate sample ID from treatment ID with '_'
bw_red$Individual_ID <- str_replace_all(string = bw_red$Individual_ID, pattern = "L", replacement = "_L")
bw_red$Individual_ID <- str_replace_all(string = bw_red$Individual_ID, pattern = "H", replacement = "_H")
bw_red <- bw_red %>%
  separate(col = Individual_ID, into = c('Individual_ID', 'envi_treatment'), sep = "_", remove = T)

# remove any coral with weight loss (likely due to breakage)
broken <- bw_red %>%
  filter(growth_g < 0) %>%
  select(Individual_ID) # isolate org ID's to remove their pairnings from larger df
bw_red <- bw_red %>%
  anti_join(broken)


# visualize high treatment vs low treatment growth (not yet normalized) by species group
bw_red %>%
  ggplot(aes(x = Individual_ID, y = growth_g, fill = envi_treatment)) +
  geom_col(position = "dodge") +
  facet_wrap(~Species_ID, scales = "free_x")

PAgrowth <- bw_red %>%
  filter(Species_ID == "Pocillopora acuta") %>%
  select(Individual_ID, growth_g, envi_treatment)
PAtest<-t.test(growth_g~envi_treatment, paired=TRUE, data=prusgrowth)
PAtest

# summarise data
bw_red %>%
  group_by(Species_ID, envi_treatment) %>%
  summarise(mean = mean(growth_g),
            sd = sd(growth_g),
            se = sd/sqrt(length(growth_g))) %>%
  ggplot(aes(x = Species_ID, y = mean, fill = envi_treatment)) +
  geom_col(position = "dodge")




############################################################
### WET WEIGHTS
############################################################

# remove organisms with breaks and/or incomplete data sets for before/after
# broken <- ww %>%
#   filter(Individual_ID %in% c("D6_H", "P9_L", "I13_L", "D9_H",
#                               "C11_H", "C8_H", "D3_H", "D4_H",
#                               "D7_L", "P10_H"))

ww_red <- ww %>%
  filter(Before_After == 'yes') %>%
  select(Individual_ID, pre_post, Weight_g, Species_ID) %>%
  pivot_wider(names_from = 'pre_post', values_from = 'Weight_g') %>%
  mutate(growth_g = post - pre) %>%
  separate(Individual_ID, sep = "_", into = c('Individual_ID', 'envi_treatment'), remove = T)

# remove organisms with breaks and/or incomplete data sets for before/after
broken <- ww_red %>%
  filter(Individual_ID %in% c("D6", "P9", "I13", "D9", "C11", "C8", "D3", "D4", "D7", "P10")) %>%
  select(Individual_ID)
ww_red <- ww_red %>%
  anti_join(broken)


# visualize high treatment vs low treatment growth (not yet normalized) by species group
ww_red %>%
  ggplot(aes(x = Individual_ID, y = growth_g, fill = envi_treatment)) +
  geom_col(position = "dodge") +
  facet_wrap(~Species_ID, scales = "free_x")


GSgrowth <- ww_red %>%
  filter(Species_ID == "Grey Sponge") %>%
  select(Individual_ID, growth_g, envi_treatment)
GStest<-t.test(growth_g~envi_treatment, paired=TRUE, data=GSgrowth)
GStest

LVgrowth <- ww_red %>%
  filter(Species_ID == "Lobophora variegata") %>%
  select(Individual_ID, growth_g, envi_treatment)
LVtest<-t.test(growth_g~envi_treatment, paired=TRUE, data=LVgrowth)
LVtest

# summarise data
ww_red %>%
  group_by(Species_ID, envi_treatment) %>%
  summarise(mean = mean(growth_g),
            sd = sd(growth_g),
            se = sd/sqrt(length(growth_g))) %>%
  ggplot(aes(x = Species_ID, y = mean, fill = envi_treatment)) +
  geom_col(position = "dodge")




############################################################
### DISPLACEMENTS
############################################################


disp_red <- disp %>%
  drop_na(c('Vol_pre_0', 'Vol_pre_1')) %>% # proxy for having pre and post measurements
  select(-c(Notes)) %>%
  mutate(Disp_0 = Vol_post_0 - Vol_pre_0,
         Disp_1 = Vol_post_1 - Vol_pre_1,
         growth_mL = Disp_1 - Disp_0) %>%
  separate(Individual_ID, sep = "_", into = c('Individual_ID', 'envi_treatment'), remove = T)

# remove organisms with breaks and/or incomplete data sets for before/after
broken <- disp_red %>%
  filter(Individual_ID %in% c("D6", "P9", "I13", "D9", "C11", "C8", "D3", "D4", "D7", "P10")) %>%
  select(Individual_ID)
disp_red <- disp_red %>%
  anti_join(broken)


# visualize high treatment vs low treatment growth (not yet normalized) by species group
disp_red %>%
  ggplot(aes(x = Individual_ID, y = growth_g, fill = envi_treatment)) +
  geom_col(position = "dodge") +
  facet_wrap(~Species_ID, scales = "free_x")


GSgrowth <- disp_red %>%
  filter(Species_ID == "Grey Sponge") %>%
  select(Individual_ID, growth_g, envi_treatment)
GStest<-t.test(growth_g~envi_treatment, paired=TRUE, data=GSgrowth)
GStest

LVgrowth <- disp_red %>%
  filter(Species_ID == "Lobophora variegata") %>%
  select(Individual_ID, growth_g, envi_treatment)
LVtest<-t.test(growth_g~envi_treatment, paired=TRUE, data=LVgrowth)
LVtest

# summarise data
disp_red %>%
  group_by(Species_ID, envi_treatment) %>%
  summarise(mean = mean(growth_g),
            sd = sd(growth_g),
            se = sd/sqrt(length(growth_g))) %>%
  ggplot(aes(x = Species_ID, y = mean, fill = envi_treatment)) +
  geom_col(position = "dodge")




