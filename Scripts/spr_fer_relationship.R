
library(tidyverse)
library(here)

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv"))) %>% rename(FE = fun_entity)
meta <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(Location == "Varari",
         CowTagID != "V13",
         CowTagID != "VSEEP") %>%
  select(CowTagID, dist_to_seep_m, meanRugosity) %>%
  arrange(dist_to_seep_m)
dist <- meta[1:19,]

sp<-ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  distinct(CowTagID, Taxa) %>%
  group_by(CowTagID) %>%
  count(Taxa) %>%
  summarise(spr = sum(n))

full<-ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(spe_fes.sgd) %>%
  distinct(CowTagID, FE) %>%
  group_by(CowTagID) %>%
  count(FE) %>%
  summarise(fer = sum(n)) %>%
  left_join(sp) %>%
  left_join(dist)
full$CowTagID <- factor(full$CowTagID, levels = dist$CowTagID)
full %>%
  ggplot(aes(x = dist_to_seep_m, y = (fer/spr))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()+
  theme(panel.grid = element_blank())

mod <- lm(spr ~ meanRugosity, data = full)
summary(mod)
spres <- residuals(mod)

mod2 <- lm(fer ~ meanRugosity, data = full)
summary(mod2)
feres <- residuals(mod2)

full <- full %>%
  cbind(spres, feres)

full %>%
  ggplot(aes(x = dist_to_seep_m, y = (feres/spres))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()+
  theme(panel.grid = element_blank())

summary(lm(spr ~ dist_to_seep_m, data = full))
summary(lm(fer ~ dist_to_seep_m, data = full))
summary(lm(spres ~ dist_to_seep_m, data = full))
summary(lm(feres ~ dist_to_seep_m, data = full))


