####  Plot relative abundance of individual functional groups across sites


library(tidyverse)
library(here)
library(PNWColors)
library(patchwork)

mypalette <- pnw_palette(name = "Bay", n = 12)
smallpalette <- pnw_palette(name = "Bay", n = 4)

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv"))) %>% rename(FE = fun_entity)
dist <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(Location == "Varari",
         CowTagID != "V13",
         CowTagID != "VSEEP") %>%
  select(CowTagID, dist_to_seep_m) %>%
  arrange(dist_to_seep_m)
dist <- dist[1:19,]


Full_data <- ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(spe_fes.sgd) %>%
  left_join(fes_traits.sgd)

ctlevels <- dist$CowTagID
Full_data$CowTagID <- factor(Full_data$CowTagID, levels = ctlevels)


myplot <- function(param){

  myfill = enquo(param)

  plota <- Full_data %>%
    ggplot(aes(x = CowTagID,
               y = pCover,
               fill = !!myfill)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top") +
    scale_fill_manual(values = mypalette)

  return(plota)

}

mysmallplot <- function(param){

  myfill = enquo(param)

  plota <- Full_data %>%
    ggplot(aes(x = CowTagID,
               y = pCover,
               fill = !!myfill)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values = smallpalette)

  return(plota)

}

pm <- myplot(Morph)
#pm2 <- myplot(Morph2)
pc <- mysmallplot(Calc)
pz <- mysmallplot(Zoox)
per <- mysmallplot(ER)
pfm <- mysmallplot(FM)

plot1 <- (pm) /
(pc + pz) /
(per + pfm)
plot1
ggsave(here("Output", "PaperFigures", "Plot_FEgroups_dist.png"), plot1, width = 6, height = 8)


FE_data <- Full_data %>%
  group_by(CowTagID, FE, Morph, Calc, Zoox, ER, FM) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup() %>%
  left_join(dist)

morphd <- FE_data %>%
  select(dist_to_seep_m, Morph, pCover) %>%
  group_by(dist_to_seep_m, Morph) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(aov(pCover ~ dist_to_seep_m*Morph, data = morphd))

calcd <- FE_data %>%
  select(dist_to_seep_m, Calc, pCover) %>%
  group_by(dist_to_seep_m, Calc) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(aov(pCover ~ dist_to_seep_m*Calc, data = calcd))

zooxd <- FE_data %>%
  select(dist_to_seep_m, Zoox, pCover) %>%
  group_by(dist_to_seep_m, Zoox) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(aov(pCover ~ dist_to_seep_m*Zoox, data = zooxd))

erd <- FE_data %>%
  select(dist_to_seep_m, ER, pCover) %>%
  group_by(dist_to_seep_m, ER) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(aov(pCover ~ dist_to_seep_m*ER, data = erd))

fmd <- FE_data %>%
  select(dist_to_seep_m, FM, pCover) %>%
  group_by(dist_to_seep_m, FM) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(aov(pCover ~ dist_to_seep_m*FM, data = fmd))


anova(aov(pCover ~ dist_to_seep_m*FE, data = Full_data %>% left_join(dist)))
Full_data %>%
  left_join(dist) %>%
  group_by(dist_to_seep_m, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  ggplot(aes(x = dist_to_seep_m, y = pCover, color = FE)) +
  geom_point()+
  geom_smooth() +
  theme_bw() +
  facet_wrap(~FE)


Full_data %>%
  filter(CowTagID == "V14" |
           CowTagID == "V20" |
           CowTagID == "V17") %>%
  group_by(CowTagID, Calc) %>%
  summarise(pCover = sum(pCover))

Full_data %>%
  filter(Calc == "Non-AC") %>%
  distinct(Taxa)
Full_data %>%
  filter(Taxa == "Galaxaura rugosa") %>%
  select(CowTagID, pCover)
Full_data %>%
  filter(CowTagID == "V20") %>%
  select(Taxa,pCover) %>%
  arrange(desc(pCover))
