####  Plot relative abundance of individual functional groups across sites

### LOAD LIBRARIES ###

library(tidyverse)
library(here)
library(RColorBrewer)
library(patchwork)
library(PNWColors)


### READ IN DATA ###

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))
dist <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(Location == "Varari",
         CowTagID != "V13",
         #CowTagID != "VSEEP"
  ) %>%
  select(CowTagID, dist_to_seep_m) %>%
  arrange(dist_to_seep_m)
dist <- dist[1:20,]
chem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Location == "Varari",
         Season == "Dry",
         Parameters == "Phosphate_umolL" | Parameters == "NN_umolL") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)


### CREATE PALETTES FOR FIGURES ###

taxonpalette <- c("#0f4a6f", "#6aa4b0", "#29a0b1",
                  # "#98d7c2",
                  "#738fa7", "#2e8bc0", "#7391c8", "#52688f")
morphpalette <- rev(pnw_palette(name = "Moth", n = 11))
calcpalette <- c("#f9eac2","#b2d2a4", "#96ad90", "#1a4314")
#symbpalette <- c("#98d7c2", "#167d7f", "#29a0b1", "#05445e")
erpalette <- c("#fbc490", "#fbaa60", "#f67b50")
#fmpalette <- c( "#deb3ad","#de847b", "#c85250")


### JOIN DATA AND CREATE PLOT FUNCTION ###

Full_data <- ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(spe_fes.sgd) %>%
  left_join(fes_traits.sgd) %>%
  left_join(alphatag) %>%
  left_join(chem)

Full_data$AlphaTag <- factor(Full_data$AlphaTag)
Full_data$Morph2 <- factor(Full_data$Morph2,
                           levels = c('Br', 'Dig', 'Fol', 'Fil', 'Stol',
                                      'Mush', 'Poly', 'Cushion', 'Mas', 'Enc', 'Sph'))


myplot <- function(param, pal){

  my_data <- Full_data %>%
    group_by(AlphaTag,{{param}}) %>%
    summarise(pCover = sum(pCover))

  myfill = enquo(param)

  plota <- my_data %>%
    ggplot(aes(x = AlphaTag,
               y = pCover,
               fill = !!myfill)) +
    geom_col(position = "stack",
             color = "white") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top") +
    scale_fill_manual(values = pal) +
    labs(x = "", y = "% Cover")

  return(plota)

}

ptplot <- function(entity, param){

  my_data <- Full_data %>%
    filter(CowTagID != "VSEEP") %>%
    group_by(AlphaTag,{{entity}}, {{param}}) %>%
    summarise(pCover = sum(pCover))

  myfacet = enquo(entity)
  independent = enquo(param)

  plota <- my_data %>%
    ggplot(aes(x = !!independent,
               y = pCover,
               fill = !!myfacet
               )) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top") +
    #scale_fill_manual(aes(values = !!myfacet)) +
    facet_wrap(myfacet, scales = "free") +
    labs(x = "", y = "% Cover")
    #ylim(min = 0, max = 100)

  return(plota)

}

### PLOT FUNCTIONS ACROSS SEEP ###
pt <- myplot(Taxon_Group, taxonpalette)
pt2 <- myplot(Taxon_Group, taxonpalette) + theme(legend.position = "none")
ptplot(Taxon_Group, Phosphate_umolL)
ptplot(Taxon_Group, NN_umolL)

pm <- myplot(Morph2, morphpalette)
pm2 <- myplot(Morph2, morphpalette) + theme(legend.position = "right")
pm3 <- myplot(Morph2, morphpalette) + theme(legend.position = "none")
ptplot(Morph2, Phosphate_umolL)
ptplot(Morph2, NN_umolL)

pc <- myplot(Calc, calcpalette)
ptplot(Calc, Phosphate_umolL)
ptplot(Calc, NN_umolL)

per <- myplot(ER, erpalette)
ptplot(ER, Phosphate_umolL)
ptplot(ER, NN_umolL)

plot1 <- (pm) /
  (per + pc) /
  (pt) +
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = 10))
plot1


 ggsave(here("Output", "PaperFigures", "Plot_FEgroups_dist.png"), plot1, width = 6, height = 10)
 ggsave(here("Output", "PaperFigures", "Plot_FEgroups_dist2.png"), plot1, width = 6, height = 6)

 ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_Taxon_dist.png"), pt2, width = 6, height = 3)
 ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm2, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_Morph_dist.png"), pm3, width = 6, height = 3)
 ggsave(here("Output", "PaperFigures", "Plot_Calc_dist.png"), pc, width = 6, height = 3.5)
 ggsave(here("Output", "PaperFigures", "Plot_ER_dist.png"), per, width = 6, height = 3.5)



 ### DISTANCE ###


 FE_data <- Full_data %>%
   group_by(CowTagID, AlphaTag, FE, Taxon_Group, Morph2, Calc, ER) %>%
   summarise(pCover = sum(pCover)) %>%
   ungroup() %>%
   left_join(dist) %>%
   left_join(chem)

 morphd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Morph2, pCover) %>%
   group_by(dist_to_seep_m, Morph2) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Morph2, data = morphd))

 taxond <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Taxon_Group, pCover) %>%
   group_by(dist_to_seep_m, Taxon_Group) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Taxon_Group, data = taxond))

 calcd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Calc, pCover) %>%
   group_by(dist_to_seep_m, Calc) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*Calc, data = calcd))


 erd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, ER, pCover) %>%
   group_by(dist_to_seep_m, ER) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ dist_to_seep_m*ER, data = erd))


 anova(lm(pCover ~ dist_to_seep_m*FE, data = Full_data %>% left_join(dist)))

 Full_data %>%
   left_join(dist) %>%
   left_join(chem) %>%
   group_by(dist_to_seep_m, FE) %>%
   summarise(pCover = sum(pCover)) %>%
   ggplot(aes(x = dist_to_seep_m, y = pCover, color = FE)) +
   geom_point()+
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
   theme_bw() +
   facet_wrap(~FE, scales = "free_y") +
   theme(legend.position = "none")



 ### PHOSPHATE ###



 FE_data <- Full_data %>%
   group_by(CowTagID, AlphaTag, FE, Taxon_Group, Morph2, Calc, ER) %>%
   summarise(pCover = sum(pCover)) %>%
   ungroup() %>%
   left_join(dist) %>%
   left_join(chem)

 morphd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Morph2, pCover) %>%
   group_by(Phosphate_umolL, Morph2) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Morph2, data = morphd))

 taxond <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Taxon_Group, pCover) %>%
   group_by(Phosphate_umolL, Taxon_Group) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Taxon_Group, data = taxond)) # taxa NS ~ phosphate

 calcd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Calc, pCover) %>%
   group_by(Phosphate_umolL, Calc) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*Calc, data = calcd))


 erd <- FE_data %>%
   select(dist_to_seep_m, NN_umolL, Phosphate_umolL, ER, pCover) %>%
   group_by(Phosphate_umolL, ER) %>%
   mutate(pCover = sum(pCover)) %>%
   distinct()
 anova(lm(pCover ~ Phosphate_umolL*ER, data = erd)) # er NS ~ phosphate


 anova(lm(pCover ~ Phosphate_umolL*FE, data = Full_data %>% left_join(chem)))


 Full_data %>%
   left_join(dist) %>%
   left_join(chem) %>%
   group_by(Phosphate_umolL, FE) %>%
   summarise(pCover = sum(pCover)) %>%
   ggplot(aes(x = Phosphate_umolL, y = pCover, color = FE)) +
   geom_point()+
   geom_smooth(method = "lm", formula = "y~x", color = "black") +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
   theme_bw() +
   facet_wrap(~FE, scales = "free") +
   theme(legend.position = "none")



### NITRATE + NITRITE ###



FE_data <- Full_data %>%
  group_by(CowTagID, AlphaTag, FE, Taxon_Group, Morph2, Calc, ER) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup() %>%
  left_join(dist) %>%
  left_join(chem)

morphd <- FE_data %>%
  select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Morph2, pCover) %>%
  group_by(NN_umolL, Morph2) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(lm(pCover ~ NN_umolL*Morph2, data = morphd))

taxond <- FE_data %>%
  select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Taxon_Group, pCover) %>%
  group_by(NN_umolL, Taxon_Group) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(lm(pCover ~ NN_umolL*Taxon_Group, data = taxond)) # taxa NS ~ phosphate

calcd <- FE_data %>%
  select(dist_to_seep_m, NN_umolL, Phosphate_umolL, Calc, pCover) %>%
  group_by(NN_umolL, Calc) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(lm(pCover ~ NN_umolL*Calc, data = calcd))


erd <- FE_data %>%
  select(dist_to_seep_m, NN_umolL, Phosphate_umolL, ER, pCover) %>%
  group_by(NN_umolL, ER) %>%
  mutate(pCover = sum(pCover)) %>%
  distinct()
anova(lm(pCover ~ NN_umolL*ER, data = erd)) # er NS ~ phosphate


anova(lm(pCover ~ NN_umolL*FE, data = Full_data %>% left_join(chem)))

Full_data %>%
  left_join(dist) %>%
  left_join(chem) %>%
  group_by(NN_umolL, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  ggplot(aes(x = NN_umolL, y = pCover, color = FE)) +
  geom_point()+
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  theme_bw() +
  facet_wrap(~FE, scales = "free") +
  theme(legend.position = "none")


