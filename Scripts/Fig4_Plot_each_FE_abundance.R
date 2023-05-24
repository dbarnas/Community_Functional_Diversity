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
fullchem <- read_csv(here("Data", "Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Location == "Varari",
         Season == "Dry")
chem <- fullchem %>%
  filter(Parameters == "Phosphate_umolL" | Parameters == "NN_umolL") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
fullchem <- fullchem %>%
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

# Full_data$NN_umolL <- factor(Full_data$NN_umolL)
# Full_data$Phosphate_umolL <- factor(Full_data$Phosphate_umolL)
Full_data$Morph2 <- factor(Full_data$Morph2,
                           levels = c('Br', 'Dig', 'Fol', 'Fil', 'Stol',
                                      'Mush', 'Poly', 'Cushion', 'Mas', 'Enc', 'Sph'))
# arrange alphatag factor for nutrient levels
AlphaOrder <- Full_data %>%
  left_join(dist) %>%
  arrange(dist_to_seep_m) %>%
  distinct(AlphaTag)
AlphaOrder <- AlphaOrder$AlphaTag
Full_data$AlphaTag <- factor(Full_data$AlphaTag, levels = AlphaOrder)

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
    summarise(pCover = sum(pCover)) %>%
    mutate(indep = round(as.numeric({{param}}),2))

  myfacet <- enquo(entity)
  x.var <- colnames(my_data[,3])

  plota <- my_data %>%
    ggplot(aes(x = indep, #!!independent,
               y = pCover,
               color = !!myfacet
               )) +
    geom_point() +
    #geom_smooth(method = "lm", formula = "y~x", color = "black") +
    #geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    facet_wrap(myfacet, scales = "free") +
    labs(x = paste(x.var), y = "% Cover")

  return(plota)

}

Param_data <- Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  mutate(entity = as.character(Morph2)) %>%
  group_by(AlphaTag,entity,NN_umolL) %>%
  summarise(pCover = sum(pCover)) %>%
  filter(entity != "Cyanobacteria" & entity != "Mush")


pval <- function(data = Full_data, entity, param, form = "poly"){ # poly or lm


  Param_data <- data %>%
    filter(CowTagID != "VSEEP") %>%
    mutate(entity = as.character({{entity}})) %>%
    group_by(AlphaTag,entity, {{param}}) %>%
    summarise(pCover = sum(pCover)) %>%
    filter(entity != "Cyanobacteria" & # single point in Taxon_Group
           entity != "Mush") %>% # single point in Morph2
    mutate(entity = if_else(entity == "Non-AC", "NAC", entity))

  distinctFT <- unique(Param_data$entity)

  p_df <- tibble(FTrait = as.character(),
                 Parameter = as.character(),
                 pvalue1 = as.numeric(),
                 pvalue2 = as.numeric(),
                 r_squared = as.numeric(),
                 adj_r_squared = as.numeric())

  for(i in 1:length(distinctFT)){

    mydata <- Param_data %>%
      filter(entity %in% distinctFT[i]) %>%
      pivot_wider(names_from = entity, values_from = pCover)

      Parameter = colnames(mydata)[2]
      FTrait = distinctFT[i]
      if(form == "poly"){
      mod <- lm(paste(FTrait, "~ poly(", Parameter, ",2)"), data = mydata)
      pvalue1 <- summary(mod)[4]$coefficients[11]
      pvalue2 <- summary(mod)[4]$coefficients[12]
      r_squared <- summary(mod)[8]$r.squared
      adj_r_squared <- summary(mod)[9]$adj.r.squared
      } else if(form == "lm"){
        mod <- lm(paste(FTrait, "~", Parameter), data = mydata)
        pvalue1 <- summary(mod)$coefficients[8]
        pvalue2 <- NA
        r_squared <- summary(mod)$r.squared
        adj_r_squared <- summary(mod)$adj.r.squared
      }

      temp <- as_tibble(cbind(FTrait, Parameter,
                              pvalue1, pvalue2, r_squared, adj_r_squared)) %>%
        mutate(pvalue1 = as.numeric(pvalue1),
               pvalue2 = as.numeric(pvalue2),
               r_squared = as.numeric(r_squared),
               adj_r_squared = as.numeric(adj_r_squared)) %>%
        mutate(FTrait = if_else(FTrait == "NAC", "Non-AC", FTrait))

      p_df <- p_df %>%
        rbind(temp)
      }
  return(p_df)
}

### PLOT FUNCTIONS ACROSS SEEP ###
### FUNCTIONAL ENTITIES
anova(lm(pCover ~ poly(Phosphate_umolL,2)*FE, data = Full_data %>%
             filter(CowTagID != "VSEEP"))) # using all abundance data at each location
### SUMMARY
# Does an individual entity change its relative abundance along the phosphate gradient? Yes!
## Chlorophyta, Fil, NC, Auto p=0.005 (no interaction...not sure what this means then)
## Cnidaria, Mas, Herm, Mix p=0.003 (no interaction...not sure what this means then)
# Does an individual entity change along the NN gradient? Yes!
## Cnidaria, Mas, Herm, Mix p=5.4e-5 (interactio bw NN and FE)
## Phaeophyta, Fol, Non-AC, Auto p=0.002 (no interaction...not sure what this means then)
## Chlorophyta, Fil, NC, Auto p=0.01 (interaction bw NN and FE)
### ANOVA
# NN:FE interaction p<0.0008, F=2.08, DF=34
# P:FE interaction p<0.002, F=1.99, DF=34

anova(lm(pCover ~ poly(Phosphate_umolL,2)*FE, data = Full_data)) # including the seep
### ANOVA
# NN:FE interaction p=0.001, F=2.03, DF=34
# P:FE no interaction


### TAXA
# stacked bar
pt <- myplot(Taxon_Group, taxonpalette)
pt2 <- myplot(Taxon_Group, taxonpalette) + theme(legend.position = "none")
pt
# regression
ppt <- ptplot(Taxon_Group, Phosphate_umolL)
npt <- ptplot(Taxon_Group, NN_umolL)
ppt + npt
# pvalues
pv1 <- pval(entity = Taxon_Group, param = Phosphate_umolL)
pv2 <- pval(entity = Taxon_Group, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*Taxon_Group, data = Full_data %>% filter(CowTagID != "VSEEP")))
anova(lm(pCover~poly(Phosphate_umolL,2)*Taxon_Group, data = Full_data %>% filter(CowTagID != "VSEEP")))
## no interaction of NN:FE or P:FE


### MORPHOLOGY
# stacked bar
pm <- myplot(Morph2, morphpalette)
pm2 <- myplot(Morph2, morphpalette) + theme(legend.position = "right")
pm3 <- myplot(Morph2, morphpalette) + theme(legend.position = "none")
# regression
ppm <- ptplot(Morph2, Phosphate_umolL) + geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black")
npm <- ptplot(Morph2, NN_umolL) + geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black")
ppm + npm
# pvalues
pv3 <- pval(entity = Morph2, param = Phosphate_umolL)
pv4 <- pval(entity = Morph2, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*Morph2, data = Full_data %>% filter(CowTagID != "VSEEP")))
anova(lm(pCover~poly(Phosphate_umolL,2)*Morph2, data = Full_data %>% filter(CowTagID != "VSEEP")))
## NN:FE p<2e-6, F=3.67, DF=18
## P:FE p<5e-7, F=3.87, DF=18
# summary
summary(lm(pCover~poly(NN_umolL,2)*Morph2, data = Full_data %>% filter(CowTagID != "VSEEP")))
summary(lm(pCover~poly(Phosphate_umolL,2)*Morph2, data = Full_data %>% filter(CowTagID != "VSEEP")))


# CALCIFICATION
# stacked bar
pc <- myplot(Calc, calcpalette)
# regression
pper <- ptplot(Calc, Phosphate_umolL)
nper <- ptplot(Calc, NN_umolL) + geom_smooth(method = "lm", formula = "y~x", color = "black")
nper2 <- ptplot(Calc, NN_umolL) + geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black")
#pper + nper
nper + nper2
# pvalues
pv5 <- pval(entity = Calc, param = Phosphate_umolL)
pv6 <- pval(entity = Calc, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*Calc, data = Full_data %>% filter(CowTagID != "VSEEP")))
anova(lm(pCover~poly(Phosphate_umolL,2)*Calc, data = Full_data %>% filter(CowTagID != "VSEEP")))
## NN:FE p<2e-6, F=3.67, DF=18
## P:FE p<5e-7, F=3.87, DF=18
# summary
summary(lm(pCover~poly(NN_umolL,2)*Calc, data = Full_data %>% filter(CowTagID != "VSEEP")))
summary(lm(pCover~poly(Phosphate_umolL,2)*Calc, data = Full_data %>% filter(CowTagID != "VSEEP")))


# ENERGETIC RESOURCE
# stacked bar
per <- myplot(ER, erpalette)
# regression
pper <- ptplot(ER, Phosphate_umolL)
nper <- ptplot(ER, NN_umolL)
pper + nper
# pvalues
pv7 <- pval(entity = ER, param = Phosphate_umolL)
pv8 <- pval(entity = ER, param = NN_umolL)
# anova
anova(lm(pCover~poly(NN_umolL,2)*ER, data = Full_data %>% filter(CowTagID != "VSEEP")))
anova(lm(pCover~poly(Phosphate_umolL,2)*ER, data = Full_data %>% filter(CowTagID != "VSEEP")))
## NN:FE p<2e-6, F=3.67, DF=18
## P:FE p<5e-7, F=3.87, DF=18
# summary
summary(lm(pCover~poly(NN_umolL,2)*ER, data = Full_data %>% filter(CowTagID != "VSEEP")))
summary(lm(pCover~poly(Phosphate_umolL,2)*ER, data = Full_data %>% filter(CowTagID != "VSEEP")))


# bind pvalue df
mypval <- rbind(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8)
mypval %>%
  filter(pvalue1 < 0.05 | pvalue2 < 0.05)



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


 # nutrients relative to distance
 summary(lm(data = Full_data %>% filter(CowTagID != "VSEEP") %>% left_join(dist), NN_umolL ~ poly(dist_to_seep_m,2)))
 Full_data %>%
   filter(CowTagID != "VSEEP") %>%
   left_join(dist) %>%
   select(AlphaTag, Phosphate_umolL, NN_umolL, dist_to_seep_m) %>%
   rename('Nitrate+Nitrite' = NN_umolL,
          'Phosphate' = Phosphate_umolL) %>%
   pivot_longer(cols = c('Phosphate', 'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
   ggplot(aes(x = dist_to_seep_m, y = Values)) +
   geom_point(aes(color = Parameters), show.legend = FALSE) +
   geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
   facet_wrap(~Parameters, scales = "free") +
   theme_bw() +
   labs(x = "Distance to seep (m)", y = "CV of Nutrient Values (umol/L)")


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


# nutrients relative to silicate and salinity
summary(lm(data = Full_data %>% filter(CowTagID != "VSEEP") %>% left_join(fullchem), NN_umolL ~ Silicate_umolL))
Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  left_join(fullchem) %>%
  select(AlphaTag, Phosphate_umolL, NN_umolL, Silicate_umolL) %>%
  rename('Nitrate+Nitrite' = NN_umolL,
         'Phosphate' = Phosphate_umolL) %>%
  pivot_longer(cols = c('Phosphate', 'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Silicate_umolL, y = Values)) +
  geom_point(aes(color = Parameters), show.legend = FALSE) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw() +
  labs(x = "CV of Silicate (umol/L)", y = "CV of Nutrient Values (umol/L)")
