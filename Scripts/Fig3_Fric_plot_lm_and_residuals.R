#### Linear Regressions of Fric and Fric residuals
### species richness, functional entity richness, and functional volume

### Created by Danielle Barnas
### Created on February 26, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(MuMIn)
library(tidytext)



###############################
# READ IN DATA
###############################
Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari", CowTagID != "V13") %>%
  #filter(CowTagID != "VSEEP") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))


### Join Sp and FE and Vol4D with metadata
resFric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13") %>%
  mutate(meanRugosity = 1-meanRugosity) %>%  # instead of lower values indicating higher structure, now high value indicate high structure
  left_join(chem)



###############################
# RESIDUAL MODELS (~ RUGOSITY)
###############################

### Calculate regressions against Rugosity

#fit model: linear relationship
# resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
# resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
resModSpRp <- lm(NbSpP ~ poly(meanRugosity,2), data=resFric) # relative species richness
resModFERp <- lm(NbFEsP ~ poly(meanRugosity,2), data=resFric) # relative entity richness
resModVol <- lm(Vol8D ~ poly(meanRugosity,2), data=resFric) # relative entity volume

#view model summary
summary(resModSpRp) # ** strong significance
summary(resModFERp) # *** strong significance
summary(resModVol) # * 0.49 weak significance

### Calculate regressions against Salinity

#fit model: linear relationship
# resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
# resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
salModSpRp <- lm(NbSpP ~ poly(Salinity,2), data=resFric) # relative species richness
salModFERp <- lm(NbFEsP ~ poly(Salinity,2), data=resFric) # relative entity richness
salModVol <- lm(Vol8D ~ poly(Salinity,2), data=resFric) # relative entity volume

#view model summary
summary(salModSpRp) # ** strong significance
summary(salModFERp) # *** strong significance
summary(salModVol) # * 0.49 weak significance





## plot richness and volume to rugosity
SuppFig1H <- resFric %>%
  rename("% SR" = NbSpP, "% FER" = NbFEsP, "% FEV" = Vol8D) %>%
  pivot_longer(cols = c('% SR', '% FER','% FEV'), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('% SR', '% FER','% FEV'))) %>%
  ggplot(aes(x = meanRugosity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  facet_wrap(~Parameters, scales = "fixed") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(x = "Mean rugosity", y = "Relative Diversity") +
  plot_annotation(tag_levels = list(c('H')))

ggsave(here("Output", "PaperFigures", "SuppFig1H_rugosity.png"), SuppFig1H, device = "png", height = 6, width = 6)



## plot richness and volume to salinity
SuppFig1Hb <- resFric %>%
  rename("% SR" = NbSpP, "% FER" = NbFEsP, "% FEV" = Vol8D) %>%
  pivot_longer(cols = c('% SR', '% FER','% FEV'), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels = c('% SR', '% FER','% FEV'))) %>%
  ggplot(aes(x = Salinity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  facet_wrap(~Parameters, scales = "fixed") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_text(size = 14)) +
  labs(x = "CV Salinity", y = "Relative Diversity") +
  plot_annotation(tag_levels = list(c('H')))

ggsave(here("Output", "PaperFigures", "SuppFig1H_salinity.png"), SuppFig1Hb, device = "png", height = 6, width = 6)



## WHEN SEEP IS REMOVED, DIST TO SEEP IS NO LONGER SIGNIFICANT FOR ALL PARAMETERS
## PHOSPHATE AND NN SIGNIFICANT FOR ALL (POLYNOMIAL): increase Rich and Vol with elevating NN and Phosphate, then rich and vol drop off
## SpR: increase with increasing ammonia and visible humidics
## Vol: Decrease with increasing M_C, Inc with Inc Ammonia, Inc with elevating Salinity, then vol drops again



###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
# Supplemental Figure 1A-F
###############################

#check residuals against other parameters
funData <- resFric %>%
  rename('Nitrate+Nitrite' = NN_umolL, 'Phosphate' = Phosphate_umolL,
         'Silicate' = Silicate_umolL, 'Temperature' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite'), names_to = "Parameters", values_to = "Values") %>%
  pivot_longer(cols = c(NbSpP, NbFEsP, Vol8D, resSpp, resFEp, resVol, meanRugosity), names_to = "Dependent", values_to = "Dep_Values") %>%
  select(CowTagID, Dependent, Dep_Values, Parameters, Values)



plotFun <- function(Y){
  data <- funData %>%
    filter(Dependent == Y) %>%
    mutate(Dependent = "DependentVar") %>%
    pivot_wider(names_from = Dependent, values_from = Dep_Values)

  plot <- data %>%
    ggplot(aes(x = Values, y = DependentVar)) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
    facet_wrap(~Parameters, scales = "free_x") +
    theme_bw() +
    labs(x = "Parameter Values", y = Y) +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
          axis.title = element_text(size = 12))

  return(plot)
}


a <- plotFun("NbSpP") + labs(y = "% SR") + theme(axis.title.x = element_blank())
b <- plotFun("NbFEsP") + labs(y = "% FER") + theme(axis.title.x = element_blank())
c <- plotFun("Vol8D") + labs(y = "% FEV")

d <- plotFun("resSpp") + labs(y = "% SR (residuals)") + theme(axis.title.x = element_blank())
e <- plotFun("resFEp") + labs(y = "% FER (residuals)") + theme(axis.title.x = element_blank())
f <- plotFun("resVol") + labs(y = "% FEV (residuals)")

rawPlots <- a/b/c +
  plot_annotation(tag_levels = list(c('A', 'B', 'C')))
resPlots <- d/e/f +
  plot_annotation(tag_levels = list(c('D', 'E', 'F')))

ggsave(here("Output", "PaperFigures", "SuppFig1ABC_regression.png"), rawPlots, device = "png", width = 7, height = 11)
ggsave(here("Output", "PaperFigures", "SuppFig1DEF_regression.png"), resPlots, device = "png", width = 7, height = 11)

# ggsave(here("Output", "PaperFigures", "LM_param_resSpp.png"), d, height = 6, width = 6)
# ggsave(here("Output", "PaperFigures", "LM_param_resFEp.png"), e, height = 6, width = 6)
# ggsave(here("Output", "PaperFigures", "LM_param_resVol.png"), f, height = 6, width = 6)




#####################################################################
### Rugosity LM (Supplemental Figure 1G)
#####################################################################

g <- plotFun("meanRugosity") + labs(y = "Mean rugosity") +
  plot_annotation(tag_levels = list(c('G')))
g
ggsave(here("Output", "PaperFigures", "SuppFig1g_regression.png"), g, device = "png", width = 6, height = 6)

# modeling
moddata <- funData %>%
  filter(Dependent == "meanRugosity") %>%
  pivot_wider(names_from = Dependent, values_from = Dep_Values) %>%
  pivot_wider(names_from = Parameters, values_from = Values)
mymod <- lm(data = moddata %>% rename(NN = 'Nitrate+Nitrite'),
            meanRugosity ~ Phosphate)
anova(mymod)
# Silicate, NN, pH (poly),

#Now let's check our assumptions
plot(mymod)
library(car)
qqp(mymod)
resid1<-residuals(mymod)
qqp(resid1, "norm") # numbers indicate outliers
#library(agricolae)
#HSD.test(mymod, "Silicate", console=TRUE)



#####################################################################
### MODELING
#####################################################################



#####################################################################
### PHOSPHATE
#####################################################################

plotfun <-function(data = resFric %>% left_join(alphatag), y){

  y <- enquo(y)

  plota <- data %>%
    ggplot(aes(x = Phosphate_umolL,
               y = !!y)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(),
          #axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1)) +
    labs(x = "")

  return(plota)

}

## Relative non-normalized richness
rug_SpR_plot <- plotfun(y = NbSpP) +
  labs(y = "% SR") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_FER_plot <- plotfun(y = NbFEsP) +
  labs(y = "% FER") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rug_Vol_plot <- plotfun(y = Vol8D) +
  labs(y = "% FEV") +
  theme(axis.title.x = element_blank()) +
  ylim(min = 0, max = 100)

rugosityplot <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot
rugosityplot
#
# ## Relative non-normalized richness
# rug_SpR_plot <- plotfun(y = NbSpP) + labs(y = "% SR", x = "")+
#   ylim(min = 0, max = 100)
# rug_SpR_plot
#
# rug_FER_plot <- plotfun(y = NbFEsP) + labs(y = "% FER", x = "")+
#   ylim(min = 0, max = 100)
# rug_FER_plot
#
# rug_Vol_plot <- plotfun(y = Vol8D) + labs(y = "% FEV", x = "")+
#   ylim(min = 0, max = 100)
# rug_Vol_plot


## Relative richness and volume residuals
rug_res_SpRp_plot <- plotfun(y = resSpp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% SR residuals",
       x = "")

rug_res_FERp_plot <- plotfun(y = resFEp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FER residuals",
       x = expression("CV Phosphate ("*mu*"mol/L)"))

rug_res_Vol_plot <- plotfun(y = resVol) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FEV residuals",
       x = "")

rugosityresplot <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
rugosityresplot

# ## Relative richness and volume residuals
# rug_res_SpRp_plot <- plotfun(y = resSpp) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(y = "% SR (Rugosity-normalized)",
#        x = "CV Phosphate (umol/L)")
# rug_res_SpRp_plot
#
# rug_res_FERp_plot <- plotfun(y = resFEp) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(y = "% FER (Rugosity-normalized)",
#        x = "CV Phosphate (umol/L)")
# rug_res_FERp_plot
#
# rug_res_Vol_plot <- plotfun(y = resVol) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(y = "% FEV (Rugosity-normalized)",
#        x = "CV Phosphate (umol/L)")
# rug_res_Vol_plot



### Save patched plots
# plot1 <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot
# plot1
# #ggsave(here("Output", "PaperFigures", "lm_relative_raw_Phos.png"), plot1, height = 6, width = 10)
#
#
# plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
# plot2
#ggsave(here("Output", "PaperFigures", "lm_relative_residuals_Phos.png"), plot2, height = 6, width = 10)

divPlots <- rugosityplot / rugosityresplot +
  plot_annotation(tag_levels = 'A')
divPlots
ggsave(here("Output", "PaperFigures", "Fig3_LM_diversity_Phosphate.png"), divPlots, height = 6, width = 10)

#####################################################################
### NITRATES + NITRITES
#####################################################################

## Raw richness residuals
rug_res_SpR_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Rugosity-normalized)",
       x = "CV of NN (umol/L)")
rug_res_SpR_plot

rug_res_FER_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resFE)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FER (Rugosity-normalized)",
       x = "CV of NN (umol/L)")
rug_res_FER_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resSpp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")
rug_res_SpRp_plot

rug_res_FERp_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resFEp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")
rug_res_FERp_plot

rug_res_Vol_plot <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = resVol)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "CV of NN (umol/L)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1
#ggsave(here("Output", "PaperFigures", "lm_rich_residuals_NN.png"), plot1, height = 6, width = 10)


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2

#ggsave(here("Output", "PaperFigures", "lm_relative_residuals_NN.png"), plot2, height = 6, width = 10)



# ratio of FE:Sp richness ~ distance to seep INCLUDING SEEP
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity

meta <- read_csv(here("Data","Full_metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
resFric <- read_csv(here("Data", "Sp_FE_Vol_res.csv")) %>%
  as_tibble() %>%
  left_join(meta) %>%
  left_join(chem) %>%
  filter(#CowTagID != "VSEEP",
    CowTagID != "V13")

resFric %>% mutate(FE_SP = NbFEs / NbSp) %>% select(CowTagID,NbFEs, NbSp, FE_SP)

pRat <- resFric %>%
  filter(CowTagID != "VSEEP") %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(panel.grid = element_blank()) +
  labs(y = "FER / SR")

pRatSeep <- resFric %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  theme(panel.grid = element_blank()) +
  labs(y = "FER / SR",
       x = expression("CV Phosphate ("*mu*"mol/L)"))


SpFERatio <- (pRat / pRatSeep) +
   plot_annotation(tag_levels = 'A')

ggsave(here("Output", "PaperFigures", "Sp_FE_Ratio.png"), SpFERatio, device = "png", width = 6, height = 6)

summary(lm(data = resFric %>% filter(CowTagID != "VSEEP"), (NbFEs / NbSp) ~ Phosphate_umolL))
summary(lm(data = resFric, (NbFEs / NbSp) ~ Phosphate_umolL))
## Can use the three values above, and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness



### checking if volume is significantly different across a nutrient cutoff point
summary(lm(data = resFric %>%
             filter(CowTagID != "VSEEP") %>%
             select(CowTagID, resSp:resVol, NN_umolL, Phosphate_umolL) %>%
             mutate(relNut = if_else(NN_umolL < 0.23, "lowN",
                                     if_else(NN_umolL > 0.48,"highN", "midN"))),
           resVol ~ relNut))
# 0.23 - 0.48
# 0.23 - 1 (including seep)
# 0.23 - 0.4


# check nutrients against silicate
summary(lm(Phosphate_umolL ~ Salinity ,data = resFric %>% filter(CowTagID != "VSEEP")))
resFric %>%
  filter(CowTagID != "VSEEP") %>%
  ggplot(aes(x = Silicate_umolL, y = Phosphate_umolL)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
