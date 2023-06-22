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
  filter(CowTagID != "VSEEP") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))


### Join Sp and FE and Vol4D with metadata
reg.Fric <- resFric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13")




###############################
# RESIDUAL MODELS (~ RUGOSITY)
###############################

### Calculate regressions again with Residuals (remove structure/substrate)

#fit model: linear relationship
resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
resModSpRp <- lm(NbSpP ~ meanRugosity, data=reg.Fric) # relative species richness
resModFERp <- lm(NbFEsP ~ meanRugosity, data=reg.Fric) # relative entity richness
resModVol <- lm(Vol8D ~ meanRugosity, data=reg.Fric) # relative entity volume

#view model summary
summary(resModSpRp) # ** strong significance
summary(resModFERp) # *** strong significance
summary(resModVol) # * 0.49 weak significance

# use above to justify use of resFric
resFric <- reg.Fric %>%
  left_join(chem)

## plot richness and volume to rugosity
resFric %>%
  pivot_longer(cols = c(NbSpP, NbFEsP, Vol8D), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = meanRugosity, y = Values)) +#, color = NN_umolL)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw()
#scale_color_continuous(low = "yellow", high = "blue")

## WHEN SEEP IS REMOVED, DIST TO SEEP IS NO LONGER SIGNIFICANT FOR ALL PARAMETERS
## PHOSPHATE AND NN SIGNIFICANT FOR ALL (POLYNOMIAL): increase Rich and Vol with elevating NN and Phosphate, then rich and vol drop off
## SpR: increase with increasing ammonia and visible humidics
## Vol: Decrease with increasing M_C, Inc with Inc Ammonia, Inc with elevating Salinity, then vol drops again



###############################
# MODEL RESIDUALS ~ SGD PARAMETERS
###############################

#check residuals against other parameters
p1param <- resFric %>%
  rename('Nitrate+Nitrite (umol/L)' = NN_umolL, 'Phosphate (umol/L)' = Phosphate_umolL,
         'Silicate (umol/L)' = Silicate_umolL, 'Temperature (C)' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite (umol/L)'), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resSpp)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free_x") +
  theme_bw() +
  labs(x = "Parameter Values", y = "Relative Species Richness (rugosity-normalized)") +
  theme(strip.background = element_rect(fill = "white"))

summary(lm(data = resFric, resSpp ~ VisibleHumidic_Like)) # *
summary(lm(data = resFric, resSpp ~ poly(Phosphate_umolL,2))) # *
summary(lm(data = resFric, resSpp~ poly(NN_umolL,2))) # *


p2param <- resFric %>%
  rename('Nitrate+Nitrite (umol/L)' = NN_umolL, 'Phosphate (umol/L)' = Phosphate_umolL,
         'Silicate (umol/L)' = Silicate_umolL, 'Temperature (C)' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite (umol/L)'), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resFEp)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free_x") +
  theme_bw() +
  labs(x = "Parameter Values", y = "Relative FE Richness (rugosity-normalized)") +
  theme(strip.background = element_rect(fill = "white"))

summary(lm(data = resFric, resFEp ~ poly(NN_umolL,2))) # *
summary(lm(data = resFric, resFEp ~ poly(Phosphate_umolL,2))) # *


p3param <- resFric %>%
  rename('Nitrate+Nitrite (umol/L)' = NN_umolL, 'Phosphate (umol/L)' = Phosphate_umolL,
         'Silicate (umol/L)' = Silicate_umolL, 'Temperature (C)' = Temperature) %>%
  select(-TA) %>%
  pivot_longer(cols = c(Salinity:'Nitrate+Nitrite (umol/L)'), names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resVol)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free_x") +
  theme_bw() +
  labs(x = "Parameter Values", y = "Relative FE Volume (rugosity-normalized)") +
  theme(strip.background = element_rect(fill = "white"))

summary(lm(data = resFric, resVol ~ poly(NN_umolL,2))) # *
summary(lm(data = resFric, resVol ~ poly(Phosphate_umolL,2))) # *
summary(lm(data = resFric, resVol ~ poly(Salinity,2))) # *

ggsave(here("Output", "PaperFigures", "LM_param_Spp.png"), p1param, height = 6, width = 6)
ggsave(here("Output", "PaperFigures", "LM_param_FEp.png"), p2param, height = 6, width = 6)
ggsave(here("Output", "PaperFigures", "LM_param_Vol.png"), p3param, height = 6, width = 6)



#####################################################################
### DISTANCE TO SEEP
#####################################################################

## Raw richness residuals
rug_res_SpR_plot <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSp)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_SpR_plot

rug_res_FER_plot <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFE)) +
  geom_point(size = 2.5) +
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
rug_res_SpRp_plot <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_SpRp_plot

rug_res_FERp_plot <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFEp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_FERp_plot

rug_res_Vol_plot <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resVol)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2


#####################################################################
### PHOSPHATE (FIGURE 2)
#####################################################################

plotfun <-function(data = resFric %>% left_join(alphatag), y){

  y <- enquo(y)

  plota <- data %>%
    ggplot(aes(x = Phosphate_umolL,
               y = !!y)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 13)) +
    theme(panel.grid = element_blank()) +
    labs(y = "",
         x = "")  +
    geom_text_repel(aes(label = AlphaTag),
                    size = 6)


  return(plota)

}

## Relative non-normalized richness
rug_SpR_plot <- plotfun(y = NbSpP) + labs(y = "% SR", x = "")+
  ylim(min = 0, max = 100)
rug_SpR_plot

rug_FER_plot <- plotfun(y = NbFEsP) + labs(y = "% FER", x = "")+
  ylim(min = 0, max = 100)
rug_FER_plot

rug_Vol_plot <- plotfun(y = Vol8D) + labs(y = "% FEV", x = "")+
  ylim(min = 0, max = 100)
rug_Vol_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- plotfun(y = resSpp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% SR (Rugosity-normalized)",
       x = "CV Phosphate (umol/L)")
rug_res_SpRp_plot

rug_res_FERp_plot <- plotfun(y = resFEp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FER (Rugosity-normalized)",
       x = "CV Phosphate (umol/L)")
rug_res_FERp_plot

rug_res_Vol_plot <- plotfun(y = resVol) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "% FEV (Rugosity-normalized)",
       x = "CV Phosphate (umol/L)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_SpR_plot + rug_FER_plot + rug_Vol_plot
plot1
#ggsave(here("Output", "PaperFigures", "lm_relative_raw_Phos.png"), plot1, height = 6, width = 10)


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2
#ggsave(here("Output", "PaperFigures", "lm_relative_residuals_Phos.png"), plot2, height = 6, width = 10)

divPlots <- plot1 / plot2 +
  plot_annotation(tag_levels = 'A') +
  theme(plot.tag = element_text(size = ))
divPlots
ggsave(here("Output", "PaperFigures", "LM_diversity_Phosphate.png"), divPlots, height = 6, width = 10)

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

nnRat <- resFric %>%
  ggplot(aes(x = NN_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FE / Sp",
       x = "CV of NN (umol/L)") +
  #ylim(min = 0, max = 1) +
  geom_text_repel(aes(label = CowTagID)) +
  labs(y = "")

pRat <- resFric %>%
  ggplot(aes(x = Phosphate_umolL,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FE / Sp",
       x = "CV of Phosphate (umol/L)") +
  #ylim(min = 0, max = 1) +
  geom_text_repel(aes(label = CowTagID)) +
  labs(y = "")


dRat <- resFric %>%
  ggplot(aes(x = dist_to_seep_m,
             y = NbFEs / NbSp)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FE / Sp",
       x = "Distance to seep (m)") +
  #ylim(min = 0, max = 1) +
  geom_text_repel(aes(label = CowTagID))

dRat + pRat + nnRat+
  plot_annotation(tag_levels = 'A')

summary(lm(data = resFric, (NbFEs / NbSp) ~ dist_to_seep_m))
summary(lm(data = resFric, (NbFEs / NbSp) ~ NN_umolL))
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
