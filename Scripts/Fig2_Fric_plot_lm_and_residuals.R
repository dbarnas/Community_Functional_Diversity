#### Linear Regressions of Fric and Fric residuals
### species richness, functional entity richness, and functional volume

### Created by Danielle Barnas
### Created on February 26, 2023


library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)

Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv")) %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         CowTagID != "V13" &
           CowTagID != "VSEEP") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal)


### Join Sp and FE and Vol4D with metadata
reg.Fric <- Fric %>%
  as_tibble() %>%
  left_join(meta) %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "V13")
  #mutate(meanRugosity = if_else(CowTagID == "VSEEP", 0.97, meanRugosity))




### Calculate regressions again with Residuals (remove structure/substrate)

#fit model: linear relationship
resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
resModSpRp <- lm(NbSpP ~ meanRugosity, data=reg.Fric) # relative species richness
resModFERp <- lm(NbFEsP ~ meanRugosity, data=reg.Fric) # relative entity richness
resModVol <- lm(Vol8D ~ meanRugosity, data=reg.Fric) # relative entity volume

#view model summary
summary(resModSpRp) # *** strong significance
summary(resModFERp) # *** significance
summary(resModVol) # *** strong significance

#calculate the standardized residuals
resSpR <- residuals(resModSpR)
resFER <- residuals(resModFER)
resSpRp <- residuals(resModSpRp)
resFERp <- residuals(resModFERp)
resVol <- residuals(resModVol)

#column bind standardized residuals back to original data frame
res_data <- reg.Fric %>%
  cbind(resSpR) %>%
  cbind(resFER) %>%
  cbind(resSpRp) %>%
  cbind(resFERp) %>%
  cbind(resVol) %>%
  left_join(chem)

## WHEN SEEP IS REMOVED, DIST TO SEEP IS NO LONGER SIGNIFICANT FOR ALL PARAMETERS
## PHOSPHATE AND NN SIGNIFICANT FOR ALL (POLYNOMIAL): increase Rich and Vol with elevating NN and Phosphate, then rich and vol drop off
## SpR: increase with increasing ammonia and visible humidics
## Vol: Decrease with increasing M_C, Inc with Inc Ammonia, Inc with elevating Salinity, then vol drops again

#check residuals against other parameters
res_data %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resSpRp)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw()
summary(lm(data = res_data, resSpRp ~ dist_to_seep_m)) # **
summary(lm(data = res_data, resSpRp ~ poly(Ammonia_umolL,2))) # *
summary(lm(data = res_data, resSpRp ~ poly(NN_umolL,2))) # *
summary(lm(data = res_data, resSpRp ~ VisibleHumidic_Like)) # *
summary(lm(data = res_data, resSpRp ~ poly(Phosphate_umolL,2))) # *




res_data %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resFERp)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw()
summary(lm(data = res_data, resFERp ~ dist_to_seep_m)) # NS
summary(lm(data = res_data, resFERp ~ poly(NN_umolL,2))) # *
summary(lm(data = res_data, resFERp ~ poly(Phosphate_umolL,2))) # *

res_data %>%
  pivot_longer(cols = Salinity:Tyrosine_Like, names_to = "Parameters", values_to = "Values") %>%
  ggplot(aes(x = Values, y = resVol)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "red") +
  facet_wrap(~Parameters, scales = "free") +
  theme_bw()
summary(lm(data = res_data, resVol ~ dist_to_seep_m)) # NS
summary(lm(data = res_data, resVol ~ Ammonia_umolL)) # *
summary(lm(data = res_data, resVol ~ M_C)) # *
summary(lm(data = res_data, resVol ~ poly(NN_umolL,2))) # *
summary(lm(data = res_data, resVol ~ poly(Phosphate_umolL,2))) # *
summary(lm(data = res_data, resVol ~ poly(Salinity,2))) # *

#plot predictor variable vs. standardized residuals

# isolate seep poing
# ptseep <- Fric %>%
#   filter(CowTagID == "VSEEP")



#####################################################################
### DISTANCE TO SEEP
#####################################################################

## Raw richness residuals
rug_res_SpR_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpR)) +
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

rug_res_FER_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFER)) +
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
  labs(y = "FER (Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_FER_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpRp)) +
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
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_SpRp_plot

rug_res_FERp_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFERp)) +
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
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_FERp_plot

rug_res_Vol_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resVol)) +
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
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "Distance to seep (m)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1
ggsave(here("Output", "PaperFigures", "lm_rich_residuals_dist.png"), plot1, height = 6, width = 10)


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2

ggsave(here("Output", "PaperFigures", "lm_relative_residuals_dist.png"), plot2, height = 6, width = 10)

#####################################################################
### PHOSPHATE
#####################################################################

## Raw richness residuals
rug_res_SpR_plot <- res_data %>%
  ggplot(aes(x = Phosphate_umolL,
             y = resSpR)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Rugosity-normalized)",
       x = "Phosphate (umol L-1)")
rug_res_SpR_plot

rug_res_FER_plot <- res_data %>%
  ggplot(aes(x = Phosphate_umolL,
             y = resFER)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FER (Rugosity-normalized)",
       x = "Phosphate (umol L-1)")
rug_res_FER_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- res_data %>%
  ggplot(aes(x = Phosphate_umolL,
             y = resSpRp)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "Phosphate (umol L-1)")
rug_res_SpRp_plot

rug_res_FERp_plot <- res_data %>%
  ggplot(aes(x = Phosphate_umolL,
             y = resFERp)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "Phosphate (umol L-1)")
rug_res_FERp_plot

rug_res_Vol_plot <- res_data %>%
  ggplot(aes(x = Phosphate_umolL,
             y = resVol)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "Phosphate (umol L-1)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1
ggsave(here("Output", "PaperFigures", "lm_rich_residuals_Phos.png"), plot1, height = 6, width = 10)


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2

ggsave(here("Output", "PaperFigures", "lm_relative_residuals_Phos.png"), plot2, height = 6, width = 10)


#####################################################################
### NITRATES + NITRITES
#####################################################################

## Raw richness residuals
rug_res_SpR_plot <- res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resSpR)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "SR (Rugosity-normalized)",
       x = "NN (umol L-1)")
rug_res_SpR_plot

rug_res_FER_plot <- res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resFER)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "FER (Rugosity-normalized)",
       x = "NN (umol L-1)")
rug_res_FER_plot

## Relative richness and volume residuals
rug_res_SpRp_plot <- res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resSpRp)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative SR (%, Rugosity-normalized)",
       x = "NN (umol L-1)")
rug_res_SpRp_plot

rug_res_FERp_plot <- res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resFERp)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FER (%, Rugosity-normalized)",
       x = "NN (umol L-1)")
rug_res_FERp_plot

rug_res_Vol_plot <- res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resVol)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "NN (umol L-1)")
rug_res_Vol_plot



### Save patched plots
plot1 <- rug_res_SpR_plot + rug_res_FER_plot
plot1
ggsave(here("Output", "PaperFigures", "lm_rich_residuals_NN.png"), plot1, height = 6, width = 10)


plot2 <- rug_res_SpRp_plot + rug_res_FERp_plot + rug_res_Vol_plot
plot2

ggsave(here("Output", "PaperFigures", "lm_relative_residuals_NN.png"), plot2, height = 6, width = 10)



# ratio of FE:Sp richness ~ distance to seep
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity
res_data %>%
  ggplot(aes(x = NN_umolL,
             y = resFER / resSpR)) +
  geom_point(size = 2.5) +
  # geom_point(data = ptseep,
  #            aes(x = dist_to_seep_m,
  #                y = resSpR),
  #            shape = 23,
  #            size = 5,
  #            color = "black",
  #            fill = "black") + # add seep point as icon
  geom_smooth(method = "lm", formula = "y~poly(x,2)", color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank()) +
  labs(y = "Relative FE Volume (%, Rugosity-normalized)",
       x = "NN (umol L-1)")


## Can use the three values above, and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness


