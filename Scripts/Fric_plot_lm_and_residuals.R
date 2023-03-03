#### Linear Regressions of Fric and Fric residuals
### species richness, functional entity richness, and functional volume

### Created by Danielle Barnas
### Created on February 26, 2023


library(tidyverse)
library(here)

Fric <- read_csv(here("Data", "Sp_FE_Vol.csv"))
Fric <- data.matrix(column_to_rownames(Fric, var = "CowTagID")) # convert Fric from tibble to matrix array
meta <- read_csv(here("Data", "Full_Metadata.csv"))



### Linear Regression Sp and FE and Vol4D ~ SGD
reg.Fric <- Fric %>%
  as_tibble() %>%
  mutate(CowTagID = relative.sgd) %>%
  relocate(CowTagID, .before = NbSp) %>%
  left_join(meta)

summary(lm(data = reg.Fric, NbSp ~ dist_to_seep_m)) # **
summary(lm(data = reg.Fric, NbFEs ~ dist_to_seep_m)) # **
summary(lm(data = reg.Fric, Vol8D ~ dist_to_seep_m)) # not significant

summary(lm(data = reg.Fric, NbSp ~ N_percent)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ N_percent)) # not significant
summary(lm(data = reg.Fric, Vol8D ~ N_percent)) # not significant

summary(lm(data = reg.Fric, NbSp ~ del15N)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ del15N)) # . 0.099 not significant
summary(lm(data = reg.Fric, Vol8D ~ del15N)) # not significant

summary(lm(data = reg.Fric, NbSp ~ meanRugosity)) # **
summary(lm(data = reg.Fric, NbFEs ~ meanRugosity)) # *
summary(lm(data = reg.Fric, Vol8D ~ meanRugosity)) # **

summary(lm(data = reg.Fric, NbSp ~ adj_CT_depth_cm)) # not significant
summary(lm(data = reg.Fric, NbFEs ~ adj_CT_depth_cm)) # not significant
summary(lm(data = reg.Fric, Vol8D ~ adj_CT_depth_cm)) # not significant




### Calculate regressions again with Residuals (remove structure/substrate)

#fit model: linear relationship
resModSpR <- lm(NbSp ~ meanRugosity, data=reg.Fric) # species richness
resModFER <- lm(NbFEs ~ meanRugosity, data=reg.Fric) # entity richness
resModVol <- lm(Vol8D ~ meanRugosity, data=reg.Fric) # entity volume

#view model summary
summary(resModSpR) # strong significance
summary(resModFER) # significance
summary(resModVol) # strong significance

#calculate the standardized residuals
resSpR <- residuals(resModSpR)
resFER <- residuals(resModFER)
resVol <- residuals(resModVol)

#column bind standardized residuals back to original data frame
res_data <- reg.Fric %>%
  cbind(resSpR) %>%
  cbind(resFER) %>%
  cbind(resVol)



#plot predictor variable vs. standardized residuals
rug_res_SpR_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resSpR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_smooth(method="lm", formula = y~poly(x,2)) +
  #geom_text_repel(aes(label = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Residuals (Species Richness ~ Rugosity)",
       x = "Distance to seep (m)")
rug_res_SpR_plot
summary(lm(data = res_data, resSpR ~ dist_to_seep_m))

rug_res_FER_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resFER)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_smooth(method="lm", formula = y~poly(x,2)) +
  #geom_text_repel(aes(label = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Residuals (Functional Entity Richness ~ Rugosity)",
       x = "Distance to seep (m)")
rug_res_FER_plot
summary(lm(data = res_data, resFER ~ dist_to_seep_m))

rug_res_Vol_plot <- res_data %>%
  ggplot(aes(x = dist_to_seep_m,
             y = resVol)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_smooth(method="lm", formula = y~poly(x,2)) +
  #geom_text_repel(aes(label = CowTagID)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Residuals (Functional Volume ~ Rugosity)",
       x = "Distance to seep (m)")
rug_res_Vol_plot
summary(lm(data = res_data, resVol ~ dist_to_seep_m))


reg.Fric %>%
  ggplot(aes(x = dist_to_seep_m, y = NbFEs/NbSp)) +
  geom_point(aes(color = dist_to_seep_m)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  geom_label(aes(label = CowTagID))
summary(lm(data = reg.Fric, (NbFEs/NbSp) ~ dist_to_seep_m)) # 0.06
summary(lm(data = reg.Fric, (NbFEs/NbSp) ~ meanRugosity)) # *


# ratio of FE:Sp richness ~ distance to seep
# value of 1 would show that each species has its own function, and any lower values show lower functional diversity


## Can use the three values above, and also community composition: either relative abundance or presence-absence
## then can do a permanova / nMDS of community comp with the volume / FErichness


