#### Model selection on linear and nonlinear regressions of Fric residuals and SGD parameters
#### Species Richness, FE Richness, FE Volume vs
#### N+N, Phosphate, Silicate, Salinity, Temperature, pH

### Created by Danielle Barnas
### Created on May 31, 2023

###############################
# LOAD LIBRARIES
###############################
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
library(tidytext)
library(AICcmodavg)
library(kableExtra)


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


# use above to justify use of resFric
resFric <- reg.Fric %>%
  left_join(chem)


#####################################################################
### MODEL SELECTION (FIGURE 3)
#####################################################################

modTable <- tibble(Y = as.character(),
                     Parameter = as.character(),
                     Reg_Type = as.character(), # linear or nonlinear
                     AICc = as.numeric(),
                     R2 = as.numeric(),
                     pVal = as.numeric())

myDep <- colnames(resFric %>% select(resSpp, resFEp, resVol))
mydata <- resFric %>%
  select(CowTagID, resSpp, resFEp, resVol, Salinity, Temperature, pH, Phosphate_umolL, Silicate_umolL, NN_umolL)

# mymods1 <- list()
# mymods2 <- list()

for(i in myDep){
  Y <- as.character(i)
  k <- which(myDep == i) # use as multiplier for list

  for(j in 5:ncol(mydata)){
    Parameter <- colnames(mydata[j])
    model1 <- lm(paste0(Y, "~", Parameter), data = mydata)
    subdata1 <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Reg_Type = "Linear",
             AICc = AICc(model1),
             R2 = summary(model1)$r.squared,
             pVal = summary(model1)$coefficients[8])
    model2 <- lm(paste0(Y, "~ poly(", Parameter, ",2)"), data = mydata)
    subdata2 <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Reg_Type = "Polynomial",
             AICc = AICc(model2),
             R2 = summary(model2)$r.squared,
             pVal = summary(model2)$coefficients[12])
    modTable <- modTable %>%
      rbind(subdata1) %>%
      rbind(subdata2)

     # m <- j-4 # use as multiplier for list
     # n <- k*m
     # mymods1[[n]] <- model1
     # mymods2[[n]] <- model2
  }
}

# View AIC and R2 values and calculate delAIC
modelTable <- modTable %>%
  group_by(Y) %>%
  arrange(AICc) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate",
                     if_else(Parameter == "Silicate_umolL", "Silicate", Parameter))))) %>%
  mutate(modelParam = paste("Y ~",Parameter)) %>%
  relocate(modelParam, .after = Y)
modelTable



# create AIC table for paper
write_csv(modelTable, here("Output", "PaperFigures","Model_Selection_Table.csv"))


# final table style
modelTable %>%
  kbl() %>% # group by multiple rows
  kable_styling(fixed_thead = TRUE)
  #collapse_rows(columns = 1:2, valign = "middle") # should but doesn't work

AICplot <- modelTable %>%
  group_by(Y) %>%
  mutate(Y = if_else(Y == "resSpp", "% SR", if_else(Y == "resFEp", "% FER", "% FEV")),
         Y = factor(Y, levels = c('% SR', '% FER', '% FEV')),
         Parameter = if_else(Parameter == "Silicate_umolL", "Silicate (umol/L)",
                             if_else(Parameter == "Phosphate_umolL", "Phosphate (umol/L)",
                                     if_else(Parameter == "NN_umolL", "N+N (umol/L)",
                                             if_else(Parameter == "Temperature", "Temperature (C)", Parameter))))) %>%
  mutate(minAICc = min(AICc),
         deltaAICc = AICc - minAICc) %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(.f = Parameter, .x = desc(deltaAICc)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  scale_fill_manual(values = c("lightgrey", "darkgrey")) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7) +
  # add mu symbol to labels
  scale_y_discrete(labels = rev(c(expression("Phosphate ("*mu*"mol/L)"),
                                  expression("N+N ("*mu*"mol/L)"),
                                  'pH',
                                  'Salinity',
                                  expression("Silicate ("*mu*"mol/L)"),
                                  'Temperature (C)')))
AICplot


ggsave(here("Output", "PaperFigures", "AIC_model_selection.png"), AICplot, width = 6, height = 6, device = "png")



#####################################################################
### MODEL SELECTION PARAMETERS (Supplemental Table 2)
#####################################################################


#mymods # list of my models in modelTable

# LM
model.names <- modelTable %>%
  unite(Y, Parameter, col = "Model.Name", remove = T, sep = " ~ ") %>%
  unite(Model.Name, Reg_Type, col = "Model.Name", sep = ", ")
model.names <- model.names$Model.Name

all.mods <- lapply(mymods, class)
check.class <- unique(all.mods)
check.class

aictab(cand.set = mymods, modnames = model.names)


# Polynomial

modelTable %>%
  mutate(llik = (2 - AICc)/2) # calculate the log-likelihood estimate
## AIC = 2K - 2ln(L)
## 2ln(L) = 2K - AIC






