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
library(PNWColors)


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
           CowTagID != "V13") %>%

  mutate(meanRugosity = 1-meanRugosity)



# use above to justify use of resFric
resFric <- reg.Fric %>%
  left_join(chem)


#####################################################################
### MODEL SELECTION (FIGURE 3)
#####################################################################

# for no interactive effect
modTable <- tibble(Y = as.character(),
                     Parameter = as.character(),
                     Reg_Type = as.character(), # linear or nonlinear
                     AICc = as.numeric(),
                     R2 = as.numeric(),
                     pVal.P = as.numeric(),
                     pVal.Rug = as.numeric())
# for rugosity as interactive effect
# modTable_rug <- tibble(Y = as.character(),
#                      Parameter = as.character(),
#                      Reg_Type = as.character(), # linear or nonlinear
#                      AICc = as.numeric(),
#                      R2 = as.numeric(),
#                      pVal = as.numeric())

myDep <- colnames(resFric %>% select(NbSpP, NbFEsP, Vol8D, resSpp, resFEp, resVol))
mydata <- resFric %>%
  select(CowTagID,NbSpP, NbFEsP, Vol8D, resSpp, resFEp, resVol, Salinity, Temperature, pH, Phosphate_umolL, Silicate_umolL, NN_umolL, meanRugosity)


for(i in myDep){
  Y <- as.character(i)
  k <- which(myDep == i) # use as multiplier for list

  for(j in 8:ncol(mydata)){
    Parameter <- colnames(mydata[j])

    # no interactive effect
    # model1 <- lm(paste0(Y, "~", Parameter), data = mydata)
    # subdata1 <- as_tibble(cbind(Y,Parameter)) %>%
    #   mutate(Reg_Type = "Linear",
    #          AICc = AICc(model1),
    #          R2 = summary(model1)$r.squared,
    #          pVal = summary(model1)$coefficients[8])

    # with rugosity as interactive effect
    model1R <- lm(paste0(Y, "~", Parameter, "+ meanRugosity"), data = mydata)
    subdata1R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter," + Rugosity")) %>%
      mutate(Reg_Type = "Linear",
             AICc = AICc(model1R),
             R2 = summary(model1R)$r.squared,
             pVal.P = summary(model1R)$coefficients[8],
             pVal.Rug = summary(model1R)$coefficients[8])

    # no interactive effect
    # model2 <- lm(paste0(Y, "~ poly(", Parameter, ",2)"), data = mydata)
    # subdata2 <- as_tibble(cbind(Y,Parameter)) %>%
    #   mutate(Reg_Type = "Polynomial",
    #          AICc = AICc(model2),
    #          R2 = summary(model2)$r.squared,
    #          pVal = summary(model2)$coefficients[12])

    # with rugosity as interactive effect
    model2R <- lm(paste0(Y, "~ poly(", Parameter, ",2)+ meanRugosity"), data = mydata)
    subdata2R <- as_tibble(cbind(Y,Parameter)) %>%
      mutate(Parameter = paste0(Parameter," + Rugosity")) %>%
      mutate(Reg_Type = "Polynomial",
             AICc = AICc(model2R),
             R2 = summary(model2R)$r.squared,
             pVal.P = summary(model2R)$coefficients[15],
             pVal.Rug = summary(model2R)$coefficients[16])


    # for rugosity + rugosity
    subdata1R <- subdata1R %>%
      mutate(pVal.P = if_else(is.na(pVal.P), summary(model1R)$coefficients[8], pVal.P),
             pVal.Rug = if_else(is.na(pVal.Rug), summary(model1R)$coefficients[8], pVal.Rug))

    subdata2R <- subdata2R %>%
      mutate(pVal.P = if_else(is.na(pVal.P), summary(model2R)$coefficients[12], pVal.P),
             pVal.Rug = if_else(is.na(pVal.Rug), summary(model2R)$coefficients[12], pVal.Rug))

    modTable <- modTable %>%
      rbind(subdata1R) %>%
      rbind(subdata2R)

    # modTable_rug <- modTable_rug %>%
    #   rbind(subdata1R) %>%
    #   rbind(subdata2R)
  }
}

# View AIC and R2 values and calculate delAIC
# modelTable <- modTable %>%
#   group_by(Y) %>%
#   arrange(AICc) %>%
#   mutate(minAIC = min(AICc)) %>%
#   mutate(delAICc = AICc - minAIC) %>%
#   select(-c(minAIC)) %>%
#   mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
#                      if_else(Parameter == "Phosphate_umolL", "Phosphate",
#                      if_else(Parameter == "Silicate_umolL", "Silicate", Parameter)))) %>%
#   mutate(modelParam = paste("Y ~",Parameter)) %>%
#   relocate(modelParam, .after = Y)
# modelTable
#
# modelTable_rug <- modTable_rug %>%
#   drop_na(pVal) %>%  # remove rugosity:rugosity
#   group_by(Y) %>%
#   arrange(AICc) %>%
#   mutate(minAIC = min(AICc)) %>%
#   mutate(delAICc = AICc - minAIC) %>%
#   select(-c(minAIC)) %>%
#   mutate(Parameter = if_else(Parameter == "NN_umolL:Rugosity", "Nitrate+Nitrite:Rugosity",
#                              if_else(Parameter == "Phosphate_umolL:Rugosity", "Phosphate:Rugosity",
#                                      if_else(Parameter == "Silicate_umolL:Rugosity", "Silicate:Rugosity", Parameter)))) %>%
#   mutate(modelParam = paste("Y ~",Parameter)) %>%
#   relocate(modelParam, .after = Y)


modelTable <- modTable %>%
  #filter(Parameter != "meanRugosity + Rugosity") %>%  # do not include alone in model test
  #drop_na(pVal) %>%  # remove rugosity:rugosity
  group_by(Y) %>%
  mutate(minAIC = min(AICc)) %>%
  mutate(delAICc = AICc - minAIC) %>%
  select(-c(minAIC)) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL", "Nitrate+Nitrite",
                     if_else(Parameter == "Phosphate_umolL", "Phosphate",
                     if_else(Parameter == "Silicate_umolL", "Silicate", Parameter)))) %>%
  mutate(Parameter = if_else(Parameter == "NN_umolL + Rugosity", "Nitrate+Nitrite + Rugosity",
                     if_else(Parameter == "Phosphate_umolL + Rugosity", "Phosphate + Rugosity",
                     if_else(Parameter == "Silicate_umolL + Rugosity", "Silicate + Rugosity",
                     if_else(Parameter == "meanRugosity + Rugosity", "Rugosity", Parameter)))))
  #mutate(modelParam = paste("Y ~",Parameter)) %>%
  #relocate(modelParam, .after = Y)

# create AIC table for paper
write_csv(modelTable, here("Output", "PaperFigures","Model_Selection_Table.csv"))
#write_csv(modelTable_rug, here("Output", "PaperFigures","Model_Selection_Table_rugosity.csv"))

modelTable <- read_csv(here("Output", "PaperFigures","Model_Selection_Table.csv"))
#modelTable_rug <- read_csv(here("Output", "PaperFigures","Model_Selection_Table_rugosity.csv"))

mypal <- pnw_palette("Starfish", n = 2)

modelTableData <- modelTable %>%
  group_by(Y) %>%
  mutate(Y = if_else(Y == "resSpp", "% SR (res)",
             if_else(Y == "resFEp", "% FER (res)",
             if_else(Y == "resVol", "% FEV (res)",
             if_else(Y == "NbSpP", "% SR",
             if_else(Y == "NbFEsP", "% FER",
             if_else(Y == "Vol8D", "% FEV", Y)))))),
         Y = factor(Y, levels = c('% SR', '% FER', '% FEV', '% SR (res)', '% FER (res)', '% FEV (res)')),
         #Parameter = if_else(Parameter == "meanRugosity", "Rugosity", Parameter),
         Reg_Type = factor(Reg_Type, levels = c("Polynomial", "Linear")),
         Parameter = factor(Parameter,
                            levels = c("Phosphate", "Nitrate+Nitrite", "pH", "Salinity", "Silicate", "Temperature", "Rugosity",
                                       "Phosphate + Rugosity", "Nitrate+Nitrite + Rugosity", "pH + Rugosity", "Salinity + Rugosity", "Silicate + Rugosity", "Temperature + Rugosity"))) %>%
  mutate(minAICc = min(AICc),
         deltaAICc = AICc - minAICc)

AICplot <- modelTableData %>%
  #filter(ModelType == "One-way") %>%
  filter(Y== "% SR" | Y == "% FER" | Y== "% FEV") %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)
#
# AICplot_two <- modelTableData %>%
#   filter(ModelType == "Two-way") %>%
#   filter(Y== "% SR" | Y == "% FER" | Y== "% FEV") %>%
#   ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
#   geom_col(position = "dodge", color = "black") +
#   facet_wrap(~Y) +
#   labs(x = expression(Delta*"AICc"),
#        y= "Parameters",
#        fill = "Regression") +
#   theme_bw() +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 10),
#         legend.position = "top",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         panel.grid = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank()) +
#   scale_fill_manual(values = mypal) +
#   geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)

AICplotres <- modelTableData %>%
  filter(Y== "% SR (res)" | Y == "% FER (res)" | Y== "% FEV (res)") %>%
  ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~Y) +
  labs(x = expression(Delta*"AICc"),
       y= "Parameters",
       fill = "Regression") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14)) +
  scale_fill_manual(values = mypal) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)
#
# AICplotres_two <- modelTableData %>%
#   filter(ModelType == "Two-way") %>%
#   filter(Y== "% SR (res)" | Y == "% FER (res)" | Y== "% FEV (res)") %>%
#   ggplot(aes(x = deltaAICc, y = fct_reorder(Parameter, desc(Parameter)), fill = Reg_Type)) +
#   geom_col(position = "dodge", color = "black") +
#   facet_wrap(~Y) +
#   labs(x = expression(Delta*"AICc"),
#        y= "Parameters",
#        fill = "Regression") +
#   theme_bw() +
#   theme(strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 10),
#         legend.position = "top",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         panel.grid = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 14)) +
#   scale_fill_manual(values = mypal) +
#   geom_vline(xintercept = 2, linetype = "dashed", size = 0.7)

AICplotAll <- AICplot / AICplotres +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = 'collect') & theme(legend.position = 'top')
AICplotAll

ggsave(here("Output", "PaperFigures", "AIC_model_selection.png"), AICplotAll, width = 8, height = 6, device = "png")




