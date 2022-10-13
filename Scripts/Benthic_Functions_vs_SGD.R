### Broad functional categorizations across SGD
### Created by Danielle Barnas
### Created on 10/8/2022

#### LIBRARIES ####
library(tidyverse)
library(here)
library(patchwork)
library(ggrepel)


#### READ IN DATA ####
taxa <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
survey <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
nutrient <- read_csv(here("Data", "Biogeochem", "Nutrient_Processed_CV.csv"))


#### CLEAN DATA ####

# group macroalgae
funTrait <- taxa %>%
  mutate(coral_algae = if_else(Taxon_Group == "Rhodophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Chlorophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Phaeophyta", true = "Macroalgae", false = Taxon_Group))))

# calculate percent cover of primary taxon groups
pCover <- funTrait %>%
  full_join(survey) %>%
  group_by(Location, CowTagID) %>%
  count(coral_algae) %>%
  pivot_wider(names_from = coral_algae, values_from = n) %>%
  mutate(total = sum(across(Abiotic:Exposed), na.rm = T)) %>% # calculate total individuals per category
  mutate(across(Abiotic:Exposed, ~ . / total * 100), .keep = c("all")) %>%  # calculate % cover for benthic categories per cowtagid
  select(-total) %>%   # remove intermediate metric
  pivot_longer(Abiotic:Exposed, names_to = "Taxa", values_to = "values") %>% # pivot to replace na with zero
  mutate(values = if_else(is.na(values), 0, values)) %>%
  pivot_wider(names_from = Taxa, values_from = values) # pivot back to wide

# binning calcifiers and fleshy+turf algae
pCover <- pCover %>%
  mutate(TotalCalc = CCA + Coral,
         TotalAlgae = Macroalgae + Turf + Abiotic)

# isolate t ornata nitrogen data
nutrient <- nutrient %>%
  select(Location, CowTagID, Season, del15N, C_N, N_percent) %>%
  distinct()

# join with biogeochemical data
pCover <- pCover %>%
  filter(CowTagID != "VSEEP", CowTagID != "CSEEP") %>%
  left_join(nutrient) %>%
  drop_na(Season)


#### PLOTTING ####

# Macroalgae vs Coral

pCover %>%
  ggplot(aes(x = TotalCalc, y = TotalAlgae)) +
  geom_point() +
  geom_text(aes(label = ifelse(TotalCalc < 5, CowTagID, NA)),
            nudge_y = -1, nudge_x = 1, size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Season)

# Remove outliers
redCover <- pCover %>%
  filter(CowTagID != "V17")

#### Comparing coral/algae ratios across sgd ####

# del15N
p1 <- redCover %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = del15N)) +
  geom_point() +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)), nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, lty = 2)+
  scale_size_binned("%N (Nutrient Loading)") +
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "del15N (Nutrient Source)")+
  annotate("text", x = 5, y = 0.5, label = "Calcifier-dominated")+
  annotate("text", x = 5, y = -0.7, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.direction = "horizontal",
        legend.position = c(.18, .1)) +
  facet_wrap(~Season, scales = "free_x")
p1
ggsave(here("Output","CoralAlgae_15N.png"), p1, width = 10, height = 10)


# percent N
p2 <- redCover %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = N_percent)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)), nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_hline(yintercept = 0, lty = 2) + #lty: line type, 2: dotted line
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "%N (Nutrient Loading)") +
  # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
  # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Season, scales = "free_x")
p2
ggsave(here("Output","CoralAlgae_Npercent.png"), p2, width = 10, height = 10)


# C:N ratio
p3 <- redCover %>%
  filter(CowTagID != "V17") %>%
  ggplot(aes(y = log((TotalCalc+1)/(TotalAlgae+1)), x = C_N))+
  geom_point() +
  geom_smooth(method = "lm", color = "black") +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  # geom_text(aes(label = ifelse(log((TotalCalc+1)/(Macroalgae+1)) > 0, CowTagID, NA)), nudge_y = 0, nudge_x = 0.1, size = 3) +
  geom_hline(yintercept = 0, lty = 2) + #lty: line type, 2: dotted line
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = "log ratio of calcifiers to fleshy algae",
       color = "",
       x = "C:N Ratio (Nutrient Loading)") +
  # annotate("text", x = 1.1, y = 0.5, label = "Calcifier-dominated")+
  # annotate("text", x = 1.1, y = -0.5, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Season, scales = "free_x")
p3
ggsave(here("Output","CoralAlgae_CN.pdf"), p3, width = 10, height = 10)


modBenthic15N<-lm(log((TotalCalc+1)/(TotalAlgae+1))~del15N*Season, data = redCover)
anova(modBenthic15N)

modBenthicNpercent<-lm(log((TotalCalc+1)/(TotalAlgae+1))~N_percent*Season, data = redCover)
anova(modBenthicNpercent)

modBenthicCN<-lm(log((TotalCalc+1)/(TotalAlgae+1))~C_N*Season, data = redCover)
anova(modBenthicCN)

#### Comparing Functional Groups across SGD ####
### Creating functions for plotting and linear regression modeling ###

# plot function
plotfun <- function(mydata, x, num, den) {

  x<-enquo(x)
  num<-enquo(num)
  den<-enquo(den)

  plot <- ggplot(data = mydata,
                 aes(y = log((!!num+1)/(!!den+1)),
                     x = !!x)) +
  geom_point() +
  geom_text_repel(aes(label = CowTagID), size = 3) +
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("#6A3937","#89A5A7"))+
  labs(y = paste("log ratio of", as_label(num), "to", as_label(den)),
       x = paste(as_label(x), "(Nutrient Source)")) +
  #annotate("text", x = 5, y = 0.5, label = "Calcifier-dominated")+
  #annotate("text", x = 5, y = -0.7, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.direction = "horizontal",
        legend.position = c(.18, .1)) +
  facet_wrap(~{{Season}})
  return(plot)
}


# evaluate pvalue of model
pval <- function(mydata, num, den, x){

  logval <- mydata %>%
    mutate(logvalue = log(({{num}} + 1)/({{den}} + 1))) %>%
    select(Location, CowTagID, Season, logvalue, {{x}})

  x <- colnames(logval[,5])
  y <- colnames(logval[,4])
  z <- colnames(logval[,3])


  modBenthic <- lm(paste(y ,"~", x ,"*", z), data = logval)
  a <- anova(modBenthic)
  p <- as_tibble(a[[5]]) %>%
    drop_na() %>%
    rename(r2_value = value) %>%
    cbind(Param = c(as.character(x),as.character(z),paste0(x,":",z)))
  return(p)
}


plotfun(mydata = redCover, x = N_percent, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = N_percent, num = TotalCalc, den = TotalAlgae)

plotfun(mydata = redCover, x = del15N, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = del15N, num = TotalCalc, den = TotalAlgae)

plotfun(mydata = redCover, x = C_N, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = C_N, num = TotalCalc, den = TotalAlgae)




