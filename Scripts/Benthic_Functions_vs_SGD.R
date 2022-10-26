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
nutrient <- read_csv(here("Data", "Biogeochem", "July2022", "Turb_NC.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))


#### CLEAN DATA ####

# group macroalgae
funTrait <- taxa %>%
  mutate(coral_algae = if_else(Taxon_Group == "Rhodophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Chlorophyta", true = "Macroalgae", false = if_else(
                                 Taxon_Group == "Phaeophyta", true = "Macroalgae", false = Taxon_Group))))
dist <- meta %>%
  select(CowTagID, dist_to_seep_m)

# calculate percent cover of primary taxon groups
pCover <- funTrait %>%
  full_join(survey) %>%
  group_by(Location, CowTagID) %>%
  count(coral_algae) %>%
  ungroup() %>%
  pivot_wider(names_from = coral_algae, values_from = n) %>%
  mutate_at(vars('Abiotic':'Cyanobacteria'), .funs = as.double) %>% # need to parse from integer to double
  mutate_at(vars('Abiotic':'Cyanobacteria'), .funs = ~if_else(is.na(.) == TRUE, 0, .)) %>% # change NAs to 0
  group_by(Location,CowTagID) %>%
  mutate(total = sum(across(Abiotic:Cyanobacteria), na.rm = T)) %>% # calculate total individuals per category
  mutate(across(Abiotic:Cyanobacteria, ~ . / total * 100), .keep = c("all")) %>%  # calculate % cover for benthic categories per cowtagid
  select(-total) # remove intermediate metric

# binning calcifiers and fleshy+turf algae
pCover <- pCover %>%
  mutate(TotalCalc = CCA + Coral,
         TotalAlgae = Macroalgae + Turf) # + Abiotic)

# join with biogeochemical data
pCover <- pCover %>%
  mutate(Location = if_else(Location == "Varari_Maya", "Varari", Location)) %>%
  left_join(nutrient) %>%
  filter(CowTagID != "VSEEP", CowTagID != "CSEEP") %>%
  drop_na(N_percent)

#write_csv(pCover, here("Data", "Biogeochem", "July2022", "Calcifier_Nutrients_June2022.csv"))

distCover <- pCover %>%
  right_join(dist)


#### PLOTTING ####

# Macroalgae vs Coral

pCover %>%
  ggplot(aes(x = TotalCalc, y = TotalAlgae)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_text(aes(label = CowTagID),
            nudge_y = -1, nudge_x = 1, size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Location)
anova(lm(data = pCover, TotalCalc ~ TotalAlgae))

pCover <- pCover %>%
  filter(CowTagID != "Sand_SGD_2",
         CowTagID != "Sand_SGD_1",
         CowTagID != "V17") # remove outliers 0 or near 0% calcifiers)


## Preview nutrient loading at survey locations
pCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = N_percent), y = N_percent)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
pCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = C_N), y = C_N)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
pCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = del15N), y = del15N)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))

## Preview with distance to seep
distCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = dist_to_seep_m), y = N_percent)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
distCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = dist_to_seep_m), y = C_N)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
distCover %>%
  ggplot(aes(x = fct_reorder(.f = CowTagID, .x = dist_to_seep_m), y = del15N)) +
  geom_point() + labs(x = "CowTagID") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))

distCover %>%
  pivot_longer(cols = c(del15N:N_percent), names_to = "Nutrients", values_to = "Values") %>%
  ggplot(aes(x = dist_to_seep_m, y = Values)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~Nutrients, scales = "free_y") +
  geom_text_repel(aes(label = CowTagID))

anova(lm(data = distCover, N_percent ~ dist_to_seep_m))
anova(lm(data = distCover, C_N ~ dist_to_seep_m))
anova(lm(data = distCover, del15N ~ dist_to_seep_m))


#### Comparing coral/algae ratios across sgd ####

# del15N
p1 <- pCover %>%
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
  annotate("text", x = 4.8, y = 0.1, label = "Calcifier-dominated")+
  annotate("text", x = 4.8, y = -0.1, label = "Fleshy algal-dominated")+
  theme_bw()+
  theme(legend.direction = "horizontal",
        legend.position = c(.18, .1))
p1
ggsave(here("Output","Calcifier_del15N.png"), p1, width = 10, height = 10)


# percent N
p2 <- pCover %>%
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
  annotate("text", x = 1.0, y = 0.1, label = "Calcifier-dominated")+
  annotate("text", x = 1.0, y = -0.1, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none")
p2
ggsave(here("Output","Calcifier_Npercent.png"), p2, width = 10, height = 10)


# C:N ratio
p3 <- pCover %>%
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
  annotate("text", x = 33, y = 0.1, label = "Calcifier-dominated")+
  annotate("text", x = 33, y = -0.1, label = "Fleshy algal-dominated")+
  theme_bw() +
  theme(legend.position = "none")
p3
ggsave(here("Output","Calcifier_CN.png"), p3, width = 10, height = 10)


modBenthic15N<-lm(log((TotalCalc+1)/(TotalAlgae+1))~del15N, data = pCover)
anova(modBenthic15N)

modBenthicNpercent<-lm(log((TotalCalc+1)/(TotalAlgae+1))~N_percent, data = pCover)
anova(modBenthicNpercent)

modBenthicCN<-lm(log((TotalCalc+1)/(TotalAlgae+1))~C_N, data = pCover)
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
        legend.position = c(.18, .1)) #+
  #facet_wrap(~{{Season}})
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


  #modBenthic <- lm(paste(y ,"~", x ,"*", z), data = logval) # with season effect
  modBenthic <- lm(paste(y ,"~", x), data = logval)
  a <- anova(modBenthic)
  p <- as_tibble(a[[5]]) %>%
    drop_na() %>%
    rename(p_value = value) %>%
    #cbind(Param = c(as.character(x),as.character(z),paste0(x,":",z))) # with season effect
    cbind(Param = as.character(x))
  return(p)
}


plotfun(mydata = redCover, x = N_percent, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = N_percent, num = TotalCalc, den = TotalAlgae)

plotfun(mydata = redCover, x = del15N, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = del15N, num = TotalCalc, den = TotalAlgae)

plotfun(mydata = redCover, x = C_N, num = TotalCalc, den = TotalAlgae)
pval(mydata = redCover, x = C_N, num = TotalCalc, den = TotalAlgae)













