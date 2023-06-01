##### PCA of richness ~ rugosity residuals with biogeochemical and physical site characteristics
##### Created by Danielle Barnas
##### Created on February 28, 2023

library(tidyverse)
library(here)
library(ggrepel)
library(PNWColors)
library(vegan)
library(pairwiseAdonis)
library(patchwork)
library(ggmap)

#set.seed(7)



# Read in data

#div <- read_csv(here("Data","Surveys","Species_Diversity.csv")) # cannot calculate species div with % cover data
#rich <- read_csv(here("Data","Surveys","Species_Richness.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv"))
comp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
taxa <- read_csv(here("Data","Surveys","Distinct_Taxa.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrient_Processed_CV.csv"))

# create color palette for plotting
mypalette <- pnw_palette("Bay", 19, type="continuous")
midpalette <- pnw_palette(name = "Bay", n = 30)
largepalette <- pnw_palette(name = "Bay", n = 100)



### set up functions
## Create function to order cowtagID's when plotting
order <- function(data = Full_data, param, by){

  param<- enquo(param)
  by<- enquo(by)

  orderParam <- data %>%
    arrange(!!param) %>%
    distinct(!!by) %>%
    mutate(CowTagID = as_factor(as.character(!!by)))

  # sort factor levels by N_percent
  ParamLevels <- paste(orderParam$CowTagID)
  # assign order to factor levels by % N
  data$CowTagID <- factor(data$CowTagID, levels = ParamLevels)
  # levels(data$CowTagID) # check

  return(data)
}

### Create function for geom_point plots with facet_wrapping
# create ggplot function
groupplot <- function(mydata = Full_data, x, y, fw) {

  x<-enquo(x)
  y<-enquo(y)

  plot <- ggplot(data = mydata,
                 aes(y = !!y,
                     x = !!x,
                     color = !!enquo(fw))) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(enquo(fw), scales = "free_y") +
    labs(y = paste(as_label(y)),
         x = paste(as_label(x))) +
    theme_bw()+
    theme(legend.direction = "horizontal",
          legend.position = "top")

  return(plot)
}






### Process and Join data

# reduce metadata, remove cabral and v13
meta <- meta %>%
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  select(Location, CowTagID, LiveCoral:Sand, dist_to_seep_m, adj_CT_depth_cm, meanRugosity) # will add nutrient data later

# calculate species richness
rich <- comp %>%
  select(Location, CowTagID, Taxa, SpeciesCounts) %>%
  mutate(Taxa = if_else(Taxa == 'Bare Rock', "Hard Substrate", Taxa),
         Taxa = if_else(Taxa == 'Rubble', "Hard Substrate", Taxa)) %>%
  filter(Taxa != "Sand" & Taxa != "Hard Substrate") %>% # actually just remove the hard substrate and sand
  group_by(Location, CowTagID, Taxa) %>%
  count(Taxa) %>%
  mutate(n = 1) %>%
  ungroup() %>%
  group_by(Location,CowTagID) %>%
  mutate(spRichness = sum(n)) %>%
  distinct(Location, CowTagID, spRichness) %>%
  ungroup()

# join dfs
Full_data <- rich %>%
  full_join(chem) %>%
  drop_na(lat) %>% # remove largely incomplete datasets
  filter(Location == "Varari",
         CowTagID != "V13") %>%
  #filter(CowTagID != "C14") %>%  # not surveyed
  left_join(meta)


comp <- comp %>%
  select(Location, CowTagID, Taxa, SpeciesCounts) %>%
  mutate(Taxa = if_else(Taxa == 'Bare Rock', "Hard Substrate", Taxa),
         Taxa = if_else(Taxa == 'Rubble', "Hard Substrate", Taxa)) %>%
  filter(Taxa != "Sand" & Taxa != "Hard Substrate") %>% # actually just remove the hard substrate and sand
  group_by(Location, CowTagID, Taxa) %>%
  mutate(SpeciesCounts = sum(SpeciesCounts)) %>%
  distinct() %>%
  ungroup() %>% # ungroup by taxa to regroup by site
  group_by(Location, CowTagID) %>%
  mutate(totalcount = sum(SpeciesCounts),
         sp_percent = SpeciesCounts / totalcount*100) %>%
  select(-totalcount) %>%
  left_join(Full_data) %>%
  ungroup() %>%
  filter(Location == "Varari",
         CowTagID != "V13") # only Varari for now and not v13

# rearrange
Full_data <- Full_data %>%
  relocate(Location, .before = CowTagID)

# only Varari for now
Full_data <- Full_data %>%
  filter(Location == "Varari",
         CowTagID != "V13")



# Normalize to substrate and rugosity

## Richness ~ Rugosity

## plot
groupplot(mydata = Full_data, x = meanRugosity, y = spRichness, fw = NA) +
  geom_text_repel(aes(label = CowTagID)) +
  scale_color_manual(values = mypalette)

# A residual is the difference between an observed value and a predicted value in a regression model. Residual = Observed value – Predicted value
# An observation has a positive residual if its value is greater than the predicted value made by the regression line and a negative residual if its value is less than the predicted value made by the regression line.

#fit model
resMod <- lm(spRichness ~ meanRugosity, data=Full_data %>% filter(CowTagID != "VSEEP"))

#view model summary
summary(resMod)

#calculate the standardized residuals
res <- residuals(resMod)

#view the standardized residuals
res

#column bind standardized residuals back to original data frame
res_data <- Full_data %>%
  filter(CowTagID != "VSEEP") %>%
  cbind(res) %>%
  relocate(res, .after = CowTagID) %>%
  rename(spR_residuals = res)



## Residuals regressions

anova(lm(data = res_data, spRichness ~ res))


res_pval <- tibble(Parameter = as.character(),
                   pvalue = as.numeric(),
                   adj_r_squared = as.numeric())


for(i in 4:ncol(res_data %>% select(-c(lat,lon)))){

  Parameter <- colnames(res_data %>% select(-c(lat,lon)))[i] # select dependent parameter
  mod <- lm(paste(Parameter, "~", "res"), data = res_data)
  pvalue <- summary(mod)[4]$coefficients[8]
  adj_r_squared <- summary(mod)[9]$adj.r.squared

  temp <- as_tibble(cbind(Parameter, pvalue, adj_r_squared)) %>%
    mutate(pvalue = as.numeric(pvalue),
           adj_r_squared = as.numeric(adj_r_squared))

  res_pval <- res_pval %>%
    rbind(temp)
}

# view table of pvalues
res_pval

# view bar graph of pvalues with 0.05 cutoff
res_pval_plot <- res_pval %>%
  ggplot(aes(x = Parameter, y = pvalue)) +
  geom_col(fill = if_else(res_pval$pvalue < 0.05, "red", "darkgrey")) +
  geom_hline(yintercept = 0.05) +
  theme_bw() +
  labs(y = "p-value: ~ Residuals") +
  theme(axis.text.x = element_text(angle = 90))

#ggsave(here("Output", "PaperFigures", "sp_res_pval_supplemental.png"), res_pval_plot, width = 6, height = 5)
res_pval_plot



# PCA with residuals


# PCA and Permanova: survey location characteristics w/o substrate

## PCA

# extract parameters for the PCA, excluding substrate
pca_fulldata3 <- res_data %>%
  relocate(spR_residuals, .after=lon) %>%
  select(spR_residuals:N_percent,dist_to_seep_m:adj_CT_depth_cm) # avoid rugosity because rugosity residuals

# Run the PCA
pca_full3 <- prcomp(pca_fulldata3, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_full3 <- as_tibble(pca_full3$x[,1:2]) # scores for each sample location
PC_loadings_full3 <- as_tibble(pca_full3$rotation) %>% # loadings for each parameter
  bind_cols(labels = rownames(pca_full3$rotation))

# Put it with all the original data
pca_fulldata_all3 <- res_data %>%
  relocate(spR_residuals, .after=lon) %>%
  bind_cols(PC_scores_full3) %>%
  select(-meanRugosity)



## scores plot

p1_full3 <- pca_fulldata_all3 %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_point(shape = 21,
             color = "black",
             aes(fill = dist_to_seep_m),
             size = 4) +
  scale_fill_gradient(low = "red", high = "yellow") +
  #scale_fill_manual(values = mypalette) +
  theme_bw()+
  theme(legend.position = "top",
        #panel.background = element_rect(fill = "azure4"),
        panel.grid = element_blank()) +
  labs(title = "Varari Scores Plot", fill = "Distance to seep (m)") +
  geom_text_repel(aes(label = CowTagID),
                  size = 4)
# ggforce::geom_mark_ellipse(
#   aes(fill = V_cluster, label = V_cluster),
#   alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
#coord_cartesian(xlim = c(-4, 7), ylim = c(-6, 4)) +
#scale_shape_manual(values = c(22,16))+
#scale_colour_hue(l = 45)+ # l: luminance/lightness 0-100
#scale_fill_hue(l = 45)
p1_full3


## loadings plot

# loadings plot
p1_full_loading3 <- PC_loadings_full3 %>%
  ggplot(aes(x = PC1, y = PC2, label = labels)) + # labels are the parameters
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = PC1*10,
                   yend = PC2*10),
               arrow = arrow(length = unit(0.1,"cm")),
               color = "grey") +
  annotate("text",
           x = PC_loadings_full3$PC1*10+0.5,
           y = PC_loadings_full3$PC2*10+0.5,
           label = PC_loadings_full3$labels,
           size = 3) +
  coord_cartesian(xlim = c(-4, 5), ylim = c(-6, 5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Varari Loadings")

p1_full_loading3


### Patch PCA plots
VarariPCA3 <- p1_full3 + p1_full_loading3 +
  patchwork::plot_annotation("Varari PCA and Loadings (w/ Residuals)",
                             theme = theme(plot.title = element_text(size = rel(1.5),
                                                                     face = "bold",
                                                                     hjust = 0.5,
                                                                     margin = margin(t = 10, b = 20,
                                                                                     unit = "pt"))))

ggsave(plot = VarariPCA3, filename = here("Output","PaperFigures","VarariPCA_meta_residuals.png"), width = 14, height = 8)
VarariPCA3


###########################################################
###########################################################
### Calculate FER
FE_data <- comp %>%
  left_join(taxa) %>%
  unite(Morphology,
        #Life_Span,
        #Max_size,
        #Growth_rate,
        Zooxanthellate,
        Calcification,
        Energetic_Resource,
        Feeding_Mode,
        sep = ",",
        remove = T,
        col = "FE") %>%
  relocate(FE, .after = Taxa) %>%
  group_by(Location,CowTagID) %>%
  count(FE) %>%
  mutate(n = 1,
         FER = sum(n)) %>%
  left_join(comp) %>%
  select(Location, CowTagID, FER,
         Salinity:N_percent,dist_to_seep_m:meanRugosity) %>%
  distinct()


# Normalize to substrate and rugosity

## Richness ~ Rugosity

## plot
groupplot(mydata = FE_data, x = meanRugosity, y = FER, fw = NA) +
  geom_text_repel(aes(label = CowTagID)) +
  scale_color_manual(values = mypalette)

# A residual is the difference between an observed value and a predicted value in a regression model. Residual = Observed value – Predicted value
# An observation has a positive residual if its value is greater than the predicted value made by the regression line and a negative residual if its value is less than the predicted value made by the regression line.

#fit model
resMod <- lm(FER ~ meanRugosity, data=FE_data %>% filter(CowTagID != "VSEEP"))

#view model summary
summary(resMod) # 0.03

#calculate the standardized residuals
res <- residuals(resMod)

#view the standardized residuals
res

#column bind standardized residuals back to original data frame
res_data <- FE_data %>%
  filter(CowTagID != "VSEEP") %>%
  cbind(res) %>%
  relocate(res, .after = CowTagID) %>%
  rename(FER_residuals = res)



## Residuals regressions

anova(lm(data = res_data, spRichness ~ res))


res_pval <- tibble(Parameter = as.character(),
                   pvalue = as.numeric(),
                   adj_r_squared = as.numeric())


for(i in 4:ncol(res_data %>% select(-c(lat,lon)))){

  Parameter <- colnames(res_data %>% select(-c(lat,lon)))[i] # select dependent parameter
  mod <- lm(paste(Parameter, "~", "res"), data = res_data)
  pvalue <- summary(mod)[4]$coefficients[8]
  adj_r_squared <- summary(mod)[9]$adj.r.squared

  temp <- as_tibble(cbind(Parameter, pvalue, adj_r_squared)) %>%
    mutate(pvalue = as.numeric(pvalue),
           adj_r_squared = as.numeric(adj_r_squared))

  res_pval <- res_pval %>%
    rbind(temp)
}

# view table of pvalues
res_pval

# view bar graph of pvalues with 0.05 cutoff
res_pval %>%
  ggplot(aes(x = Parameter, y = pvalue)) +
  geom_col(fill = if_else(res_pval$pvalue < 0.05, "red", "darkgrey")) +
  geom_hline(yintercept = 0.05) +
  theme_bw() +
  labs(y = "p-value: ~ Residuals") +
  theme(axis.text.x = element_text(angle = 90))




# PCA with residuals


# PCA and Permanova: survey location characteristics w/o substrate

## PCA

# extract parameters for the PCA, excluding substrate
pca_fulldata3 <- res_data %>%
  relocate(spR_residuals, .after=lon) %>%
  select(spR_residuals:N_percent,dist_to_seep_m:adj_CT_depth_cm) # avoid rugosity because rugosity residuals

# Run the PCA
pca_full3 <- prcomp(pca_fulldata3, scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores_full3 <- as_tibble(pca_full3$x[,1:2]) # scores for each sample location
PC_loadings_full3 <- as_tibble(pca_full3$rotation) %>% # loadings for each parameter
  bind_cols(labels = rownames(pca_full3$rotation))

# Put it with all the original data
pca_fulldata_all3 <- res_data %>%
  relocate(spR_residuals, .after=lon) %>%
  bind_cols(PC_scores_full3) %>%
  select(-meanRugosity)



## scores plot

p1_full3 <- pca_fulldata_all3 %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_point(shape = 21,
             color = "black",
             aes(fill = dist_to_seep_m),
             size = 4) +
  scale_fill_gradient(low = "red", high = "yellow") +
  #scale_fill_manual(values = mypalette) +
  theme_bw()+
  theme(legend.position = "top",
        #panel.background = element_rect(fill = "azure4"),
        panel.grid = element_blank()) +
  labs(title = "Varari Scores Plot", fill = "Distance to seep (m)") +
  geom_text_repel(aes(label = CowTagID),
                  size = 4)
# ggforce::geom_mark_ellipse(
#   aes(fill = V_cluster, label = V_cluster),
#   alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
#coord_cartesian(xlim = c(-4, 7), ylim = c(-6, 4)) +
#scale_shape_manual(values = c(22,16))+
#scale_colour_hue(l = 45)+ # l: luminance/lightness 0-100
#scale_fill_hue(l = 45)
p1_full3


## loadings plot

# loadings plot
p1_full_loading3 <- PC_loadings_full3 %>%
  ggplot(aes(x = PC1, y = PC2, label = labels)) + # labels are the parameters
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = PC1*10,
                   yend = PC2*10),
               arrow = arrow(length = unit(0.1,"cm")),
               color = "grey") +
  annotate("text",
           x = PC_loadings_full3$PC1*10+0.5,
           y = PC_loadings_full3$PC2*10+0.5,
           label = PC_loadings_full3$labels,
           size = 3) +
  coord_cartesian(xlim = c(-4, 5), ylim = c(-6, 5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Varari Loadings")

p1_full_loading3


### Patch PCA plots
VarariPCA3 <- p1_full3 + p1_full_loading3 +
  patchwork::plot_annotation("Varari PCA and Loadings (w/ Residuals)",
                             theme = theme(plot.title = element_text(size = rel(1.5),
                                                                     face = "bold",
                                                                     hjust = 0.5,
                                                                     margin = margin(t = 10, b = 20,
                                                                                     unit = "pt"))))

