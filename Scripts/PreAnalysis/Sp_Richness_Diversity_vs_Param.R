### CORRELATING SPECIES RICHNESS AND DIVERSITY TO SGD PARAMETERS
### Created by Danielle Barnas
### Created on 9/30/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)


#### READ IN DATA ####
rich <- read_csv(here("Data", "Surveys", "Species_Richness.csv"))
div <- read_csv(here("Data", "Surveys", "Species_Diversity.csv"))
envi <- read_csv(here("Data", "Biogeochem", "AugNutrient_Processed_CV.csv"))


#### JOIN DF AND CLEAN ####
Full_data <- full_join(rich, div) %>%
  drop_na(Location) %>%
  full_join(envi)


#### REGRESSION ####

# create function to calculate r2 value for species richness or diversity vs each parameter
rsq_fun <- function(data = Full_data, x, y, metric) {
  fun_data <- data %>%
    filter(Location != "Cabral") %>% # for now
    filter(CowTagID != 'VSEEP' & CowTagID != 'CSEEP') # remove seeps to not skew

  Vmodel <- lm(data = as_tibble(fun_data[fun_data$Location == 'Varari',]), paste0({{y}}, "~", {{x}}))
  #Cmodel <- lm(data = as_tibble(fun_data[fun_data$Location == 'Cabral',]), paste0({{y}}, "~", {{x}}))

  Vsum <- summary(Vmodel)$r.squared
  Vsum <- as_tibble(Vsum) %>%
    #mutate(Location = "Varari") %>%
    mutate(Parameter = x,
           Metric = metric)

  # Csum <- summary(Cmodel)$r.squared
  # Csum <- as_tibble(Csum) %>%
  #   mutate(Location = "Cabral")


  # sums <- rbind(Vsum, Csum) %>%
  #   mutate(Parameter = y)

  return(Vsum)
}


# set empty df
rsq_rich <- tibble(Parameter = as.character(),
                   value = as.numeric(),
                   Metric = as.character())
rsq_div <- tibble(Parameter = as.character(),
                  value = as.numeric(),
                  Metric = as.character())

# loop through parameters
# richness
for(i in 8:ncol(Full_data)){

  xparam = colnames(Full_data[,i])
  yparam = "spR"

  sums <- rsq_fun(data = Full_data, x = xparam, y = yparam, metric = "Richness")

  rsq_rich <- rsq_rich %>% full_join(sums)
}

# diversity
for(i in 8:ncol(Full_data)){

  xparam = colnames(Full_data[,i])
  yparam = "ShannonDivSpecies"

  sums <- rsq_fun(data = Full_data, x = xparam, y = yparam, metric = "Diversity")

  rsq_div <- rsq_div %>% full_join(sums)
}

# join richness and diversity r2 dfs
rsq_data <- full_join(rsq_rich, rsq_div)

#### GRAPH CORRELATIONS ####

# Species Richness
NamesR <- colnames(Full_data)

# long form dataframe
long_data <- Full_data %>%
  pivot_longer(cols = adj_CT_depth_cm:N_percent, names_to = "Parameter", values_to = "Values")


# set empty list for figures
r <- list()

for(i in 8:ncol(Full_data)){

  # select parameter; as.name calls the identifier name, not a character string
  param <- as.name(NamesR[i])

  # select parameter
  param_data <- long_data %>%
    filter(Location == "Varari",
           Parameter == param,
           CowTagID != "VSEEP")
  r2data <- param_data %>%
    left_join(rsq_data) %>%
    rename(r2 = value)


  # graph relationships
  p <- param_data %>%
    ggplot(aes(x = Values,
               y = spR,
               label = CowTagID),
           color = "black") +
    geom_point() +
    geom_text_repel(size = 2, color = "black") +
    geom_smooth(method = "lm", color = "black") +
    labs(x = "{.x}",
         y = "Species Richness",
         title = paste("Species Richness vs", NamesR[i]),
         caption = paste("R-squared:",round(r2data[r2data$Metric == "Richness",]$r2, 2)
                         #"\n",
                         #"Cabral R-squared (Max val):", round(rsq_data[rsq_data$Location == "Cabral" & rsq_data$MaxMin == "Richness",]$value, 2)
         )) +
    theme_bw()


  q <- param_data %>%
    ggplot(aes(x = Values,
               y = ShannonDivSpecies,
               label = CowTagID),
           color = "black") +
    geom_point() +
    geom_text_repel(size = 2, color = "black") +
    geom_smooth(method = "lm", color = "black") +
    labs(x = "{.x}",
         y = "Species Diversity (H')",
         title = paste("Species Diversity vs", NamesR[i]),
         caption = paste("R-squared:",round(r2data[r2data$Metric == "Diversity",]$r2, 2)
           )) +
    theme_bw()

  # save plot in list
  r[[i]] <- p + q

}


# Save all plots in a single pdf
pdf(here("Output","Species_vs_Param.pdf"), onefile = TRUE)
for (i in 8:length(r)) {
  tplot <- r[[i]]
  print(tplot)
}
dev.off()


