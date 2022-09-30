### CREATE REGRESSION PLOTS FOR PARAMETERS RELATIVE TO DEPTH OF CT LOGGERS ###
### Created by Danielle Barnas
### Created on 9/29/2022


#### LOAD LIBRARIES ####
library(tidyverse)
library(lubridate)
library(ggmap)
library(PNWColors)
library(here)
library(curl) # pull data from url
library(tidytext) # arrange data in a facet_wrap
library(ggrepel) # labeling
library(ggpmisc) # add R2 to figures

#### READ IN DATA ####
## biogeochem data processed in Nutrient_Processing script
Full_data<-read_csv(here("Data","Biogeochem","AugNutrient_Processed_CV.csv"))
MaxMin_data<-read_csv(here("Data","Biogeochem","AugNutrient_Processed_MaxMin.csv"))


#### CHECK CORRELATIONS ####
## ~ silicate
Full_data %>%
  filter(CowTagID != "VSEEP" &
         CowTagID != "CSEEP") %>%
  pivot_longer(cols = c(adj_CT_depth_cm ,Salinity:Phosphate_umolL, NN_umolL:N_percent),
               names_to = "Parameters",
               values_to = "Values") %>%
  ggplot(aes(x = Silicate_umolL,
             y = Values,
             color = Location),
         color = "black") +
  geom_point() +
  geom_smooth(method = "lm", color = "black", aes(group = Location)) +
  facet_wrap(~Parameters, scales = "free") +
  labs(title = "All parameters ~ Silicate",
       x = "CV of Silicate (umolL)",
       y = "Value Coefficient of Variation") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))



## ~ depth
Full_data %>%
  filter(CowTagID != "VSEEP" &
         CowTagID != "CSEEP") %>%
  pivot_longer(cols = c(Salinity:N_percent),
               names_to = "Parameters",
               values_to = "Values") %>%
  ggplot(aes(x = adj_CT_depth_cm,
             y = Values,
             color = Location
             ),
  color = "black") +
  geom_point() +
  geom_smooth(method = "lm", color = "black", aes(group = Location)) +
  facet_wrap(~Parameters, scales = "free") +
  labs(title = "All parameters ~ Depth",
       x = "Depth to CT logger (cm)",
       y = "Value Coefficient of Variation") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))


# create function to calculate rsq value for each parameter in
# relation to CT depth (aka water sample depth)
rsq_fun <- function(data = Full_data, x = "adj_CT_depth_cm", y) {
  fun_data <- data %>%
    filter(CowTagID != 'VSEEP' & CowTagID != 'CSEEP') # remove seeps to not skew

  Vmodel <- lm(data = as_tibble(fun_data[fun_data$Location == 'Varari',]), paste0({{y}}, "~", {{x}}))
  Cmodel <- lm(data = as_tibble(fun_data[fun_data$Location == 'Cabral',]), paste0({{y}}, "~", {{x}}))

  Vsum <- summary(Vmodel)$r.squared
  Vsum <- as_tibble(Vsum) %>%
    mutate(Location = "Varari")

  Csum <- summary(Cmodel)$r.squared
  Csum <- as_tibble(Csum) %>%
    mutate(Location = "Cabral")


  sums <- rbind(Vsum, Csum) %>%
    mutate(Parameter = y)

  return(sums)
}



# set empty df
rsq_data <- tibble(Location = as.character(),
                   Parameter = as.character(),
                   value = as.numeric())
# empty list for figures
q <- list()

# loop through parameters
for(i in 6:ncol(Full_data)){

  yparam = colnames(Full_data[,i])
  xparam = "adj_CT_depth_cm"

  sums <- rsq_fun(data = Full_data, x = xparam, y = yparam)

  rsq_data <- rsq_data %>% full_join(sums)
}


  # full data to long form and join with rsq df
  reduc_data <- Full_data %>%
    filter(CowTagID != "VSEEP" &
             CowTagID != "CSEEP")
  Names <- colnames(Full_data)

for(i in 6:ncol(reduc_data)){

  # select parameter; as.name calls the identifier name, not a character string
  param <- as.name(Names[i])

  # select Rsquared based on param
  param_rsq_data <- rsq_data %>%
    filter(Parameter == param)

  # graph relationship
  p <- reduc_data %>%
    ggplot(aes(x = adj_CT_depth_cm,
               y = reduc_data[[Names[i]]],
               label = CowTagID,
               color = Location),
           color = "black") +
    geom_point() +
    #geom_text_repel(size = 2, color = "black") +
    geom_smooth(method = "lm", color = "black", aes(group = Location)) +
    labs(x = "Depth to CT logger (cm)",
         y = "Coefficient of Variation",
         title = paste("CV of", Names[i], "vs Depth (cm)"),
         caption = paste("Varari R-squared:",round(param_rsq_data[param_rsq_data$Location == "Varari",]$value, 2), "\n",
                         "Cabral R-squared:", round(param_rsq_data[param_rsq_data$Location == "Cabral",]$value, 2)
         )) +
    theme_bw()

  # save plot in list
  q[[i]] <- p

}


# Save all plots in a single pdf
pdf(here("Output","Param_vs_Depth.pdf"), onefile = TRUE)
for (i in 6:length(q)) {
  tplot <- q[[i]]
  print(tplot)
}
dev.off()






##### USING MAXIMUM AND MINIMUM VALUES ####

# create function to calculate rsq value for each parameter in
# relation to CT depth (aka water sample depth)
rsq_fun_mm <- function(data = MaxMin_data, x = "adj_CT_depth_cm", y) {
  ## maximum values
  maxfun_data <- data %>%
    filter(CowTagID != 'VSEEP' & CowTagID != 'CSEEP') %>%  # remove seeps to not skew
    filter(MaxMin == "Maximum")

  maxVmodel <- lm(data = as_tibble(maxfun_data[maxfun_data$Location == 'Varari',]), paste0({{y}}, "~", {{x}}))
  maxCmodel <- lm(data = as_tibble(maxfun_data[maxfun_data$Location == 'Cabral',]), paste0({{y}}, "~", {{x}}))

  maxVsum <- summary(maxVmodel)$r.squared
  maxVsum <- as_tibble(maxVsum) %>%
    mutate(Location = "Varari") %>%
    mutate(MaxMin = "Maximum")

  maxCsum <- summary(maxCmodel)$r.squared
  maxCsum <- as_tibble(maxCsum) %>%
    mutate(Location = "Cabral") %>%
    mutate(MaxMin = "Maximum")

  maxsums <- rbind(maxVsum, maxCsum) %>%
    mutate(Parameter = y)

  ## minimum values
  minfun_data <- data %>%
    filter(CowTagID != 'VSEEP' & CowTagID != 'CSEEP') %>%  # remove seeps to not skew
    filter(MaxMin == "Minimum")

  minVmodel <- lm(data = as_tibble(minfun_data[minfun_data$Location == 'Varari',]), paste0({{y}}, "~", {{x}}))
  minCmodel <- lm(data = as_tibble(minfun_data[minfun_data$Location == 'Cabral',]), paste0({{y}}, "~", {{x}}))

  minVsum <- summary(minVmodel)$r.squared
  minVsum <- as_tibble(minVsum) %>%
    mutate(Location = "Varari") %>%
    mutate(MaxMin = "Minimum")

  minCsum <- summary(minCmodel)$r.squared
  minCsum <- as_tibble(minCsum) %>%
    mutate(Location = "Cabral") %>%
    mutate(MaxMin = "Minimum")

  minsums <- rbind(minVsum, minCsum) %>%
    mutate(Parameter = y)

  sums <- rbind(maxsums, minsums)


  return(sums)
}


# set empty df
rsq_data_mm <- tibble(Location = as.character(),
                   Parameter = as.character(),
                   value = as.numeric(),
                   MaxMin = as.character())
# empty list for figures
r <- list()

# loop through parameters
for(i in 7:ncol(MaxMin_data)){

  yparam = colnames(MaxMin_data[,i])
  xparam = "adj_CT_depth_cm"

  sums <- rsq_fun_mm(data = MaxMin_data, x = xparam, y = yparam)

  rsq_data_mm <- rsq_data_mm %>% full_join(sums)
}


# full data to long form and join with rsq df
reduc_data_mm <- MaxMin_data %>%
  filter(CowTagID != "VSEEP" &
           CowTagID != "CSEEP")
Names_mm <- colnames(MaxMin_data)

for(i in 7:ncol(reduc_data_mm)){

  # select parameter; as.name calls the identifier name, not a character string
  param <- as.name(Names_mm[i])

  # select Rsquared based on param
  param_rsq_data <- rsq_data_mm %>%
    filter(Parameter == param)

  # graph relationship
  p <- reduc_data_mm %>%
    ggplot(aes(x = adj_CT_depth_cm,
               y = reduc_data_mm[[Names_mm[i]]],
               label = CowTagID,
               color = Location),
           color = "black") +
    geom_point() +
    #geom_text_repel(size = 2, color = "black") +
    geom_smooth(method = "lm", color = "black", aes(group = Location)) +
    labs(x = "Depth to CT logger (cm)",
         y = "Coefficient of Variation",
         title = paste("CV of", Names_mm[i], "vs Depth (cm)"),
         caption = paste("Varari R-squared (Max val):",round(param_rsq_data[param_rsq_data$Location == "Varari" & param_rsq_data$MaxMin == "Maximum",]$value, 2),
                         "Varari R-squared (Min val):",round(param_rsq_data[param_rsq_data$Location == "Varari" & param_rsq_data$MaxMin == "Minimum",]$value, 2), "\n",
                         "Cabral R-squared (Max val):", round(param_rsq_data[param_rsq_data$Location == "Cabral" & param_rsq_data$MaxMin == "Maximum",]$value, 2),
                         "Cabral R-squared (Min val):", round(param_rsq_data[param_rsq_data$Location == "Cabral" & param_rsq_data$MaxMin == "Minimum",]$value, 2)
         )) +
    facet_wrap(~MaxMin) +
    theme_bw()

  # save plot in list
  r[[i]] <- p

}


# Save all plots in a single pdf
pdf(here("Output","MaxMin_Param_vs_Depth.pdf"), onefile = TRUE)
for (i in 6:length(r)) {
  tplot <- r[[i]]
  print(tplot)
}
dev.off()


