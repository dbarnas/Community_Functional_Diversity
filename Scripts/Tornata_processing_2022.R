### July 2022 T. ornata nitrogen laoding cleaning script
### Created by Danielle Barnas
### Created on 10/13/2022

#### LOAD LIBRARIES ####
library(tidyverse)
library(here)


#### BRING IN DATA ####
raw <- read_csv(here("Data","Biogeochem","July2022","RawBarnasTornata.csv"))
ref <- read_csv(here("Data","Biogeochem","Nutrient_Processed_MaxMin.csv"))


turb <- raw %>%
  mutate(C_N = C_ug / N_ug,
         N_percent = (N_ug/1000) / Weight_mg * 100) %>%
  select(-c(Weight_mg, C_ug, N_ug))
turb
ref %>%  # to check if new values are in similar range
  select(colnames(turb))

write_csv(turb, here("Data","Biogeochem","July2022","Turb_NC.csv"))
