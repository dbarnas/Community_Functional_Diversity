

#### LIBRARIES ####
library(tidyverse)
library(here)
library(tidytext)



#### READ IN DATA ####
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))

### Both Survey locations
taxa_a <- read_csv(here("Data", "Surveys", "Distinct_Taxa.csv"))

a <- survey %>%
  distinct(Taxa)

### Join with functional trait dataframe
taxa_a <- taxa_a %>%
  full_join(a)

# view new taxa added
taxa_a %>% anti_join(a)

### Write over csv
write_csv(taxa_a, here("Data", "Surveys", "Distinct_Taxa.csv"))


### Varari only
taxa_b <- read_csv(here("Data", "Surveys", "Distinct_Varari_Taxa.csv"))

b <- survey %>%
  filter(Location != "Cabral") %>%
  distinct(Taxa)

### Join with functional trait dataframe
taxa_b <- taxa_b %>%
  full_join(b)

# view new taxa added
taxa_b %>% anti_join(b)

### Write over csv
write_csv(taxa_b, here("Data", "Surveys", "Distinct_Varari_Taxa.csv"))
