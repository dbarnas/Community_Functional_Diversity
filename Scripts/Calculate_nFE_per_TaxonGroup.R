#### Calculate total FE per taxonomic group

### Created by Danielle Barnas
### Created on March 6, 2023


library(tidyverse)
library(here)
library(PNWColors)
library(patchwork)

mypalette <- pnw_palette(name = "Bay", n = 12)
smallpalette <- pnw_palette(name = "Bay", n = 4)

ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))
spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv")))
fes_taxa <- read_csv(here("Data", "Surveys", "Distinct_Taxa.csv"))
fes_taxa <- fes_taxa %>%
  select(Taxa, Taxon_Group)
dist <- read_csv(here("Data","Full_metadata.csv")) %>%
  filter(Location == "Varari",
         CowTagID != "V13",
         CowTagID != "VSEEP") %>%
  select(CowTagID, dist_to_seep_m) %>%
  arrange(dist_to_seep_m)
dist <- dist[1:19,]


Full_data <- ab.sgd %>%
  pivot_longer(cols = 2:ncol(ab.sgd), names_to = "Taxa", values_to = "pCover") %>%
  filter(pCover > 0) %>%
  left_join(fes_taxa) %>%
  left_join(spe_fes.sgd) %>%
  left_join(fes_traits.sgd)

ctlevels <- dist$CowTagID
Full_data$CowTagID <- factor(Full_data$CowTagID, levels = ctlevels)

Full_data %>%
  distinct(Taxon_Group, FE) %>%
  dplyr::count(Taxon_Group)
