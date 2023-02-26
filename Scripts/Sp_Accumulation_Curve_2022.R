### Species Accumulation Curve from June 2022 ###
### Created by Danielle Barnas
### Created on November 13, 2021

library(tidyverse)
library(here)


### READ IN DATA
survey<-read_csv(here("Data","Surveys","Species_Composition_2022.csv"))

### ANALYSIS

## Total quads by species/taxa counts
## species accumulation curve shows how many new species you get by doing more surveys

# Remove abiotic from "Taxa"
abiotic <- survey %>%
  filter(Taxa %in% c("Sand", "Rubble", "Bare Rock", "Bare Rock - exposed"))

Vquad <- survey %>%
  anti_join(abiotic) %>%
  filter(Location == "Varari") %>%
  filter(CowTagID != "VSEEP" & CowTagID != "V13") %>%
         #| Location == "Varari_Maya") %>% # remove maya sites while not used for analysis
  select(CowTagID, Taxa)

# remove any duplicate species from dataframe to only display one of each across full survey
Vquad_order <- Vquad[order(Vquad[,'CowTagID'],Vquad[,'Taxa']),]
Vquad_order <- Vquad_order[!duplicated(Vquad_order$Taxa),]

# rejoin with df containing all quads
Vquad <- Vquad %>%
  select(CowTagID) %>%
  left_join(Vquad_order) %>%
  distinct() %>%
  arrange(CowTagID)

# create column "n" containing a value of 1 next to each species
Vquad <- Vquad %>%
  group_by(CowTagID) %>%
  dplyr::count(Taxa) %>%
  ungroup()

# replace values of 1 with values of 0 for all NA species (indicating no new species for a later survey)
Vquad <- Vquad %>%
  mutate(n = if_else(is.na(Taxa) == TRUE, true = (n = 0), false = (n = 1)))

# add total species per quadrat survey
Vquad <- Vquad %>%
  group_by(CowTagID) %>%
  summarise(sp.sum = sum(n)) %>%
  ungroup() %>%
  arrange(desc(sp.sum))


# add column to enumerate each quadrat survey
Vquad$row.num=seq.int(nrow(Vquad)) # add column of row numbers aka quad number

# sequentially add species richness as number of surveys increases
Vquad<-Vquad %>%
  mutate(sp_accumulation = cumsum(sp.sum))

plot1 <- Vquad %>%
  ggplot(aes(y = sp_accumulation, x = row.num)) +
  geom_point() +
  labs(x = "Total Quadrat Surveys",
       y = "Species Richness",
       title = "Varari Species Accumulation Curve",
       subtitle = "June 2022") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "#e8b4b8"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_smooth(color = "black")
plot1

# ggsave(here("Output", "Species_Accumulation_Varari.png"), plot1, height = 5, width = 5, device = "png")
# ggsave(here("Output", "Species_Accumulation_Varari.pdf"), plot1, height = 5, width = 5, device = "pdf")


#############################################################################################
#############################################################################################

Cquad <- survey %>%
  anti_join(abiotic) %>%
  filter(Location == "Cabral") %>%
  select(CowTagID, Taxa)

# remove any duplicate species from dataframe to only display one of each across full survey
Cquad_order <- Cquad[order(Cquad[,'CowTagID'],Cquad[,'Taxa']),]
Cquad_order <- Cquad_order[!duplicated(Cquad_order$Taxa),]

# rejoin with df containing all quads
Cquad <- Cquad %>%
  select(CowTagID) %>%
  left_join(Cquad_order) %>%
  distinct() %>%
  arrange(CowTagID)

# create column "n" containing a value of 1 next to each species
Cquad <- Cquad %>%
  group_by(CowTagID) %>%
  count(Taxa) %>%
  ungroup()

# replace values of 1 with values of 0 for all NA species (indicating no new species for a later survey)
Cquad <- Cquad %>%
  mutate(n = if_else(is.na(Taxa) == TRUE, true = (n = 0), false = (n = 1)))

# add total species per quadrat survey
Cquad <- Cquad %>%
  group_by(CowTagID) %>%
  summarise(sp.sum = sum(n)) %>%
  ungroup() %>%
  arrange(desc(sp.sum))


# add column to enumerate each quadrat survey
Cquad$row.num=seq.int(nrow(Cquad)) # add column of row numbers aka quad number

# sequentially add species richness as number of surveys increases
Cquad<-Cquad %>%
  mutate(sp_accumulation = cumsum(sp.sum))

plot2 <- Cquad %>%
  ggplot(aes(y = sp_accumulation, x = row.num)) +
  geom_point() +
  labs(x = "Total Quadrat Surveys",
       y = "Species Richness",
       title = "Cabral Species Accumulation Curve",
       subtitle = "June 2022") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "azure4"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_smooth(color = "black")
plot2

# ggsave(here("Output", "Species_Accumulation_Cabral.png"), plot2, height = 10, width = 10, device = "png")
# ggsave(here("Output", "Species_Accumulation_Cabral.pdf"), plot2, height = 10, width = 10, device = "pdf")


