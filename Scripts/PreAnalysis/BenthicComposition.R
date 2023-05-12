### Summer 2022 Benthic Community Composition Surveys

#### LIBRARIES ####
library(tidyverse)
library(here)

library(patchwork)
library(PNWColors)

library(ggmap)
library(viridis)
library(maptools)
#library(kriging)
#library(ggnewscale)
#library(wql)
#library(glue)
#library(gridExtra)


#### READ IN DATA ####
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv")) %>% filter(Location == "Varari")
#richness <- read_csv(here("Data", "Surveys", "Species_Richness.csv"))
#diversity <- read_csv(here("Data", "Surveys", "Species_Diversity.csv"))
func <- read_csv(here("Data", "Surveys", "Distinct_Taxa.csv"))

meta <- read_csv(here("Data", "Full_Metadata.csv"))
#AugChem <- read_csv(here("Data","Biogeochem","AugNutrient_Processed_CV.csv")) #%>% filter(Location == "Varari")
Turb22 <- read_csv(here("Data","Biogeochem","July2022", "Turb_NC.csv"))

# create color palette for plotting
mypalette <- pnw_palette(name="Bay", n=12)


#### CLEAN AND COMPILE DATA ####

### relative abundance
# species
Sp.abundance <- survey %>%
  # filter(Taxa != 'Bare Rock',
  #        Taxa != 'Sand',
  #        Taxa != 'Rubble') %>%
  group_by(Location,CowTagID) %>% # group by species at each cowtag location
  mutate(total = sum(SpeciesCounts)) %>%
  mutate(pCoverSpecies = SpeciesCounts / total * 100) %>% # calculate percent cover of each species by cowtag location
  select(Location, CowTagID, Taxa, SpeciesCounts, pCoverSpecies)
# genera
G.abundance <- survey %>%
  left_join(func) %>%
  select(Location, CowTagID, Taxa, Genus, SpeciesCounts) %>%
  mutate(Genus = replace_na(Genus, "Abiotic")) %>% # include hard substrate as a category
  group_by(Location,CowTagID,Genus) %>%
  mutate(GenusCounts = sum(SpeciesCounts)) %>% # calculate total genus counts per cowtag location
  ungroup() %>%
  group_by(Location,CowTagID) %>%
  mutate(total = sum(SpeciesCounts),  # total counts
         pCoverGenus = GenusCounts / total * 100) %>%
  select(Location, CowTagID, Genus, GenusCounts, pCoverGenus) %>%
  distinct()
# broader taxon groups
T.abundance <- survey %>%
  left_join(func) %>%
  select(Location, CowTagID, Taxa, SpeciesCounts, Taxon_Group) %>%
  mutate(Taxon_Group = if_else(Taxon_Group == "Sand", "Abiotic", Taxon_Group),
         Taxon_Group = if_else(Taxon_Group == "Hard Substrate", "Abiotic", Taxon_Group)) %>%
  group_by(Location,CowTagID,Taxon_Group) %>%
  mutate(TaxonCounts = sum(SpeciesCounts)) %>% # calculate total taxon group counts per cowtag location
  ungroup() %>%
  group_by(Location,CowTagID) %>%
  mutate(total = sum(SpeciesCounts),  # total counts
         pCoverTaxon = TaxonCounts / total * 100) %>%
  select(Location, CowTagID, Taxon_Group, TaxonCounts, pCoverTaxon) %>%
  distinct()

# species richness
richness <- survey %>%
  filter(Taxa != 'Bare Rock',
         Taxa != 'Sand',
         Taxa != 'Rubble') %>%
  group_by(Location,CowTagID) %>%
  count(Taxa) %>%  # counts every time a distinct species appears at a cowtag location
  mutate(n=1) %>%  # make sure all species counts are 1 per species
  mutate(spRich = sum(n)) %>% # calculate species richness
  select(Location, CowTagID, spRich) %>%
  distinct() %>% # remove redundant rows
  left_join(meta)  # join all metadata


# order CowTagID's by distance to seep (m)
orderdist <- richness %>%
  mutate_all( ~ if_else(dist_to_seep_m < 0, -1*dist_to_seep_m + 200, dist_to_seep_m)) %>%  # add large distance to 13 to make sure it comes out at far end of ordering (becuase 13 is upstream of seep)
  arrange(dist_to_seep_m) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID)))
distLevels <- paste(orderdist$CowTagID)


# assign order to factor levels
richness$CowTagID <- factor(richness$CowTagID, levels = distLevels)
Sp.abundance$CowTagID <- factor(Sp.abundance$CowTagID, levels = distLevels)
G.abundance$CowTagID <- factor(G.abundance$CowTagID, levels = distLevels)
T.abundance$CowTagID <- factor(T.abundance$CowTagID, levels = distLevels)
#levels(richness$CowTagID) # check



#### VISUALIZATIONS ####


# SPECIES RICHNESS BARPLOT
richness %>%
  ggplot(aes(x = CowTagID, y = spRich)) +
  geom_col(fill = "palevioletred4", color = "black") +
  labs(y = "Speices Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))



# SPECIES ABUNDANCE STACKED BARPLOT
Sp.abundance %>%
  ggplot(aes(x = CowTagID, y = pCoverSpecies, fill = Taxa)) +
  geom_col(position = "stack")

G.abundance %>%
  ggplot(aes(x = CowTagID, y = pCoverGenus, fill = Genus)) +
  geom_col(position = "stack")

T.abundance %>%
  ggplot(aes(x = CowTagID, y = pCoverTaxon, fill = Taxon_Group)) +
  geom_col(position = "stack")




facet_taxon_abundance <- T.abundance %>%
  select(Location, CowTagID, Taxon_Group, pCoverTaxon) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = pCoverTaxon,
             fill = Taxon_Group)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 8),
        axis.text.y = element_text(size = 12)
        ) +
  labs(x = "Survey Location",
       y = "% Benthic Cover",
       fill = "Taxon Groups") +
  scale_y_continuous(breaks = c(0, 25, 50, 100)) +
  scale_fill_manual(values = mypalette) +
  facet_wrap(~ Taxon_Group)
facet_taxon_abundance

ggsave(here("Output", "PaperFigures", "facet_taxon_abundance.png"), facet_taxon_abundance, device = "png")


stacked_taxon_abundance <- T.abundance %>%
  select(Location, CowTagID, Taxon_Group, pCoverTaxon) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = pCoverTaxon,
             fill = Taxon_Group)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 8),
        axis.text.y = element_text(size = 12)
  ) +
  labs(x = "Survey Location",
       y = "% Benthic Cover",
       fill = "Taxon Groups") +
  scale_y_continuous(breaks = c(0, 25, 50, 100)) +
  scale_fill_manual(values = mypalette)
stacked_taxon_abundance

ggsave(here("Output", "PaperFigures", "stacked_taxon_abundance.png"), stacked_taxon_abundance, device = "png")

stacked_taxon_bio_abundance <- T.abundance %>%
  select(Location, CowTagID, Taxon_Group, pCoverTaxon) %>%
  distinct() %>%
  filter(Taxon_Group != "Abiotic") %>%
  group_by(CowTagID) %>%
  mutate(totalsum = sum(pCoverTaxon),
         newpercent = pCoverTaxon / totalsum * 100) %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = newpercent,
             fill = Taxon_Group)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x.bottom = element_text(size = 8),
        axis.text.y = element_text(size = 12)
  ) +
  labs(x = "Survey Location",
       y = "% Benthic Cover",
       fill = "Taxon Groups") +
  scale_y_continuous(breaks = c(0, 25, 50, 100)) +
  scale_fill_manual(values = mypalette)
stacked_taxon_bio_abundance

ggsave(here("Output", "stacked_taxon_bio_abundance.png"), stacked_taxon_bio_abundance, device = "png")


####
