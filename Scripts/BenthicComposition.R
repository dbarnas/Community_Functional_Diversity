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
survey <- read_csv(here("Data","Surveys","Species_Composition_2022.csv")) #%>% filter(Location == "Varari")
richness <- read_csv(here("Data", "Surveys", "Species_Richness.csv"))
diversity <- read_csv(here("Data", "Surveys", "Species_Diversity.csv"))
func <- read_csv(here("Data", "Surveys", "Distinct_Taxa.csv"))

meta <- read_csv(here("Data", "Surveys", "Survey_Metadata.csv"))
depth <- read_csv(here("Data","Adj_Sandwich_Depth.csv"))
dist <- read_csv(here("Data", "Plate_Distance_to_Seep.csv"))
AugChem <- read_csv(here("Data","Biogeochem","AugNutrient_Processed_CV.csv")) #%>% filter(Location == "Varari")
Turb22 <- read_csv(here("Data","Biogeochem","July2022", "Turb_NC.csv"))

# create color palette for plotting
mypalette <- pnw_palette(name="Bay", n=11)


#### CLEAN AND COMPILE DATA ####
# remove unnecessary columns
depth <- depth %>%
  select(Location, CowTagID, adj_CT_depth_cm)
dist <- dist %>%
  select(CowTagID, lat, lon, dist_to_seep_m) %>%
  mutate(dist_to_seep_m = if_else(CowTagID == "V13", (-1 * dist_to_seep_m), dist_to_seep_m))

Full_meta <- meta %>%
  select(Location, CowTagID, Chain1:Chain3) %>%
  full_join(depth) %>%
  full_join(dist) %>%
  distinct() %>%
  group_by(CowTagID) %>%
  mutate(meanChain = sum(Chain1, Chain2, Chain3, na.rm = T)/3) %>% #calculate mean chain length across benthos
  drop_na() %>% # remove surveys with no rugosity measurements
  mutate(meanRugosity = meanChain / 2.03) %>%  # chain length = 2.03m
  select(-c(meanChain, Chain1:Chain3))

# save full metadata together as csv
write_csv(Full_meta, here("Data","Full_Metadata.csv"))

species <- full_join(richness, diversity)

Full_species <- Full_meta %>%
  left_join(species) %>%
  left_join(Turb22)

Full_survey <- survey %>%
  select(Location, CowTagID, Taxa, SpeciesCounts) %>%
  full_join(func) %>%
  select(-c(Link)) %>%
  distinct() %>%
  select(Location:Calcification)


# associate silicate CV order to Top Plate ID order
# orderSilicate <- AugChem %>%
#   select(Silicate_umolL, Location, CowTagID) %>%
#   distinct() %>%
#   arrange(Silicate_umolL) %>%
#   mutate(CowTagID = as_factor(as.character(CowTagID))) # as_factor creates levels based on current position

# order CowTagIDs by Nitrogen loading
orderNpercent <- Full_species %>%
  arrange(N_percent) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID)))
NLevels <- paste(orderNpercent$CowTagID) # sort factor levels by N
Full_species$CowTagID <- factor(Full_species$CowTagID, levels = NLevels) # assign order to factor levels by N
levels(Full_species$CowTagID) # check

# order CowTagID's by distance to seep (m)
orderdist <- Full_species %>%
  arrange(dist_to_seep_m) %>%
  mutate(CowTagID = as_factor(as.character(CowTagID)))
distLevels <- paste(orderdist$CowTagID)


# assign order to factor levels
Full_species$CowTagID <- factor(Full_species$CowTagID, levels = NLevels)
Full_species$CowTagID <- factor(Full_species$CowTagID, levels = distLevels)
levels(Full_species$CowTagID) # check

### Processing

# percent cover of species
percent <- Full_survey %>%
  select(Location, CowTagID, Taxa, SpeciesCounts) %>%
  group_by(Location, CowTagID) %>%
  mutate(TotalCounts = sum(SpeciesCounts)) %>%
  ungroup() %>%
  mutate(PercentTaxa = SpeciesCounts / TotalCounts * 100) %>%
  select(-c(SpeciesCounts, TotalCounts))

# percent cover of genera
percent <- Full_survey %>%
  select(Location, CowTagID, Taxon_Group) %>%
  # mutate(Genus = if_else(Genus %in% c('coral1', 'coral2', 'coral3', 'coral4', 'coral5',
  #                                     'coral6', 'coral7', 'coral8', 'coral9', 'coral10',
  #                                     'coral11', 'coral12', 'coral13', 'coral14', 'brown algae1'),
  #                        "Unidentified",
  #                        Genus)) %>%
  group_by(Location, CowTagID) %>%
  count(name = 'TaxonCounts', Taxon_Group) %>%
  mutate(TotalTaxon = sum(TaxonCounts)) %>%
  mutate(PercentTaxon = TaxonCounts / TotalTaxon * 100) %>%
  ungroup() %>%
  select(-c(TaxonCounts, TotalTaxon)) %>%
  right_join(percent)

# percent cover of broad taxon groups
percent <- taxon %>%
  right_join(survey) %>%
  select(CowTagID, Taxon_Group) %>%
  group_by(CowTagID) %>%
  count(name = 'TaxonCounts', Taxon_Group) %>%
  mutate(TotalTaxon = sum(TaxonCounts)) %>%
  mutate(PercentTaxon = TaxonCounts / TotalTaxon * 100) %>%
  ungroup() %>%
  select(-c(TaxonCounts, TotalTaxon)) %>%
  right_join(percent)

# total species richness
richness <- richness %>%
  left_join(AugChem) # join with site data

# total genus richness
richness <- genera %>%
  right_join(survey) %>%
  select(CowTagID, Genus) %>%
  distinct() %>%
  group_by(CowTagID) %>%
  count(name = 'Count', Genus) %>%
  mutate(GenusRichness = sum(Count)) %>%
  distinct(CowTagID, GenusRichness) %>%
  right_join(richness)




### Plotting

# anti_join unidentified species
SpNotIn <- percent %>%
  filter(Taxa %in% c('coral1', 'coral2', 'coral3', 'coral4', 'coral5',
                     'coral6', 'coral7', 'coral8', 'coral9', 'coral10',
                     'coral11', 'coral12', 'coral13', 'coral14', 'brown algae1',
                     'Halimeda2', 'Montipora', 'Porites', 'Psammocora')) #, # remove unidentified taxa
                     #'Bare Rock', 'Bare Rock - exposed', 'Rubble', 'Sand' # remove abiotic substrates


Genus_perc_plot <- percent %>%
  drop_na(Silicate_umolL) %>%
  select(Location, CowTagID, Genus, PercentGenus, Silicate_umolL) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = PercentGenus,
             fill = Genus)) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )


percent %>%
  drop_na(Silicate_umolL) %>%
  select(Location, CowTagID, Taxon_Group, PercentTaxon, Silicate_umolL) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = PercentTaxon,
             fill = Taxon_Group)) +
  geom_col() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

taxon_perc_plot <- percent %>%
  drop_na(Silicate_umolL) %>%
  select(Location, CowTagID, Taxon_Group, PercentTaxon, Silicate_umolL) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = PercentTaxon,
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
taxon_perc_plot

ggsave(here("Output", "Taxon_perc_cover.pdf"), taxon_perc_plot, device = "pdf")


taxon_perc_plot_stacked <- percent %>%
  drop_na(Silicate_umolL) %>%
  select(Location, CowTagID, Taxon_Group, PercentTaxon, Silicate_umolL) %>%
  distinct() %>%
  ggplot(aes(x = CowTagID, # ordered in ascending Silicate_umolL
             y = PercentTaxon,
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
taxon_perc_plot_stacked

ggsave(here("Output", "Taxon_perc_cover_stacked.pdf"), taxon_perc_plot_stacked, device = "pdf")

taxon_perc_plot_bioticstacked <- percent %>%
  drop_na(Silicate_umolL) %>%
  select(Location, CowTagID, Taxon_Group, PercentTaxon, Silicate_umolL) %>%
  distinct() %>%
  filter(Taxon_Group != "Abiotic") %>%
  group_by(CowTagID) %>%
  mutate(totalsum = sum(PercentTaxon),
         newpercent = PercentTaxon / totalsum * 100) %>%
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
taxon_perc_plot_bioticstacked

ggsave(here("Output", "Taxon_perc_cover_bioticstacked.pdf"), taxon_perc_plot_bioticstacked, device = "pdf")


####
