#### Begin processing community composition Dominant Species surveys from Varari in August 2021
#### Created by Danielle Barnas

rm(list=ls())

### INSTALL LIBRARIES AS NEEDED
ifelse('tidyverse' %in% rownames(installed.packages()),"",install.packages('tidyverse'))
ifelse('here' %in% rownames(installed.packages()),"",install.packages('here'))
ifelse('geosphere' %in% rownames(installed.packages()),"",install.packages('geosphere'))
ifelse('RColorBrewer' %in% rownames(installed.packages()),"",install.packages('RColorBrewer'))
ifelse('patchwork' %in% rownames(installed.packages()),"",install.packages('patchwork'))
ifelse('ggmap' %in% rownames(installed.packages()),"",install.packages('ggmap'))
ifelse('maps' %in% rownames(installed.packages()),"",install.packages('maps'))
ifelse('maptools' %in% rownames(installed.packages()),"",install.packages('maptools'))
ifelse('ggsn' %in% rownames(installed.packages()),"",install.packages('ggsn'))
ifelse('ggrepel' %in% rownames(installed.packages()),"",install.packages('ggrepel'))
ifelse('scatterpie' %in% rownames(installed.packages()),"",install.packages('scatterpie'))
ifelse('curl' %in% rownames(installed.packages()),"",install.packages('curl'))
ifelse('vegan' %in% rownames(installed.packages()),"",install.packages('vegan'))
ifelse('pairwiseAdonis' %in% rownames(installed.packages()),"",devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"))


### LIBRARIES
library(tidyverse)
library(here)
library(geosphere)

library(RColorBrewer)
library(patchwork)
library(ggmap)
library(maps)
library(maptools)
library(ggsn)
library(ggrepel)
library(scatterpie)
library(curl)

library(vegan)
library(pairwiseAdonis)

API<-names(read_table(here("Data","API.txt")))
register_google(key = API) ### uses my API in separate txt file

### READ IN DATA
survey<-read_csv(here("Data","Surveys","Varari_CC_Survey_reduced_taxon.csv"))
traits<-read_csv(here("Data","Functional_Traits.csv"))
turb<-read_csv(here("Data","Biogeochem","Turb_NC.csv"))
AllChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))



### PERCENT COVER OF DOMINANT SUBSTRATE

## Calculate percentages by quad

sp_count<-survey %>%
  pivot_longer(cols = DS.1:DS.25, names_to = 'quadID', values_to = 'Taxa') %>%
  group_by(Top_Plate_ID,plate.num, Taxa) %>%
  count(Taxa) %>%
  rename(spCount = n) %>%
  mutate(percentCover = spCount/100 * 100) # percent cover

## Even out surveys only observed for 0.5m2 instead of 1m2

halfsurvey<-sp_count %>%  # two sandwiches only surveyed 0.5m2
  filter(Top_Plate_ID == 'V10' | Top_Plate_ID == 'V8') %>%
  mutate(percentCover = percentCover * 2)
# rebind to OG dataframe
sp_count<-sp_count %>%
  filter(Top_Plate_ID != 'V10' & Top_Plate_ID !='V8') %>%
  rbind(halfsurvey)
# remove unnecessary df from environment
rm(list='halfsurvey')



### JOIN DATA (COUNTS, BIOGEOCHEM, TRAIT CATEGORIES)

full_data <- sp_count %>%
  left_join(AllChemData, by = "Top_Plate_ID") %>%
  left_join(turb, by = 'CowTagID') %>%
  left_join(traits) %>%
  distinct()


### CALCULATE DISTANCES TO SGD TREE

# isolate seep lat and lon at Varari
seepData <- AllChemData %>%
  filter(Plate_Seep == 'Seep' | Plate_Seep == 'Spring',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon, Top_Plate_ID, Plate_Seep)

# isolate single numeric value for lat and lon
seepLat <- as.numeric(seepData$lat[1])
seepLon <- as.numeric(seepData$lon[1])

# select distinct points for each plate location to calculate distances
distData <- AllChemData %>%
  filter(Plate_Seep == 'Plate',
         Location == 'Varari') %>%
  select(CowTagID, lat, lon, Top_Plate_ID, Plate_Seep) %>%
  distinct() %>%
  mutate(lat_seep = seepLat,
         lon_seep = seepLon) %>%
  # find Haversine distance
  mutate(dist_to_seep_m = distHaversine(cbind(lon_seep, lat_seep), cbind(lon, lat))) %>%
  # group by sample Site Number
  group_by(CowTagID) %>%
  # choose only minimum distances
  slice(which.min(dist_to_seep_m)) %>%
  select(-c(lat_seep, lon_seep))

# isolate V13, which is in ambient upcurrent of SGD
V13dist <- distData %>%
  filter(CowTagID == 'V13') %>%
  mutate(dist_to_seep_m = -dist_to_seep_m)
# remove V13 from distData to readd new value
distData <- distData %>%
  filter(CowTagID != 'V13') %>%
  rbind(V13dist)

# Relevel based on relative distance to SGTree
full_data <- full_data %>%
  left_join(distData)

# associate distance order to Top Plate ID order
orderPlates <- full_data %>%
  ungroup() %>%
  select(dist_to_seep_m,Top_Plate_ID) %>%
  distinct() %>%
  arrange(dist_to_seep_m) %>%
  # as_factor creates levels based on current position
  mutate(Top_Plate_ID = as_factor(as.character(Top_Plate_ID)))

# sort factor levels by distance
# distLevels <- paste(sort(as.numeric(levels(full_data$dist_to_seep_m))))
distLevels <- paste(levels(orderPlates$Top_Plate_ID))

# assign order to factor levels by distance
full_data$Top_Plate_ID <- factor(full_data$Top_Plate_ID, levels = distLevels)
levels(full_data$Top_Plate_ID) # check

# Relevel based on relative distance to SGTree
# sp_count$Top_Plate_ID<-factor(x = sp_count$Top_Plate_ID,
#                                 levels = c('V14','V20','V17','V10','V9','V18',
#                                            'V8','V7','V15','V6','V4','V16','V3','V2','V1','V13')) # relevel top plate id's for plotting


### SUMMARISE and GRAPH

# Total species counts across all surveys
# bar graph with abundance labels
full_data %>%
  select(Top_Plate_ID,Taxa,spCount) %>%
  distinct() %>%
  group_by(Taxa) %>%
  filter(Taxa != 'Sand/Bare Rock') %>%
  summarise(representation = sum(spCount)) %>%
  ggplot(aes(x = reorder(Taxa,representation), y = representation))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  geom_text(aes(label = reorder(representation,Taxa)),
            nudge_y = 10)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Taxonomic Group",
       y = "Total Species Count")

# Nutrient content across plates
# bar graph colored by tide and faceted by day/night
# x axis is plate location ordered by distance to seep. Graph is distinct to specified nutrient
full_data %>%
  select(Top_Plate_ID,Silicate_umolL,Tide,Day_Night) %>%
  ggplot(aes(x = Top_Plate_ID, y = Silicate_umolL,
             fill = Tide))+
  geom_col(position = 'dodge')+
  facet_wrap(~Day_Night, scales = 'fixed')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90))+
  labs(x = "Plate ID",
       y = "Nutrient Concentration")

# species richness across plates
# bar graph with total richness above columb
# x axis is plate location ordered by distance to seep
sp_richness<-full_data %>%
  select(Top_Plate_ID,Taxa) %>%
  distinct() %>%
  group_by(Top_Plate_ID) %>%
  count(Taxa) %>%
  rename(Sp_Richness = n) %>%
  summarise(Sp_Richness = sum(Sp_Richness))
p1<-sp_richness %>%
  ggplot(aes(x = Top_Plate_ID, y = Sp_Richness))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = Sp_Richness),
            nudge_y = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Plate ID",
       y = "Species Richness")

# functional richness
# bar plot colored and grouped by functional groups
# y axis shows functional traits ordered by groups
# total observations of traits shown on bar
mypalette<-brewer.pal(5,"BuPu")

FR <- full_data %>%
  pivot_longer(cols = c(Morphology:lifespan),
               names_to = "Functional_group",
               values_to = "Functional_trait") %>%
  drop_na(Taxa) %>%
  drop_na(Functional_trait) %>%
  filter(Functional_group != "lifespan" & Functional_group != "Growth_Rate") %>%
  group_by(Functional_group) %>%
  count(Functional_trait) %>%
  ggplot(aes(x = fct_reorder(Functional_trait, desc(Functional_group)),
             y = n,
             fill = Functional_group)) +
  geom_col() +
  coord_flip() +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(aes(label = n),
            nudge_y = 3) +
  labs(x = "Functional Traits",
       y = "Trait Richness",
       fill = "Functional Groups") +
  scale_fill_manual(values = mypalette)




###################################
###### nMDS AND perMANOVA
###################################

# perMANOVA: is there a difference in response variables among groups?

# 3 steps in using classification functions:
# 1. Get the coefficients and constant of the classification equation for each group
# 2. Get a classification score for each observation for each group
#     a. Calculated by using the actual values for each variable to solve the classification for that group
# 3. Classify each observation into the group it matches most closely
#     a. This may or may not be the actual group from which the observation came
#     b. You have to know what group the observation goes into.
#        You have to know the classifications already, and then you can see how good the DFA is at fitting the data

### nMDS plot


# select numerical and single category data only
sp_wide <- full_data %>%
  select(-c(spCount)) %>%
  pivot_wider(values_from = percentCover, names_from = Taxa) %>%
  replace_na(list('Algal Turf' = 0, 'Crustose Corallines' = 0, 'Dictyosphaeria versluysii' = 0, 'Dictyota' = 0,
                  'Galaxaura' = 0, 'Leptastrea transversa'= 0, 'Porites rus' = 0, 'Sand/Bare Rock' = 0,
                  'Turbinaria ornata' = 0, 'Amansia rhodantha' = 0, 'Halimeda' = 0,
                  'Montipora' = 0, 'Sponge' = 0, 'Valonia' = 0, 'Padina' = 0, 'Peyssonnelia' = 0,
                  'Lobed Porites' = 0, 'Pocillopora acuta' = 0, 'Soft Coral' = 0, 'Sargassum' = 0,
                  'Caulerpa' = 0, 'Filamentous Green Algae' = 0, 'Seastar' = 0, 'Anemone' = 0,
                  'Cyphastrea microphthalma' = 0, 'Pavona varians' = 0, 'Ascidian' = 0, 'Pocillopora eydouxi' = 0))
# change code to is.na = TRUE



mds_data<-sp_wide[-c(1:2)] # numerical data only

#create the ordination output using bray curtis dissimilarity matrix
# numerical data only
ord<-metaMDS(mds_data,k=2, distance='bray')

#let's look at the stress with k=2 dimensions. Is it < 0.3?
ord$stress

# Let's look at the stress plot
stressplot(ord)
# looks like a good fit, want to minimize scatter

# basic plot
ordiplot(ord) # dots represent sites (sandwich locations) and + represents substrates

# add text
ordiplot(ord, type = 'text')

# add on plate ID's again
mds_data<-mds_data %>%
  cbind(sp_wide[c(1:2)])

### Plot with ggplot
ord_data<-data.frame(ord$points) %>% # bind with environmental data
  cbind(sp_wide[c(1:2)]) %>%
  left_join(turb, by = "Top_Plate_ID") %>%
  left_join(morph, by = "Top_Plate_ID") %>%
  left_join(sp_richness, by = "Top_Plate_ID")
p4<-ggplot(ord_data, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(size = del15N, color = FT_Richness))+ # size = deltaN
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel(aes(label = Top_Plate_ID)) +
  labs(title = "Morphological Trait Richness")
p4b<-ggplot(ord_data, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(size = del15N, color = Sp_Richness))+ # size = deltaN
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel(aes(label = Top_Plate_ID)) +
  labs(#title = "Species Richness",
       x = "nMDS 1", y = "nMDS 2",
       color = "Species Richness",
       size = "delta Nitrogen-15")

#p4 + p4b

###  MAPPING

### read in API key for maps
API<-names(read_table(here("Data","API.txt")))
register_google(key = API) ### use your own API

gps_reduced<-gps %>%
  select(Top_Plate_ID,lat,lon)
lon = gps_reduced$lon
lat = gps_reduced$lat

ord_data<-ord_data %>%
  left_join(gps_reduced, by = "Top_Plate_ID")

myMap <-get_map(location = c(lon = -149.8999, lat = -17.5405), # set center point
                maptype = 'satellite',
                zoom = 19) # zooms enough to see all 16 points

# Plot site location with nutrient gradient on map
p5<-ggmap(myMap) +
  geom_point(data = ord_data,
             aes(x = lon, y = lat, color = del15N,
                 size = 5)) +
  scale_color_gradient(low = "yellow", high = "blue")+
  geom_label_repel(data = ord_data,
                  aes(label = Top_Plate_ID), # text label
                  #segment.color = 'grey80',
                  colour = "black")+ # text color
  labs(x="", y="",
       color="delta Nitrogen-15")

#p4 / p5

    # ## subset by relative distance to sgd
# highsgd<-mds_data %>%
#   filter(Top_Plate_ID == 'V14' | Top_Plate_ID == 'V20' | Top_Plate_ID == 'V17') %>%
#   mutate(rel_dist='adjacent')
# medsgd<-mds_data %>%
#   filter(Top_Plate_ID == 'V10' | Top_Plate_ID == 'V9' | Top_Plate_ID == 'V18' | Top_Plate_ID == 'V8') %>%
#   mutate(rel_dist='near')
# modsgd<-mds_data %>%
#   filter(Top_Plate_ID == 'V7' | Top_Plate_ID == 'V15' | Top_Plate_ID == 'V6') %>%
#   mutate(rel_dist='moderate')
# lowsgd<-mds_data %>%
#   filter(Top_Plate_ID == 'V4' | Top_Plate_ID == 'V16' | Top_Plate_ID == 'V3' | Top_Plate_ID == 'V2' |
#            Top_Plate_ID == 'V1' | Top_Plate_ID == 'V13') %>%
#   mutate(rel_dist='far')
#
# mds_data<-rbind(highsgd,medsgd,modsgd,lowsgd)
#

#We could formally test whether sites are different using a perMANOVA.
#This requires the adonis command in the vegan package

#numeric data only
permanovamodel<-adonis(mds_data[,1:28]~rel_dist, mds_data, permutations = 999,
                       method="bray")
permanovamodel
# just barely significant between four groups
# p = 0.046

#If we are to trust the results of the permanova, then we have to assume that the dispersion among
#data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(mds_data[,1:28])
betadisper(disper, mds_data$rel_dist)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

#An option for doing post-hoc pairwise comparisons in R
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pairwise.adonis(mds_data[,1:28], mds_data$rel_dist, perm=999)
# p < 0.05 for adjacent vs near and adjacent vs far

#Get coefficients to see which species are most important in explaining site differences:
permanovamodel$coefficients


