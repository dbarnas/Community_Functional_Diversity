### Select random points within Varari and Cabral sites


## load library#####
library(here)
library(tidyverse)
library(curl)
library(ggmap)
library(maptools)
library(kriging)
library(sf)
library(patchwork)

#
# library(patchwork)
# library(viridis)
# library(ggnewscale)
# library(wql)
# library(glue)
# library(gridExtra)



## Read in data from github repository url
AllChemData<-read_csv(curl('https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC.csv'))
#View(AllChemData)

# mean lat and long for the maps
MeanGPS<-AllChemData %>%
  filter(Location != "Offshore")%>%
  group_by(Location) %>%
  summarise(lon = median(lon, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))

SiteGPS<-AllChemData %>%
  filter(Location != "Offshore")%>%
  group_by(Location, CowTagID) %>%
  summarise(lon = mean(lon, na.rm = TRUE),
            lat = mean(lat, na.rm = TRUE))
VGPS <- SiteGPS %>%
  filter(Location == "Varari", CowTagID != "Varari_Well")
CGPS <- SiteGPS %>%
  filter(Location == "Cabral", CowTagID != "CSPRING")

# Base Maps
# API<-names(read_table(here("Data","API.txt")))
# register_google(key = API) ### use your own API

# Varari
#VarariBaseMap<-get_map(MeanGPS %>% filter(Location == "Varari") %>% select(lon,lat), maptype = 'satellite', zVaraoom = 18)

# Cabral
#CabralBaseMap<-get_map(MeanGPS %>% filter(Location == "Cabral") %>% select(lon,lat), maptype = 'satellite', zoom = 18)


### Bring in the polygons for the sites
#Varari
V_kml <- getKMLcoordinates(kmlfile=here("Data","Polygons","Varari_Polygon.kml"), ignoreAltitude=T)
#Cabral
C_kml <- getKMLcoordinates(kmlfile=here("Data","Polygons","Cabral_polygon.kml"), ignoreAltitude=T)

## RANDOM VARARI POINTS

# R will pick the 7th random set of points every time (gives me static "random" points,
# so I don't get different points every time I run)
# set.seed(7) for Varari, but set.seed(2) for Cabral to have different, more even distribution
set.seed(7)

V_kml_sf<-sf::st_read(here("Data","Polygons","Varari_Polygon.kml")) #read in kml polygon file as a shape file (sf)
myrandomV<-sf::st_sample(V_kml_sf, size=80, type="random") #sample random points from polygon shape file
myrandomV

#lon_list<-lapply(myrandom, '[[', 1)  # This returns a list with only the first element
Vlon_list<-unlist(lapply(myrandomV, '[[', 1)) # This returns a vector with the third element
Vlon_list

Vlat_list<-unlist(lapply(myrandomV, '[[', 2)) # This returns a vector with the third element
Vlat_list

myrandomVpts<-cbind(Vlon_list,Vlat_list) # merge lat and lon lists
myrandomVpts<-data.frame(myrandomVpts) %>% # parse to dataframe and rename column headings
  rename(lon="Vlon_list", lat="Vlat_list")

## RANDOM CABRAL POINTS
set.seed(2)

C_kml_sf<-sf::st_read(here("Data","Polygons","Cabral_Polygon.kml")) #read in kml polygon file as a shape file (sf)
myrandomC<-sf::st_sample(C_kml_sf, size=80, type="random") #sample random points from polygon shape file
myrandomC

#lon_list<-lapply(myrandom, '[[', 1)  # This returns a list with only the third element
Clon_list<-unlist(lapply(myrandomC, '[[', 1)) # This returns a vector with the third element
Clon_list

Clat_list<-unlist(lapply(myrandomC, '[[', 2)) # This returns a vector with the third element
Clat_list

myrandomCpts<-cbind(Clon_list,Clat_list) # merge lat and lon lists
myrandomCpts<-data.frame(myrandomCpts) %>% # parse to dataframe and rename column headings
  rename(lon="Clon_list", lat="Clat_list")

## Join points of both locations
Vpts<-myrandomVpts %>%
  mutate(Location = "Varari") %>% # add location column
  arrange(lat) # organize data from south to north
Cpts<-myrandomCpts %>%
  mutate(Location = "Cabral") %>%
  arrange(lat)
Rpts<-Vpts %>%
  full_join(Cpts)
Rpts <- Rpts %>%
  mutate(CowTagID = paste0("random_",rownames(Rpts)))
Allpts <- Rpts %>%
  full_join(VGPS) %>%
  full_join(CGPS)

#write_csv(Allpts, here("Data","Randomized_Locations.csv"))
## Import csv into Google Earth to visualize


# Maps with random points

Allpts <- read_csv(here("Data","Randomized_Locations.csv"))

# Varari
VarariRandomMap<-get_map(MeanGPS %>% filter(Location == "Varari") %>% select(lon,lat) %>% mutate(lon = lon + 0.00005, lat = lat + 0.0001), #Allpts %>% filter(Location == "Varari") %>% select(lon,lat),
                         maptype = 'satellite', zoom = 19)
VSeepPt<-VGPS %>% filter(CowTagID == "VSEEP") # isolate Varari seep
VMap<-ggmap(VarariRandomMap) +
  geom_point(data=VGPS, aes(x=lon, y=lat),
             size=4, shape = 23, alpha = 0.8, fill = "brown2", color = "black") + # selected 20 points
  geom_point(data=Vpts, aes(x=lon, y=lat),
             size=4, shape = 21, alpha = 0.8, fill = "white", color = "black") + # randomized points
  geom_point(data=VSeepPt, aes(x=lon, y=lat),
             size=4, shape = 8, color = "yellow") + # seepage point
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude")
VMap

ggsave(here("Output","Random_Varari_map.pdf"),VMap, "pdf", width=6.88, height=4.35, units="in")
ggsave(here("Output","Random_Varari_map.png"),VMap, "png", width=6.88, height=4.35, units="in")

# Cabral
CabralRandomMap<-get_map(MeanGPS %>% filter(Location == "Cabral") %>% select(lon,lat) %>% mutate(lon = lon + 0.00005, lat = lat - 0.00008),
                                 maptype = 'satellite', zoom = 18)
CSeepPt<-CGPS %>% filter(CowTagID == "CSEEP") # isolate Varari seep
CMap<-ggmap(CabralRandomMap) +
  geom_point(data=CGPS, aes(x=lon, y=lat),
             size=4, shape = 23, alpha = 0.8, fill = "brown2", color = "black") + # selected 20 points
  geom_point(data=Cpts, aes(x=lon, y=lat),
             size=4, shape = 21, alpha = 0.8, fill = "white", color = "black") + # randomized points
  geom_point(data=CSeepPt, aes(x=lon, y=lat),
             size=4, shape = 8, color = "yellow") + # seepage point
  theme(legend.position = "none") +
  labs(x = "Longitude", y = "Latitude")
CMap

ggsave(here("Output","Random_Cabral_maps.pdf"), CMap, "pdf", width=6.88, height=4.35, units="in")
ggsave(here("Output","Random_Cabral_maps.png"), CMap, "png", width=6.88, height=4.35, units="in")

## Patchwork Maps

SurveyMaps<-VMap + CMap +
  plot_layout(nrow = 1) +
  patchwork::plot_annotation(tag_levels = 'A',
                             theme = theme(plot.title = element_text(size = rel(1.5),
                                                                     face = "bold",
                                                                     hjust = 0.5,
                                                                     margin = margin(t = 10, b = 20,
                                                                                     unit = "pt"))))
SurveyMaps

ggsave(plot = SurveyMaps, filename = here("Output","Random_Survey_Maps.pdf"), width = 14, height = 10)
ggsave(plot = SurveyMaps, filename = here("Output","Random_Survey_Maps.png"), width = 14, height = 10)



