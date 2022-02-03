### Select random points within Varari and Cabral sites


## load library#####
library(here)
library(tidyverse)
library(curl)
library(ggmap)
library(maptools)
library(kriging)
library(sf)

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

# Base Maps
API<-names(read_table(here("Data","API.txt")))
register_google(key = API) ### use your own API

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
write_csv(Rpts, here("Data","Randomized_Locations.csv"))
## Import csv into Google Earth to visualize


# Maps with random points
# Varari
VarariRandomMap<-get_map(myrandomVpts %>% select(lon,lat), 
                         maptype = 'satellite', zoom = 18)
VMap<-ggmap(VarariRandomMap) +
  geom_point(data=myrandomVpts, aes(x=lon, y=lat, colour="red"), 
             size=2, alpha=0.8) +
  theme(legend.position = "none")
ggsave(here("Output","Random_Varari_map.pdf"),VMap, "pdf", width=6.88, height=4.35, units="in")

# Cabral
CabralRandomMap<-get_map(myrandomCpts %>% select(lon,lat), 
                         maptype = 'satellite', zoom = 18)
CMap<-ggmap(CabralRandomMap) +
  geom_point(data=myrandomCpts, aes(x=lon, y=lat, colour="red"), 
             size=2, alpha=0.8) +
  theme(legend.position = "none")
ggsave(here("Output","Random_Cabral_maps.pdf"), CMap, width=6.88, height=4.35, units="in")




