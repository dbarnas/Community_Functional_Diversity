##### Map of Moorea and survey locations #####
##### Created by: Danielle Barnas #####
##### Created on: 10/22/2021 #####
##### Edited on: 2/25/2023 #####

##### LOAD LIBRARY #####
library(here)
library(tidyverse)
library(ggmap)
library(maptools)
library(ggrepel)

##### READ IN DATA #####
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))
meta <- read_csv(here("Data", "Full_Metadata.csv")) %>%
  left_join(alphatag) %>%
  filter(Location == "Varari",
         CowTagID != "V13")


# isolate seep point for mapping
seeppt <- meta %>%
  filter(CowTagID == "VSEEP")

# mean lat and long for the maps
LocationGPS <- meta %>%
  group_by(Location) %>% # varari vs cabral
  summarise(lon = median(lon + 0.00013, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))


# Varari
VarariBaseMap<-get_map(LocationGPS %>% filter(Location == "Varari") %>% select(lon,lat),
                       maptype = 'satellite',
                       zoom = 19)

# base map
# Varari
VmapSites <- ggmap(VarariBaseMap) +
  labs(x = "Longitude", y = "Latitude") +  #label x and y axes
  geom_point(data = meta,
             aes(x = lon, y = lat),
             size = 8,
             shape = 22,
             fill = "white",
             color = "black") +
  geom_text(data = meta,
             aes(x = lon, y = lat,
                 label = AlphaTag),
             size = 4) +
  # add the seep point separately
  geom_label(data = seeppt,
             aes(x = lon, y = lat + 0.00003),
             label = "Seep\nA",
             fill = "white") +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 16))

VmapSites

ggsave(here("Output", "PaperFigures", "Varari_map_sites.png"), VmapSites, height = 6, width = 6, device = "png")


# Moorea
MooreaMap<-get_map('Moorea', maptype = 'satellite', zoom = 12)

MooreaMapPlot <- ggmap(MooreaMap) + # base map

  geom_rect(data=LocationGPS, aes(xmin=lon[1] - 0.006, xmax=lon[1] + 0.006, ymin=lat[1] - 0.01, ymax=lat[1] + 0.01), color="white", alpha=0, size = 2) + # Cabral square

  #geom_rect(data=LocationGPS, aes(xmin=lon[2] - 0.006, xmax=lon[2] + 0.006, ymin=lat[2] - 0.01, ymax=lat[2] + 0.01), color="white", alpha=0, size = 2) + # Varari square

  geom_point(data=LocationGPS, aes(x = lon, y = lat, label = Location),
             shape = 18, color = "white", fill = "white", size = 8) + # adds symbol at center of Location

  labs(x = "Longitude", y = "Latitude") +

  geom_text(data = LocationGPS, aes(label = Location), color = "white", hjust = -1, size = 10) + # adds Location names to the right of the boxes

  geom_segment(x = LocationGPS$lon[1] + 0.006, y = LocationGPS$lat[1], xend = LocationGPS$lon[1] + 0.023, yend = LocationGPS$lat[1], color = "white", size = 2)  # adds horizontal line from edge of box to Location name

#geom_segment(x = LocationGPS$lon[2] + 0.006, y = LocationGPS$lat[2], xend = LocationGPS$lon[2] + 0.023, yend = LocationGPS$lat[2], color = "white", size = 2)

MooreaMapPlot

ggsave(here("Output", "PaperFigures", "Moorea_Map.png"), MooreaMapPlot,height = 10, width = 10)

