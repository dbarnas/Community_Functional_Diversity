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
meta <- read_csv(here("Data", "Full_Metadata.csv"))


# mean lat and long for the maps
LocationGPS <- meta %>%
  #filter(CowTagID != "V13") %>%
  group_by(Location) %>% # varari vs cabral
  summarise(lon = median(lon + 0.00013, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE))


# Varari
VarariBaseMap<-get_map(LocationGPS %>% filter(Location == "Varari") %>% select(lon,lat), maptype = 'satellite', zoom = 19)

# base map
# Varari
VmapSites <- ggmap(VarariBaseMap) +
  # geom_point(data = meta,
  #            aes(x = lon, y = lat),
  #            color = "white",
  #            size = 2) +
  labs(x = "Longitude", y = "Latitude",  #label x and y axes
       title = "Varari Sample Locations") +
  geom_label_repel(data = meta,
             aes(x = lon, y = lat,
                 label = CowTagID),
             size = 3,
             alpha = 0.9,
             box.padding = 0.05) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

VmapSites

ggsave(here("Output", "PaperFigures", "Varari_map_sites.png"), VmapSites, height = 5, width = 5, device = "png")


# Moorea
MooreaMap<-get_map('Moorea', maptype = 'satellite', zoom = 12)

MooreaMapPlot <- ggmap(MooreaMap) + # base map

  geom_rect(data=LocationGPS, aes(xmin=lon[1] - 0.006, xmax=lon[1] + 0.006, ymin=lat[1] - 0.01, ymax=lat[1] + 0.01), color="white", alpha=0, size = 2) + # Cabral square

  geom_rect(data=LocationGPS, aes(xmin=lon[2] - 0.006, xmax=lon[2] + 0.006, ymin=lat[2] - 0.01, ymax=lat[2] + 0.01), color="white", alpha=0, size = 2) + # Varari square

  geom_point(data=LocationGPS, aes(x = lon, y = lat, label = Location),
             shape = 18, color = "white", fill = "white", size = 8) + # adds symbol at center of Location

  labs(x = "Longitude", y = "Latitude") +

  geom_text(data = LocationGPS, aes(label = Location), color = "white", hjust = -1, size = 10) + # adds Location names to the right of the boxes

  geom_segment(x = LocationGPS$lon[1] + 0.006, y = LocationGPS$lat[1], xend = LocationGPS$lon[1] + 0.023, yend = LocationGPS$lat[1], color = "white", size = 2) + # adds horizontal line from edge of box to Location name

  geom_segment(x = LocationGPS$lon[2] + 0.006, y = LocationGPS$lat[2], xend = LocationGPS$lon[2] + 0.023, yend = LocationGPS$lat[2], color = "white", size = 2)

MooreaMapPlot

ggsave(here("Output", "PaperFigures", "Moorea_Map.png"), MooreaMapPlot,height = 10, width = 10)

