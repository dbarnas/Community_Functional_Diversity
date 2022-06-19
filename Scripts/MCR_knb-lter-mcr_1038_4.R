# Package ID: knb-lter-mcr.1038.4 Cataloging System:https://pasta.lternet.edu.
# Data set title: MCR LTER: Coral Reef: Long-term Community Dynamics: Backreef (Lagoon) Corals Annual Survey, ongoing since 2005.
# Data set creator:    - Moorea Coral Reef LTER
# Data set creator:  Peter Edmunds - Moorea Coral Reef LTER
# Contact:    - Information Manager LTER Network Office  - tech-support@lternet.edu
# Contact:    - Information Manager Moorea Coral Reef LTER  - mcrlter@msi.ucsb.edu
# Metadata Link: https://portal.lternet.edu/nis/metadataviewer?packageid=knb-lter-mcr.1038.4
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@Virginia.edu
#
#install package tidyverse if not already installed
if(!require(tidyverse)){ install.packages("tidyverse") }
library("tidyverse")
infile1 <- trimws("https://pasta.lternet.edu/package/data/eml/knb-lter-mcr/1038/4/b7ca625d17a5630441af8c229cc0ab92")
infile1 <-sub("^https","http",infile1)
# This creates a tibble named: dt1
dt1 <-read_delim(infile1
                 ,delim=","
                 ,skip=1
                 , col_names=c(
                   "year",
                   "site",
                   "patch",
                   "quad",
                   "coral",
                   "macroalgae",
                   "cb",
                   "sand",
                   "turf",
                   "millepora",
                   "sub_squares",
                   "note"   ),
                 col_types=list(
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_number() ,
                   col_number() ,
                   col_number() ,
                   col_number() ,
                   col_number() ,
                   col_number() ,
                   col_number() ,
                   col_character()),
                 na=c(" ",".","NA")  )


# Convert Missing Values to NA for individual vectors
dt1$coral <- ifelse((trimws(as.character(dt1$coral))==trimws("na")),NA,dt1$coral)
suppressWarnings(dt1$coral <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$coral))==as.character(as.numeric("na"))),NA,dt1$coral))
dt1$coral <- ifelse((trimws(as.character(dt1$coral))==trimws("BW")),NA,dt1$coral)
suppressWarnings(dt1$coral <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$coral))==as.character(as.numeric("BW"))),NA,dt1$coral))
dt1$macroalgae <- ifelse((trimws(as.character(dt1$macroalgae))==trimws("na")),NA,dt1$macroalgae)
suppressWarnings(dt1$macroalgae <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$macroalgae))==as.character(as.numeric("na"))),NA,dt1$macroalgae))
dt1$macroalgae <- ifelse((trimws(as.character(dt1$macroalgae))==trimws("BW")),NA,dt1$macroalgae)
suppressWarnings(dt1$macroalgae <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$macroalgae))==as.character(as.numeric("BW"))),NA,dt1$macroalgae))
dt1$cb <- ifelse((trimws(as.character(dt1$cb))==trimws("na")),NA,dt1$cb)
suppressWarnings(dt1$cb <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$cb))==as.character(as.numeric("na"))),NA,dt1$cb))
dt1$cb <- ifelse((trimws(as.character(dt1$cb))==trimws("BW")),NA,dt1$cb)
suppressWarnings(dt1$cb <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$cb))==as.character(as.numeric("BW"))),NA,dt1$cb))
dt1$sand <- ifelse((trimws(as.character(dt1$sand))==trimws("na")),NA,dt1$sand)
suppressWarnings(dt1$sand <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$sand))==as.character(as.numeric("na"))),NA,dt1$sand))
dt1$sand <- ifelse((trimws(as.character(dt1$sand))==trimws("BW")),NA,dt1$sand)
suppressWarnings(dt1$sand <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$sand))==as.character(as.numeric("BW"))),NA,dt1$sand))
dt1$turf <- ifelse((trimws(as.character(dt1$turf))==trimws("na")),NA,dt1$turf)
suppressWarnings(dt1$turf <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$turf))==as.character(as.numeric("na"))),NA,dt1$turf))
dt1$turf <- ifelse((trimws(as.character(dt1$turf))==trimws("BW")),NA,dt1$turf)
suppressWarnings(dt1$turf <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$turf))==as.character(as.numeric("BW"))),NA,dt1$turf))
dt1$millepora <- ifelse((trimws(as.character(dt1$millepora))==trimws("na")),NA,dt1$millepora)
suppressWarnings(dt1$millepora <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$millepora))==as.character(as.numeric("na"))),NA,dt1$millepora))
dt1$millepora <- ifelse((trimws(as.character(dt1$millepora))==trimws("BW")),NA,dt1$millepora)
suppressWarnings(dt1$millepora <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt1$millepora))==as.character(as.numeric("BW"))),NA,dt1$millepora))
dt1$note <- ifelse((trimws(as.character(dt1$note))==trimws("na")),NA,dt1$note)
suppressWarnings(dt1$note <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt1$note))==as.character(as.numeric("na"))),NA,dt1$note))


# Observed issues when reading the data. An empty list is good!
problems(dt1)
# Here is the structure of the input data tibble:
glimpse(dt1)
# And some statistical summaries of the data
summary(dt1)
# Get more details on character variables

summary(as.factor(dt1$site))
summary(as.factor(dt1$patch))
summary(as.factor(dt1$quad))
summary(as.factor(dt1$note))
infile2 <- trimws("https://pasta.lternet.edu/package/data/eml/knb-lter-mcr/1038/4/fafae7af67a337e366ee58db25938f48")
infile2 <-sub("^https","http",infile2)
# This creates a tibble named: dt2
dt2 <-read_delim(infile2
                 ,delim=","
                 ,skip=1
                 , col_names=c(
                   "year",
                   "site",
                   "patch",
                   "quad",
                   "benthic_category",
                   "percent_cover"   ),
                 col_types=list(
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_number() ),
                 na=c(" ",".","NA")  )


# Convert Missing Values to NA for individual vectors
dt2$percent_cover <- ifelse((trimws(as.character(dt2$percent_cover))==trimws("na")),NA,dt2$percent_cover)
suppressWarnings(dt2$percent_cover <- ifelse(!is.na(as.numeric("na")) & (trimws(as.character(dt2$percent_cover))==as.character(as.numeric("na"))),NA,dt2$percent_cover))
dt2$percent_cover <- ifelse((trimws(as.character(dt2$percent_cover))==trimws("BW")),NA,dt2$percent_cover)
suppressWarnings(dt2$percent_cover <- ifelse(!is.na(as.numeric("BW")) & (trimws(as.character(dt2$percent_cover))==as.character(as.numeric("BW"))),NA,dt2$percent_cover))


# Observed issues when reading the data. An empty list is good!
problems(dt2)
# Here is the structure of the input data tibble:
glimpse(dt2)
# And some statistical summaries of the data
summary(dt2)
# Get more details on character variables

summary(as.factor(dt2$site))
summary(as.factor(dt2$patch))
summary(as.factor(dt2$quad))
summary(as.factor(dt2$benthic_category))


####################################################################################
####################################################################################

dt2 %>%
  group_by(year,site,benthic_category) %>%
  drop_na(percent_cover) %>%
  filter(benthic_category != "BW") %>%
  summarise(mean_pc = mean(percent_cover, na.rm = T)) %>%
  ggplot(aes(x = benthic_category,
             y = mean_pc,
             fill = site)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~year)


## by 2013 coral and turf  are dominant
