# Package ID: knb-lter-mcr.8.33 Cataloging System:https://pasta.edirepository.org.
# Data set title: MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Benthic Algae and Other Community Components, ongoing since 2005.
# Data set creator:    - Moorea Coral Reef LTER
# Data set creator:  Robert Carpenter - Moorea Coral Reef LTER
# Metadata Provider:    - Moorea Coral Reef LTER
# Contact:    - Information Manager Moorea Coral Reef LTER  - mcrlter@msi.ucsb.edu
# Stylesheet for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@Virginia.edu
#
#install package tidyverse if not already installed
if(!require(tidyverse)){ install.packages("tidyverse") }
library("tidyverse")
infile1 <- trimws("https://pasta.lternet.edu/package/data/eml/knb-lter-mcr/8/33/54d54c25616a48b9ec684118df9d6fca")
infile1 <-sub("^https","http",infile1)
# This creates a tibble named: dt1
dt1 <-read_delim(infile1
                 ,delim=","
                 ,skip=1
                 ,quote='"'
                 , col_names=c(
                   "Year",
                   "Date",
                   "Location",
                   "Site",
                   "Habitat",
                   "Transect",
                   "Quadrat",
                   "Taxonomy_Substrate_Functional_Group",
                   "Percent_Cover"   ),
                 col_types=list(
                   col_character(),
                   col_date("%Y-%m-%d"),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_character(),
                   col_number() ),
                 na=c(" ",".","NA")  )


# Convert Missing Values to NA for individual vectors
dt1$Percent_Cover <- ifelse((trimws(as.character(dt1$Percent_Cover))==trimws("-1")),NA,dt1$Percent_Cover)
suppressWarnings(dt1$Percent_Cover <- ifelse(!is.na(as.numeric("-1")) & (trimws(as.character(dt1$Percent_Cover))==as.character(as.numeric("-1"))),NA,dt1$Percent_Cover))


# Observed issues when reading the data. An empty list is good!
problems(dt1)
# Here is the structure of the input data tibble:
glimpse(dt1)
# And some statistical summaries of the data
summary(dt1)
# Get more details on character variables

summary(as.factor(dt1$Location))
summary(as.factor(dt1$Site))
summary(as.factor(dt1$Habitat))
summary(as.factor(dt1$Transect))
summary(as.factor(dt1$Quadrat))
summary(as.factor(dt1$Taxonomy_Substrate_Functional_Group))




