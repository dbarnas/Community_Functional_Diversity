### PROCESSING DEPTH TO PLATES AND CTs
### Created by Danielle Barnas
### Created on 9/28/2022


#### LOAD LIBRARIES ####
library(tidyverse)
library(here)


#### READ IN DATA ####
depth <- read_csv(here("Data", "Sandwich_Depth.csv"))
Vsled <- read_csv(here("Data","Vsled_WL_06032022.csv"))

#### PROCESS RAW DEPTH DATA ####

# round all depth times to nearest 2 minutes
depth <- depth %>%
  mutate(Time = if_else(CowTagID == 'V11', hms("09:28:00"), hms(Time))) %>% # add temporary median time stamp for V11
  mutate(TimeB = if_else(minute(Time) %% 2 == 0, hms(Time), hms(Time) + minutes(1))) %>% # add 1 minute to odd time stamps to match sled times
  mutate(TimeB = if_else(CowTagID == 'C12', hms("12:00:00"), hms(TimeB))) %>% # modify from 11:60:00 to 12:00:00
  mutate(TimeB = if_else(CowTagID == 'V20', hms("07:50:00"), hms(TimeB))) %>% # sled does not have data until 2 minutes after V20 depth recording
  unite(Date, TimeB, col = 'DateTime', sep = " ", remove = T) %>% # unite date and time
  select(-Time) %>%
  mutate(DateTime = mdy_hms(DateTime))

Vsled <- Vsled %>%
  unite(Date, Time, col = 'DateTime', sep = " ", remove = T) %>% # unite date and time
  mutate(DateTime = ymd_hms(DateTime))

# join dataframes based on depth df times
fullDepth <- left_join(depth, Vsled) %>%
  filter(Site == "Varari") %>%  # only use Varari data for now
  select(CowTagID, Dist_CT_cm, Depth) %>%  # only keep relevant data
  mutate(Depth = Depth * 100) # calculate Depth in cm


#### CALCULATE TIDAL DIFFERENCES ####

# identify first value as having zero tidal difference
firstval <- tibble(Tidal_diff = as.numeric(0))

# parse tidal differences vector into tibble
diffs <- as_tibble(diff(fullDepth$Depth, lag = 1)) %>%
  rename(Tidal_diff = value)

# rbind tidal differences, with firstval on top
diffs <- rbind(firstval, diffs)

# bind with full dataframe
fullDepth <- fullDepth %>%
  cbind(diffs) %>%
  mutate(adj_CT_depth_cm = Dist_CT_cm - Tidal_diff) # adjusted CT depths with tidal differences subtracted

# write csv
write_csv(fullDepth, here("Data", "Adj_Sandwich_Depth.csv"))






