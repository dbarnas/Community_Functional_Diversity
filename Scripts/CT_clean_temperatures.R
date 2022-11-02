### CLEAN CT FILES FROM BOTH DEPLOYMENT SEASONS ###
### Created by Danielle Barnas
### Created on 11/2/2022


library(tidyverse)
library(here)
library(lubridate)

#list all csv file names in the CT folder

ct1 <- read_csv(here("Data", "Biogeochem", "CT_files", "Full_CT_07232021_final.csv"))
ct2 <- read_csv(here("Data", "Biogeochem", "CT_files", "Full_CT_08112021.csv"))
ct3 <- read_csv(here("Data", "Biogeochem", "CT_files", "Full_CT_20220318.csv"))
ct4 <- read_csv(here("Data", "Biogeochem", "CT_files", "Full_CT_20220323_b.csv"))
ct5 <- read_csv(here("Data", "Biogeochem", "CT_files", "Full_CT_20220323.csv"))
meta <- read_csv(here("Data", "Biogeochem", "CT_files", "CTLoggerIDMetaData_Plates.csv"))


## join and clean CT temperatures
meta <- meta %>%
  select(Location = Site, CowTagID, LoggerID, date = Date_launched) %>%
  mutate(date = as.character(mdy(date)))

ct1 <- ct1 %>%
  select(date, LoggerID, TempInSitu, Salinity_psu) %>%
  mutate(date = mdy_hm(date)) %>%
  rename(datetime = date) %>%
  separate(col = datetime, into = c("date", "time"), sep = " ", remove = F) %>%
  left_join(meta) %>%
  group_by(LoggerID) %>%
  fill(Location, CowTagID, .direction = c("down"))

ct2 <- ct2 %>%
  select(date, LoggerID, TempInSitu, Salinity_psu) %>%
  mutate(date = ymd_hms(date)) %>%
  rename(datetime = date) %>%
  separate(col = datetime, into = c("date", "time"), sep = " ", remove = F) %>%
  left_join(meta) %>%
  group_by(LoggerID) %>%
  fill(Location, CowTagID, .direction = c("down"))

ct3 <- ct3 %>%
  select(date, LoggerID, TempInSitu, Salinity_psu) %>%
  mutate(date = ymd_hms(date)) %>%
  rename(datetime = date) %>%
  separate(col = datetime, into = c("date", "time"), sep = " ", remove = F) %>%
  left_join(meta) %>%
  group_by(LoggerID) %>%
  fill(Location, CowTagID, .direction = c("down"))

ct4 <- ct4 %>%
  select(date, LoggerID, TempInSitu, Salinity_psu) %>%
  mutate(date = ymd_hms(date)) %>%
  rename(datetime = date) %>%
  separate(col = datetime, into = c("date", "time"), sep = " ", remove = F) %>%
  left_join(meta) %>%
  group_by(LoggerID) %>%
  fill(Location, CowTagID, .direction = c("down"))

ct5 <- ct5 %>%
  select(date, LoggerID, TempInSitu, Salinity_psu) %>%
  mutate(date = ymd_hms(date)) %>%
  rename(datetime = date) %>%
  separate(col = datetime, into = c("date", "time"), sep = " ", remove = F) %>%
  left_join(meta) %>%
  group_by(LoggerID) %>%
  fill(Location, CowTagID, .direction = c("down"))

# join
ct_full <- ct1 %>%
  full_join(ct2) %>%
  full_join(ct3) %>%
  full_join(ct4) %>%
  full_join(ct5) %>%
  select(Location, CowTagID, LoggerID, datetime, TempInSitu, Salinity_psu)

# remove outlier temperatures
rm1 <- ct_full %>% # temp jumps up
  filter(LoggerID == "350",
         datetime > "2022-03-22 00:00:00")
rm2 <- ct_full %>% # temp jumps up
  filter(LoggerID == "340",
         datetime > "2022-03-22 00:00:00")
ct_full <- ct_full %>%
  anti_join(rm1) %>%
  anti_join(rm2)

# plot
ct_full %>%
  ggplot(aes(x = datetime, y = TempInSitu)) +
  geom_point() +
  facet_wrap(~LoggerID, scales = "free_x")


ct_full <- ct_full %>%
  group_by(Location, CowTagID) %>%
  summarise(meanTemp = mean(TempInSitu, na.rm = TRUE),
            sdTemp = sd(TempInSitu, na.rm = TRUE),
            maxTemp = max(TempInSitu, na.rm = TRUE),
            minTemp = min(TempInSitu, na.rm = TRUE)) %>%
  mutate(rangeTemp = maxTemp - minTemp,
         CVTemp = sdTemp / meanTemp * 100)

write_csv(ct_full, here("Data","Biogeochem","CT_files","cleaned_CT.csv"))
