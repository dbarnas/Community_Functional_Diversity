
library(tidyverse)
library(here)
rm(list=ls())

## read in metadata sheets
biometa<-read_csv(here("Data","RespoFiles","BioMetadata.csv"))
respometa<-read_csv(here("Data","RespoFiles","RespoMetadata.csv"))
sp<-read_csv(here("Data","RespoFiles","SpeciesMetadata.csv"))

## Assign species ID's from species metadata based on assemblage ID's from bio metadata
for(i in 1:nrow(sp)){
  Srow<-sp$Assemblage.ID[i] # locate ith row of species metadata and assign character value
  for(j in 1:nrow(biometa)){
    if(str_detect(string = biometa$AssemblageID[j],pattern = Srow) == T){ # if character is in assemblage ID string in jth row of bio metadata
      if(is.na(biometa$SpeciesID[j]) == F){
        biometa$SpeciesID[j] = paste0(biometa$SpeciesID[j],",",sp$SpeciesID[i]) # add species ID character to associated assemblage rows
      } else {
        biometa$SpeciesID[j] = paste0(sp$SpeciesID[i])
      }}}}
View(biometa)
#write_csv(biometa,here("Data","RespoFiles","BioMetadata.csv"))

## Create the FileID column in RespoMetadata based on sample and block id's
respometa<-respometa %>%
  unite(col="FileID",sep="_",SampleID,treatment.block,block,remove=F) %>%
  mutate(FileID = paste0(FileID,"_O2.csv"))
View(respometa)
#write_csv(respometa,here("Data","RespoFiles","RespoMetadata.csv"))
