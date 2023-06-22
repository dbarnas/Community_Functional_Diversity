

############################
### LOAD LIBRARIES
############################

library(tidyverse)
library(here)
library(PNWColors)


############################
### READ IN DATA
############################

traits <- read_csv(here("Data", "Surveys","Distinct_Taxa.csv"))
myspecies <- read_csv(here("Data", "Surveys", "Species_Composition_2022.csv"))
alphatag <- read_csv(here("Data", "CowTag_to_AlphaTag.csv"))

#meta <- read_csv(here("Data", "Full_Metadata.csv"))
chem <- read_csv(here("Data","Biogeochem", "Nutrients_Processed_All.csv"))


############################
### CLEAN DATA
############################

chem <- chem %>%
  filter(Season == "Dry") %>%
  filter(Location == "Varari",
         #CowTagID != "VSEEP" &
         CowTagID != "V13") %>%
  select(CowTagID, Parameters, CVSeasonal) %>%
  pivot_wider(names_from = Parameters, values_from = CVSeasonal) %>%
  select(CowTagID, Phosphate_umolL) %>%
  left_join(alphatag)

AlphaOrder <- chem %>%
  arrange(Phosphate_umolL)
AlphaOrder <- AlphaOrder$AlphaTag

myspecies <- myspecies %>%
  filter(Location == "Varari", # only analyze varari for now
         CowTagID != "V13") %>%
  left_join(chem) %>%
  mutate(AlphaTag = factor(AlphaTag, levels = AlphaOrder))

traits <- traits %>%
  select(Taxa, Taxon_Group, Morph2, Calc, ER)


############################
### ANALYZE SPECIES PRESENCE
############################
mypal <- pnw_palette("Bay", n = 13)

# calculate percent cover of species including substrate points
# to get actual abundance of species along gradient
topSp <- myspecies %>%
  select(CowTagID:SpeciesCounts, AlphaTag) %>% # remove unnecessary columns
  select(-Date) %>%
  group_by(CowTagID) %>%
  mutate(totalCounts = sum(SpeciesCounts)) %>%
  group_by(CowTagID, Taxa) %>%
  mutate(SpeciesCounts = sum(SpeciesCounts)) %>% # add all of same species together
  distinct() %>% # remove duplicates
  ungroup() %>%
  mutate(pCover = SpeciesCounts / totalCounts*100) %>%
  distinct() %>%
  #filter(pCover > 10) %>%  # select top abundant species
  filter(Taxa != "Sand" &
           Taxa != "Bare Rock" &
           Taxa != "Rubble")

# fixed y axes
# topSp %>%
#   ggplot(aes(x = AlphaTag, y = pCover, fill = Taxa)) +
#   geom_col(color = "black", position = "dodge") +
#   theme_bw() + theme(panel.grid.major.x = element_blank()) +
#   #scale_fill_manual(values = mypal) +
#   facet_wrap(~Taxa)

# put in zero values
zerotopSp <- topSp %>%
  pivot_wider(names_from = Taxa, values_from = pCover) %>%
  mutate_at(.vars = 5:ncol(.), .funs= ~if_else(is.na(.),0,.)) %>%
  pivot_longer(cols = 5:ncol(.), names_to = "Taxa", values_to = "pCover")

# relative abundance across sites, zeros included
Supp.Abundance.Plot <- zerotopSp %>%
  # adds line break into string for graphing
  separate(Taxa, into = c("t1","t2"), sep = " ",) %>%
  mutate(t3 = "\n",
         t2 = if_else(is.na(t2),"",t2)) %>%
  unite(t1,t3,t2, col = "Taxa", sep = "") %>%
  ggplot(aes(x = AlphaTag, y = pCover)) +
  geom_col(color = "black", fill = "grey",
           position = "dodge", show.legend = FALSE) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 8)) +
  labs(x = "Survey Location",
       y = "Abundance (%)") +
  scale_y_continuous(n.breaks = 3) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 5)
Supp.Abundance.Plot
ggsave(here("Output","PaperFigures","Supp_Sp_Abundance.png"), Supp.Abundance.Plot, device = "png", height = 11, width = 9)

# how many sites were species found?
Supp.Abundance.Plot.Red <- zerotopSp %>%
  filter(pCover > 0) %>%
  # adds line break into string for graphing
  separate(Taxa, into = c("t1","t2"), sep = " ",) %>%
  mutate(t3 = "\n",
         t2 = if_else(is.na(t2),"",t2)) %>%
  unite(t1,t3,t2, col = "Taxa", sep = "") %>%
  ggplot(aes(x = AlphaTag, y = pCover, fill = Taxa)) +
  geom_col(color = "black", fill = "grey",
           position = "dodge", show.legend = FALSE) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 8)) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = "Survey Location",
       y = "Abundance (%)") +
  facet_wrap(~Taxa, scales = "free")
Supp.Abundance.Plot.Red
ggsave(here("Output","PaperFigures","Supp_Red_Sp_Abundance.png"), Supp.Abundance.Plot.Red, device = "png", height = 10, width = 10)



# species richness
zerotopSp %>%
  filter(pCover > 0) %>%
  distinct(CowTagID,Taxa) %>%
  count(CowTagID) # counts taxa in cowtagid

# abundance across sites (occurrences)
High.occur <- zerotopSp %>%
  filter(pCover > 0) %>%
  distinct(CowTagID,Taxa) %>%
  count(Taxa) %>%  # counts cowtagid's where taxa are found
  arrange(desc(n))
#View(High.occur)

##### highest occurring species
High.sp <- (High.occur %>%
  filter(n >=15))$Taxa
zerotopSp %>%
  filter(Taxa %in% High.sp) %>%
  ggplot(aes(x = AlphaTag, y = pCover, fill = Taxa)) +
  geom_col(color = "black", position = "dodge", show.legend = FALSE) +
  theme_bw() + theme(panel.grid.major.x = element_blank()) +
  #scale_fill_manual(values = mypal) +
  facet_wrap(~Taxa, scales = "free")

View(topSp %>%
  arrange(desc(pCover)) %>%
    group_by(Taxa) %>%
    filter(pCover == max(pCover)))

### calculate total percent of species representing majority of plots
High.occur %>%
  filter(n <= 3) %>%
  count() %>%
  mutate(n / 51*100)

Low.occur <- (High.occur %>% filter(n <= 3))$Taxa
View(topSp %>% filter(Taxa %in% Low.occur))





############################
### ANALYZE FUNCTIONAL ENTITY PRESENCE
############################
mypal <- pnw_palette("Bay", n = 13)

# calculate percent cover of species including substrate points
# to get actual abundance of species along gradient
topFE <- myspecies %>%
  select(CowTagID:SpeciesCounts, AlphaTag) %>% # remove unnecessary columns
  select(-Date) %>%
  left_join(traits) %>%
  unite(col = "FE", Taxon_Group, Morph2, Calc, ER, sep = ",", remove = F) %>%
  group_by(CowTagID) %>%
  mutate(totalCounts = sum(SpeciesCounts)) %>%
  filter(Taxa != "Sand" &
           Taxa != "Bare Rock" &
           Taxa != "Rubble") %>%
  select(-Taxa) %>%
  group_by(CowTagID, FE) %>%
  mutate(SpeciesCounts = sum(SpeciesCounts)) %>% # add all of same species together
  distinct() %>% # remove duplicates
  ungroup() %>%
  mutate(pCover = SpeciesCounts / totalCounts*100) %>%
  distinct()
  #filter(pCover > 10) %>%  # select top abundant species


# fixed y axes
# topSp %>%
#   ggplot(aes(x = AlphaTag, y = pCover, fill = Taxa)) +
#   geom_col(color = "black", position = "dodge") +
#   theme_bw() + theme(panel.grid.major.x = element_blank()) +
#   #scale_fill_manual(values = mypal) +
#   facet_wrap(~Taxa)

# put in zero values
fetraits <- topFE %>%
  select(FE:ER) %>%
  distinct()
zerotopFE <- topFE %>%
  select(-c(Taxon_Group:ER)) %>%
  pivot_wider(names_from = FE, values_from = pCover) %>%
  mutate_at(.vars = 5:ncol(.), .funs= ~if_else(is.na(.),0,.)) %>%
  pivot_longer(cols = 5:ncol(.), names_to = "FE", values_to = "pCover") %>%
  left_join(fetraits)

# relative abundance across sites, zeros included
Supp.Abundance.Plot.FE <- zerotopFE %>%
  # adds line break into string for graphing
  # separate(Taxa, into = c("t1","t2"), sep = " ",) %>%
  # mutate(t3 = "\n",
  #        t2 = if_else(is.na(t2),"",t2)) %>%
  # unite(t1,t3,t2, col = "Taxa", sep = "") %>%
  ggplot(aes(x = AlphaTag, y = pCover)) +
  geom_col(color = "black", fill = "grey",
           position = "dodge", show.legend = FALSE) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 8)) +
  labs(x = "Survey Location",
       y = "Abundance (%)") +
  scale_y_continuous(n.breaks = 3) +
  facet_wrap(~FE, scales = "free_y")
Supp.Abundance.Plot.FE
ggsave(here("Output","PaperFigures","Supp_FE_Abundance.png"), Supp.Abundance.Plot.FE, device = "png", height = 11, width = 9)

# how many sites were species found?
Supp.Abundance.Plot.FE.Red <- zerotopFE %>%
  filter(pCover > 0) %>%
  # adds line break into string for graphing
  # separate(Taxa, into = c("t1","t2"), sep = " ",) %>%
  # mutate(t3 = "\n",
  #        t2 = if_else(is.na(t2),"",t2)) %>%
  # unite(t1,t3,t2, col = "Taxa", sep = "") %>%
  ggplot(aes(x = AlphaTag, y = pCover)) +
  geom_col(color = "black", fill = "grey",
           position = "dodge", show.legend = FALSE) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 8)) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = "Survey Location",
       y = "Abundance (%)") +
  facet_wrap(~FE, scales = "free")
Supp.Abundance.Plot.FE.Red
ggsave(here("Output","PaperFigures","Supp_Red_FE_Abundance.png"), Supp.Abundance.Plot.FE.Red, device = "png", height = 10, width = 10)



# FE richness
zerotopFE %>%
  filter(pCover > 0) %>%
  distinct(CowTagID,FE) %>%
  count(CowTagID) # counts FE in cowtagid

# abundance across sites (occurrences)
High.occur <- zerotopFE %>%
  filter(pCover > 0) %>%
  distinct(CowTagID,FE) %>%
  count(FE) %>%  # counts cowtagid's where FE are found
  arrange(desc(n))
#View(High.occur)

##### highest occurring species
High.FE <- (High.occur %>%
              filter(n >=15))$FE
zerotopFE %>%
  filter(FE %in% High.FE) %>%
  ggplot(aes(x = AlphaTag, y = pCover, fill = FE)) +
  geom_col(color = "black", position = "dodge", show.legend = FALSE) +
  theme_bw() + theme(panel.grid.major.x = element_blank()) +
  #scale_fill_manual(values = mypal) +
  facet_wrap(~FE, scales = "free")

View(topFE %>%
       arrange(desc(pCover)) %>%
       group_by(FE) %>%
       filter(pCover == max(pCover)))

### calculate total percent of species representing majority of plots
High.occur %>%
  filter(n <= 3) %>%
  count() %>%
  mutate(n / 51*100)

Low.occur <- (High.occur %>% filter(n <= 3))$FE
View(topFE %>% filter(FE %in% Low.occur))


## get driving species for FE
View(topSp %>%
  select(-c(SpeciesCounts,totalCounts)) %>%
  left_join(traits) %>%
  left_join(fetraits) %>%
  drop_na(FE) %>%
  filter(FE %in% High.FE) %>%
  distinct() %>%
  relocate(pCover, .after = FE) %>%
  arrange(FE, desc(pCover)))


### get proportions of functional traits within FE
fetraits %>%
  count(Taxon_Group) %>%
  mutate(FEprop = n / 22*100)

fetraits %>%
  count(Morph2) %>%
  mutate(FEprop = n / 22*100)

fetraits %>%
  count(Calc) %>%
  mutate(FEprop = n / 22*100)

fetraits %>%
  count(ER) %>%
  mutate(FEprop = n / 22*100)

