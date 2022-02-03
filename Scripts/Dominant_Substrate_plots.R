#### Begin processing community composition Dominant Substrate surveys from Varari in August 2021


### LIBRARIES
library(tidyverse)
library(here)

library(ggmap)
library(maps)
library(maptools)
library(ggsn)
library(ggrepel)
library(scatterpie)


### READ IN DATA
survey<-read_csv(here("Data","Varari_Sub_Survey.csv"))
gps<-read_csv(here("Data","SandwichLocations_swapped_18and17.csv"))
nutrient<-read_csv(here("Data","Nutrient_Snap_May2021.csv"))

### PERCENT COVER OF DOMINANT SUBSTRATE
sp_count<-survey %>% 
  pivot_longer(cols = Sub.1:Sub.25, names_to = 'quadID', values_to = 'substrate') %>% 
  group_by(Top_Plate_ID,plate.num, substrate) %>% 
  count(substrate) %>% 
  rename(spCount = n) %>% 
  mutate(percentCover = spCount/100 * 100) # percent cover

halfsurvey<-sp_count %>%  # two sandwiches only surveyed 0.5m2
  filter(Top_Plate_ID == 'V10' | Top_Plate_ID == 'V8') %>% 
  mutate(percentCover = percentCover * 2)

sp_count<-sp_count %>%
  filter(Top_Plate_ID != 'V10' & Top_Plate_ID !='V8') %>% 
  rbind(halfsurvey) %>% 
  filter(substrate != 'photo?') # temporarily remove unidentified substrate

### Relevel based on relative distance to SGTree
sp_count$Top_Plate_ID<-factor(x = sp_count$Top_Plate_ID, levels = c('V14','V20','V17','V10','V9','V18',
                                             'V8','V7','V15','V6','V4','V16','V3',
                                             'V2','V1','V13')) # relevel top plate id's for plotting


### GRAPHING 

## Stacked bar plot of % cover
sp_count %>%
  ggplot(aes(x = Top_Plate_ID, y = percentCover, fill = substrate))+
  geom_col()+
  theme_bw()

## Line plot
sp_count %>% 
  ggplot(aes(x = Top_Plate_ID, y = percentCover, group = substrate, color = substrate))+
  geom_line()+
  theme_bw()

# join nutrient data
sp_count <-sp_count %>% 
  left_join(nutrient, by = 'Top_Plate_ID')

# join gps locations
gps_count<-sp_count %>%
  left_join(gps, by = 'Top_Plate_ID')
 


###  MAPPING

### read in API key for maps
API<-names(read_table(here("Data","API.txt")))
register_google(key = API) ### use your own API


lon = gps_count$lon
lat = gps_count$lat

myMap <-get_map(location = c(lon = -149.8999, lat = -17.5405), # set center point
                maptype = 'satellite',
                zoom = 19) # zooms enough to see all 16 points

# Plot site location
ggmap(myMap) +
  geom_point(data = gps_count,
             aes(x = lon, y = lat)) +
  geom_label(data = gps_count,
                   mapping = aes(x = lon, y = lat, 
                                 label = Top_Plate_ID))


### nutrients gradient on map

ggmap(myMap)+
  geom_point(data = gps_count, mapping = aes(x=lon, y=lat, color = Silicate), alpha = .60)+
  scale_color_gradient(low = "yellow", high = "blue")


### pie charts on map ###
### nutrients pie charts on map
pie.cat.nut<-names(gps_count[16:19])
gps_count[16:19] <- lapply(gps_count[16:19], as.numeric)
#gps_count$radius <- 5

ggmap(myMap) +
  geom_scatterpie(data = gps_count,
                  aes(x = lon, y = lat, group = Top_Plate_ID),
                      #r = radius),
                  cols = pie.cat.nut,
                  color = NA) # remove black outline within pie charts

### substrate pie charts on map
wide_gps<-gps_count %>% 
  select(Top_Plate_ID,plate.num,substrate,percentCover,Location,lat,lon) %>% 
  pivot_wider(names_from = substrate, values_from = percentCover)

pie.cat.sub<-names(wide_gps[6:14])
wide_gps[6:14] <- lapply(wide_gps[6:14], as.numeric)

wide_gps<-wide_gps %>% 
  replace_na(list('bare rock' = 0, 'CCA' = 0, 'sand' = 0, 'turf' = 0,
                  'macroalgae' = 0, 'stony coral'= 0, 'soft coral' = 0,
                  'anemone' = 0, 'sponge' = 0))

#wide_gps$radius<-6*abs(rnorm(nrow(wide_gps)))
ggmap(myMap) +
  geom_scatterpie(data = wide_gps,
                  aes(x = lon, y = lat, group = Top_Plate_ID),
                  cols = pie.cat.sub,
                 # r = radius,
                  color = NA) # remove black outline within pie charts

###################################
###### nMDS AND perMANOVA
###################################

# perMANOVA: is there a difference in response variables among groups?

# 3 steps in using classification functions:
# 1. Get the coefficients and constant of the classification equation for each group
# 2. Get a classification score for each observation for each group 
#     a. Calculated by using the actual values for each variable to solve the classification for that group
# 3. Classify each observation into the group it matches most closely 
#     a. This may or may not be the actual group from which the observation came 
#     b. You have to know what group the observation goes into. 
#        You have to know the classifications already, and then you can see how good the DFA is at fitting the data

### nMDS plot
library(vegan)

# select numerical and single category data only 
mds_data<-wide_gps[-c(2:5)]

#create the ordination output using bray curtis dissimilarity matrix
# numerical data only
ord<-metaMDS(mds_data[,-1],k=2, distance='bray') 

#let's look at the stress with k=2 dimensions. Is it < 0.3? 
ord$stress

# Let's look at the stress plot
stressplot(ord)
# looks like a good fit, want to minimize scatter

# basic plot
ordiplot(ord) # dots represent sites (sandwich locations) and + represents substrates

# add text
ordiplot(ord, type = 'text')

## subset by relative distance to sgd
highsgd<-mds_data %>% 
  filter(Top_Plate_ID == 'V14' | Top_Plate_ID == 'V20' | Top_Plate_ID == 'V17') %>% 
  mutate(rel_dist='adjacent')
medsgd<-mds_data %>% 
  filter(Top_Plate_ID == 'V10' | Top_Plate_ID == 'V9' | Top_Plate_ID == 'V18' | Top_Plate_ID == 'V8') %>% 
  mutate(rel_dist='near')
modsgd<-mds_data %>% 
  filter(Top_Plate_ID == 'V7' | Top_Plate_ID == 'V15' | Top_Plate_ID == 'V6') %>% 
  mutate(rel_dist='moderate')
lowsgd<-mds_data %>% 
  filter(Top_Plate_ID == 'V4' | Top_Plate_ID == 'V16' | Top_Plate_ID == 'V3' | Top_Plate_ID == 'V2' | 
         Top_Plate_ID == 'V1' | Top_Plate_ID == 'V13') %>% 
  mutate(rel_dist='far')

mds_data<-rbind(highsgd,medsgd,modsgd,lowsgd)

# let's make a better plot
# play with the x and y lim to get a graph that best shows differences
plot(1, # creates empty plot frame
     type='n', 
     xlim=c(-1,1), ylim=c(-0.75,0.75), 
     xlab='nMDS1', 
     ylab='nMDS2',
     xaxt='n', yaxt='n')



# build points on frame
points(ord$points[mds_data$rel_dist=='adjacent',1],ord$points[mds_data$rel_dist=='adjacent',2], 
       pch=19, col='purple', cex=2)

points(ord$points[mds_data$rel_dist=='near',1],ord$points[mds_data$rel_dist=='near',2], 
       pch=19, col='darkgreen', cex=2)

points(ord$points[mds_data$rel_dist=='moderate',1],ord$points[mds_data$rel_dist=='moderate',2], 
       pch=19, col='orange', cex=2)

points(ord$points[mds_data$rel_dist=='far',1],ord$points[mds_data$rel_dist=='far',2], 
       pch=19, col='blue', cex=2)


# let's add a circle around all points by groups
ordiellipse(ord, 
            groups=mds_data$rel_dist, 
            label=F, 
            kind='ehull', 
            #kind = 'sd', # to make the circles standard deviation, comment out the above two
            border='white', 
            col=c('purple','blue','orange','darkgreen'), # adjacent, far, near, moderate
            lwd=2, 
            draw ='polygon')

# can add or remove labels 
#add a legend with stress
legend('topleft', legend = paste('2D stress = ', round(ord$stress,2)), bty='n')
#add a Site legend
legend('bottomleft',legend=c('Adjacent to SGD','Near SGD','Moderate Dist to SGD','Far from SGD'),
       col=c('purple','darkgreen','orange','blue'), pch=19, bty='n')


#We could formally test whether sites are different using a perMANOVA. 
#This requires the adonis command in the vegan package

#numeric data only
permanovamodel<-adonis(mds_data[,2:10]~rel_dist, mds_data, permutations = 999, 
                       method="bray")
permanovamodel
# just barely significant between four groups
# p = 0.046

#If we are to trust the results of the permanova, then we have to assume that the dispersion among
#data is the same in each group. We can test with assumption with a PermDisp test:
disper<-vegdist(mds_data[,2:10])
betadisper(disper, mds_data$rel_dist)
#Look at the Average distance to median...these numbers should be reasonably similar
#A rule of thumb is that one number should not be twice as high as any other

#An option for doing post-hoc pairwise comparisons in R
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis(mds_data[,2:10], mds_data$rel_dist, perm=999)
# p < 0.05 for adjacent vs near and adjacent vs far

#Get coefficients to see which species are most important in explaining site differences:
#permanovamodel$coefficients



######################################
### Other plots
######################################

# nutrient data
sp_longer<-sp_count %>% 
  pivot_longer(cols = c(Phosphate,Silicate,Nitrite,Ammonia), 
               names_to = 'Values',
               values_to = 'Nutrients')
sp_longer %>% 
  ggplot(aes(x=Top_Plate_ID,  y = Nutrients)) +
  geom_col(position = 'dodge') +
  facet_wrap(~Values, scales = 'free') +
  theme_bw()






