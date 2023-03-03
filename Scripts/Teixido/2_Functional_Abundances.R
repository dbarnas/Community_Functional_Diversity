####################################################################
# Teixido et al. submitted. Functional biodiversity along a natural CO2 gradient. Nature Communications
#
# Script written by: Valeriano Parravacini, Nuria Teixido, Sebastien Villeguer
#
# Code to calculate:
#1) Functional entities with abundances
#2) Functional redundancy defined as Ripley's K
#3) Functional vulnerability
#4) Relative abundance of functional trait categories among pH zones
#
####################################################################

# set working directory

setwd("XXXX/")



############################### FUNCTIONAL SPACE


#load libraries
library('FD')
library('tripack')
library('geometry')
library('matrixStats')


#load data
ab <- read.csv2("Data/Teixido/Data_Abundance.csv", sep=";", dec=",", row.names=1)
#fes_traits <- read.csv2("Data/Teixido/Data_FEs.csv", sep=";", dec=",", row.names=1)
fes_traits <- read_csv("Data/Teixido/Data_FEs.csv")
spe_fes <- read.csv2("Data/Teixido/Data_Species_FEs.csv", sep=";", dec=",")
sites <- read.table("Data/Teixido/Data_Sites.txt", sep="\t", header=T, row.names=1)

condition <- c("Ambient", "Low" , "Extreme Low")

fd.coord<- read.csv2 ("Data/Teixido/FE_4D_coord.csv", sep=",", dec=",", row.names=1)


fd.coord<- as.matrix(fd.coord)


################################## Data manipulation and arrangements

ab.conditions <- lapply(condition, function(x) {

  quad <- rownames(sites[sites$pH.conditions == x,])

  colSums(ab[rownames(ab) %in% quad,])

})#eo lapply


ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

ab.conditions <- ab.conditions/2400*100 #number of quadrats 24 per condition and expressed as %

################################# compute abundance of FEs for the three conditions

fes <- levels(as_factor(spe_fes$FE))

ab.fe.conditions <- lapply(condition, function (z) {

                       abund.fes <-  sapply(fes, function (x) {

                                            spec <- as.character(spe_fes[which(spe_fes$FE == x),]$Species)

                                            sum(ab.conditions[z,spec])


                                     })#eo sapply

                       abund.fes

})#eo lapply

names(ab.fe.conditions) = condition

ab.fe.conditions <- do.call(rbind, ab.fe.conditions)

######################


# Figure 2. Overall distribution of FE abundance across the functional space


#define number of axes

n_axes = 4

labels_fig_cv <- c("NbSp", "NbFE", paste0("Vol. ",n_axes,"D"))


tiff(filename="Figure_2.tif", height=6, width=15, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(mfrow=c(1,3))


plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Ambient pH", col.main="#3A5FCD")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#3A5FCD" , bg="#3A5FCD70", cex=sqrt(ab.fe.conditions["Ambient",]/100*30)*3)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Low pH", col.main="#FFA500")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#FFA500" , bg="#FFA50070", cex=sqrt(ab.fe.conditions["Low",]/100*30)*3)

plot(fd.coord[,1], fd.coord[,2], xlab="PCoA1", ylab="PCoA2", type="n", main="Extreme low pH", col.main="#CD2626")
points(fd.coord[,1], fd.coord[,2], pch=21, col="#CD2626" , bg="#CD262670", cex=sqrt(ab.fe.conditions["Extreme Low",]/100*30)*3)


dev.off()



################################# compute Redundancy defined as Ripley's K


#get data

ab <- read.csv2("Data/Teixido/Data_Abundance.csv", sep=";", dec=",", row.names=1)

sites <- read.table("Data/Teixido/Data_Sites.txt", sep="\t", header=T, row.names=1)

condition <- c("Ambient", "Low" , "Extreme Low")

spe_fes <- read.csv2("Data/Teixido/Data_Species_FEs.csv", sep=";", dec=",")

spe_fes$Species = as.character(spe_fes$Species)

# presence/absence k = 1

d <- as.matrix(dist(fd.coord))

k <- max(d)/100*10 # k = 1%


Redundancy <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)] # list species in condition

  fes_cond <- unique(as.character(spe_fes[spe_fes$Species %in% species, "FE"])) # list unique FE from species list

  fe_red <- sapply(fes_cond, function(j) {

    d_red <- d[rownames(d) == j,] # for each row (mean to have 1 row per FE?)

    d_red <- d_red[names(d_red) %in% fes_cond] # only select for which have an FE present

    d_red <- d_red[d_red <= k] # assigns same value to all?

    FE_red <- names(d_red)

    if(length(FE_red) < 2) {

      sp_FE_red <- as.character(spe_fes[which(spe_fes$FE == FE_red),]$Species)

      NbSpec <- length(sp_FE_red)

      Ab <- sum(ab.conditions[x, sp_FE_red ])/2400*100

    } else {

      sp_FE_red <- as.character(spe_fes[which(spe_fes$FE %in% FE_red),]$Species)

      NbSpec <- length(sp_FE_red)

      Ab <- sum(ab.conditions[x, sp_FE_red ])/2400*100

    }#eo ifelse
    c(NbSpec, Ab)

  })#eo lapply

  rownames(fe_red) = c("NbSpec", "Ab")

  return(fe_red)

})#eo lapply

names(Redundancy) = condition


###Supplementary Figure 4. Functional redundancy among pH zones.

tiff(filename="Figure_S4.tif", height=25, width=20, units="cm", compression = c("lzw"), res=300, pointsize=16)

cols <- c("#CD2626", "#FFA500", "#3A5FCD")
names(cols) <- c("Extreme Low", "Low", "Ambient")


par(mfrow=c(3,2))


for (i in condition) {


  dat = Redundancy[[i]]

  dat_s <- dat[,order(dat[1,], decreasing = T)]

  plot(dat_s[1,], xlab="Rank of Functional Entity", ylab="# species", type="n", main=i, xlim=c(0,68))

  lines(dat_s[1,], col=cols[i], lwd=3)

  dat_ab <- dat[,order(dat[2,], decreasing = T)]

  plot(dat_ab[2,], xlab="Rank of Functional Entity", ylab="Abundance (%)", type="n", main=i, xlim=c(0,68))

  lines(dat_ab[2,], col=cols[i], lwd=3)

}

dev.off()


###### Vulnerability of FEs: The FEs with only 1 species

library('tidyverse')


tspe_fes<-t(spe_fes)
tab.conditions<-t(ab.conditions)

#names(tab.conditions)<-c("Species")
tab.conditions<- rownames_to_column(as.data.frame(tab.conditions),var = "Species")
vu.fe.condition<- full_join(tab.conditions, spe_fes)

fes.unique <- vu.fe.condition %>%
  group_by(FE) %>%
  filter(n()==1)

plot_data=fes.unique %>%
  group_by(FE) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate(Ambient_bin = ifelse(Ambient >0, 1, 0),
         Low_bin = ifelse(Low >0, 1, 0),
         ExLow_bin = ifelse(`Extreme Low` >0, 1, 0)) %>%
  group_by(FE) %>%
  summarise(count_am = (sum(Ambient_bin)/dim(fes.unique)[1])*100,
            count_low = (sum(Low_bin)/dim(fes.unique)[1])*100,
            count_exlow = (sum(ExLow_bin)/dim(fes.unique)[1])*100) %>%
  ungroup() %>%
  gather(key = condition, value = FE, count_am,count_low,count_exlow)

mylevels = c("count_am", "count_low", "count_exlow")
plot_data
# Supplementary Figure 5. Vulnerability of FEs among pH zones.

vu.fe.plot<-ggplot(plot_data, aes(x = factor(condition, levels=mylevels), y=FE, fill=condition))+
  theme_classic()+
  geom_col(width=0.5)+
  scale_fill_manual(values=c("#3A5FCD", "#CD2626", "#FFA500"))+
  theme(legend.position="none",
        axis.text.x=element_text(size=14, colour="black"),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size = 18, colour="black"),
        axis.line= element_line(size=1))+
  labs(x = "")+
  annotate("text", label = c("n = 49", "n = 40", "n = 20"),
           x=c(1.0, 2.0, 3.0),
           y= c(80, 70, 40), size=4)+
  scale_y_continuous(name="Vulnerability of FE (%)", limits=c(0, 100), breaks=seq (0,100,20), expand = c(0,0))+
  scale_x_discrete(labels=c("count_am"="Ambient pH","count_low"="Low pH","count_exlow"= "Extreme Low pH"  ))

vu.fe.plot


ggsave("Figure_S5.tif", plot=vu.fe.plot, device="tiff", scale=1, height=10, width=11, units=c("cm"), compression = c("lzw"), dpi = 300, pointsize=8)


#######################################################


# computing relative abundances of trait categories among pH conditions

# list to store results
ab.trait_condition<-list()


# ordering FEs in trait table as in abundance table
fes_traits<-fes_traits[ colnames(ab.fe.conditions), ]


# for each trait, computing total relative abundance for each modality among FEs among pH zones
for (t in colnames(fes_traits) )
{
  # levels of trait t
  levels_t<-as.character( sort( unique(fes_traits[,t]) ) )

  # empty vectors to store results
  trait_t<-c()
  condition_t<-c()
  rel_ab_t<-c()

  # computing relative abundance of each trait modality in each condition
  for (i in condition)
    for (j in levels_t)
      {
        trait_t<-c(trait_t, j)
        condition_t<-c(condition_t, i)
        rel_ab_t<-c(rel_ab_t, sum( ab.fe.conditions[ i , which(fes_traits[,t]==j) ] )  )
  }# end of i,j

  # setting correcdt order for levels of conditions
  condition_t <- factor(condition_t, levels = condition )

  # storing results in a dataframe
  ab.trait_condition[[t]]<-data.frame(  trait_val= trait_t,  condition= condition_t, rel_ab_FE=rel_ab_t )

}# end of t


# Figure 3. Change in relative abundance of functional trait categories along the pH gradient.
#ggplot does not like list of dataframes; thus we code 15 individual plots and later we grouped together

# define color palette for 13, 7, 6, 5, 3, 2 trait categories

cc13<-c("#a6cee3","#1f78b4","#b2df8a", "#33a02c", "#fb9a99", "#e31a1c","#fdbf6f","#ff7f00", "#cab2d6",  "#6a3d9a", "#ffff99","#b15928", "#1c1b1b")
cc7<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#ffb24f")
cc6<-c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c")
cc5<-c("#e31a1c", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
cc3<-c("#ffb24f",  "#1f78b4","#33a02c")
cc2<-c( "#a6cee3","#1f78b4")


#plot1: Morphological.form, 13 trait values, cc13
plot1 <-ggplot(data=ab.trait_condition$Morphological.form, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  guides(fill=guide_legend(ncol=2))+
  labs( fill="Traits")+
  ggtitle("Morphological form")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc13)

plot1

#plot2: Solitary.Colonial, 3 trait values, cc3
plot2 <-ggplot(data=ab.trait_condition$Solitary.Colonial, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Solitary-Colonial")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc3)

plot2

#plot 3 Max.Longevity, 7 trait values, cc7
plot3 <-ggplot(data=ab.trait_condition$Max.Longevity, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Max Longevity")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc7)

plot3

#plot 4 Height, 5 trait values, cc5

plot4 <-ggplot(data=ab.trait_condition$Height, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Height")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc5)

plot4


#plot 5 Width, 6 trait values, cc6

plot5 <-ggplot(data=ab.trait_condition$Width, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Width")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc6)

plot5

#plot 6 Epibiosis, 3 trait values, cc3
plot6 <-ggplot(data=ab.trait_condition$Epibiosis, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Epibiosis")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc3)

plot6

#plot7: Energetic.resource, 3 trait values, cc3
plot7 <-ggplot(data=ab.trait_condition$Energetic.resource, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(y="Relative Abundance of Functional Categories (%)", fill="Traits")+
  ggtitle("Energetic resource")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc3)

plot7


#plot 8 Major.photosynthetic.pigments, 7 trait values, cc7
plot8 <-ggplot(data=ab.trait_condition$Major.photosynthetic.pigments, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Photosynthetic pigments")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc7)

plot8

#plot 9 Feeding, 5 trait values, cc5

plot9 <-ggplot(data=ab.trait_condition$Feeding, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Feeding ")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc5)

plot9

#plot 10 Age.reproductive.maturity, 6 trait values, cc6

plot10 <-ggplot(data=ab.trait_condition$Age.reproductive.maturity, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(fill="Traits")+
  ggtitle("Age reproductive maturity")+
  scale_y_continuous( expand = c(0,0))+
  scale_fill_manual(values = cc6)

plot10


#plot11: Asexual.Reproduction, 2 trait values, cc2

plot11<-ggplot(data=ab.trait_condition$Asexual.Reproduction, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Asexual Reproduction")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc2)
plot11

#plot12: Growth.rates, 5 trait values, cc5

plot12<-ggplot(data=ab.trait_condition$Growth.rates, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        legend.title=element_blank())+
  labs( fill="Traits")+
  ggtitle("Growth rates")+
  scale_y_continuous( expand = c(0,0))+
  #scale_x_discrete(labels=c("Ambient"="Ambient pH", "Low"="Low pH", "Extreme Low"= "Extreme Low pH"))+
  scale_fill_manual(values = cc5)
plot12

#plot13: Calcification, 5 trait values, cc5

plot13<-ggplot(data=ab.trait_condition$Calcification, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Calcification")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc5)
plot13


#plot14: Chemical.defenses, 2 trait values, cc2

plot14<-ggplot(data=ab.trait_condition$Chemical.defenses, aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Traits")+
  ggtitle("Chemical defenses")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient ", "Low"="Low ", "Extreme Low"= "Extreme Low "))+
  scale_fill_manual(values = cc2)
plot14

#plot15: Mobility, 2 trait values, cc2

plot15<-ggplot(data=ab.trait_condition$Mobility , aes( fill=trait_val, y=rel_ab_FE, x=condition ))+
  geom_bar(stat='identity', width=0.5 )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.x=element_text(angle=60,hjust=1, size=12),
        plot.title = element_text(size = 12, face = "plain", hjust=0.5),
        axis.title.x = element_text(size=14),
        legend.title=element_blank())+
  labs(x="pH conditions",  fill="Trait Value")+
  ggtitle("Mobility")+
  scale_y_continuous( expand = c(0,0))+
  scale_x_discrete(labels=c("Ambient"="Ambient", "Low"="Low", "Extreme Low"= "Extreme Low"))+
  scale_fill_manual(values = cc2)
plot15



library ('patchwork')



all.plots.vert<-(plot1|plot2|plot3 ) / (plot4|plot5| plot6 )/(plot7| plot8 |plot9)/
  (plot10| plot11|plot12) / (plot13| plot14| plot15)
all.plots.vert

ggsave("Figure_3.tif", plot= all.plots.vert, device="tiff", height=25, width=20, units="cm", dpi=300)





###########################################################################################################################
