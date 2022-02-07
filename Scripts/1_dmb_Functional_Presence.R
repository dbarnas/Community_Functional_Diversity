####################################################################
#### Adaptation by Danielle Barnas
#### Adaptation created on 2022-2-03
#### Code accompanying:
#
# Teixido et al. submitted. Functional biodiversity loss along a natural CO2 gradient. Nature Communications
#
# Script written by: Valeriano Parravacini, Nuria Teixido, Sebastien Villeguer
#
# Functions required and written by Sebastien Villeguer:
# quality_funct_space.R: This function is an updated version of the Appendix S1 associated to Maire et al. 2015 (Global Ecology and Biogeography)
# intersect.R: function to compute convex hull
#
# Code to calculate:
#1) Number of species, number of FEs and functional volume filled by each assemblage
#2) Null model of functional richness
#3) Intersection of the three functional volumes
#4) Taxonomic and functional beta diversity
#5) Distribution of functional trait categories across the functional space
#
####################################################################




############################### FUNCTIONAL SPACE


#### Load required packages
library('tidyverse')
library('here')
library('FD') # includes vegan and geometry
library('tripack')
library('geometry')
library('matrixStats')


############################### DATA
# Load FEs data
fes <- read_csv(here("Data","Teixido","Data_FEs.csv")) #, sep=";", dec=",", row.names=1)

#Load Species and FE data
spe_fes <- read.csv2(here("Data","Teixido","Data_Species_FEs.csv"), sep=";", dec=",") #, row.names=1)

#Load Abundance data per quadrat
ab <- read.csv2(here("Data","Teixido","Data_Abundance.csv"), sep=";", dec=",") %>% rename(Quadrats = X) #, row.names=1)

#Load sites and quadrats

#sites <- read_csv(here("Data","Teixido","Data_Sites.csv")) # brings in column 1 as column
sites <- read.table(here("Data","Teixido","Data_Sites.txt"), sep="\t", header=T, row.names=1) # brings in column 1 as row headings

# Defining pH conditions
condition <- c("Ambient", "Low" , "Extreme Low")

############################### FUNCTIONAL SPACE
# computing  multidimensional functional spaces (2 to 14 D)
#load additional functions

source(here("Scripts","Teixido","quality_funct_space.R"))

#receive error with character columns. will remove to test
fes <- fes %>%
  select_if(is.numeric)

qfs <- quality_funct_space(fes, traits_weights=NULL, nbdim=14, metric="Gower", dendro=FALSE, plot="quality_funct_space")

# quality of spaces (low meanSD = high quality)
round( qfs$meanSD , 4)

# keeping coordiantes on the 4 dimensions, meanSD<0.004
fd.coord <- qfs$details_funct_space$mat_coord[,1:4]

write.csv(fd.coord, here("Data","Teixido","FE_4D_coord.csv")) #to use it for further analyses

#see variance explained by the PCoA axes
gower<-qfs$details_funct_space$mat_dissim

fit <- cmdscale(gower,eig=TRUE, k=4) # PCoA

# variance explained by the axes
cumsum(fit$eig[fit$eig>=0]) / sum(fit$eig[fit$eig>0])




#################################  FUNCTIONAL RICHNESS

# Data manipulation and arrangements

ab.conditions <- lapply(condition, function(x) { ### returns all zeros - need to reassess

  quad <- rownames(sites[(sites$pH_conditions == x),])

  colSums(ab[rownames(ab) %in% quad,])

})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions) # rbind all lists of ab.conditions into table

rownames(ab.conditions) = condition # rename rownames as conditions: Ambient, Low, Extreme Low

#### Calculate convex hull

Fric <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]

  ch <- convhulln(m, options = "FA")

  chg <- convhulln(fd.coord, options = "FA")

  c(length(species), length(species)/72*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)

})#eo lapply

names(Fric) = condition

# Fric contains the number of species(NbSp) and FEs (NbFEs), relative percentages (NbSpP,NbFEsP ) , and the volume among the 3 pH zones
Fric <- do.call(rbind, Fric)

colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol8D")



############################################# plot convex hull
cols <- c("#CD2626", "#FFA500", "#3A5FCD")
colstr <- c("#CD262670", "#FFA50070", "#3A5FCD70")

names(cols) <- c("Extreme Low", "Low", "Ambient")
names(colstr) <- c("Extreme Low", "Low", "Ambient")


# Figure 1. Species and functional diversity changes among pH zones.
# All volumes in distinct plots

n_axes = 4
labels_fig_cv <- c("Sp", "FE", paste0("Vol. ",n_axes,"D"))

tiff(filename="Figure_1.tif", height=10, width=11, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(mfrow = c(2,3))

Fric <- Fric[,c(2,4,5)]

colnames(Fric) = labels_fig_cv

for (i in condition) {

  midpoints <-  barplot(Fric[i,], ylab="Relative richness (%)", col=cols[i], ylim=c(0,105), main=paste0(i, " pH"), col.main=cols[i], cex.main=1.2 )

  if(Fric[3,3] > 1) {

    lab = round(Fric[i,],0)

  } else {

    lab = round(Fric[i,],2)

  }#eo ifelse

  text(midpoints, Fric[i,]+8, labels=lab, col=cols[i], cex=1.2, font=2)

}#eo for barplot

for (i in condition) {

  species <- colnames(ab.conditions)[which(ab.conditions[i,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]

  tr <-tri.mesh(m[,1],m[,2])
  ch <- convex.hull(tr)

  plot(fd.coord[,1], fd.coord[,2], xlab = "PCoA 1", ylab = "PCoA 2", type="n")

  polygon(ch, col=colstr[i], border=cols[i])
  points(m[,1:2], pch = 16, col=cols[i])

}#eo for convex

dev.off()

# Supplementary Figure 1. Intersection of the three functional volumes among pH zones.
# all volumes in 1 figure,

tiff(filename ="Figure_S1.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)

par(mfrow = c(1,1))
plot(fd.coord[,1], fd.coord[,2], xlab = "PCoA 1", ylab = "PCoA 2", type="n")

condition = c("Ambient", "Low", "Extreme Low")



for (i in condition) {

  species <- colnames(ab.conditions)[which(ab.conditions[i,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond, ]

  tr <-tri.mesh(m[,1],m[,2])
  ch <- convex.hull(tr)


  polygon(ch, col=colstr[i], border=cols[i])
  points(m[,1:2], pch = 16, col=cols[i])

  legend("topleft",legend=c("Extreme Low pH", "Low pH", "Ambient pH"),col=cols ,bty="n",pch=rep(16, length(lab)))


}

dev.off()



############### null model of Functional richness among pH zones

n_perm = 100

spe_fes_r = spe_fes

Fric_perm <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]


  perm <- sapply((1:n_perm), function (z) {


    spe_fes_r$FE <- sample(spe_fes$FE)

    fes_cond <- spe_fes_r[rownames(spe_fes_r) %in% species, ]

    m <- fd.coord[rownames(fd.coord) %in% fes_cond,]

    ch <- convhulln(m, options = "FA")

    chg <- convhulln(fd.coord, options = "FA")

    c(length(species), length(species)/72*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)



  })#eo sapply

  rownames(perm) <- c("NbSp", "NbSpP", "NbFE", "NbFEP", "Vol")


  perm

})#eo lapply

names(Fric_perm) = condition



Fric_perm_Q <- lapply(Fric_perm, function (x) {

  rowQuantiles(x, probs=c(0.05, 0.95))

})#eo lapply



Fric = as.data.frame(Fric)

Fric$lowerFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,1] })#eo sapply
Fric$upperFE <- sapply(condition, function (x) { Fric_perm_Q[[x]][3,2] })#eo sapply
Fric$lowerVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,1] })#eo sapply
Fric$upperVol <- sapply(condition, function (x) { Fric_perm_Q[[x]][5,2] })#eo sapply
Fric$cond <- condition
condition <- factor(condition, levels = c("Ambient", "Low", "Extreme Low"))

Fric$cond <- as.factor(condition)
levels(Fric$cond)
colnames(Fric) <- c("NbSp", "NbFE", "Vol8D", "lowerFE", "upperFE", "lowerVol", "upperVol", "cond")


#Plot the null model
#Supplementary Figure 2. Null model of functional richness (functional volume) among pH zones.


tiff(filename="Figure_S2.tif", height=10, width=10, units="cm", compression = c("lzw"), res=300, pointsize=8)


plot(Vol8D ~ cond, data=Fric, border="white", xlab="", ylab="Relative Richness (%)", ylim=c(0,100))
points(Vol8D ~ cond, data=Fric, pch=16, col=cols[c(3,2,1)], cex=2)

Ambient <- rbind(Fric[1,], Fric[1,])
Ambient$Vol8D <- c(Ambient$lowerVol[1], Ambient$upperVol[1])
lines(Vol8D ~ cond, data=Ambient, lwd=3, col=cols["Ambient"])

Low <- rbind(Fric[2,], Fric[2,])
Low$Vol8D <- c(Low$lowerVol[1], Low$upperVol[1])
lines(Vol8D ~ cond, data=Low, lwd=3, col=cols["Low"])

ELow<- rbind(Fric[3,], Fric[3,])
ELow$Vol8D <- c(ELow$lowerVol[1], ELow$upperVol[1])
lines(Vol8D ~ cond, data=ELow, lwd=3, col=cols["Extreme Low"])

dev.off()


################################################ Convex Hull Intersect

#load intersect function to compute convex hull (vertices + volume) of two set of points and their intersection



source(here("Scripts","Teixido","intersect.R"))


mat_int <- Fric <- lapply(condition, function (x) {

  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]

  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]

  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]

  return(m)

})#eo lapply

names(mat_int) = condition

###############intersect Ambient with Low

Amb_int_Low <- CHVintersect(mat_int[["Ambient"]],mat_int[["Low"]])

#pergcentage of the Low volume within Ambient
Amb_int_Low$vol[3]/Amb_int_Low$vol[2]




###############intersect Ambient with Extreme Low

Amb_int_Ex_Low <- CHVintersect(mat_int[["Ambient"]],mat_int[["Extreme Low"]])

#pergcentage of the Extreme Low volume within Ambient
Amb_int_Ex_Low$vol[3]/Amb_int_Ex_Low$vol[2]



###############intersect Low with Extreme Low

Low_int_Ex_Low <- CHVintersect(mat_int[["Low"]],mat_int[["Extreme Low"]])

#pergcentage of the Extreme Low volume within Ambient
Low_int_Ex_Low$vol[3]/Low_int_Ex_Low$vol[2]





########################################### BETA DIVERSITY

library('betapart')


###### taxonomic (Jaccard)
ab.conditions[which(ab.conditions>0)] = 1
bata.taxo <- beta.pair(ab.conditions, index.family="jaccard")


###### functional (Jaccard like)

# Compute abundances of FEs for the three conditions

# Load again the spe_fes matrix, 2 column variables

spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=",")

fes <- levels(spe_fes$FE)
ab.fe.conditions <- lapply(condition, function (z) {
  abund.fes <-  sapply(fes, function (x) {
    spec <- as.character(spe_fes[which(spe_fes$FE == x),]$Species)
    sum(ab.conditions[z,spec])
  })#eo sapply
  abund.fes
})#eo lapply

names(ab.fe.conditions) = condition

ab.fe.conditions <- do.call(rbind, ab.fe.conditions)
ab.fe.conditions[which(ab.fe.conditions>0)] = 1

# true functional beta
beta.fun <- functional.beta.pair(ab.fe.conditions, fd.coord, index.family="jaccard")

#######Plot categories of the 15 functional traits across the functional space

#get data to plot the traits
fes <- read.csv2("Data_FEs.csv", sep=";", dec=",", row.names=1)

spe_fes <- read.csv2("Data_Species_FEs.csv", sep=";", dec=",", row.names=1)


###### Supplementary Figure 6. Distribution of functional trait categories across the functional space

tiff(filename="Figure_S6.tif", height=20, width=30, units="cm", compression = c("lzw"), res=300, pointsize=10)


ftr <- colnames(fes)

par(mfrow=c(3,5))



for (i in ftr) {

  lab <- as.factor(sort(unique(fes[,i])))

  plot(fd.coord[,1], fd.coord[,2], pch=16, cex=1.2, col = as.numeric(fes[,i]), xlim = c(-0.3, 0.7),  main=gsub("."," ", i, fixed = T), xlab="PCoA 1", ylab="PCoA 2")
  legend(x=0.45, y=0.35, legend=lab, pch=rep(16, length(lab)), col=as.numeric(as.factor(lab)), bty = "n")


}

dev.off()









