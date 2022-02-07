####################################################################
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
#
# Sensitivity analyses based on a decreased number of categories for each trait
#
####################################################################
# set working directory

setwd("/XXXX/")

#define number of axes




############################### FUNCTIONAL SPACE


#load libraries
library('FD')
library('tripack')
library('geometry')
library('matrixStats')


#load data
fes <- read.csv2("Data_FEs_sensitivity.csv", sep=";", dec=",")

fes <- fes[!duplicated(fes), ]

fes_name <- fes[,1]

fes <- fes[,2:dim(fes)[2]]

rownames(fes) <- fes_name 

spe_fes <- read.csv2("Data_Species_FEs_sensitivity.csv", sep=";", dec=",", row.names=1)

#functional space
#gower matrix

d <- gowdis(fes)

fit <- cmdscale(d,eig=TRUE, k=4)


#see variance explained by the axes
cumsum(fit$eig[fit$eig>=0]) / sum(fit$eig[fit$eig>0])

#extract coordinates
fd.coord <- fit$points


################################# functional richness based on a decreased number of categories

ab <- read.csv2("Data_Abundance.csv", sep=";", dec=",", row.names=1)

spe_fes <- read.csv2("Data_Species_FEs_sensitivity.csv", sep=";", dec=",", row.names=1)

sites <- read.table("Data_Sites.txt", sep="\t", header=T, row.names=1)

condition <- c("Ambient", "Low" , "Extreme Low")

ab.conditions <- lapply(condition, function(x) {
  
            quad <- rownames(sites[sites$pH.conditions == x,])
            
            colSums(ab[rownames(ab) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition
  
# convex hull

Fric <- lapply(condition, function (x) {
  
        species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
        fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
        
        m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
        
        ch <- convhulln(m, options = "FA")
        
        chg <- convhulln(fd.coord, options = "FA")
        
        c(length(species), length(species)/72*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
        
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol8D")



############################################# plot convex hull
cols <- c("#CD2626", "#FFA500", "#3A5FCD")
colstr <- c("#CD262670", "#FFA50070", "#3A5FCD70")

names(cols) <- c("Extreme Low", "Low", "Ambient")
names(colstr) <- c("Extreme Low", "Low", "Ambient")


# Supplementary Figure 3. Sensitivity analyses based on a decreased number of categories for each trait. 

n_axes = 4

labels_fig_cv <- c("Sp", "FE", paste0("Vol. ",n_axes,"D"))

tiff(filename="Figure_S3.tif", height=10, width=11, units="cm", compression = c("lzw"), res=300, pointsize=8)

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

# functional richness from north and south venting sites and decreased number of categories 


condition <- levels(sites$Site.and.pH)

condition = c(condition[1],condition[2],condition[5],condition[6],condition[3],condition[4])


ab.conditions <- lapply(condition, function(x) {
  
  quad <- rownames(sites[sites$Site.and.pH == x,])
  
  colSums(ab[rownames(ab) %in% quad,])
  
})#eo lapply

ab.conditions <- do.call(rbind, ab.conditions)

rownames(ab.conditions) = condition

# convex hull

Fric <- lapply(condition, function (x) {
  
  species <- colnames(ab.conditions)[which(ab.conditions[x,] > 0)]
  
  fes_cond <- spe_fes[rownames(spe_fes) %in% species, ]
  
  m <- fd.coord[rownames(fd.coord) %in% fes_cond,]
  
  ch <- convhulln(m, options = "FA")
  
  chg <- convhulln(fd.coord, options = "FA")
  
  c(length(species), length(species)/72*100, dim(m)[1], dim(m)[1]/dim(fd.coord)[1]*100, ch$vol/chg$vol*100)
  
})#eo lapply

names(Fric) = condition
Fric <- do.call(rbind, Fric)
colnames(Fric) <- c("NbSp", "NbSpP", "NbFEs","NbFEsP", "Vol8D")



############################################# plot convex hull
#Supplementary Figure 7. Analysis based on data from north and south venting sampling sites and a decreased number of categories for each trait. 

cols <- c("#CD2626", "#CD2626", "#FFA500", "#FFA500", "#3A5FCD", "#3A5FCD")
colstr <- c("#CD262670", "#CD262670", "#FFA50070", "#FFA50070", "#3A5FCD70","#3A5FCD70")

names(cols) <- c("Extreme Low_N", "Extreme Low_S", "Low_N", "Low_S", "Ambient_N", "Ambient_S")
names(colstr) <-  c("Extreme Low_N", "Extreme Low_S", "Low_N", "Low_S", "Ambient_N", "Ambient_S")

# plot all volumes in distinct plots


tiff(filename="Figure_S7.tif", height=10, width=20, units="cm", compression = c("lzw"), res=300, pointsize=8)


par(mfrow = c(2,6))

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

