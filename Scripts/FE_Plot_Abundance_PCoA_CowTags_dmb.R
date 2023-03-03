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


############################### FUNCTIONAL SPACE


#load libraries
library(tidyverse)
library('FD')
library('tripack')
library('geometry')
library('matrixStats')
library(PNWColors)
library(here)


mypalette <- pnw_palette(name = "Bay", n = 19)
alphapalette <- c(paste0(mypalette[1],"70"), paste0(mypalette[2],"70"), paste0(mypalette[3],"70"),
                  paste0(mypalette[4],"70"), paste0(mypalette[4],"70"), paste0(mypalette[6],"70"),
                  paste0(mypalette[7],"70"), paste0(mypalette[8],"70"), paste0(mypalette[9],"70"),
                  paste0(mypalette[10],"70"), paste0(mypalette[11],"70"), paste0(mypalette[12],"70"),
                  paste0(mypalette[13],"70"), paste0(mypalette[14],"70"), paste0(mypalette[15],"70"),
                  paste0(mypalette[16],"70"), paste0(mypalette[17],"70"), paste0(mypalette[18],"70"),
                  paste0(mypalette[19],"70"))

#load data
ab.sgd <- read_csv(here("Data", "Species_Abundances_wide.csv"))
ab.sgd <- as.data.frame(column_to_rownames(ab.sgd, var = 'CowTagID')) # move tag names to rownames and make data.frame class

fes_traits.sgd <- read_csv(here("Data", "Distinct_FE.csv"))

spe_fes.sgd <- as.data.frame(read_csv(here("Data", "Species_FE.csv"))) %>% rename(FE = fun_entity)

#sites <- read.table("Data/Teixido/Data_Sites.txt", sep="\t", header=T, row.names=1)

# condition is just the cowtagids
condition <- rownames(ab.sgd)

dist <- read_csv(here("Data", "Full_metadata.csv")) %>% filter(Location == "Varari", CowTagID != "V13") %>% select(CowTagID, dist_to_seep_m)
dist <- dist[1:19,]

fd.coord.sgd <- read_csv(here("Data", "FE_4D_coord_dmb.csv"))
fd.coord.sgd <- as.data.frame(column_to_rownames(fd.coord.sgd, var = '...1'))
#fd.coord<- read.csv2 ("Data/Teixido/FE_4D_coord.csv", sep=",", dec=",", row.names=1)


fd.coord.sgd <- as.matrix(fd.coord.sgd)
#fd.coord <- as.matrix(fd.coord)





################################## Data manipulation and arrangements

ab.conditions.sgd <- ab.sgd

################################# compute abundance of FEs for the three conditions

fes.sgd <- levels(as_factor(spe_fes.sgd$FE))

ab.conditions.sgd <- rownames_to_column(ab.conditions.sgd, var = "CowTagID")
ab.conditions.sgd2 <- ab.conditions.sgd %>%
  pivot_longer(names_to = "Taxa", values_to = "pCover", cols = 2:ncol(ab.conditions.sgd)) %>%
  left_join(spe_fes.sgd) %>%
  group_by(CowTagID, FE) %>%
  summarise(pCover = sum(pCover)) %>%
  ungroup()
ab.conditions.sgd <- ab.conditions.sgd2 %>%
  pivot_wider(names_from = FE, values_from = pCover)


######################




# Figure 2. Overall distribution of FE abundance across the functional space


#define number of axes

# n_axes = 4
#
# labels_fig_cv <- c("NbSp", "NbFE", paste0("Vol. ",n_axes,"D"))


CTlevels <- dist %>%
  arrange(desc(dist_to_seep_m)) %>%
  filter(CowTagID != "V13",
         CowTagID != "VSEEP") %>%
  select(CowTagID)
CTlevels <- CTlevels$CowTagID

## relative abundance in ggplot
fig2.fd.sgd <- rownames_to_column(as.data.frame(fd.coord.sgd), var = "FE") %>%
  full_join(ab.conditions.sgd2) %>%
  left_join(dist)
fig2.fd.sgd$CowTagID <- factor(fig2.fd.sgd$CowTagID, levels = CTlevels)
fig2 <- fig2.fd.sgd %>%
  arrange(CowTagID) %>%
  filter(pCover > 0) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(size = pCover,
                 color = CowTagID)) +
  labs(x = "PCoA1", y = "PCoA2") +
  facet_wrap(~CowTagID) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_color_manual(values = alphapalette)
fig2
ggsave(here("Output", "PaperFigures", "FE_abundance_faceted.png"),fig2, width = 7, height = 6)

#View(fig2.fd.sgd)



################################# compute Redundancy defined as Ripley's K


#review data

ab.sgd
spe_fes.sgd
spe_fes.sgd$Taxa = as.character(spe_fes.sgd$Taxa)

# presence/absence k = 1

d <- as.matrix(dist(fd.coord.sgd)) # gets distance matrix for all FE pair combinations

k <- max(d)/100*10 # k = 1%

# count number of species in redundant FE groups and % of community accounted for by that FE

Redundancy <- lapply(condition, function (x) {

  species <- colnames(ab.sgd)[which(ab.sgd[x,] > 0)] # list species in condition

  fes_cond <- unique(as.character(spe_fes.sgd[spe_fes.sgd$Taxa %in% species, "FE"])) # list unique FE from species list

  fe_red <- sapply(fes_cond, function(j) {

    d_red <- d[rownames(d) == j,] # for each row where each row is an FE

    d_red <- d_red[names(d_red) %in% fes_cond] # only select for which have an FE present

    d_red <- d_red[d_red <= k] # looks for FEs with value less than or equal to k

    FE_red <- names(d_red)

    if(length(FE_red) < 2) {


      sp_FE_red <- as.character(spe_fes.sgd[which(spe_fes.sgd$FE == FE_red),]$Taxa)

      NbSpec <- length(sp_FE_red)

      Ab <- sum(ab.sgd[x, sp_FE_red ]) #/2400*100

    } else {

      sp_FE_red <- as.character(spe_fes.sgd[which(spe_fes.sgd$FE %in% FE_red),]$Taxa)

      NbSpec <- length(sp_FE_red)

      Ab <- sum(ab.sgd[x, sp_FE_red ]) #/2400*100

    }#eo ifelse
    c(NbSpec, Ab)

  })#eo lapply

  rownames(fe_red) = c("NbSpec", "Ab")

  return(fe_red)

})#eo lapply

names(Redundancy) = condition







###Supplementary Figure 4. Functional redundancy among pH zones.

tiff(filename="Figure_S4_dmb.tif", height=25, width=20, units="cm", compression = c("lzw"), res=300, pointsize=16)

#cols <- c("#CD2626", "#FFA500", "#3A5FCD")
#names(cols) <- c("Extreme Low", "Low", "Ambient")
cols <- mypalette
names(cols) <- CTlevels

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

# t function transforms
tspe_fes<-t(spe_fes.sgd)
tab.conditions<-t(ab.sgd)


# cowtags wide format, species as rows
#names(tab.conditions)<-c("Species")

tab.conditions<- rownames_to_column(as.data.frame(tab.conditions),var = "Taxa")
vu.fe.condition<- full_join(tab.conditions, spe_fes.sgd)

fes.unique <- vu.fe.condition %>%
  group_by(FE) %>%
  filter(n()==1) # filter by which FE are only observed once

plot_data <- fes.unique %>%
  group_by(FE) %>%
  filter(n()==1) %>%
  ungroup() %>%
  mutate_at(vars(V4:V7), ~if_else(. > 0, 1, 0)) %>% # presence absence as 1's and 0's
  group_by(FE) %>%
  summarise_at(vars(V4:V7), .funs = ~sum(.)/dim(fes.unique)[1]*100) %>%
  ungroup() %>%
  gather(key = condition, value = FE, V4:V7)

mylevels = CTlevels
plot_data
# Supplementary Figure 5. Vulnerability of FEs among pH zones.

vu.fe.plot<-ggplot(plot_data, aes(x = factor(condition, levels=mylevels), y=FE, fill=condition))+
  theme_classic()+
  geom_col(width=0.5)+
  #scale_fill_manual(values=c("#3A5FCD", "#CD2626", "#FFA500"))+
  theme(legend.position="none",
        axis.text.x=element_text(size=14, colour="black"),
        axis.text.y=element_text(size=14, colour="black"),
        axis.title.y=element_text(size = 18, colour="black"),
        axis.line= element_line(size=1))+
  labs(x = "")+
  annotate("text", label = c("n = 49", "n = 40", "n = 20"),
           x=c(1.0, 2.0, 3.0),
           y= c(80, 70, 40), size=4)+
  scale_y_continuous(name="Vulnerability of FE (%)", limits=c(0, 100), breaks=seq (0,100,20), expand = c(0,0))#+
  #scale_x_discrete(labels=c("count_am"="Ambient pH","count_low"="Low pH","count_exlow"= "Extreme Low pH"  ))

vu.fe.plot


#ggsave("Figure_S5.tif", plot=vu.fe.plot, device="tiff", scale=1, height=10, width=11, units=c("cm"), compression = c("lzw"), dpi = 300, pointsize=8)








