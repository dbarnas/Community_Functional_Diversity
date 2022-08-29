#### Species Diversity across SGD zones ####
### Created by Danielle Barnas
### Created on Aug 29, 2022

#### Load Libraries ####
library(tidyverse)
library(here)
library(stats) # lm() and anova()
library(emmeans)
library(agricolae) # HSD.test()

#### Read in Data ####
spcomp <- read_csv(here("Data","Surveys","Species_Composition_2022.csv"))
cluster <- read_csv(here("Data", "Surveys", "Cluster_metadata.csv"))
turbcluster <- read_csv(here("Data", "Surveys", "Cluster_turb_metadata.csv"))

#### CLUSTERS FROM FULL BIOGEOCHEM ####

### Clean data ###
spcluster <- spcomp %>%
  filter(Location == "Varari") %>%
  select(-c(PhotoNum, Notes)) %>%
  left_join(cluster)
head(spcluster)


### calculate species richness ###
richness <- spcluster %>%
  select(Location, CowTagID, Taxa, V_cluster) %>%
  group_by(Location, CowTagID, V_cluster) %>%
  count(name = 'spR', Taxa) %>%
  mutate(spR = sum(spR)) %>%
  ungroup() %>%
  select(-Taxa) %>%
  distinct() %>%
  mutate(V_cluster = as.factor(V_cluster))
head(richness)


### anova models ###

model1 <- aov(data = richness, spR ~ V_cluster)
anova(model1)
summary(model1)


#If I  wanted to do Tukey post-hoc tests
TukeyHSD(model1)
# 3 is distinct from 1 & 2
# 1 and 2 are not distinct from each other

# associate significant difference indicator with cluster ID
richness <- richness %>%
  mutate(sigdif = ifelse(V_cluster == 3, "a", "b"))


### make a graph of those data ###
# gather summary data using emmeans
graphdata1<-as.data.frame(emmeans(model1, ~ V_cluster))
graphdata1<-graphdata1 %>%
  mutate(sigdif = ifelse(V_cluster == 3, "a", "b")) # as above, associate significant difference indicator with cluster ID
graphdata1

# graph
ggplot(graphdata1,
       aes(x=V_cluster,
           y=emmean,
           fill = V_cluster)) +
           #group=factorV_cluster)) +
  geom_bar(stat="identity",
           position="dodge",
           size=0.6) + #determines the bar width
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE),
                stat="identity",
                position=position_dodge(width=0.9),
                width=0.1) + #adds error bars
  geom_text(aes(label = sigdif),
                vjust = c(-4, -3, -2)) +
  labs(x="Varari Cluster",
       y="Species Richness") + #labels the x and y axes
  scale_fill_manual(values=c("red2", "dodgerblue2", "green4")) #fill colors for the bars





#### CLUSTERS FROM TURBINARIA CLUSTERS ####

### Clean data ###
spclusterT <- spcomp %>%
  filter(Location == "Varari") %>%
  select(-c(PhotoNum, Notes)) %>%
  left_join(turbcluster)
head(spclusterT)

### calculate species richness ###
richnessT <- spclusterT %>%
  rename(V_cluster = V_cluster_turb) %>%
  select(Location, CowTagID, Taxa, V_cluster) %>%
  group_by(Location, CowTagID, V_cluster) %>%
  count(name = 'spR', Taxa) %>%
  mutate(spR = sum(spR)) %>%
  ungroup() %>%
  select(-Taxa) %>%
  distinct() %>%
  mutate(V_cluster = as.factor(V_cluster))
head(richnessT)


### anova models ###

model2 <- aov(data = richnessT, spR ~ V_cluster)
anova(model2)
summary(model2)


#If I  wanted to do Tukey post-hoc tests
TukeyHSD(model2)
# 3 is distinct from 1
# 1 and 2 are not distinct from each other
# 2 and 3 are not distinct from each other

### associate significant difference indicator with cluster ID
richnessT <- richnessT %>%
  mutate(sigdif = ifelse(V_cluster == 3, "a", "b")) %>% # 3 gets a, all else gets b
  mutate(sigdif = ifelse(V_cluster == 2, "ab", sigdif)) # replace b with ab for 2, all else maintain a or b



### make a graph of those data ###
# gather summary data using emmeans
graphdata2<-as.data.frame(emmeans(model2, ~ V_cluster))
graphdata2<-graphdata2 %>%
  mutate(sigdif = ifelse(V_cluster == 3, "a", "b")) %>% # 3 gets a, all else gets b
  mutate(sigdif = ifelse(V_cluster == 2, "ab", sigdif)) # replace b with ab for 2, all else maintain a or b
graphdata2


# graph
ggplot(graphdata2,
       aes(x=V_cluster,
           y=emmean,
           fill = V_cluster)) +
  geom_bar(stat="identity",
           position="dodge",
           size=0.6) + #determines the bar width
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_errorbar(aes(ymax=emmean+SE, ymin=emmean-SE),
                stat="identity",
                position=position_dodge(width=0.9),
                width=0.1) + #adds error bars
  geom_text(aes(label = sigdif),
            vjust = c(-2, -3, -5)) +
  labs(x="Varari Cluster",
       y="Species Richness") + #labels the x and y axes
  scale_fill_manual(values=c("red2", "dodgerblue2", "green4")) #fill colors for the bars



