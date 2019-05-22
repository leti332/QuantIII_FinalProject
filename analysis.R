## Analysis
# DJ Nilson 5/21/2019
# data file headers: imgloc.cy3,Nucleus_puncta,Cyto_puncta,sum(nucShadow),sum(cytShadow)

# Libraries

library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(ggthemes)

# Input data

cells <- fread('/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data.csv',
              header=FALSE,
              sep=",")
colnames(cells) <- c("Name","NucleusP","CytoP","NucleusA","CytoA")

cells.data <- cells[,2:5]

boxplot(cells.data[,1:2])
boxplot(log10(cells.data[,3:4]))

# ggplot boxplots
Gcells.data.puncta  <-gather(cells.data[,1:2],Puncta,Counts)

p <- ggplot(Gcells.data.puncta,aes(Puncta,Counts)) + labs(title="Puncta relative to POI")
p.fin <- p + geom_boxplot() + scale_y_continuous(trans='log10') + theme_gray()
ggsave(p.fin, filename = "Puncta_v_POI.png", width = 5, height=5)

Gcells.data.area  <-gather(cells.data[,3:4],POI,Area)
q <- ggplot(Gcells.data.area,aes(POI,Area)) + labs(title="Area Comparison")
q.fin <- q+ geom_boxplot() + scale_y_continuous(trans='log10') + theme_gray()
ggsave(q.fin, filename = "Area_v_POI.png", width = 5, height=5)

# Removing Errant Data

cells.data<-subset(cells.data,NucleusP>10 & cells.data$NucleusA>10 & cells.data$CytoP>10 & cells.data$NucleusA>10)

Gcells.data.puncta  <-gather(cells.data[,1:2],Puncta,Counts)

b <- ggplot(Gcells.data.puncta,aes(Puncta,Counts)) + labs(title="Puncta relative to POI",subtitle="Noise removed")
b.fin <- b + geom_boxplot() + scale_y_continuous(trans='log10') + theme_gray()
ggsave(b.fin, filename = "Puncta_v_POI_RemovedNoise.png", width = 5, height=5)

Gcells.data.area  <-gather(cells.data[,3:4],POI,Area)
c <- ggplot(Gcells.data.area,aes(POI,Area)) + labs(title="Area Comparison",subtitle="Noise removed")
c.fin <- c + geom_boxplot() + scale_y_continuous(trans='log10') + theme_gray()
ggsave(c.fin, filename = "Area_v_POI_RemovedNoise.png", width = 5, height=5)

# Normalizing Puncta with Respect to Area
Normcells.data <-cells.data[,1:2]


for (i in dim(cells.data)[1]){
  Normcells.data$NucleusP[i]<-cells.data$NucleusP[i]/cells.data$NucleusA[i]
  Normcells.data$CytoP[i]<-cells.data$CytoP[i]/cells.data$CytoA[i]
}

GNormcells.data.puncta  <-gather(Normcells.data[,1:2],Puncta,Counts) 

z <- ggplot(GNormcells.data.puncta,aes(Puncta,Counts)) + labs(title="Puncta normalized to POI Area")
z.fin <- z + geom_boxplot() + theme_gray()
ggsave(z.fin, filename = "NormalizedPuncta_v_POI_RemovedNoise.png", width = 5, height=5)

# Scatter Plot of Puncta v. Area separated by POI
Gcells.data.puncta  <-gather(cells.data[,1:2],Puncta,Counts)
Gcells.data.area  <-gather(cells.data[,3:4],POI,Area)
Gcells.data.scatter <- cbind(Gcells.data.puncta,Gcells.data.area$Area)
colnames(Gcells.data.scatter)[3] <- "Area"

m <- ggplot(Gcells.data.scatter, aes(x=log10(Area), y=Counts,color=Puncta)) +
  geom_point(size=2, shape=23) +labs(title="Puncta vs. Area",subtitle="Noise removed") + theme_gray()
m
ggsave(m, filename = "Scatter_Counts_v_Area.png", width = 5, height=5)

