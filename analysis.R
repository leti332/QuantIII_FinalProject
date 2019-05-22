## Analysis
# DJ Nilson 5/21/2019
# data file headers: imgloc.cy3,Nucleus_puncta,Cyto_puncta,sum(nucShadow),sum(cytShadow)

# Libraries

library(data.table)
library(ggplot2)
library(gridExtra)

# Input data

cells <- fread('/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data.csv',
              header=FALSE,
              sep=",")
colnames(cells) <- c("Name","NucleusP","CytoP","NucleusA","CytoA")

cells.data <- cells[,2:5]

boxplot(cells.data[,1:2])
boxplot(log10(cells.data[,3:4]))
summary(cells.data)
