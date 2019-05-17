# Script is intended to take a TIFF of multiple imbedded images,
# a z-stack, and create a maximum intensity image.
# Attempts of increasing speed via parallelization has been unsuccessful
# https://support.bioconductor.org/p/80566/ - helpful thread on Image Normalization
# 5/15/2019
# Daniel Nilson

#library("tiff")
#library("imager")
library(tictoc)
library("EBImage")
library("FISHalyseR")




tic("Total")
tic("Loading Image")
## import 
# img = readImage(f)
imgloc.cy3 = "/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif"
imgloc.dapi = "/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/w1_HelaKyoto_Gapdh_2597_p01_dapi__Cell_CP_6.tif"

img.cy3 = readImage("/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif")

display(normalize(img.cy3), method="browser")
#img.cy3 <- tiff::readTIFF(imgloc.cy3, all = TRUE,native=FALSE) # Yields smaller ram usage
#img.dapi <- tiff::readTIFF(imgloc.dapi, all = TRUE,native=FALSE) # Yields smaller ram usage
toc()

tic("Plotting Image")
# Test plot
image(img.cy3[[32]],col=gray((0:32)/32),axes=FALSE)
image(img.dapi[[32]],col=gray((0:32)/32),axes=FALSE)
toc()

#Max Intensity
#Simple Approach
tic("Z-stack Simple")

MaxZstack.cy3 <- apply(simplify2array(img.cy3) ,c(1,2),max)
MaxZstack.dapi <- apply(simplify2array(img.dapi) ,c(1,2),max)


toc()


tic("Plotting ZStack")
image(MaxZstack.cy3,col=gray((0:32)/32),axes=FALSE)
image(MaxZstack.dapi,col=gray((0:32)/32),axes=FALSE)
toc()
toc()



## Practicing Particle Detection
# Max Entropy Threshold
img<-MaxZstack.cy3
t = calculateMaxEntropy(img)
img[img<t] <- 0
img[img>=t] <- 1

image(img,col=gray((0:32)/32),axes=FALSE)

# Otsu Threshold
img<-MaxZstack.cy3
t = calculateThreshold(img)
img[img<t] <- 0
img[img>=t] <- 1
image(img,col=gray((0:32)/32),axes=FALSE)