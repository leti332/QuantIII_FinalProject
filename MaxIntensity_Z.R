# Script is intended to take a TIFF of multiple imbedded images,
# a z-stack, and create a maximum intensity image.
# Attempts of increasing speed via parallelization has been unsuccessful
# 5/13/2019
# Daniel Nilson

library("tiff")
library("imager")
library(tictoc)


## If "tiff" package does not compile, additional C libraries are necessary
# ergo: sudo apt-get install libfftw3-dev libtiff5-dev
# Package names proably vary with OS
#Tiff: https://cran.r-project.org/web/packages/ijtiff/vignettes/the-imagej-problem.html


tic("Total")
tic("Loading Image")
# import 
imgloc = "/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/GAPDH/images/w1_HelaKyoto_Gapdh_2597_p03_cy5.TIF"
img2 <- tiff::readTIFF(imgloc, all = TRUE,native=TRUE) # Yields smaller ram usage
toc()

tic("Plotting Image")
# Test plot
image(img2[[32]],col=gray((0:32)/32),axes=FALSE)
toc()

#Max Intensity
#Simple Approach
tic("Z-stack Simple")

MaxZstack <- apply(simplify2array(img2) ,c(1,2),max)
toc()


tic("Plotting ZStack")
image(MaxZstack,col=gray((0:32)/32),axes=FALSE)
toc()
toc()







