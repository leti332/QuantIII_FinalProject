# Script is intended to take a TIFF of multiple imbedded images,
# a z-stack, and create a maximum intensity image.
# Attempts of increasing speed via parallelization has been unsuccessful
# 5/13/2019
# Daniel Nilson

library("tiff")
library("imager")
# library(plyr)
library(tictoc)
# library(parallel)
# library(future.apply)
# library(foreach)
# library(doParallel)



## If "tiff" package does not compile, additional C libraries are necessary
# ergo: sudo apt-get install libfftw3-dev libtiff5-dev
# Package names proably vary with OS
#Tiff: https://cran.r-project.org/web/packages/ijtiff/vignettes/the-imagej-problem.html

# # Creating 30 R instances over each CPU
# plan(multicore)
# registerDoParallel(makeCluster(10))

tic("Total")
tic("Loading Image")
# import 
imgloc = "/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/GAPDH/images/w1_HelaKyoto_Gapdh_2597_p03_cy5.TIF"
# img1 <- tiff::readTIFF(imgloc, all = TRUE,native=FALSE)
img2 <- tiff::readTIFF(imgloc, all = TRUE,native=TRUE) # Yields smaller ram usage
#str(img)  
str(img2) # Smaller in the memory!
toc()

tic("Plotting Image")
# Test plot
image(img2[[32]],col=gray((0:32)/32),axes=FALSE)
toc()

#Max Intensity
#Simple Approach
tic("Z-stack Simple")
#ListtoArray <- simplify2array(img2) # Makes an array with a Z/3rd dimension
MaxZstack <- apply(simplify2array(img2) ,c(1,2),max)
toc()

# #Max Intensity
# #Multi Core
# tic("Z-stack MC")
# ListtoArray <- simplify2array(img2) # Makes an array with a Z/3rd dimension
# MaxZstack <- future_apply(ListtoArray,c(1,2),
#                           max,future.scheduling=1.0,
#                           future.chunk.size=NULL)
# #MaxZstack <- future_apply(simplify2array(img2) ,c(1,2),max)
# toc()
# tic()
# MaxZStackMaybe <- foreach(n = 1:2048, m = 1:2048, .combine = rbind) %dopar% max(ListtoArray[n,m,])
# toc()




tic("Plotting ZStack")
image(MaxZstack,col=gray((0:32)/32),axes=FALSE)
toc()
toc()







