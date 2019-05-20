setwd("~/Desktop/QSRBIII (G)/QSRBIII_project/Particle_Detection")
getwd()

# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()

# BiocManager::install("EBImage")
# EBImage will be used for loading the image into an editable form in R.
library("EBImage")

# BiocManager::install("MaxContrastProjection")
# MaxContrastProjection has functions that will help stack the tif file.
library("MaxContrastProjection")


# BiocManager::install("FISHalyseR")
# Hopefully FISHalyseR will be able to do our particle detection.
library(FISHalyseR)

# install.packages("devtools")
# devtools::install_github("federicomarini/flowcatchR")
# library(flowcatchR)


# input file names here
fileCy = "w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_21 copy.tif"
fileDapi = "w1_HelaKyoto_Gapdh_2597_p01_dapi__Cell_CP_21 copy.tif"

## Load a sample image from EBImage.
# f = system.file("images", file, package="EBImage")

## read a cyan-stain image "fileCy" into variable name cyImg
## read a dapi-stain image "fileDapi" into variable dapiImg
cyImg = readImage(fileCy)
class(cyImg)

frames = numberOfFrames(cyImg)

# display(cyImg, method="raster", all=TRUE)
cyImgTrunc <- cyImg[,,(frames/2 -7):(frames/2 +7)]
display(cyImgTrunc, method = "raster", all=TRUE)
hist(cyImgTrunc)


dapiImg = readImage(fileDapi)
class(dapiImg)
frames2 = numberOfFrames(dapiImg)
dapiImgTrunc <- dapiImg[,,(frames2/2 -7):(frames2/2 +7)]
display(dapiImgTrunc, method="raster", all=TRUE)

# max_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "max")
# display(normalize(max_intensity_proj2))
# min_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "min")
# display(normalize(min_intensity_proj2))
# mean_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "mean")
# display(normalize(mean_intensity_proj2))
# median_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "median")
# display(normalize(median_intensity_proj2))
# sd_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "sd")
# display(normalize(sd_intensity_proj2))
# sum_intensity_proj2 = intensityProjection(
#   imageStack = dapiImgTrunc, projType = "sum")
# display(normalize(sum_intensity_proj2))

## Testing Otsu threshold as described in FISHalyseR

#img = normalize(max_intensity_proj2)
#t = calculateThreshold(img)
#img[img<t] <- 0
#img[img>=t] <- 1
#display(img)
# rm(img)
##

### Thresholding ###
# try to mask the nucleus and then do analysis

## Otsu Global Threshold
# threshold = otsu(normalize(dapiImgTrunc))
# threshold
# otsuThresh = combine( mapply(function(frame, th) frame > th, 
#                              getFrames(normalize(dapiImgTrunc)), 
#                              threshold, SIMPLIFY=FALSE) )
# display(otsuThresh, method = "raster", all=TRUE)
##

## Adaptive Threshold
# I think this uses a Laplacian filter
# disc = makeBrush(31, "disc")
# disc = disc / sum(disc)
# offset = 0.05
# background = filter2( normalize(dapiImgTrunc), disc )
# adapThresh = normalize(dapiImgTrunc) > background + offset
# display(adapThresh, method = ("raster"),all=TRUE)
##

## Adaptive Threshold
###### USE THIS ONE!
# The other methods give Image objects with storage mode Boolean

# This uses a linear filter "thresh"
adapThresh2 <- thresh(normalize(dapiImgTrunc), w=15, h=15, offset=0.05)
display(adapThresh2, method = "raster", all=TRUE )
max_intensity_proj = intensityProjection(imageStack = adapThresh2, projType = "max")
display(normalize(max_intensity_proj))

nucShadow <- fillHull(normalize(max_intensity_proj))
display(nucShadow)
##
print("stop")





# max_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "max")
# display(normalize(max_intensity_proj))
# min_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "min")
# display(normalize(min_intensity_proj))
# mean_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "mean")
# display(normalize(mean_intensity_proj))
# median_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "median")
# display(normalize(median_intensity_proj))
# sd_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "sd")
# display(normalize(sd_intensity_proj))
# sum_intensity_proj = intensityProjection(imageStack = cyImgTrunc, projType = "sum")
# display(normalize(sum_intensity_proj))
# 
