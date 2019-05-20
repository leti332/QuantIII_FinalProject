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
display(normalize(cyImgTrunc), method = "raster", all=TRUE)
hist(cyImgTrunc)


dapiImg = readImage(fileDapi)
frames2 = numberOfFrames(dapiImg)
dapiImgTrunc <- dapiImg[,,(frames/2 -7):(frames/2 +7)]

## Adaptive Threshold for nucleus ##
###### USE THIS ONE!
# The other methods give Image objects with storage mode Boolean
# I guess I could have used the other methods.

# This uses a linear filter "thresh"
adapThresh2 <- thresh(normalize(dapiImgTrunc), w=15, h=15, offset=0.05)
display(adapThresh2, method = "raster", all=TRUE )
max_intensity_proj_thresh2 = intensityProjection(imageStack = adapThresh2, projType = "max")
display(normalize(max_intensity_proj_thresh2))

nucShadow <- fillHull(normalize(max_intensity_proj_thresh2))

# Sanity check
display(nucShadow)
class(nucShadow)
imageData(nucShadow)


## Otsu Global Threshold for cytoplasm ##
threshold = otsu(normalize(cyImgTrunc))
threshold
otsuThresh = combine( mapply(function(frame, th) frame > th,
                             getFrames(normalize(cyImgTrunc)),
                             threshold, SIMPLIFY=FALSE) )
display(otsuThresh, method = "raster", all=TRUE)
##########
nmask = watershed( distmap(otsuThresh), 3 )     # Creating a cytoplasm outline
display(bwlabel(nmask), all=TRUE)

cytShadow <- fillHull(nmask)                    # Filling in holes in the cytoplasm
display(cytShadow)
##########################


# testFilter <- getFrame(cyImgTrunc, 1)
# display(normalize(testFilter))
# 
# testFilter[nucShadow==1] <-0
# display(normalize(testFilter))

## Masking the nucleus in cyImgTrunc
cyImgFilter1 <- cyImgTrunc
cyImgFilter1[nucShadow == 1] <- 0
display(normalize(cyImgFilter1), method = "raster", all=TRUE)

## Masking the cytoplasm in cyImgTrunc
cyImgFilter2 <- cyImgFilter1
cyImgFilter2[cytShadow == 0] <- 0
display(normalize(cyImgFilter2), method = "raster", all=TRUE)

##################################################
## A cytoplasm outline has been made AHAHAHAHA  ##
##################################################

max_contrast = contrastProjection(imageStack = cyImgFilter2, w_x = 1, w_y = 1,
                                  smoothing = 0, brushShape = "box")
display(normalize(max_contrast), method = "raster")
text(x = 20, y = 20, label = "Contrast", adj = c(0,1), col = "blue", cex = 2)

max_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "max")
display(normalize(max_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Max", adj = c(0,1), col = "blue", cex = 2)

min_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "min")
display(normalize(min_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Min", adj = c(0,1), col = "blue", cex = 2)

mean_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "mean")
display(normalize(mean_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Mean", adj = c(0,1), col = "blue", cex = 2)

median_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "median")
display(normalize(median_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Med", adj = c(0,1), col = "blue", cex = 2)

sd_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "sd")
display(normalize(sd_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Sd", adj = c(0,1), col = "blue", cex = 2)

sum_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "sum")
display(normalize(sum_intensity_proj), method = "raster")
text(x = 20, y = 20, label = "Sum", adj = c(0,1), col = "blue", cex = 2)


a <- normalize(max_contrast^10 - min_intensity_proj^10)
display(a)


img = (a)
t = calculateThreshold(img)
img[img<t] <- 0
img[img>=t] <- 1
display(img)


##
print("stop")

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

## Otsu Global Threshold ##
# threshold = otsu(normalize(dapiImgTrunc))
# threshold
# otsuThresh = combine( mapply(function(frame, th) frame > th, 
#                              getFrames(normalize(dapiImgTrunc)), 
#                              threshold, SIMPLIFY=FALSE) )
# display(otsuThresh, method = "raster", all=TRUE)
###########################

## Adaptive Threshold ##
# I think this uses a Laplacian filter
# disc = makeBrush(31, "disc")
# disc = disc / sum(disc)
# offset = 0.05
# background = filter2( normalize(dapiImgTrunc), disc )
# adapThresh = normalize(dapiImgTrunc) > background + offset
# display(adapThresh, method = ("raster"),all=TRUE)
###########################






