# Merged Code and Methods from EBI_Load_MaxIntensity_Z & Particle\ Detect\ V2
# 

library("EBImage")
library(FISHalyseR)
library("MaxContrastProjection")
#library("spatialfil")



# input file names here
setwd("/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/")

imgloc.cy3 = "w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif"
imgloc.dapi = "w1_HelaKyoto_Gapdh_2597_p01_dapi__Cell_CP_6.tif"
x<-40
y<-50

## Load a sample image from EBImage.
# f = system.file("images", file, package="EBImage")

## read a cyan-stain image "fileCy" into variable name cyImg
## read a dapi-stain image "fileDapi" into variable dapiImg
cyImg = readImage(imgloc.cy3)

frames = numberOfFrames(cyImg)

# display(cyImg, method="raster", all=TRUE)
#cyImgTrunc <- cyImg[,,(frames/2 -7):(frames/2 +7)]
cyImgTrunc <- cyImg[,,x:y]
#display(normalize(cyImgTrunc), method = "raster", all=TRUE)
hist(cyImgTrunc)


dapiImg = readImage(imgloc.dapi)
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
imageData(nucShadow)


## Otsu Global Threshold for cytoplasm ##
threshold = otsu(normalize(cyImgTrunc))
#threshold
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


## Making the nucleus and cytoplasm masks from cyImgTrunc
cyImgFilter1 <- cyImgTrunc
NucleusFilter <- cyImgTrunc
NucleusFilter[nucShadow == 0] <- 0
cyImgFilter1[nucShadow == 1] <- 0
#display(normalize(cyImgFilter1), method = "raster", all=TRUE)
#display(normalize(NucleusFilter), method = "raster", all=TRUE)

## Masking the Extracellular in cyImgTrunc
cyImgFilter2 <- cyImgFilter1
cyImgFilter2[cytShadow == 0] <- 0
CytoplasmImage<-cyImgFilter2
#display(normalize(cyImgFilter2), method = "raster", all=TRUE)

##################################################
## A cytoplasm outline has been made AHAHAHAHA  ##
##################################################

## Cytoplasm

max_contrast = contrastProjection(imageStack = cyImgFilter2, w_x = 2, w_y = 2,
                                  smoothing = 0, brushShape = "box")

min_intensity_proj = intensityProjection(imageStack = cyImgFilter2, projType = "min")



a <- normalize(max_contrast^10 - min_intensity_proj^10)
img = (a)
t = calculateThreshold(normalize(img))
img[img<t] <- 0
img[img>=t] <- 1
display(img)
img_2<-analyseParticles(img, 20, 2,0) # Very simple "Clean up"
display(img_2)




## Nucleus

max_contrast_N = contrastProjection(imageStack = NucleusFilter, w_x = 2, w_y = 2,
                                  smoothing = 0, brushShape = "box")

min_intensity_proj_N = intensityProjection(imageStack = NucleusFilter, projType = "min")
a_N <- normalize(max_contrast_N^8 - min_intensity_proj_N^8)
img_N = (a_N)
t_N = calculateThreshold(normalize(img_N))
img_N[img_N<t] <- 0
img_N[img_N>=t] <- 1
display(img_N)

img_N_2<-analyseParticles(img_N, 20, 2,0) # Very simple "Clean up"
display(img_N_2)

# Counting

img_N_2.counted <- bwlabel(img_N_2) # Labels connected (connected sets) objects in a binary image. 
img_2.counted <- bwlabel(img_2)

Nucleus_puncta <- max(img_N_2.counted)
Cyto_puncta <- max(img_2.counted)

