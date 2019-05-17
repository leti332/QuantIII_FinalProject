setwd("~/Desktop/QSRBIII (G)/QSRBIII_project/Testing EIB")
getwd()

# BiocManager::install("EBImage")
# BiocManager::install("MaxContrastProjection")
library("EBImage")

# input your file names here
fileCy = "w1_HelaKyoto_Gapdh_2597_p09_cy3__Cell_CP_9 copy.tif"
fileDapi = "w1_HelaKyoto_Gapdh_2597_p09_dapi__Cell_CP_9 copy.tif"

## Load a sample image from EBImage.
# f = system.file("images", file, package="EBImage")

## read a cyan-stain image "fileCy" into variable name cyImg
## read a dapi-stain image "fileDapi" into variable dapiImg
cyImg = readImage(fileCy)
display(cyImg, method="raster", all=TRUE)

dapiImg = readImage(fileDapi)
display(dapiImg, method="raster", all=TRUE)

## The next two lines of code normalize all the images in img ##
# imgNorm <- normalize(img)
# display(imgNorm, method = "raster", all = TRUE)



### A problem when recording 3D fluorescent microscopy images is 
### how to properly present these results in 2D. Maximum intensity 
### projections are a popular method to determine the focal plane 
### of each pixel in the image. The problem with this approach, 
### however, is that out-of-focus elements will still be visible, 
### making edges and fine structures difficult to detect. 

### This package aims to resolve this problem by using the contrast 
### around a given pixel to determine the focal plane, 
### allowing for a much cleaner structure detection than 
### would be otherwise possible. For convenience, this package also 
### contains functions to perform various other types of projections, 
### including a maximum intensity projection.

library("MaxContrastProjection")

## Generating one image from a max intensity projection
## Then, normalizing intensity values so image can be visualized
## Let 1 be cyan channel and 2 be dapi channel

### Should we consider using the function contrastProjection? ###
maxInt1 <- intensityProjection(imageStack = cyImg, projType = "max")
maxIntNorm1 <- normalize(maxInt1)
display(maxIntNorm1)

maxInt2 <- intensityProjection(imageStack = dapiImg, projType = "max")
maxIntNorm2 <- normalize(maxInt2)
display(maxIntNorm2)


## giving normalized images blue and green color
# cyImgBlue = channel(maxIntNorm1, 'asblue')
# dapiImgGr = channel(maxIntNorm2, 'asgreen')

## Using one function to generate bicolored image
stainCells = rgbImage(green = maxIntNorm1, blue = maxIntNorm2)
display(stainCells)
