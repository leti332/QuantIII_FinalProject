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
library("MaxContrastProjection")
library("FISHalyseR")
library("spatialfil")



##



tic("Total")

tic("Loading Image")



############

## import ##

############



# img = readImage(f)

setwd("/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/")

imgloc.cy3 = "w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif"

imgloc.dapi = "w1_HelaKyoto_Gapdh_2597_p01_dapi__Cell_CP_6.tif"



img.cy3 = readImage("w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif")

display(normalize(img.cy3), method="browser") # Normalizing enables proper function of image processing tools.

#img.cy3<-normalize(img.cy3) #this is bad

#img.cy3<-normalize(img.cy3) #this is bad



#img.cy3 <- tiff::readTIFF(imgloc.cy3, all = TRUE,native=FALSE) # Yields smaller ram usage

#img.dapi <- tiff::readTIFF(imgloc.dapi, all = TRUE,native=FALSE) # Yields smaller ram usage



##





toc()



#tic("Plotting Image")

# Test plot

#image(img.cy3[[32]],col=gray((0:32)/32),axes=FALSE)

#image(img.dapi[[32]],col=gray((0:32)/32),axes=FALSE)

#toc()



#################

#### Z-Stack ####

#################





# Max Intensity

# Simple Approach

tic("Z-stack Simple")

#MaxZstack.cy3 <- apply(simplify2array(img.cy3) ,c(1,2),max)
x=40
y=50

MaxZstack.cy3 <- apply(simplify2array(img.cy3[,,x:y]^2),c(1,2),max)
MaxZstack.dapi <- apply(simplify2array(img.dapi) ,c(1,2),max)


# Different Kinds of Projects
max_contrast_large = contrastProjection(imageStack = img.cy3[,,x:y],w_x = 15, w_y = 15, smoothing = 15,brushShape = "box")
max_contrast_small = contrastProjection(imageStack = img.cy3[,,x:y],w_x = 3, w_y = 3, smoothing = 2,brushShape = "box")
max_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "max")
min_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "min")
mean_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y],projType = "mean")
median_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y],projType = "median")
sd_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "sd")
sum_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "sum")

#MaxZstack.cy3 <- apply(simplify2array(img.cy3) ,c(1,2),max)

x=40

y=50



MaxZstack.cy3 <- apply(simplify2array(img.cy3[,,x:y]^2),c(1,2),max)
MaxZstack.dapi <- apply(simplify2array(img.dapi) ,c(1,2),max)





# Different Kinds of Projects

max_contrast_large = contrastProjection(imageStack = img.cy3[,,x:y],w_x = 15, w_y = 15, smoothing = 15,brushShape = "box")

max_contrast_small = contrastProjection(imageStack = img.cy3[,,x:y],w_x = 3, w_y = 3, smoothing = 2,brushShape = "box")

max_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "max")

min_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "min")

mean_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y],projType = "mean")

median_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y],projType = "median")

sd_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "sd")

sum_intensity_proj = intensityProjection(imageStack = img.cy3[,,x:y], projType = "sum")

sd_maxcontrast <-normalize(max_contrast_large+sd_intensity_proj)

mean_maxcontrast <-normalize(max_contrast_large-median_intensity_proj)

toc()

#Plotting the Projections



#Plotting the Projections

tic("Plotting ZStack")

image(normalize(MaxZstack.cy3),col=gray((0:32)/32),axes=FALSE) + title(main = "Simple-Max Intensity")

image(MaxZstack.dapi,col=gray((0:32)/32),axes=FALSE) + title(main = "Simple-Max Intensity")



image(normalize(max_contrast_large+max_contrast_small),col=gray((0:32)/32),axes=FALSE)  + title(main = "teste")

image(normalize(max_contrast_large),col=gray((0:32)/32),axes=FALSE)  + title(main = "max_contrast_large")

image(normalize(max_contrast_small),col=gray((0:32)/32),axes=FALSE)  + title(main = "max_contrast_small")

image(normalize(max_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "max_intensity_proj")

image(normalize(min_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "min_intensity_proj")

image(normalize(mean_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "mean_intensity_proj")

image(normalize(median_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "median_intensity_proj")

image(normalize(sd_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "sd_intensity_proj")

image(normalize(sum_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "sum_intensity_proj")



image(normalize(sd_maxcontrast),col=gray((0:32)/32),axes=FALSE)  + title(main = "sd_maxcontrast")

image(normalize(mean_maxcontrast),col=gray((0:32)/32),axes=FALSE)  + title(main = "mean_maxcontrast")



toc()

toc()


####
# Making the Image that the Particle Analyzer actually likes

sd_intensity_proj_LP = applyFilter(sd_intensity_proj, kernel = convKernel(sigma = 1.4, k = "LoG"))
image(sd_intensity_proj_LP,col=gray((0:32)/32),axes=FALSE)
StageTwo<-(max_contrast_small+1*sd_intensity_proj_LP)^1
image(StageTwo,col=gray((0:32)/32),axes=FALSE)
image(MaxZstack.cy3,col=gray((0:32)/32),axes=FALSE)


# Experimentation

# filter= "sobel"
# 
# sd_intensity_proj_LP = applyFilter(sd_intensity_proj, kernel = convKernel(sigma = 1.4, k =filter))
# image(sd_intensity_proj_LP^1,col=gray((0:32)/32),axes=FALSE)
# image(normalize(sd_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "sd_intensity_proj")
# 
# 
# max_contrast_small_LP = applyFilter(max_contrast_small^1, kernel = convKernel(sigma = 2, k =filter))
# image(max_contrast_small_LP*100,col=gray((0:32)/32),axes=FALSE)
# image(normalize(max_contrast_small),col=gray((0:32)/32),axes=FALSE)  + title(main = "max_contrast_small")

####




####

# Making the Image that the Particle Analyzer actually likes



sd_intensity_proj_LP = applyFilter(sd_intensity_proj, kernel = convKernel(sigma = 1.4, k = "LoG"))
image(sd_intensity_proj_LP,col=gray((0:32)/32),axes=FALSE)
StageTwo<-(max_contrast_small+1*sd_intensity_proj_LP)^1
image(StageTwo,col=gray((0:32)/32),axes=FALSE)
image(MaxZstack.cy3,col=gray((0:32)/32),axes=FALSE)





# Experimentation

# filter= "sobel"
# 
# sd_intensity_proj_LP = applyFilter(sd_intensity_proj, kernel = convKernel(sigma = 1.4, k =filter))
# image(sd_intensity_proj_LP^1,col=gray((0:32)/32),axes=FALSE)
# image(normalize(sd_intensity_proj),col=gray((0:32)/32),axes=FALSE)  + title(main = "sd_intensity_proj")
# 
# 
# max_contrast_small_LP = applyFilter(max_contrast_small^1, kernel = convKernel(sigma = 2, k =filter))
# image(max_contrast_small_LP*100,col=gray((0:32)/32),axes=FALSE)
# image(normalize(max_contrast_small),col=gray((0:32)/32),axes=FALSE)  + title(main = "max_contrast_small")



####





## Practicing Particle Detection

# Max Entropy Threshold

img<-normalize(StageTwo)


img<-normalize(StageTwo)


t = calculateMaxEntropy(img)

img[img<t] <- 0

img[img>=t] <- 1



image(img,col=gray((0:32)/32),axes=FALSE)

img2<-analyseParticles(img, 20, 1,0) # Very simple "Clean up"

par(mfrow=c(1,3))
image(MaxZstack.cy3,col=gray((0:32)/32),axes=FALSE)
image(StageTwo,col=gray((0:32)/32),axes=FALSE)
image(img2,col=gray((0:32)/32),axes=FALSE)



# Otsu Threshold
img<-normalize(max_contrast_small_LP)

img2<-analyseParticles(img, 20, 1,0) # Very simple "Clean up"



par(mfrow=c(1,3))

image(MaxZstack.cy3,col=gray((0:32)/32),axes=FALSE)

image(StageTwo,col=gray((0:32)/32),axes=FALSE)

image(img2,col=gray((0:32)/32),axes=FALSE)







# Otsu Threshold

img<-normalize(max_contrast_small_LP)

t = calculateThreshold(img)

img[img<t] <- 0

img[img>=t] <- 1

image(img,col=gray((0:32)/32),axes=FALSE)