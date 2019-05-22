#!/usr/bin/env Rscript
#!/usr/bin/Rscript
options(warn=-1)

# Merged Code and Methods from EBI_Load_MaxIntensity_Z & Particle\ Detect\ V2
# Image is loaded. Then the image is cut apart from two masks: a nucleus mask and cytoplasm mask
# Cytoplasm mask is made from the DAPI image
# Nucleus mask is simply the inverse of the Cytoplasm mask (extracellular removed)
#suppressMessages(library(base64))
suppressMessages(library("EBImage"))
suppressMessages(library(FISHalyseR))
suppressMessages(library("MaxContrastProjection"))
#library("spatialfil")
suppressMessages(library("optparse"))

# Command Line

option_list = list(
  make_option(c("-c", "--cy"), type="character", default=NULL, 
              help="tiff cy3 file", metavar="character"),
  make_option(c("-d", "--dapi"), type="character", default=NULL, 
              help="tiff dapi file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$cy)){
  print_help(opt_parser)
  stop("Need corresponding Cy3 and DAPI TIFF files.", call.=FALSE)
}



# # input file names here
# setwd("/home/daniel/Documents/Code/R/Einstein/QuantIII_FinalProject/data_simulation/cropped_img/")
# 
# imgloc.cy3 = "w1_HelaKyoto_Gapdh_2597_p01_cy3__Cell_CP_6.tif"
# imgloc.dapi = "w1_HelaKyoto_Gapdh_2597_p01_dapi__Cell_CP_6.tif"
imgloc.cy3 = opt$cy
imgloc.dapi = opt$dapi
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
#hist(cyImgTrunc)


dapiImg = readImage(imgloc.dapi)
frames2 = numberOfFrames(dapiImg)
dapiImgTrunc <- dapiImg[,,(frames/2 -7):(frames/2 +7)]

######### Adaptive Threshold for nucleus ##
# USE THIS ONE!
# The other methods give Image objects with storage mode Boolean
# I guess I could have used the other methods.

threshold2 = otsu(normalize(dapiImgTrunc))
otsuThresh2 = combine( mapply(function(frame, th) frame > th,
                              getFrames(normalize(dapiImgTrunc)),
                              threshold2, SIMPLIFY=FALSE) )
#display(otsuThresh2, method = "raster", all=TRUE)

##########

nmask2 = watershed( distmap(otsuThresh2), 3 )   # Creating a nucleus outline

#display(colorLabels(nmask2), all=TRUE)

nucShadow <- intensityProjection(imageStack = fillHull(nmask2),projType = "max")  # Filling in holes in the mask and stacking 15 masks)


#display(nucShadow)



## Otsu Global Threshold for cytoplasm ##
threshold = otsu(normalize(cyImgTrunc))
#threshold
otsuThresh = combine( mapply(function(frame, th) frame > th,
                             getFrames(normalize(cyImgTrunc)),
                             threshold, SIMPLIFY=FALSE) )
#display(otsuThresh, method = "raster", all=TRUE)
##########
nmask = watershed( distmap(otsuThresh), 3 )     # Creating a cytoplasm outline
#display(bwlabel(nmask), all=TRUE)

cytShadow <- fillHull(nmask)                    # Filling in holes in the cytoplasm
#display(cytShadow)
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
## A cytoplasm outline has been made   ##
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
#display(img)
img_2<-analyseParticles(img, 20, 2,0) # Very simple "Clean up"
#display(img_2)




## Nucleus

max_contrast_N = contrastProjection(imageStack = NucleusFilter, w_x = 2, w_y = 2,
                                  smoothing = 0, brushShape = "box")

min_intensity_proj_N = intensityProjection(imageStack = NucleusFilter, projType = "min")
a_N <- normalize(max_contrast_N^8 - min_intensity_proj_N^8)
img_N = (a_N)
t_N = calculateThreshold(normalize(img_N))
img_N[img_N<t_N] <- 0
img_N[img_N>=t_N] <- 1
#display(img_N)

img_N_2<-analyseParticles(img_N, 20, 2,0) # Very simple "Clean up"
#display(img_N_2)

# Counting

img_N_2.counted <- bwlabel(img_N_2) # Labels connected (connected sets) objects in a binary image. 
img_2.counted <- bwlabel(img_2)

Nucleus_puncta <- max(img_N_2.counted)
Cyto_puncta <- max(img_2.counted)

##### Area Calculations
## Area counter of cytoplasm

#i = 1
#areaCountCyt = 0
#for (i in 1:length(cytShadow)){
#  if (cytShadow[i] == 1){
#    areaCountCyt = areaCountCyt + 1
#  }
 # i+1
  #print(areaCountCyt)
  #print (i)
#}

## Area counter of nucleus. areaCountNuc outputs number of pixels comprising nucleus

#i = 1
#areaCountNuc = 0
#for (i in 1:length(nucShadow)){
#  if (nucShadow[i] == 1){
#    areaCountNuc = areaCountNuc + 1
#  }
#  i+1
#  print(areaCountNuc)
#  print (i)
#}

#areaCountNuc


# Output to Terminal


output<- c(imgloc.cy3,Nucleus_puncta,Cyto_puncta,sum(nucShadow),sum(cytShadow))
cat(paste(shQuote(output, type="cmd"), collapse=", "),"\n")
#print("\n",quote=FALSE)

