#ASSESSMENT-1
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview","spatstat"))
#Loading the libraries
library(terra)
library(sf)
library(spatstat)
library(mapview)
#Read in the .csv file 
meles<- read.csv("Melesmeles.csv")
#Check the first  six rows of the data 
head(meles)
#Checking the class of the loaded .csv file
class(meles)
#Mapping the data
meles<-meles[!is.na(meles$Latitude),]
meles<-meles[!is.na(meles$Longitude),]

#remove all points with uncertainty > 100m
meles[meles$Coordinate.uncertainty_m<=1000,]

#creating crs object
meles.latlong<-data.frame(x=meles$Longitude,y=meles$Latitude)

#Using coordinates to create spatial object points
meles.sp<-vect(meles.latlong,geom=c("x","y"))

#check that the points now have our desired crs. 
crs(meles.sp)<-"epsg:4326"

#Plot the points on the map
plot(meles.sp)
mapview(st_as_sf(meles.sp))

#Read in the raster data
LCM=rast("LCMUK.tif")

#First set the extent to something workable
studyExtent<-c(-4.2,-2.7,56.5,57.5) #list coordinates in the order: min x, max x, min y, max y

#now crop our points to this area
C<-crop(meles.sp,studyExtent)

melesFin<-project(C,crs(LCM))
melesCoords<-crds(melesFin)


x.min <- min(melesCoords[,1]) - 5000
x.max <- max(melesCoords[,1]) + 5000
y.min <- min(melesCoords[,2]) - 5000
y.max <- max(melesCoords[,2]) + 5000

extent.new <- ext(x.min, x.max, y.min, y.max)


LCM <- crop(LCM$LCMUK_1, extent.new)

plot(LCM)

plot(melesFin,add=TRUE, col='blue')
