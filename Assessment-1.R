#ASSESSMENT-1
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview","spatstat"))
install.packages("mlr")
install.packages("rnaturalearth")
#Loading the libraries
library(terra)
library(sf)
library(spatstat)
library(mapview)
library(dismo)
library(rnaturalearth)
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

