#ASSESSMENT-1
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview","spatstat"))
#Loading the libraries
library(terra)
library(sf)
library(spatstat)
#Read in the .csv file 
meles<- read.csv("Melesmeles.csv")
#Check the first  six rows of the data 
head(meles)
#Checking the class of the loaded .csv file
class(meles)
#Mapping the data
meles<-meles[!is.na(meles$Latitude),]

#creating crs object
meles.latlong<-data.frame(x=meles$Longitude,y=meles$Latitude)

#Creating spatial point objects from coordinates

meles.sp<-vect(meles.latlong,geom=c("x","y"))

#Setting the crs. 
crs(meles.sp)<-"epsg:4326"

#Plot the points on the map
plot(meles.sp)
