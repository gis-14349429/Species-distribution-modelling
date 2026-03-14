#ASSESSMENT-1
#Set working directory 
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview","spatstat","dismo","glmnet","maxnet"))
install.packages("mlr")
install.packages("rnaturalearth")
#Loading the libraries
library(terra)
library(sf)
library(spatstat)
library(mapview)
library(dismo)
library(rnaturalearth)
library(glmnet)
library(maxnet)
library(mlr)
#Read in the .csv file 
meles<- read.csv("Melesmeles.csv")
#Check the first  six rows of the data 
head(meles)
#Checking the class of the loaded .csv file
class(meles)
#Selecting all rows with complete values for latitude
meles<-meles[!is.na(meles$Latitude),]

#remove all points with uncertainty > 1000m
meles <- meles[meles$Coordinate.uncertainty_m<=1000, ]

#creating CRS object
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

#Setting the extent in the order: min x, max x, min y, max y
studyExtent<-c(-4.2, -2.7, 56.5, 57.5)
meles.crop<-crop(meles.sp, studyExtent)

#Projecting badger points to same CRS as LCM (EPSG:27700)
meles.proj<-project(meles.crop, crs(LCM))

#Setting crop extent with 5 km buffer
melesCoords<-crds(meles.proj)
x.min<-min(melesCoords[,1]) - 5000
x.max<-max(melesCoords[,1]) + 5000
y.min<-min(melesCoords[,2]) - 5000
y.max<-max(melesCoords[,2]) + 5000
extent.new<-ext(x.min,x.max,y.min,y.max)
LCM.crop <- crop(LCM$LCMUK_1, extent.new)
plot(LCM.crop)

#Reclassification of LCM to binary braodleaf woodland layer 
LCM.factor<-as.factor(LCM.crop)
levels(LCM.factor)
#create an vector object called reclass
reclass <- c(0,1,rep(0,19))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix<- cbind(levels(LCM.factor)[[1]],reclass)

RCmatrix<-RCmatrix[,2:3]

#apply function to make sure new columns are numeric 
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
broadleaf<-classify(LCM.crop, RCmatrix)

#Generating pseudo-absence points
set.seed(11)
back.xy<-spatSample(LCM.crop, size=1000, as.points=TRUE)
plot(LCM.crop)
plot(meles.proj,add=T)
plot(back.xy,add=TRUE,col='red')

#Creating a single dataframe to hold both presence and absence data
Abs<-data.frame(crds(back.xy),Pres=0)
head(Abs)
Pres<-data.frame(crds(meles.proj),Pres=1)
head(Pres)
melesData<-rbind(Pres,Abs)
head(melesData)
#Add coordinates to the dataframe using sf
melesSF=st_as_sf(melesData,coords=c("x","y"),crs="EPSG:27700")

#Creating a function to automate buffer calculations
landBuffer<-function(speciesData, r){
  #Buffer each point
  melesBuffer<-st_buffer(speciesData, dist=r)
  #Crop the woodland layer to the buffer extent
  bufferlandcover<-crop(broadleaf, melesBuffer)
  #Extracting the raster values
  masklandcover<-extract(bufferlandcover,melesBuffer,fun="sum")
  #Get woodland area
  landcoverArea<-masklandcover$LCMUK_1*625
  #Convert to percentage cover but convert to numeric object
  percentcover<-landcoverArea/as.numeric(st_area(melesBuffer))*100
  #Return the result
  return(percentcover)
}
#Loop over 100-2000 meter in 100 m steps
radii<-seq(100,2000,by=100)
resList=list()

for(i in radii){
  res.i=landBuffer(speciesData=melesSF,r=i)
  res.i
  resList[[i/100]]=res.i
  print(i)
  
}
#Combining results into data frame()
resFin=do.call("cbind",resList)

#convert to data frame
glmData=data.frame(resFin)
colnames(glmData)=paste("radius",radii,sep="")
glmData$Pres<-melesData$Pres

#make an empty data frame to then add a sequence of results from GLM
glmRes=data.frame(radius=NA,loglikelihood=NA)

#for loop to iterate over radius values and run a general linear model with glm()
for(i in radii){
  
  #build the model formula
  n.i=paste0("Pres~","radius",i,sep ="")
  
  glm.i=glm(formula(n.i),family = "binomial",data = glmData)
  
  ll.i=as.numeric(logLik(glm.i))
  
  glmRes=rbind(glmRes,c(i,ll.i))
  
}

