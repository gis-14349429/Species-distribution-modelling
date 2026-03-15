#ASSESSMENT-1
#Set working directory 
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview","spatstat","dismo","glmnet","maxnet","ggplot2","precrec","ranger"))
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
library(ggplot2)
library(precrec)
library(ranger)
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

#Finding optimum buffer size

#remove the NAs in the first row
glmRes=glmRes[!is.na(glmRes),]

#Subset the data frame to just the row containing the max log likelihood value

opt<-glmRes[which.max(glmRes$loglikelihood),]

#print
opt

#Plotting the loglikelyhood
plot(glmRes$radius, glmRes$loglikelihood, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "Buffer", ylab = "Log-Likelihood")

#Store optimum scale 
opt_scale_broadleaf<-opt$radius

#############################################################################

##### Model testing for Species distribution (GLM & Maxnet) ####

scot=st_read('scotSamp.shp')
#Loading in the land cover map 
LCM_sdm = rast("LCMUK.tif")
#Crop to the extent of the study area 
LCM_sdm = crop(LCM,st_buffer(scot, dist= 1000))
#Aggregate LCM raster
LCM_sdm = aggregate(LCM_sdm$LCMUK_1,fact=4,fun="modal")

meles_sdm <- vect(data.frame(x=meles$Longitude, y=meles$Latitude),geom=c("x","y"))
#Defining the project of the badger points
crs(meles_sdm)<-"epsg:4326"

meles.proj1=st_as_sf(meles_sdm,coords=c("x","y"),crs="epsg:4326")
#project badger data
meles.proj1 = st_transform(meles.proj1,crs(LCM))
#now crop the points to the study area
melesFin = meles.proj1[scot,]
#finally, mask the LCM to this boundary
LCM_sdm = crop(LCM_sdm,scot,mask=TRUE)

##### Set up covariates (environmental and land cover variable) ##### 
#Access levels of the raster
LCM_factor = as.factor(LCM_sdm$LCMUK_1)

#create an vector object called reclass
reclass_sdm = c(0,1,rep(0,20))

# combine with the LCM categories into a matrix of old and new values.
RCmatrix_sdm = cbind(levels(LCM_factor)[[1]],reclass_sdm)

RCmatrix_sdm = RCmatrix_sdm[,2:3]

#apply function to make sure new columns are numeric
RCmatrix_sdm = apply(RCmatrix_sdm,2,FUN=as.numeric)

#Assign new values to LCM with reclassification matrix
broadleaf_sdm = classify(LCM_sdm, RCmatrix_sdm)

####Creating spatial weights matrix#####
#Number of pixels needed to cover the optimum radius
nPix <- round(opt_scale_broadleaf/res(LCM_sdm)[1])
nPix <- (nPix*2)+1
#Build weights matrix
weightsMatrix = matrix(1:nPix^2,nrow=nPix,ncol=nPix)
#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)
focalCell=weightsMatrix[x,y]
indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)

#compute distances
distances=list()

for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM_sdm)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM_sdm)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}

weightsMatrix[] = unlist(distances)
weightsMatrix[weightsMatrix > opt_scale_broadleaf] <- NA
weightsMatrix[!is.na(weightsMatrix)] <- 1/length(weightsMatrix[!is.na(weightsMatrix)])

#Sum neighbourhood values from surrounding cells
lcm_wood=focal(broadleaf_sdm,w=weightsMatrix,fun="sum")

plot(lcm_wood)

## Creating a layer for urban class at 2300 meter ##
#create an vector object called reclassUrban
reclassUrban = c(rep(0,19),1,1)

# combine with the LCM categories into a matrix of old and new values.
RCmatrixUrban= cbind(levels(LCM_factor)[[1]],reclassUrban)

RCmatrixUrban=RCmatrixUrban[,2:3]

#apply function to make sure new columns are numeric
RCmatrixUrban=apply(RCmatrixUrban,2,FUN=as.numeric)
#Use the reclassify() function to asssign new values to LCM with our reclassification matrix

urban = classify(LCM_sdm, RCmatrixUrban)

#neighbourhood weights matrix to sum all available resources for each cell

#get number of picels needed to cover the 2300 metre radius for the urban class

nPixUrban=round(2300/res(LCM_sdm)[1])

#double this number (one for rows and one for columns) and add one

nPixUrban=(nPixUrban*2)+1

#buiild weights matrix
weightsMatrixUrban=matrix(1:nPixUrban^2,nrow=nPixUrban,ncol=nPixUrban)

#get focal cell 
x=ceiling(ncol(weightsMatrixUrban)/2)
y=ceiling(nrow(weightsMatrixUrban)/2)


focalCell=weightsMatrixUrban[x,y]


indFocal=which(weightsMatrixUrban==focalCell,arr.ind = TRUE)

#compute distances
distancesUrban=list()

for(i in 1:nPixUrban^2){
  ind.i=which(weightsMatrixUrban==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM_sdm)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM_sdm)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distancesUrban[[i]]=dist.i
  
}


#add distance values to the weights matrix
weightsMatrixUrban[]=unlist(distancesUrban)

#set cells outside search radius to NA
weightsMatrixUrban[weightsMatrixUrban>2300]=NA

#normalise the weights matrix by dividing all cell values by the number of cells. 
weightsMatrixUrban[!is.na(weightsMatrixUrban)]=1/length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

#sum urban class from all surrounding cells

lcm_urban_2300=focal(urban,w=weightsMatrixUrban,fun="sum")

#Loading a Digital elevation model
dem = rast('demScotland.tif')
#Resampling the DEM to match its resolution, origin and extent to another raster
dem = terra::resample(dem,lcm_wood)
#Stack the covariate layers together
allEnv = c(lcm_wood,lcm_urban_2300,dem)
names(allEnv) = c("broadleaf","urban","elev")

#Generating background points
set.seed(11)
#sample background - one point for every cell (9775)
back_sdm = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back_sdm =back_sdm[!is.na(back_sdm$broadleaf),]

back_sdm = st_as_sf(back_sdm,crs="EPSG:27700")

##Extraction of data to points##
#Get environmental covariates at presence locations
eP_sdm = terra::extract(allEnv, melesFin)

#bind together the presence data using cbind() which binds together objects by column (i.e. with different columns but the same number of rows)
Pres.cov=st_as_sf(cbind(eP_sdm, melesFin))

Pres.cov$Pres=1


#Remove the first column which is just an ID field.
Pres.cov=Pres.cov[,-1]

#get coordinates for spatial cross-validation later
coordsPres=st_coordinates(Pres.cov)

#drop geometry column using st_drop_geometry()
Back.cov=st_as_sf(data.frame(back_sdm,Pres=0))


#get coordinates of background points for cross validation later
coordsBack=st_coordinates(back_sdm)

#combine
coords=data.frame(rbind(coordsPres,coordsBack))

#assign column names
colnames(coords)=c("x","y")

#Combine all data

#combine presence and background
all.cov = rbind(Pres.cov,Back.cov)
#add coordinates
all.cov = cbind(all.cov,coords)
#remove any NAs
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)

