#ASSESSMENT-1
setwd("C:/Users/User/OneDrive - The University of Manchester/GEOG71922/Assessment-1/spatialecology-assessment1")
#Install packages
install.packages(c("terra","sf","mapview"))
#Loading the libraries
library(terra)
library(sf)
#Read in the csv file 
meles <- read.csv("Melesmeles.csv")
#Check the first  six rows of the data 
head(meles)
class(meles)
