## Mini project Advanced GI Science for Earth and Environment
## Period 4, Academic Year 2015-2016
## Case 1 - Random Animal Trajectories
##M ark ten Vregelaar – Stijn Wijdeven – Joris Wever – Erwin van den Berg – Bob Houtkooper
## 29-02-2016

## Import modules
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(spacetime)
library(lattice)
library(maptools)
library(plyr)
library(raster)
## Source functions
source('Functions/create_random_traj.R')
source('Functions/maxspeed2009.R')
source('Functions/getbeardata.R')
source('Functions/LineDensityRaster.R')
source('Functions/EncounterProbabilityRaster.R')
## Create directories
data_dir = 'Data'
output_dir = 'Output'
functions_dir = 'Functions'

dir_list = list(data_dir, output_dir, functions_dir)
for (i in dir_list)
{if (!file.exists(i)){
  dir.create(i)
} 
}

## Load the data
load("Data/mating2009.Rdata")
DEM<- readGDAL("Data/DEM.tif") 
data2007 <- read.table("Data/August2007.txt", header = T, sep=",")
data2007$GMT_date <- NULL
data2007$LMT_date <- as.POSIXct(data2007$LMT_date, 
                                format='%d-%m-%Y %H:%M:%S')


Koski2007 <- subset(data2007, PubName == "Koski (2310)")
Koski2007 <- Koski2007[order(Koski2007$LMT_date),]      # order on data_time
tstartkoski <- min(Koski2007$LMT_date)
Koski2007<-getbeardata(data2007,Pubname = "Koski (2310)",tstart = tstartkoski )

mating2009$LMT_date <- as.POSIXct(mating2009$LMT_date, 
                                  format='%d-%m-%Y %H:%M:%S')

Krut2009 <- subset(mating2009, PubName == "Krut (2937)")
tstartkrut <- min(Krut2009$LMT_date)
Krut2009<-getbeardata(mating2009,Pubname = "Krut (2937)",tstart = tstartkrut )


mspeedKrut <- maxspeed2009('Krut') 

mypoints<-Krut2009
#mypoints<-subset(mypoints,tspannum>450 & tspannum<700)
mypoints<-subset(mypoints,tspannum>690 & tspannum<700)
## Compute random trajectories
V<-mspeedKrut*1000+100
Rtrajectories<-list()

for (i in seq(1:10)){
  
  Rtrajectory<-create_random_traj(mypoints,V,4)
  Rtrajectories<-c(Rtrajectories,Rtrajectory)
  
  
}


plot(mypoints,col='blue',xlim=c(bbox(mypoints)[1][1]-(V/5),bbox(mypoints)[,2][1]+(V/5)),axes=T,
     ylim=c(bbox(mypoints)[2][1]-(V/5),bbox(mypoints)[,2][2]+(V/5)),xlab='X (meters)', ylab= 'Y (meters)',
     main='Random trajectory between known points')
lapply(Rtrajectories,function(x) lines(x[1]@coords,add=T,col='red'))





















#plot(mypoints,add=T,col='blue')
## Create line density surface of trajectories 
system.time(LineDensity <- LineDensityRaster(Rtrajectories,500))


plot(LineDensity, col=colorRampPalette(c("white", "orangered", "black"))(101))
lapply(Rtrajectories,function(x) lines(x[1]@coords,add=T,col=rgb(0,1,0,alpha=0.1)))
plot(mypoints,add=T,pch=19,col='white')
lapply(Rtrajectories,function(x) lines(x[1]@coords,add=T,col=rgb(0,1,0,alpha=0.025)))






# system.time(Probability <- EncounterProbabilityRaster(Rtrajectories1, Rtrajectories2, 50))
# plot(Probability, col=colorRampPalette(c("lightyellow", "orangered", "black"))(101))
# 
# lapply(Rtrajectories_A,function(x) lines(x[1]@coords,add=T,col=rgb(0,1,0,alpha=0.1)))
# lapply(Rtrajectories_B,function(x) lines(x[1]@coords,add=T,col=rgb(1,0,0,alpha=0.1)))
# 
# plot(mypoints_A,add=T,pch=19,col='white')
# plot(mypoints_B,add=T,pch=19,col='yellow')





## PROBABILITY
load('data/all1equi.Rdata')
load('data/allgood1.Rdata')
load('data/all2equi.Rdata')
load('data/allgood2.Rdata')


# comparison 1: use n of 3 or 4
n1 <- 1

equi1_extent <- extent(all1_equi[[n1]])
good1_extent <- extent(allgood1[[n1]])
output_extent <- extent( min(equi1_extent[1], good1_extent[1]) ,
                         max(equi1_extent[2], good1_extent[2]) ,
                         min(equi1_extent[3], good1_extent[3]) ,
                         max(equi1_extent[4], good1_extent[4]))

# plot(all1_equi[[n1]]@coords, xlim=c(output_extent[1], output_extent[2]), ylim=c(output_extent[3], output_extent[4]))
# plot(allgood1[[n1]]@coords, xlim=c(output_extent[1], output_extent[2]), ylim=c(output_extent[3], output_extent[4]))




# comparison 2: use n of 2 or 3
n2 <- 2

equi2_extent <- extent(all2_equi[[n2]])
good2_extent <- extent(allgood2[[n2]])
output_extent <- extent( min(equi2_extent[1], good2_extent[1]) ,
                         max(equi2_extent[2], good2_extent[2]) ,
                         min(equi2_extent[3], good2_extent[3]) ,
                         max(equi2_extent[4], good2_extent[4]))


# get probability map
Rtrajectories<-list()

# A
for (i in seq(1:200)){
  Rtrajectory<-create_random_traj(allgood2[[n2]],V,1)
  Rtrajectories<-c(Rtrajectories,Rtrajectory)
}
Rtrajectories_A <- Rtrajectories
Rtrajectories<-list()

# B
for (i in seq(1:200)){
  Rtrajectory<-create_random_traj(all2_equi[[n2]],V,1)
  Rtrajectories<-c(Rtrajectories,Rtrajectory)
}
Rtrajectories_B <- Rtrajectories

ProbabilityRaster <- EncounterProbabilityRaster(Rtrajectories_A, Rtrajectories_B, 50)


# visualize the probability raster
t <- ProbabilityRaster
t[t < 0.075] <- NA

plot(t, col=colorRampPalette(c("snow", "orangered", "black"))(101), 
     ylab="Y (meters)",
     xlab="X (meters)",
     main="Probability (%) of an encounter between Koski and Grivla", 
     sub=mtext("between 2007-08-11 21:00 and 2007-08-11 22:00"))

plot(allgood2[[2]],add=T,pch=19,col='royalblue', cex=1.25)
plot(all2_equi[[2]],add=T,pch=19,col='limegreen', cex=1.25)

legend( "bottomright", 
        legend=c("Koski","Grivla"),
        col=c("royalblue","limegreen"),lwd = 1, lty=c(NA,NA), 
        pch=c(19,19) )

# visualize the random trajectories

plot(t, col=colorRampPalette(c("white"))(101), 
     ylab="Y (meters)",
     xlab="X (meters)",
     main="Random trajectories between known points of Koski and Grivla")
plot(allgood2[[2]],add=T,pch=19,col='royalblue', cex=1.25)
plot(all2_equi[[2]],add=T,pch=19,col='limegreen', cex=1.25)

lapply(Rtrajectories1,function(x) lines(x[1]@coords,add=T,col=rgb(0,1,0,alpha=0.045)))
lapply(Rtrajectories2,function(x) lines(x[1]@coords,add=T,col=rgb(0,0,1,alpha=0.045)))
legend( "bottomright", 
        legend=c("Koski","Grivla"),
        col=c("royalblue","limegreen"),lwd = 1, lty=c(NA,NA), 
        pch=c(19,19) )
