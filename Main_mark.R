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

data2007 <- subset(data2007, Locale_E > 1400000)

data2007$GMT_date <- NULL

data2007$LMT_date <- as.POSIXct(data2007$LMT_date, 
                                
                                format='%d-%m-%Y %H:%M:%S')

mating2009$LMT_date <- as.POSIXct(mating2009$LMT_date, 
                                  
                                  format='%d-%m-%Y %H:%M:%S')




# 2007 bears

Koski2007 <- subset(data2007, PubName == "Koski (2310)")

Koski2007 <- Koski2007[order(Koski2007$LMT_date),]      # order on data_time

Grivla2007 <- subset(data2007, PubName == "Grivla (2911)")

Grivla2007 <- Grivla2007[order(Grivla2007$LMT_date),]   # order on data_time




tstart <- min(min(Koski2007$LMT_date),min(Grivla2007$LMT_date))

Koski2007<-getbeardata(data2007,Pubname = "Koski (2310)",tstart = tstart)

Grivla2007<-getbeardata(data2007,Pubname = "Grivla (2911)",tstart = tstart)

mspeedKoski <- maxspeed2009(data2007,'Koski (2310)') 

mspeedGrivla <- maxspeed2009(data2007,'Grivla (2911)') 

mspeed<-max(mspeedKoski,mspeedGrivla)




#2009 bear

Krut2009 <- subset(mating2009, PubName == "Krut (2937)")

tstartkrut <- min(Krut2009$LMT_date)

Krut2009<-getbeardata(mating2009,Pubname = "Krut (2937)",tstart = tstartkrut )







mspeedKrut <- maxspeed2009(mating2009,'Krut (2937)') 




mypoints<-Krut2009

#mypoints<-subset(mypoints,tspannum>450 & tspannum<700)# mooi lijntje

mypoints<-subset(mypoints,tspannum>680 & tspannum<700)




# define 2 bears for analysis

bear1<-Koski2007[200:1423,]

bear2<-Grivla2007[200:1423,]




#3D plot od trajextories 2 bears




xscale <- c(min(bear2@coords[,1],bear1@coords[,1]),max(bear2@coords[,1],bear1@coords[,1]))

yscale <- c(min(bear2@coords[,2],bear1@coords[,2]),max(bear2@coords[,2],bear1@coords[,2]))

zscale <- c(min(bear2@data$tspannum,bear1@data$tspannum),max(bear2@data$tspannum,bear1@data$tspannum))




open3d(windowRect=c(100,100,700,700))

plot3d(bear1@coords[,1], bear1@coords[,2], bear1@data$tspannum, type="l", col="red", lwd=1.5,
       
       xlab="East", ylab="North", zlab="time", xlim=xscale, ylim=yscale,zlim=zscale)

lines3d(bear2@coords[,1], bear2@coords[,2], bear2@data$tspannum, col="blue", lwd= 1.5)




# 2D INTERSECTION OF TRAJECTORIES (using spatial objects and methods from sp & rgeos libraries)

lbear2 <- Lines(Line(cbind(bear2@coords[,1], bear2@coords[,2])),"1")

lbear2 <- SpatialLines(list(lbear2), proj4string = CRS("+init=epsg:2400"))

lbear2 <- SpatialLinesDataFrame(lbear2, data=data.frame(ID="1",name="Grivla"), match.ID = T)




lbear1<- Lines(Line(cbind(bear1@coords[,1], bear1@coords[,2])),"1")

lbear1 <- SpatialLines(list(lbear1), proj4string = CRS("+init=epsg:2400"))

lbear1 <- SpatialLinesDataFrame(lbear1, data=data.frame(ID="1",name="Koski"), match.ID = T)




# spatial intersection, no time yet

scoin <- gIntersection(lbear1, lbear2)




spplot(lbear2, zcol="name",col.regions="blue", lwd=1, colorkey=FALSE,
       
       xlim=xscale, ylim=yscale, sp.layout=list(list("sp.lines", lbear1, col="red"),
                                                
                                                list("sp.points", scoin, pch=19, col="black")))

########## finding points close trajectory other bear in both space and time

OnLineSegment <- function(x0,y0, x1,y1, x2,y2){
  
  # Returns TRUE if (x0,y0) is on the line segment between (x1,y1) and (x2,y2)
  
  # Thanks to http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
  
  if (x0 < min(x1,x2) || y0 < min(y1,y2) || x0 > max(x1,x2) || y0 > max(y1,y2))
    
    return(F)
  
  else
    
    dist <- ifelse((x2 == x1 && y2 == y1), sqrt((x1-x0)^2 + (y1-y0)^2),
                   
                   abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / sqrt((x2-x1)^2 + (y2-y1)^2))
  
  dist <= 0.1
  
} 




FindTime <- function(x0,y0,index,dframe){
  
  # Based on (x0,y0) Linearly interpolates time between location[index-1] and location[index]
  
  dist1 <- sqrt((dframe@coords[,1][index-1] - x0)^2 + (dframe@coords[,2][index-1] - y0)^2)
  
  dist2 <- sqrt((dframe@coords[,1][index] - x0)^2 + (dframe@coords[,2][index] - y0)^2)
  
  dist3 <- dist1 + dist2
  
  as.numeric(dframe@data$tspannum[index]*dist1/dist3 + dframe@data$tspannum[index-1]*dist2/dist3)
  
}







# apply the functions to the trajectories of both bears

tdataframe <- data.frame(tbear2=numeric(nrow(scoin@coords)), tbear1=numeric(nrow(scoin@coords)))

for (j in 1:nrow(scoin@coords)){
  
  for (ind in 2:nrow(bear1))
    
    if (OnLineSegment(scoin@coords[j,1],scoin@coords[j,2],bear1@coords[,1][ind-1], 
                      
                      bear1@coords[,2][ind-1],bear1@coords[,1][ind], 
                      
                      bear1@coords[,2][ind])){
      
      tdataframe$tbear1[j] <- FindTime(scoin@coords[j,1],scoin@coords[j,2],ind,bear1)
      
      break # once a match is found you are done
      
    }
  
  for (ind in  2:nrow(bear2))
    
    if (OnLineSegment(scoin@coords[j,1],scoin@coords[j,2],bear2@coords[,1][ind-1], 
                      
                      bear2@coords[,2][ind-1],bear2@coords[,1][ind], 
                      
                      bear2@coords[,2][ind])){
      
      tdataframe$tbear2[j] <- FindTime(scoin@coords[j,1],scoin@coords[j,2],ind,bear2)
      
      break # once a match is found you are done
      
    }
  
}

delta_t <- 1.0

stClose <- scoin[which(abs(tdataframe$tbear2 - tdataframe$tbear1) < delta_t, arr.ind=T),]

stTimes <- apply(tdataframe[which(abs(tdataframe$tbear2 - tdataframe$tbear1) < delta_t, arr.ind=T),],1,mean)




SpatLagOvl <- function(x0, y0, line, dframe, delta_s){
  
  pnt <- SpatialPoints(cbind(x0, y0), CRS("+init=epsg:2400"))
  
  buf <- gBuffer(pnt, width=delta_s)
  
  if(gDisjoint(buf,line))
    
    return(F)
  
  else{
    
    # find time point on line
    
    sgm <- gIntersection(buf, line)
    
    if (nrow(sgm@lines[[1]]@Lines[[1]]@coords)!=2)
      
    {
      
      ip <- as.integer((0.5+nrow(sgm@lines[[1]]@Lines[[1]]@coords)/2))
      
      # use centre point
      
      xp <- sgm@lines[[1]]@Lines[[1]]@coords[ip,1]
      
      yp <- sgm@lines[[1]]@Lines[[1]]@coords[ip,2]
      
      tp <- dframe@data$tspannum[which(abs(dframe@coords[,1] - xp) < 0.01 & abs(dframe@coords[,2] - yp) < 0.01)][1]
      
    }
    
    else{
      
      xp <- 0.5*(sgm@lines[[1]]@Lines[[1]]@coords[1,1]+sgm@lines[[1]]@Lines[[1]]@coords[2,1])
      
      yp <- 0.5*(sgm@lines[[1]]@Lines[[1]]@coords[1,2]+sgm@lines[[1]]@Lines[[1]]@coords[2,2])
      
      for (ind in 2:nrow(dframe))
        
        if (OnLineSegment(xp,yp,dframe@coords[,1][ind-1],dframe@coords[,2][ind-1],
                          
                          dframe@coords[,1][ind], dframe@coords[,2][ind]))
          
          tp <- FindTime(xp,yp,ind,dframe)
      
    }
    
    return(tp)
    
  }
  
}




delta_s <- 200.0

near_1 <- numeric()

for (ii in 1: nrow(bear2)){
  
  t1 <- SpatLagOvl(bear2@coords[,1][ii], bear2@coords[,2][ii],
                   
                   lbear1, bear1, delta_s)
  
  if(t1!= F)
    
    near_1 <- rbind(near_1,c(bear2@coords[,1][ii],
                             
                             bear2@coords[,2][ii],bear2@data$tspannum[ii],t1))
  
}




near_2 <- numeric()

for (ii in 1: nrow(bear1)){
  
  t2 <- SpatLagOvl(bear1@coords[,1][ii], bear1@coords[,2][ii],
                   
                   lbear2, bear2, delta_s)
  
  if(t2 != F)
    
    near_2 <- rbind(near_2,c(bear1@coords[,1][ii],
                             
                             bear1@coords[,2][ii],bear1@data$tspannum[ii],t2))
  
}

colnames(near_1) <- c("x","y","tbear2","tbear1")

colnames(near_2) <- c("x","y","tbear1","tbear2")

near_1 <- data.frame(near_1)

near_2 <- data.frame(near_2)

stnear1 <- subset(near_1, abs(tbear2 -tbear1) < delta_t)

stnear2 <- subset(near_2, abs(tbear1 -tbear2) < delta_t)







# plot point close to other lines other bear. for both bears

open3d(windowRect=c(100,100,700,700))

plot3d(bear1@coords[,1], bear1@coords[,2], bear1@data$tspannum, type="l", col="red", lwd=1.5,
       
       xlab="East", ylab="North", zlab="time", xlim=xscale, ylim=yscale,zlim=zscale)

lines3d(bear2@coords[,1], bear2@coords[,2], bear2@data$tspannum, col="blue", lwd= 1.5)

points3d(stnear1[,1], stnear1[,2], stnear1[,3], col="green", size = 10)

points3d(stnear2[,1], stnear2[,2], stnear2[,3], col="yellow", size = 10)







# find points corresponding to buffer match




allgood1<-c()




for (i in seq(1:length(stnear2[,1]))){
  
  goodpoint12<-subset(bear1,tspannum==stnear2['tbear1'][i,] )
  
  loc1<-which(bear1@data$tspannum == stnear2['tbear1'][i,])
  
  goodpoint11<-bear1[(loc1-1),]
  
  goodpoint13<-bear1[loc1+1,]
  
  goodpoints<- spRbind(goodpoint11,goodpoint12)
  
  bear1points<- spRbind(goodpoints,goodpoint13)
  
  allgood1<-c(allgood1,bear1points)
  
}




allgood2<-c()




for (i in seq(1:length(stnear1[,1]))){
  
  goodpoint22<-subset(bear2,tspannum==stnear1['tbear2'][i,] )
  
  loc1<-which(bear2@data$tspannum == stnear1['tbear2'][i,])
  
  goodpoint21<-bear2[(loc1-1),]
  
  goodpoint23<-bear2[loc1+1,]
  
  goodpoints<- spRbind(goodpoint21,goodpoint22)
  
  bear2points<- spRbind(goodpoints,goodpoint23)
  
  allgood2<-c(allgood2,bear2points)
  
}




#Find equivelant point in other trajectory

timebuf<-0.6

all1_equi<-c()

for (i in seq(1:length(allgood1))){
  
  
  
  time<-allgood1[[i]]@data[2,]
  
  
  
  
  allgood1_equi<- subset(bear2, tspannum>(time-timebuf) & tspannum<(time+timebuf))
  
  
  
  
  all1_equi<-c(all1_equi,allgood1_equi)
  
}







all2_equi<-c()

for (i in seq(1:length(allgood2))){
  
  
  
  time<-allgood2[[i]]@data[2,]
  
  
  
  allgood2_equi<- subset(bear1, tspannum>(time-timebuf) & tspannum<(time+timebuf))
  
  
  
  all2_equi<-c(all2_equi,allgood2_equi)
  
}




goodpoint1<-allgood2[[2]]

goodpoint2<-all2_equi[[2]]




plot(goodpoint1,xlim=c(bbox(goodpoint1)[1][1]-750,bbox(goodpoint1)[,2][1]+750),axes=T,
     
     ylim=c(bbox(goodpoint1)[2][1]-750,bbox(goodpoint1)[,2][2]+750))

plot(goodpoint2,add=T,col='red')







## Compute random trajectories for both bears

V=mspeed*1000-1000

Rtrajectories1<-list()




for (i in seq(1:1)){
  
  
  
  Rtrajectory1<-create_random_traj(bear1,V,2)
  
  Rtrajectories1<-c(Rtrajectories1,Rtrajectory1)
  
}




Rtrajectories2<-list()




for (i in seq(1:1)){
  
  
  
  Rtrajectory2<-create_random_traj(bear2,V,2)
  
  Rtrajectories2<-c(Rtrajectories2,Rtrajectory2)
  
  
  
  
  
}




# plot last random created trejectory of both bears

open3d(windowRect=c(100,100,700,700))

plot3d(bear1@coords[,1], bear1@coords[,2], bear1$tspannum, type="p", col="blue", size=2,
       
       xlab="East", ylab="North", zlab="time", xlim=xscale, ylim=yscale,zlim=zscale)

lines3d(Rtrajectory1@coords[,1], Rtrajectory1@coords[,2], Rtrajectory1$tspannum, col='blue', lwd= 0.3)




plot3d(bear2@coords[,1], bear2@coords[,2], bear2$tspannum, type="p", col="red", size=2,
       
       xlab="East", ylab="North", zlab="time", xlim=xscale, ylim=yscale,zlim=zscale,add=T)

lines3d(Rtrajectory2@coords[,1], Rtrajectory2@coords[,2], col='red', 
        
        Rtrajectory2$tspannum, lwd= 0.3)




# plott all trajectories in 2D

plot(bear1,col='yellow',xlim=c(bbox(bear1)[1][1]-(V/5),bbox(bear1)[,2][1]+(V/5)),axes=T,
     
     ylim=c(bbox(bear1)[2][1]-(V/5),bbox(bear1)[,2][2]+(V/5)),xlab='X (meters)', ylab= 'Y (meters)',
     
     main='Random trajectories between known points',pch=8)

lapply(Rtrajectories1,function(x) lines(x[1]@coords,col='red',lwd=0.1))

plot(bear1,col='yellow',add=T,pch=19)




plot(bear2,col='green',xlim=c(bbox(bear2)[1][1]-(V/5),bbox(bear2)[,2][1]+(V/5)),axes=T,
     
     ylim=c(bbox(bear2)[2][1]-(V/5),bbox(bear2)[,2][2]+(V/5)),xlab='X (meters)', ylab= 'Y (meters)',
     
     main='Random trajectories between known points',pch=8)

lapply(Rtrajectories2,function(x) lines(x[1]@coords,col='blue',lwd=0.1))

plot(bear2,col='green',add=T,pch=19)




## Create line density surface of trajectories 

system.time(LineDensity <- LineDensityRaster(Rtrajectories,1000))




## space time event

# bear trajectories need to be subsetted to timne interval first: before the RTG;create_random_traj()




t1=0.0; x1=-1.0; y1=-1.0; 

t2=1.0; x2=-1.5; y2=-1.0; 

v1=4.0;




t3=0.0; x3=1.0; y3=1.0;

t4=2.0; x4=1.0; y4=1.0;

v2=2.26




t1=0.0

t2=0.0

#xyt points for visulization space time event

xyt <- expand.grid(x=seq(1440000,1450000,length=100), y=seq(6811000,6812000,length=100),
                   
                   t=(seq(200, 300,length=100)))




# plot space time events

plot3d(bear1@coords[,1], bear1@coords[,2], bear1@data$tspannum, type="l", col="red", lwd=1.5,
       
       xlab="East", ylab="North", zlab="time", xlim=xscale, ylim=yscale,zlim=zscale,add=T)

lines3d(bear2@coords[,1], bear2@coords[,2], bear2@data$tspannum, col="blue", lwd= 1.5)

points3d(xyt, col="orange")







# create poltgon space time event in 2D

polymatrix<-matrix(c(1440000,6811000,1440000,6812000,1450000,6812000,1450000,6811000,1440000,6811000),
                   
                   nrow = 5,ncol = 2,byrow = T)

polyevent<-list(Polygons(list(Polygon(polymatrix)),ID = 'coords'))

polyevent<- SpatialPolygons(polyevent)




# find intersections between measured points and space time event. should be 0

intersectb1<- gIntersection(bear1tg,polyevent) 

intersectb2<- gIntersection(bear2tg,polyevent) 




#find intersections and calucalte number of intersects for random trajectories of bear 

nintersectsB1<-0




for(i in seq(1:length(Rtrajectories1)) )
  
{
  
  result<-gIntersection(Rtrajectories1[[i]],polyevent)
  
  
  
  
  if (length(result)>0){
    
    nintersectsB1<-nintersectsB1+1
    
  }
  
}




nintersectsB2<-0




for(i in seq(1:length(Rtrajectories2)) )
  
{
  
  result<-gIntersection(Rtrajectories2[[i]],polyevent)
  
  
  
  if (length(result)>0){
    
    nintersectsB2<-nintersectsB2+1
    
  }
  
  
  
}












## Evaluate results




## Visualize results



