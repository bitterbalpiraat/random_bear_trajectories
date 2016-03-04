create_random_traj <- function(xyt,V,tinterval){
  #' Standardise values
  #'      Arg:
  #'          xyt(DF): Spatialpointdataframe with spatiotemporal points: x coordinate, y coordinate and colomn with time. 
  #'          V(int) = maximum Speed
  #'          tinterval = number of points to add between measured points
  #'      Return:
  #'        	A random trajectory(SpatialPointsDataFrame) 
  
  tinterval<-tinterval+1
  allpoints<-xyt[length(xyt),]
  xy<-c()
  time<-c()
  tvector<- sample(seq(1:(tinterval-1)))
  
  for(j in seq(1:(length(xyt)-1))){
    startpoint<-xyt[j,]
    endpoint<- xyt[j+1,]
    newpoints<-xyt[j,]
    npoint<-startpoint
    last_t<-0
    t<-(xyt[j+1,]@data$tspannum-xyt[j,]@data$tspannum)/tinterval
    
    for(i in tvector){
      if (i<last_t)
      { endpoint<-npoint
      }
      if (i>last_t)
      { startpoint<-npoint
      }
      buffer1<-gBuffer(startpoint,width=(V*(t*i)),quadsegs=7)
      buffer2<-gBuffer(endpoint,width=(V*(t*(tinterval-i))),quadsegs=7)
      PpA<-gIntersection(buffer1,buffer2) 
      npoint<-spsample(PpA,1,type = 'random')
      npoint<-SpatialPointsDataFrame(npoint,data=as.data.frame(xyt@data$tspannum[j]+(i*t)))
      names(npoint)<-'tspannum'
      newpoints<-spRbind(newpoints,npoint)
      last_t<-i
    }
    xy<-rbind(xy,newpoints@coords)
    time <- c(time,newpoints@data$tspannum)
  }
  
  allpoints<-SpatialPointsDataFrame(xy[order(time),],data=as.data.frame(c(time[order(time)])))
  names(allpoints)<- 'tspannum'
  allpoints<-spRbind(allpoints,xyt[length(xyt),])
  
  return(allpoints)
}
