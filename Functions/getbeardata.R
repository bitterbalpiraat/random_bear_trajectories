getbeardata <- function(data,Pubname,tstart){
  #' Standardise values
  #'      Arg:
  #'          data(dataframe): data
  #'          bearname(str): name of the bear
  #'      Return:
  #'        	point tjacetory bear
  
  # SUBSETS and ordering
  
  data$LMT_date <- as.POSIXct(data$LMT_date, 
                              format='%d-%m-%Y %H:%M:%S')
  beardata <- subset(data, PubName == Pubname)
  beardata <- beardata[order(beardata$LMT_date),]      # order on data_time
  
  beardata$tspan <- difftime(beardata$LMT_date, tstart, units="hours")
  
  mypoints <-  SpatialPointsDataFrame(cbind(beardata$Locale_E, beardata$Locale_N), data =beardata,
                                      proj4string=CRS("+init=epsg:2400"))
  mypoints@data$tspannum<- as.numeric(mypoints@data$tspan)
  mypoints<-SpatialPointsDataFrame(mypoints@coords,as.data.frame(mypoints@data$tspannum))
  names(mypoints)<-'tspannum'
  return(mypoints)
}
