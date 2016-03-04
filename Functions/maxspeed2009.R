maxspeed2009 <- function(data,bearname){
  #' Standardise values
  #'      Arg:
  #'          bearname(str): name of the bear
  #'      Return:
  #'        	max speed of bear
  
  bear <- subset(data, PubName==bearname)
  
  n <- nrow(bear)
  
  tdif <- difftime(bear$LMT_date[2:n], bear$LMT_date[1:(n-1)],
                   units="hours")
  tdif <- as.numeric(tdif)
  speed <- 0.001 *sqrt((bear$Locale_N[2:n] - 
                          bear$Locale_N[1:(n-1)])^2 + 
                         (bear$Locale_E[2:n] - bear$Locale_E[1:(n-1)])^2)/tdif
  mspeed<-max(speed)
  return(mspeed)
}








