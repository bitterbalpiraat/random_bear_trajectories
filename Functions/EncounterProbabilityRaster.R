EncounterProbabilityRaster <- function(listOfSpatialPointsDF_A, listOfSpatialPointsDF_B, cellsize){
  #' Standardise values
  #'      Arg:
  #'          ListOfSpatialPointsDF_A(list): A list of SpatialPointsDataframes containing the coordinates of the points, thus the paths, between given points of track A.
  #'          ListOfSpatialPointsDF_B(list): A list of SpatialPointsDataframes containing the coordinates of the points, thus the paths, between given points of track B.
  #'          cellsize(number): The cellsize (in meters) of the output.
  #'      Return:
  #'        	A raster containing the propability of an encounter between the two given input's.
  
  
  
  
  ## 1 - convert the spPointsDF's to a spLines
  pathList_A <- c()
  pathList_B <- c()
  
  #       Track A
  for(i in listOfSpatialPointsDF_A){
    line_obj <- Lines(Line(i@coords), "1")
    line_sp <- SpatialLines(list(line_obj), proj4string=CRS("+init=epsg:2400"))
    pathList_A <- c(pathList_A, line_sp)
  }
  max_A <- length(pathList_A)
  
  #       Track B
  for(i in listOfSpatialPointsDF_B){
    line_obj <- Lines(Line(i@coords), "1")
    line_sp <- SpatialLines(list(line_obj), proj4string=CRS("+init=epsg:2400"))
    pathList_B <- c(pathList_B, line_sp)
  }
  max_B <- length(pathList_B)
  
  
  
  
  
  
  
  
  ## 2 - specify output extent based on input extents
  maxExtent_A <- vector()
  maxExtent_B <- vector()
  outRasterExtent <- vector()
  
  #       Track A
  for( i in 1:length(pathList_A)){
    if(i == 1){
      maxExtent_A['xmin'] <- xmin(pathList_A[[i]])
      maxExtent_A['xmax'] <- xmax(pathList_A[[i]])
      maxExtent_A['ymin'] <- ymin(pathList_A[[i]])
      maxExtent_A['ymax'] <- ymax(pathList_A[[i]])
    } else {
      if(extent(pathList_A[[i]])[1] < maxExtent_A['xmin']){
        maxExtent_A['xmin'] <- extent(pathList_A[[i]])[1] }
      if(extent(pathList_A[[i]])[2] > maxExtent_A['xmax']){
        maxExtent_A['xmax'] <- extent(pathList_A[[i]])[2] }
      if(extent(pathList_A[[i]])[3] < maxExtent_A['ymin']){
        maxExtent_A['ymin'] <- extent(pathList_A[[i]])[3] }
      if(extent(pathList_A[[i]])[4] > maxExtent_A['ymax']){
        maxExtent_A['ymax'] <- extent(pathList_A[[i]])[4] }
    }
  }
  maxExtent_A <- extent(maxExtent_A['xmin'],
                        maxExtent_A['xmax'],
                        maxExtent_A['ymin'],
                        maxExtent_A['ymax'])
  
  #       Track B
  for( i in 1:length(pathList_B)){
    if(i == 1){
      maxExtent_B['xmin'] <- xmin(pathList_B[[i]])
      maxExtent_B['xmax'] <- xmax(pathList_B[[i]])
      maxExtent_B['ymin'] <- ymin(pathList_B[[i]])
      maxExtent_B['ymax'] <- ymax(pathList_B[[i]])
    } else {
      if(extent(pathList_B[[i]])[1] < maxExtent_B['xmin']){
        maxExtent_B['xmin'] <- extent(pathList_B[[i]])[1] }
      if(extent(pathList_B[[i]])[2] > maxExtent_B['xmax']){
        maxExtent_B['xmax'] <- extent(pathList_B[[i]])[2] }
      if(extent(pathList_B[[i]])[3] < maxExtent_B['ymin']){
        maxExtent_B['ymin'] <- extent(pathList_B[[i]])[3] }
      if(extent(pathList_B[[i]])[4] > maxExtent_B['ymax']){
        maxExtent_B['ymax'] <- extent(pathList_B[[i]])[4] }
    }
  }
  maxExtent_B <- extent(maxExtent_B['xmin'],
                        maxExtent_B['xmax'],
                        maxExtent_B['ymin'],
                        maxExtent_B['ymax'])
  
  #       Output Extent
  outRasterExtent['xmin'] <- ifelse(maxExtent_A[1] < maxExtent_B[1], maxExtent_A[1], maxExtent_B[1])
  outRasterExtent['xmax'] <- ifelse(maxExtent_A[2] > maxExtent_B[2], maxExtent_A[2], maxExtent_B[2])
  outRasterExtent['ymin'] <- ifelse(maxExtent_A[3] < maxExtent_B[3], maxExtent_A[3], maxExtent_B[3])
  outRasterExtent['ymax'] <- ifelse(maxExtent_A[4] > maxExtent_B[4], maxExtent_A[4], maxExtent_B[4])
  
  outRasterExtent <- extent(outRasterExtent['xmin'],
                            outRasterExtent['xmax'],
                            outRasterExtent['ymin'],
                            outRasterExtent['ymax'])
  
  
  
  
  
  
  
  
  
  # 3 - LineDensity rasters
  LineDensity_A <- raster(outRasterExtent, crs = projection(pathList_A[1]), vals=0, 
                          ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                          nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
  LineDensity_B <- raster(outRasterExtent, crs = projection(pathList_B[1]), vals=0, 
                          ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                          nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
  
  #     create raster track A
  for(i in pathList_A){
    #       create single raster per line
    singleRaster <- raster(outRasterExtent, crs = projection(i), vals=0,
                           ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                           nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
    singleRaster <- rasterize(i, singleRaster, fun='count', background=0)
    #       add the single raster to the final raster
    LineDensity_A <- LineDensity_A + singleRaster
  }
  
  #     create raster track B
  for(i in pathList_B){
    #       create single raster per line
    singleRaster <- raster(outRasterExtent, crs = projection(i), vals=0,
                           ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                           nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
    singleRaster <- rasterize(i, singleRaster, fun='count', background=0)
    #       add the single raster to the final raster
    LineDensity_B <- LineDensity_B + singleRaster
  }
  
  
  
  
  
  
  
  # 4 - ProbabilityRaster
  ProbabilityRaster <- (LineDensity_A * LineDensity_B) / (max_A * max_B) * 100
  
  
  return(ProbabilityRaster)
}