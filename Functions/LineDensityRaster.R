LineDensityRaster <- function(listOfSpatialPointsDF,cellsize){
  #' Standardise values
  #'      Arg:
  #'          ListOfSpatialPointsDF(list): A list of SpatialPointsDataframes containing the coordinates of the points, thus the paths, between given points.
  #'          cellsize(number): The cellsize (in meters) of the output.
  #'      Return:
  #'        	A raster containing the amount of lines per cell(Raster) 
  
  
  
  
  # convert the spPointsDF to a spLines
  pathList <- c()
  iterations <- 0
  
  for(i in listOfSpatialPointsDF){
    iterations <- iterations + 1 
    
    line_obj <- Lines(Line(i@coords), "1")
    line_sp <- SpatialLines(list(line_obj), proj4string=CRS("+init=epsg:2400"))
    pathList <- c(pathList, line_sp)
  }
  
  
  # specify output extent based on input extents
  outRasterExtent <- vector()
  
  for( i in 1:length(pathList)){
    if(i == 1){
      outRasterExtent['xmin'] <- xmin(pathList[[i]])
      outRasterExtent['xmax'] <- xmax(pathList[[i]])
      outRasterExtent['ymin'] <- ymin(pathList[[i]])
      outRasterExtent['ymax'] <- ymax(pathList[[i]])
    } else {
      if(extent(pathList[[i]])[1] < outRasterExtent['xmin']){
        outRasterExtent['xmin'] <- extent(pathList[[i]])[1] }
      if(extent(pathList[[i]])[2] > outRasterExtent['xmax']){
        outRasterExtent['xmax'] <- extent(pathList[[i]])[2] }
      if(extent(pathList[[i]])[3] < outRasterExtent['ymin']){
        outRasterExtent['ymin'] <- extent(pathList[[i]])[3] }
      if(extent(pathList[[i]])[4] > outRasterExtent['ymax']){
        outRasterExtent['ymax'] <- extent(pathList[[i]])[4] }
    }
  }
  outRasterExtent <- extent(outRasterExtent['xmin'],
                            outRasterExtent['xmax'],
                            outRasterExtent['ymin'],
                            outRasterExtent['ymax'])
  
  # empty Line Density Raster
  LineDensity <- raster(outRasterExtent, crs = projection(pathList[1]), vals=0, 
                        ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                        nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
  
  # create raster
  for(i in pathList){
    
    # create single raster per line
    singleRaster <- raster(outRasterExtent, crs = projection(i), vals=0,
                           ncols=(outRasterExtent[2]/cellsize) - (outRasterExtent[1]/cellsize), #xmax - xmin
                           nrows=(outRasterExtent[4]/cellsize) - (outRasterExtent[3]/cellsize)) #ymax - ymin
    
    singleRaster <- rasterize(i, singleRaster, fun='count', background=0)
    
    # add the single raster to the final raster
    LineDensity <- LineDensity + singleRaster
  }
  
  
  return(LineDensity)
}