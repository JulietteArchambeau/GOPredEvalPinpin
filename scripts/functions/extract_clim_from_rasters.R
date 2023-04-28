# =======================================================
# Function to extract the climatic variables from rasters
# =======================================================

extract_clim_from_rasters <- function(x,
                             gcm,
                             period,
                             ssp,
                             path="data/ClimaticData/ClimateDTRasters/",
                             pop_coord){
  
  if(period=="1901-1950"){
    path <- paste0(path,"1km_1901-1950_Extent-JulietteA/",x,".tif")} else {
    path <- paste0(path,"1km_",gcm,"_",period,"_",ssp,"_Extent-JulietteA/",x,".tif")
    }
  
  # For the variable Eref, there are some NAs for some populations
  # So we start with a simple extraction (extraction of the cell value if the point falls within the cell)
  # and if there are some NAs, we perform a bilinear extraction
  # bilinear interpolation = the returned values are interpolated from the values of the four nearest raster cells.
  # And we replace the NAs by the extracted values obtained with the bilinear method
  # if there are still NAs, we extract cell values using a buffer of 2km
  # If the distance between the sampling point and the center of a cell is less than or equal to the buffer, the cell is included. 
  # We take the mean (fun=mean) of the cell values obtained with the extraction done with the buffer method
  # And we replace the NAs by the values obtained with the buffer method (with a buffer of 2km)
  # If there are still NAs, we do the same with a buffer of 3km
  
  ext_simple <- raster::raster(here(path)) %>% raster::extract(pop_coord)
  i <- is.na(ext_simple)
  
  if(length(which(i==TRUE)) !=0){
    
    ext_bilinear <- raster::raster(here(path)) %>% raster::extract(pop_coord,method="bilinear")
    ext_simple[i] <- ext_bilinear[i]
    i <- is.na(ext_simple)
    
    if(length(which(i==TRUE)) !=0){
      
      ext_buff <- raster::raster(here(path)) %>% raster::extract(pop_coord,buffer=2000,fun=mean)
      ext_simple[i] <- ext_buff[i]
      i <- is.na(ext_simple)
      
      if(length(which(i==TRUE)) !=0){
        
        ext_buff <- raster::raster(here(path)) %>% raster::extract(pop_coord,buffer=3000,fun=mean)
        ext_simple[i] <- ext_buff[i]
        i <- is.na(ext_simple)
        
        
      }
      
    }
    
  }
  
  return(ext_simple)
  
}
