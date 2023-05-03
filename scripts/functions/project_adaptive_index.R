###################################
# function checking_NAs_matching 
###################################

checking_NAs_matching <- function(x){
  
  # We extract the values from the raster stack  
  vals <- lapply(names(x),function(var) raster::getValues(x[[var]])) %>% 
    setNames(clim_var) %>% 
    as.data.frame() 
  
  # checking that the the function `getValues` returns ALL raster values (not only the values with no NAs)
  # ncell(ref_rasts)==nrow(vals)
  
  # replace no NAs values by 1
  vals[!is.na(vals)] <- 1
  
  # calculate the sum of the rows (ie the values of each climatic variable at each pixel)
  vals <- vals %>% mutate(sum= rowSums(.))
  
  # we want that all the rows sum to 6
  return(length(table(vals$sum)))
  
}


###################################
# function project_adaptive_index 
###################################

# Function to project the adaptive index across the landscape

# Arguments
# =========
# clim_var: selected climatic variables
# clim_ref_adj: TRUE or FALSE, specify whether the point estimate climatic data used to scale the rasters should be adjusted or not for elevation
# ref_period: the reference period used to calculate the adaptive index, can be 1901-1950 or 1960-1991
# range_buffer: a mask can be provided to project the adaptive index on a specific area
# `NULL` => no mask 
# The default mask is based on the EUFORGEN distribution + 10 km around the NFI plots
# method: `loadings` or `predict` 
# K: number of RDA axes used for the projection of the adaptive index 

project_adaptive_index <- function(clim_var,
                                   snp_set,
                                   clim_ref_adj = FALSE, 
                                   K,
                                   range_buffer = shapefile(here('data/Mapping/PinpinDistriEUforgen_NFIplotsBuffer10km.shp')),
                                   ref_period = "ref_1901_1950",
                                   method="loadings"){
  
  # Load point estimate climatic data of the reference period
  if(clim_ref_adj==TRUE){adj <- "ADJ"} else {adj <- "noADJ"}  
  clim_ref_pe <- readRDS(here(paste0("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationPointEstimates_ReferencePeriods_",adj,".rds")))[[ref_period]]
  
  # Extract scaling parameters, i.e. mean and variance  
  scale_params <- lapply(clim_var, function(x){
    
    vec_var <- clim_ref_pe$ref_means[,x] %>% pull()
    
    list(mean = mean(vec_var),
         sd = sd(vec_var))
    
  }) %>% setNames(clim_var)
  
  
  # Reference period rasters
  path <- here(paste0("data/ClimaticData/ClimateDTRasters/1km_",clim_ref_pe$range[[1]],"-",clim_ref_pe$range[[2]],"_Extent-JulietteA/"))
  tif_paths <- lapply(clim_var, function(x) paste0(path,"/",x,".tif"))
  ref_rasts <- raster::stack(tif_paths)
  
  # checking that the CRS is the same for the buffer and the rasters
  if(identical(crs(range_buffer),crs(ref_rasts))==FALSE){stop(paste0("CRS of the range buffer is not the same as the raster CRS."))}  
  
  # checking that the rasters have no NAs
  if(!checking_NAs_matching(ref_rasts)==1){stop(paste0("Presence of NAs in the rasters of the reference period!"))}
  
  # extracting coordinates and climatic values
  vals <- as.data.frame(rasterToPoints(ref_rasts))
  # nrow(var_env_proj_pres) # not the same row number as getValues!
  # I do not understand what `rasterToPoints` exactly does, i.e. why it does not have the same number of rows than `getValues`
  
  # mean-centering the climatic variables
  for(i in clim_var){vals[,i] <- (vals[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd}
  
  # we will store the rasters with the adaptive index in a list, each element correspond to a different RDA axis
  ai_proj <- list()
  
  # predicting pixels genetic component based on RDA axes
  if(method == "loadings"){
    for(i in 1:K){
      ai_rast <- rasterFromXYZ(data.frame(vals[,c("x","y")], Z = as.vector(apply(vals[,clim_var], 1, function(x) sum( x * snp_set$rda_model$CCA$biplot[,i])))), crs = crs(ref_rasts))
      names(ai_rast) <- paste0("RDA", as.character(i))
      ai_proj[[i]] <- ai_rast
      names(ai_proj)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Prediction with RDA model and linear combinations
  if(method == "predict"){ 
    pred <- predict(snp_set$rda_model, vals[,clim_var], type = "lc")
    for(i in 1:K){
      ai_rast <- rasterFromXYZ(data.frame(vals[,c("x","y")], Z = as.vector(pred[,i])), crs = crs(ref_rasts))
      names(ai_rast) <- paste0("RDA", as.character(i))
      ai_proj[[i]] <- ai_rast
      names(ai_proj)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Mask with the range if supplied
  if(!is.null(range_buffer)){
    ai_proj <- lapply(ai_proj, function(x) mask(x, range_buffer))
  }
}

