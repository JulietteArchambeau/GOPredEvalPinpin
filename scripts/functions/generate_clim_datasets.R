# Function to generate two sets of present and future climatic variables
# scaled with the past climatic data

# Past climatic data is point estimate adjusted or not for elevation

generate_clim_datasets <- function(clim_var, clim_past = NULL, clim_fut = NULL, clim_past_adj = FALSE, ref_period = "ref_1901_1950"){

if(clim_past_adj==TRUE){adj <- "ADJ"} else {adj <- "noADJ"}  

  # Load the past climatic data
if(is.null(clim_past)){
clim_past <- readRDS(here(paste0("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationPointEstimates_ReferencePeriods_",adj,".rds")))[[ref_period]]$ref_means %>%
  dplyr::select(pop,any_of(clim_var)) # we keep only the climatic variables of interest
}
  
  # We extract the scaling parameters, i.e. the mean and the variance  
scale_params <- lapply(clim_var, function(x){
    
    vec_var <- clim_past[,x] %>% pull()
    
    list(mean = mean(vec_var),
         sd = sd(vec_var))
    
  }) %>% setNames(clim_var)

# We scale the past climatic variables
clim_past <- clim_past %>% 
  dplyr::mutate(across(where(is.numeric), ~ (. - mean(.)) / sd(.)))
  
# we load the future climatic data of the climatic variables of interest
if(is.null(clim_fut)){
list_clim_fut <- readRDS(here("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationValuesExtractedFromRasters_FiveGCMs_2041-2070_SSP370.rds"))
}


list_clim_fut <- lapply(list_clim_fut, function(clim_fut){
  
clim_fut <- clim_fut %>% dplyr::select(pop,gcm,any_of(clim_var))

# We scale the future climatic variables with the scaling parameters of the climatic values from the period 1901-1950
for(i in clim_var){
  clim_fut[,i] <- (clim_fut[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
}

return(clim_fut)

})



clim_df <- list(clim_past = clim_past,
                clim_fut = list_clim_fut)
  
}
