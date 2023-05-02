

# Function to scale the climatic datasets to be used for the genomic offset predictions

# Function arguments:

# clim_var: set of selected climatic variables
# clim_ref: a dataset of averaged climatic data over the reference period
# clim_pred: a dataset or a list of datasets (eg one for each GCM) averaged over the period for which we want to calculate the genomic offset
# clim_adj: whether the point estimate climatic data of the reference period should be adjusted for elevation
# ref_period: specify the reference period
# id_spatial_points: name of the column with the ID of the spatial points

generate_scaled_clim_datasets <- function(clim_var,
                                          clim_ref = NULL,
                                          clim_pred = NULL,
                                          clim_ref_adj = FALSE, 
                                          ref_period = "ref_1901_1950"){

if(clim_ref_adj==TRUE){adj <- "ADJ"} else {adj <- "noADJ"}  

# Load climatic data of the reference period
if(is.null(clim_ref)){
clim_ref <- readRDS(here(paste0("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationPointEstimates_ReferencePeriods_",adj,".rds")))[[ref_period]]$ref_means
}
  
  
# Extract scaling parameters, i.e. mean and variance  
scale_params <- lapply(clim_var, function(x){
    
    vec_var <- clim_ref[,x] %>% pull()
    
    list(mean = mean(vec_var),
         sd = sd(vec_var))
    
  }) %>% setNames(clim_var)


# Scale climatic data across the reference period
clim_ref <- clim_ref %>% 
  dplyr::select(any_of(c(colnames(clim_ref)[1],clim_var))) %>% 
  dplyr::mutate(across(where(is.numeric), ~ (. - mean(.)) / sd(.)))
  

# Load climatic data for the period for which we want the make the genomic offset predictions
if(is.null(clim_pred)){
clim_pred <- readRDS(here("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationValuesExtractedFromRasters_FiveGCMs_2041-2070_SSP370.rds"))
}

# If those climatic data are a list (eg a list of future climatic data, each element corresponding to predictions from a different GCM)
if(!is_tibble(clim_pred)){
clim_pred <- lapply(clim_pred, function(x){
  
x <- x %>% dplyr::select(pop,gcm,any_of(clim_var))

# We scale the future climatic variables with the scaling parameters of the climatic values from the period 1901-1950
for(i in clim_var){
  x[,i] <- (x[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
}

return(x)

})} else{ # If those climatic data are not a list (eg climatic data of the NFI plots for the survey periods specific to each plot)
 
clim_pred <- clim_pred %>% dplyr::select(plotcode,any_of(clim_var)) 

for(i in clim_var){
  clim_pred[,i] <- (clim_pred[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
}

}


clim_df <- list(clim_ref = clim_ref,
                clim_pred = clim_pred)
  
}
