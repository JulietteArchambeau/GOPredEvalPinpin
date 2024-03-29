##########################################
# function generate_scaled_nfi_clim_datasets 
##########################################

# Function to scale the climatic variables at the location of the NFI plots using the scaling parameters used to estimate the genomic offset 
  # i.e., the mean and variance of the climatic variables at the location of the populations under the reference period 1901-1950

# Function arguments:

# clim_var: set of selected climatic variables
# clim_ref: NFI climatic data under the reference period (1901-1950)
# clim_pred: NFI climatic data under the survey period

generate_scaled_nfi_clim_datasets <- function(clim_var,
                                              clim_ref = NULL,
                                              clim_pred = NULL){
  

# We load the climatic data used to estimate the genomic offset
  # = climatic data at the location of the populations for the 1901-1950 reference period (not adjusted)
 
pop_ref <- readRDS(here(paste0("data/ClimaticData/MaritimePinePops/ClimatePopulationLocationPointEstimates_ReferencePeriods_noADJ.rds")))[["ref_1901_1950"]]$ref_means
  
# We extract the scaling parameters, i.e. mean and variance
scale_params <- lapply(clim_var, function(x){
  
  vec_var <- pop_ref[,x] %>% pull()
    
  list(mean = mean(vec_var),
       sd = sd(vec_var))
    
  }) %>% setNames(clim_var)
  
  
# Scale climatic data across the reference period
  clim_ref <- clim_ref %>% dplyr::select(plotcode,longitude,latitude,any_of(clim_var))
    
  for(i in clim_var){
    clim_ref[,i] <- (clim_ref[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
    }
  
# Scale climatic data across the NFI period
  clim_pred <- clim_pred %>% dplyr::select(plotcode,longitude,latitude,any_of(clim_var))
  
  for(i in clim_var){
    clim_pred[,i] <- (clim_pred[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
  }

  
  clim_df <- list(clim_ref = clim_ref,
                  clim_pred = clim_pred)
  
}


