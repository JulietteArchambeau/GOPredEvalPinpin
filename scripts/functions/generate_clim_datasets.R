# Function to generate two sets of present and future climatic variables
# scaled with the past climatic data


generate_clim_datasets <- function(clim_var, clim_past = NULL, clim_fut = NULL){
  
# Load the past climatic data
if(is.null(clim_past)){
clim_past <- read_csv(here("data/DryadRepo/PopulationCoordinatesPastClimateInformation.csv"),
                        show_col_types = FALSE) %>%
  dplyr::select(pop,clim_var$variables) # we keep only the climatic variables of interest
}
  
  
# We extract the scaling parameters, i.e. the mean and the variance  
scale_params <- lapply(clim_var$variables, function(x){
    
    vec_var <- clim_past[,x] %>% pull()
    
    list(mean = mean(vec_var),
         sd = sd(vec_var))
    
  }) %>% setNames(clim_var$variables)

# We scale the past climatic variables
clim_past <- clim_past %>% 
  dplyr::mutate(across(where(is.numeric), ~ (. - mean(.)) / sd(.)))
  
# we load the future climatic data of the climatic variables of interest
if(is.null(clim_fut)){
clim_fut <- read_csv(here("data/DryadRepo/PopulationCoordinatesFutureClimateInformation.csv"),
                       show_col_types = FALSE) %>% 
  dplyr::select(pop,clim_var$variables) # we keep only the climatic variables of interest
}


# We scale the future climatic variables with the scaling parameters of the climatic values from the period 1901-1950
for(i in clim_var$variables){
    clim_fut[,i] <- (clim_fut[,i] - scale_params[[i]]$mean) / scale_params[[i]]$sd
  }

clim_df <- list(clim_past = clim_past,
                clim_fut = clim_fut)
  
}
