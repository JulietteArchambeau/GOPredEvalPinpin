############################
# calc_avg_clim_var function
############################

# Function to calculate the average of the climatic variables at the location of some spatial points (e.g. populations, NFI plots)

library(magrittr)
library(SPEI)

# ===========================================================
# Window function (needed for the calc_avg_clim_var function) 
# ===========================================================

# Function from the dismo R package
window <- function(x)  { 
  lng <- length(x)
  x <- c(x,  x[1:3])
  m <- matrix(ncol=3, nrow=lng)
  for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
  apply(m, MARGIN=1, FUN=sum)}




# Arguments of the calc_avg_clim_var function
# ===========================================

# clim_df = a dataframe containing:
# the ID of the spatial points
# the longitude, latitude, and elevation of the spatial points
# the annual values of the climatic variables at the location of the spatial points

# ref_period = the climatic variables will be calculated for this period, which has to be c(minimum year,maximum year)

# id_spatial_points = name of the column with the ID of the spatial points


calc_avg_clim_var <- function(clim_df,ref_period,id_spatial_points = "pop"){
  
  # Selecting the years of interest
  if(ref_period[[1]]==2041){
    
    clim_df <- clim_df %>%
      dplyr::filter(year=="2041-2070") %>% # we remove future climatic data (which are not noted as a unique year)
      dplyr::select(-year) %>%
      dplyr::rename(id=all_of(id_spatial_points)) # rename the column with ID of the spatial points to 'id'
    
  } else {
    
    clim_df <- clim_df %>%
      dplyr::filter(!year %in% c("2041-2070","baseline")) %>% # we remove future climatic data (which are not noted as a unique year)
      dplyr::mutate(year=as.numeric(year)) %>% # year column as numeric so that we can remove years after a given date
      { if(length(ref_period)==1) dplyr::filter(.,year %in% ref_period[[1]]) else # we keep the years of the period of interest
        dplyr::filter(.,ref_period[[1]]<=year& year<=ref_period[[2]]) } %>% 
      dplyr::select(-year) %>%
      dplyr::rename(id=all_of(id_spatial_points)) # rename the column with ID of the spatial points to 'id'
  }
  
  
  
  # Calculating the mean of tmn, tmx and prc over the period considered
  tab <- clim_df %>% 
    group_by(id) %>% 
    summarise_at(vars(contains("tmn"),contains("tmx"),contains("prc")),"mean")
  
  tavg <- (tab %>% dplyr::select(contains("tmn")) + tab %>% dplyr::select(contains("tmx")))  / 2
  tmx <- tab %>% dplyr::select(contains("tmx"))
  tmn <- tab %>% dplyr::select(contains("tmn"))
  prc <- tab %>% dplyr::select(contains("prc"))
  wet <- t(apply(prc, 1, window))
  tmp <- t(apply(tavg, 1, window)) / 3
  
  
  tab$bio1 <- apply(tavg ,MARGIN=1, FUN=mean)
  tab$bio2 <- apply(tmx-tmn , MARGIN=1, FUN=mean)
  tab$bio4 <- 100 * apply(tavg, 1, sd)
  tab$bio5 <- apply(tmx, 1, max)
  tab$bio6 <- apply(tmn, 1, min)
  tab$bio7 <- tab$bio5 - tab$bio6
  tab$bio3 <- 100 * tab$bio2 / tab$bio7
  tab$bio12 <- apply(prc, MARGIN=1, FUN= sum)
  tab$bio13 <- apply(prc,1,max)
  tab$bio14 <- apply(prc,1,min)
  tab$bio15 <- apply(prc+1, 1, raster::cv) # the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
  tab$bio16 <- apply(wet,1,max)
  tab$bio17 <- apply(wet,1,min)
  
  
  if (all(is.na(wet))) {
    tab$bio8 <- NA		
    tab$bio9 <- NA		
  } else {
    wetqrt <- cbind(1:nrow(wet), as.integer(apply(wet, 1, which.max)))
    tab$bio8 <- tmp[wetqrt]
    
    dryqrt <- cbind(1:nrow(wet), as.integer(apply(wet, 1, which.min)))
    tab$bio9 <- tmp[dryqrt]
  }
  
  
  
  tab$bio10 <- apply(tmp, 1, max)
  tab$bio11 <- apply(tmp, 1, min) 
  
  if (all(is.na(tmp))) {
    tab$bio18 <- NA		
    tab$bio19 <- NA
  } else {
    hot <- cbind(1:nrow(tmp), as.integer(apply(tmp, 1, which.max)))
    tab$bio18 <- wet[hot]
    
    cold <- cbind(1:nrow(tmp), as.integer(apply(tmp, 1, which.min)))
    tab$bio19 <- wet[cold]
  }
  
  tab$AHM <- (tab$bio1+10)/(tab$bio12*1e-3)
  tab$MWMT <- apply(tavg,MARGIN=1, FUN=max)
  tab$SP <- apply(tab %>% dplyr::select(prc05,prc06,prc07,prc08,prc09),MARGIN=1,FUN=sum)
  tab$SHM <- (tab$MWMT)/((tab$SP+0.1)*1e-3)
  
  
  tab <- clim_df %>%
    dplyr::select(id,contains("ude"),elevation) %>% 
    distinct() %>% 
    right_join(tab,by="id") # adding the latitude and longitude of the spatial points
  
  
  tab$Eref <- lapply(c("tmn","tmx","prc"), function(x) {
    
    tab %>% 
      dplyr::select(id,contains(x)) %>% 
      pivot_longer(cols=contains(x),names_to="variable",values_to = "mean_ref") %>% 
      mutate(month=str_sub(variable,4,-1),
             variable=str_sub(variable,1,3)) 
    
  }) %>% 
    bind_rows() %>% 
    pivot_wider(names_from="variable", values_from=mean_ref) %>% 
    left_join(tab[,c("id","latitude")],by="id") %>% 
    group_by(id) %>% 
    group_split() %>% 
    purrr::map(\(x){ 
      
      hargreaves(Tmin = x$tmn, Tmax = x$tmx, lat = unique(x$latitude), Pre=x$prc, verbose=F) %>% sum() 
      
    }) %>% unlist()
  
  tab <- tab %>% dplyr::rename_with(~paste0(id_spatial_points),id)
  
  return(tab)
  
  
}

