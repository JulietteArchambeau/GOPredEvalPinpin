# Function to 

# Arguments

# go_limits: limits of the GO values for the legend scale across the different plots (ie to help comparison across plots)

# type: mapping genomic predictions for:
          # populations under future climates => type = "pop"
          # populations under current climates at the location of the common gardens => type = "CG"
          # NFI plots 


library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

make_go_map <- function(dfcoord, 
                        snp_set, 
                        point_size=2, 
                        x_limits = c(-10, 12),
                        y_limits = c(33, 50),
                        go_limits=NULL,
                        ggtitle=NULL,
                        type="pop", # either "NFI", "CG" or "pop"
                        gcm="GFDL-ESM4",
                        cg_name=NULL,
                        cg_coord=NULL){
  
  if(type=="NFI"){
    point_go <- dfcoord %>% mutate(GO=snp_set$go_nfi)
    plot_title <- paste0(snp_set$set_name)
  } else if(type=="pop"){
    point_go <- dfcoord %>% mutate(GO=snp_set$go[[gcm]])  
    plot_title <- paste0(snp_set$set_name," - ",gcm)
  } else if(type=="CG"){
    point_go <- dfcoord %>% 
      left_join(snp_set$go_cg[,c("pop",cg_name)], by="pop") %>% 
      dplyr::rename(GO=all_of(cg_name))
    plot_title <- paste0(str_to_title(cg_name), " - ",snp_set$set_name)
  }
  
  if(is.null(ggtitle)){ggtitle <- plot_title}
  
  p <-  ggplot() + 
    geom_sf(data = world, fill="gray98") + 
    theme_bw() +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits) + 
    geom_point(data=point_go, aes(x=longitude,y=latitude,color=GO), size=point_size) + 
    xlab("") + ylab("") +
    ggtitle(ggtitle) +
    scale_color_gradientn(name = "GO", colours = rev(rainbow(5)), limits=go_limits)
  
  if(type=="CG"){
    p <- p + geom_point(data=filter(cg_coord, cg == cg_name),
                        aes(x=longitude,y=latitude), 
                        size=5, color="black", shape=8)}
  
  return(p)
}
