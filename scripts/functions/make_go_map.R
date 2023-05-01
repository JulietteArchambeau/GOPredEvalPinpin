library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

make_go_map <- function(dfcoord, 
                        snp_set, 
                        point_size=2, 
                        x_limits = c(-10, 12),
                        y_limits = c(33, 50),
                        gcm="GFDL-ESM4",
                        ggtitle=NULL,
                        go_limits=NULL,
                        CG=F,
                        cg_name=NULL,
                        cg_coord=NULL){
  
  if(CG==FALSE){
    point_go <- dfcoord %>% mutate(GO=snp_set$go[[gcm]])  
    plot_title <- paste0(snp_set$set_name," - ",gcm)
  } else if(CG==TRUE){
    point_go <- dfcoord %>% 
      left_join(snp_set$go[,c("pop",cg_name)], by="pop") %>% 
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
  
  if(CG==T){
    p <- p + geom_point(data=cg_coord[cg_coord$site==cg_name,],
                        aes(x=longitude,y=latitude), 
                        size=5, color="black", shape=8)}
  
  return(p)
}