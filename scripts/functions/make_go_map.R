# Arguments

# go_limits: limits of the GO values for the legend scale across the different plots (ie to help comparison across plots)

# type: mapping genomic predictions for:
          # populations under future climates => type = "pop"
          # populations under current climates at the location of the common gardens => type = "CG"
          # NFI plots 


library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

make_go_map <- function(df, 
                        point_size=2, 
                        x_limits = c(-10, 12),
                        y_limits = c(33, 50),
                        legend_position = "right",
                        legend_box_background = "gray",# "gray80"
                        axis_text_size = 14,
                        go_limits=NULL,
                        type="pop", # either "NFI", "CG" or "pop"
                        cg_coord=NULL,
                        plot_title){

  
  p <-  ggplot() + 
    geom_sf(data = world, fill="gray98") + 
    theme_bw() +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits) + 
    geom_point(data=df, aes(x=longitude, y=latitude, color=GO), size= point_size) + 
    xlab("") + ylab("") +
    ggtitle(plot_title) +
    theme(legend.position = legend_position,
          axis.text = element_text(size=axis_text_size),
          legend.box.background = element_rect(colour = legend_box_background, linewidth=0.6, fill="white"),
          )  +
    scale_color_gradientn(name = "Genomic offset", colours = rev(rainbow(5)), limits=go_limits)
  
  if(type=="CG"){
    p <- p + geom_point(data=cg_coord,
                        aes(x=longitude,y=latitude), 
                        size=5, color="black", shape=8)}
  
  return(p)
}

# 
# make_go_map <- function(dfcoord, 
#                         snp_set, 
#                         point_size=2, 
#                         x_limits = c(-10, 12),
#                         y_limits = c(33, 50),
#                         legend_position = "right",
#                         legend_box_background = "gray",# "gray80"
#                         axis_text_size = 14,
#                         go_limits=NULL,
#                         ggtitle=NULL,
#                         type="pop", # either "NFI", "CG" or "pop"
#                         gcm="GFDL-ESM4",
#                         cg_name=NULL,
#                         cg_coord=NULL){
#   
#   if(type=="NFI"){
#     point_go <- dfcoord %>% mutate(GO=snp_set$go_nfi)
#     plot_title <- paste0(snp_set$set_name)
#   } else if(type=="pop"){
#     point_go <- dfcoord %>% mutate(GO=snp_set$go[[gcm]])  
#     plot_title <- paste0(snp_set$set_name," - ",gcm)
#   } else if(type=="CG"){
#     point_go <- dfcoord %>% 
#       left_join(snp_set$go_cg[,c("pop",cg_name)], by="pop") %>% 
#       dplyr::rename(GO=all_of(cg_name))
#     plot_title <- paste0(str_to_title(cg_name), " - ",snp_set$set_name)
#   }
#   
#   if(is.null(ggtitle)){ggtitle <- plot_title}
#   
#   p <-  ggplot() + 
#     geom_sf(data = world, fill="gray98") + 
#     theme_bw() +
#     scale_x_continuous(limits = x_limits) +
#     scale_y_continuous(limits = y_limits) + 
#     geom_point(data=point_go, aes(x=longitude,y=latitude,color=GO), size=point_size) + 
#     xlab("") + ylab("") +
#     ggtitle(ggtitle) +
#     theme(legend.position = legend_position,
#           axis.text = element_text(size=axis_text_size),
#           legend.box.background = element_rect(colour = legend_box_background, linewidth=0.6, fill="white"),
#     )  +
#     scale_color_gradientn(name = "Genomic offset", colours = rev(rainbow(5)), limits=go_limits)
#   
#   if(type=="CG"){
#     p <- p + geom_point(data=filter(cg_coord, cg == cg_name),
#                         aes(x=longitude,y=latitude), 
#                         size=5, color="black", shape=8)}
#   
#   return(p)
# }
