
library(rnaturalearth)


make_high_go_pop_maps <- function(pop_coord,
                                  list_go = NULL,
                                  point_size=2, 
                                  x_limits = c(-10, 12),
                                  y_limits = c(33, 50),
                                  ggtitle=NULL,
                                  nb_id_pop = 5){
  
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  colfunc <- colorRampPalette(c("darkolivegreen2","gold","red"))
  
  df_sum_go <- lapply(names(list_go), function(gcm){
    
    pop_coord %>% 
      mutate(GO=list_go[[gcm]]) %>% 
      arrange(GO) %>% 
      mutate(GO = c(rep(1,nb_id_pop),rep(0,(nrow(pop_coord)-nb_id_pop))))
    
  }) %>% 
    setNames(names(list_go)) %>% 
    bind_rows(.id="GCM") %>% 
    pivot_wider(names_from = "GCM", values_from = "GO") %>% 
    mutate(sum_go = (`GFDL-ESM4` + `IPSL-CM6A-LR` + `MPI-ESM1-2-HR` + `MRI-ESM2-0` + `UKESM1-0-LL`) %>% as.factor) %>% 
    arrange(sum_go)

  
  ggmap <- ggplot() + 
    geom_sf(data = world, fill="gray98") + 
    theme_bw() +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits) + 
    geom_point(data=df_sum_go, aes(x=longitude,y=latitude,color=sum_go), size=3) + 
    xlab("") + ylab("") +
    ggtitle(ggtitle) +
    theme(legend.position = c(0.8,0.2),
          legend.box.background = element_rect(colour = "gray50"))  +
    scale_colour_manual(values=colfunc(6), 
                        name = "Number of selections")
  
  
  list(df_sum_go,ggmap)
}
