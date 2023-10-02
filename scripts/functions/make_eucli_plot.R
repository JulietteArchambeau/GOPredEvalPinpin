# Functions to plot the relationship between the Euclidean distance and the genomic offset
# ========================================================================================

make_eucli_plot <- function(X,Y,
                            colors, color_names, 
                            plot_title=NULL, 
                            xlim=NULL,ylim=NULL,
                            cex=1.2,pch=19,
                            cex_legend=0.8,
                            ylab,xlab ="Euclidean distance",
                            legend_position="topleft") {
  plot(X, Y, 
       xlab = xlab,  ylab =ylab, 
       cex = cex , pch = pch, 
       col = colors,
       xlim = xlim, ylim = ylim,
       main = plot_title)
  legend(legend_position, 
         legend=unique(color_names),
         col=unique(colors),
         bty="n", pch=pch,
         cex=cex_legend, ncol=1,
         title="Main gene pools")
  recordPlot()
}



# With ggplots2
# =============
make_ggscatterplot <- function(x,y,title, max_go, range_eucli){
  
  gps <- readRDS(here("data/GenomicData/MainGenePoolPopulations.rds")) %>%  arrange(pop)
  
  
  df <- data.frame(x, y,
                   GP = factor(gps$main_gp_pop, levels=unique(gps$main_gp_pop)))
  
  df %>% ggplot(aes(x = x, y = y)) +
    geom_smooth(method = "lm", se=FALSE, color="gray40", formula = y ~ x) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.2, vjust = 1.5,
             label=TeX(paste0("$R^2$ = ",format(summary(lm(y ~ x, df))$r.squared, digits = 3)))) +
    # geom_text(x=-Inf, y=Inf,
    #           hjust = -0.2, vjust = 1.5,
    #          label=TeX(paste0("$R^2$ = ",format(summary(lm(y ~ x, df))$r.squared, digits = 3)))) +
    geom_point(aes(colour=GP)) +
    xlab("Euclidean distance") +
    ylab("Genomic offset") +
    coord_cartesian(xlim = c(range_eucli[[1]]-0.04, range_eucli[[2]]+0.04), 
                    ylim = c(0, max_go)) +
    ggtitle(title) +
    scale_color_manual(values=unique(gps$color_main_gp_pop),
                       name="Gene pools") +
    theme_bw() +
    theme(plot.title = element_text(size = 11),
          legend.text = element_text(size=11),
          legend.title = element_text(size=12))
  
}
