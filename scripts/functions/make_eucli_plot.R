# Function to plot the relationship between the Euclidean distance and the genomic offset

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
         bty="n", pch=pch, cex=cex_legend, ncol=1,
         title="Main gene pools")
  recordPlot()
}
