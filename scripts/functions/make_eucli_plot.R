# Function to plot the relationship between the Euclidean distance and the genomic offset

make_eucli_plot <- function(X,Y,colors, color_names, plot_title=NULL, xlim=NULL,ylim=NULL,cex=1.2,pch=19,ylab,legend_position="topleft") {
  plot(X, Y, 
       xlab ="Euclidean distance",  ylab =ylab, 
       cex = cex , pch = pch, 
       col = colors,
       xlim = xlim, ylim = ylim,
       main = plot_title)
  legend(legend_position, 
         legend=unique(color_names),
         col=unique(colors),
         bty="n", pch=pch, cex=cex, ncol=1,
         title="Main gene pools")
  recordPlot()
}