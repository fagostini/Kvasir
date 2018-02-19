#' Plot annotation track function
#' 
#' The function extracts the annotated features from the TxDb object and distributes them on different plotting layers to prevent overlaps. 
#' @param data A List object (the output of createAnnoTrack() ).
#' @param reg.range A vector of two integers, which define start and end for the plot left- and right-most borders.
#' @param gene.height A positive numeric value used as a multiplying factor for the gene height; the default is 1.
#' @param colorVar A character string indicating the colour of the gene lines/arrows; the default is grey.
#' @param clean A logical indicating whether to plot the frame around the annotation track.
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import ggthemes
#' @export
#' @examples
#' pl = plotanno(anno.pl, range(region.pl[, list(start, end)]))

plotAnno <- function(data=list(), reg.range=NULL, gene.height = 1, colorVar="grey30", clean=FALSE){
   plot <- ggplot() +
      geom_segment(data=data[["Arrows"]][strand=="+"], mapping=aes(x=start, xend=end, y=Y, yend=Y), colour=colorVar, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
      geom_segment(data=data[["Arrows"]][strand=="-"], mapping=aes(x=end, xend=start, y=Y, yend=Y), colour=colorVar, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
      geom_segment(data=data[["Genes"]], mapping=aes(x=start, xend=start, y=Y-(0.05*gene.height)/2, yend=Y+(0.05*gene.height)/2), colour=colorVar, lineend="square") +
      geom_segment(data=data[["Genes"]], mapping=aes(x=end,   xend=end,   y=Y-(0.05*gene.height)/2, yend=Y+(0.05*gene.height)/2), colour=colorVar, lineend="square") + 
      geom_rect(data=data[["Exons"]], mapping=aes(xmin=start, xmax=end, ymin=Y-(0.05*gene.height), ymax=Y+(0.05*gene.height), fill=as.factor(1))) +
      coord_cartesian(xlim=reg.range) +
      scale_fill_tableau() + scale_colour_tableau() +
      theme_few() + theme(axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) + guides(fill=FALSE, colour=FALSE) 
   if( clean )
      plot = plot + theme_blank()
   print(plot)
   return(plot)
}
