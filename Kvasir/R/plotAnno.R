#' Plot annotation track function
#' 
#' The function extracts the annotated features from the TxDb object and distributes them on different plotting layers to prevent overlaps. 
#' @param data A List object (the output of createAnnoTrack() ).
#' @param reg.range A vector of two integers, which define start and end for the plot left- and right-most borders.
#' @param gene.height A positive numeric value used as a multiplying factor for the gene height; the default is 1.
#' @param colorVar A character string indicating the colour of the gene lines/arrows; the default is grey.
#' @param addID A logical indicating whether to add a label to each gene. IMPORTANT: the data object must have a 'gene_name' or 'gene_id' column. If both are present, the former has priority.
#' @param clean A logical indicating whether to plot the frame around the annotation track.
#' @param ... Other arguments passed on to ‘theme_bw()’.
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import ggthemes
#' @import ggrepel
#' @export

plotAnno <- function(data=list(), reg.range=NULL, gene.height=1, colorVar="grey30", addID=TRUE, clean=FALSE, ...){
   plot <- ggplot() +
      geom_segment(data=data[["Arrows"]][strand=="+"], mapping=aes(x=start, xend=end, y=Y, yend=Y), colour=colorVar, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
      geom_segment(data=data[["Arrows"]][strand=="-"], mapping=aes(x=end, xend=start, y=Y, yend=Y), colour=colorVar, arrow = arrow(length = unit(0.1, "cm"), type="closed")) +
      geom_segment(data=data[["Genes"]], mapping=aes(x=start, xend=start, y=Y-(0.05*gene.height)/2, yend=Y+(0.05*gene.height)/2), colour=colorVar, lineend="square") +
      geom_segment(data=data[["Genes"]], mapping=aes(x=end,   xend=end,   y=Y-(0.05*gene.height)/2, yend=Y+(0.05*gene.height)/2), colour=colorVar, lineend="square") + 
      geom_rect(data=data[["Exons"]], mapping=aes(xmin=start, xmax=end, ymin=Y-(0.05*gene.height), ymax=Y+(0.05*gene.height), fill=as.factor(1))) +
      coord_cartesian(xlim=reg.range, ylim=c(data[["Genes"]][, min(Y)-(0.15*gene.height)/2]-0.1, data[["Genes"]][, max(Y)]+0.1)) +
      scale_fill_tableau() + scale_colour_tableau() +
      theme_few() + theme(axis.title.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank()) + guides(fill=FALSE, colour=FALSE)
   if( addID==TRUE & any(grepl("gene_name|gene_id", colnames(data[["Genes"]]))) )
      if( "gene_name"%in%colnames(data[["Genes"]]) ){
         plot = plot + geom_label_repel(data=data[["Genes"]], mapping=aes(x=(start+end)/2, y=Y-(0.2*gene.height)/2, label=gene_name), vjust=1, size=3)
      }else{
         plot = plot + geom_label_repel(data=data[["Genes"]], mapping=aes(x=(start+end)/2, y=Y-(0.2*gene.height)/2, label=gene_id), vjust=1, size=3)
      }
   if( clean==TRUE ){
      theme_blank <- function(...) {
         ret <- theme_bw(...)
         ret$line <- element_blank()
         ret$rect <- element_blank()
         ret$strip.text <- element_blank()
         ret$axis.text <- element_blank()
         ret$plot.title <- element_blank()
         ret$axis.title <- element_blank()
         # ret$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")
         ret
      }
      plot = plot + theme_blank(...)
   }
   return(plot)
}
