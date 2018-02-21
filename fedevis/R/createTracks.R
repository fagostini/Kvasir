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
#' @import cowplot
#' @export
#' @examples
#' pl = plotanno(anno.pl, range(region.pl[, list(start, end)]))

createTracks <- function(x, track.names=NULL, region, use.smooth=TRUE, win.smooth=1001){
   
   bed = checkList(x)

   if( bed %in% "GRanges" ){
      bed = x
   }else{
      bed = lapply(x, rtracklayer::import)
   }

   message("Subset the GRanges...")
   bed = lapply(bed,
                function(x)
                   subsetByOverlaps(x, resize(region, width(region)+2.5e3, fix="center")))
   for( i in seq_along(bed) ){
      seqinfo(bed[[i]]) = seqinfo(TxDb)[seqnames(seqinfo(bed[[i]])),] 
      bed[[i]] = keepSeqlevels(bed[[i]], value=seqlevels(region))
   }

   message("Calculate coverage...")
   bed = lapply(bed,
                function(x){
                   temp = ifelse(use.smooth, runmean(coverage(x), win.smooth), coverage(x))
                   data.table(as.data.frame(temp[region]))
   })

   for( i in seq_along(bed) ){
      bed[[i]][, `:=`(index = start(region):end(region), group=NULL, group_name=NULL)]
      setkey(bed[[i]], index)
   }
   
   bed = lapply(bed, 
      function(x){
        x = x[index>=start(region) & index<=end(region)]
        x[, i.next :=  data.table::shift(value, type="lag")]
        x[, i.prev :=  data.table::shift(value, type="lead")]
        x = x[value!=i.next | value!=i.prev | is.na(i.next) | is.na(i.prev)]
        x[, `:=`(i.next=NULL, i.prev=NULL)]
     })

   # message("Add pseudo-count...")
   # reg = data.table(as.data.frame(coverage(resize(region, width(region)+2.5e3, fix="center"))))
   # reg[, `:=`(index = seq_along(value), group=NULL, group_name=NULL)]
   # setkey(reg, index)
   
   # bed = lapply(bed,
   #              function(x)
   #                 x[reg,][, value := sum(value, i.value, na.rm=TRUE), by="index"][, i.value := NULL])
   
   # message("Apply smooth function...")
   # bed = lapply(bed, 
   #              function(x){
   #                 if( use.smooth ){
   #                    x[, smooth := rollmean(value-1, k = win.smooth, fill = NA)]
   #                 }else{
   #                    x[, smooth := value-1]
   #                 }
   #                 x = x[index>=start(region) & index<=end(region)]
                   
   #                 x[, i.next :=  data.table::shift(smooth, type="lag")]
   #                 x[, i.prev :=  data.table::shift(smooth, type="lead")]
   #                 x = x[smooth!=i.next | smooth!=i.prev | is.na(i.next) | is.na(i.prev)]
   #                 x[, `:=`(i.next=NULL, i.prev=NULL)]
   #              })
   # message("Remove pseudo-count...")
   
   if( is.null(track.names) & is.null(names(bed)) ) 
      track.names = rep("", length(bed))
   names(bed) = track.names
   
   message("Create plot tracks...")
   scale_max = max(sapply(bed, function(x) max(x[, value])))
   ggl = lapply(seq_along(bed),
                function(x)
                   ggplot(bed[[x]], aes(x=index, y=value)) +
                   geom_area(alpha=0.8, aes(fill=as.factor(1))) +
                   scale_y_continuous("Average coverage", limits = c(0, scale_max+(scale_max/10))) +
                   coord_cartesian(xlim=range(as.data.table(region)[, list(start, end)])) +
                   scale_colour_solarized("red") +
                   annotate("text", x = Inf, y = Inf, label = names(bed)[[x]], hjust = 2, vjust = 3, parse = FALSE) + 
                   theme_few() + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) + guides(fill=FALSE, colour=FALSE)
                )
   
   # print(plot_grid(plotlist=ggl, nrow=length(ggl), axis="lr", align="v", rel_heights=rep(1,length(ggl))))
   
   message("Done")
   return(ggl)
}

# ggl = createTracks(files, region, use.smooth=FALSE)
