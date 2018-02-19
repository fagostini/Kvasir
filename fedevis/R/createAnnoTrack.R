#' Create annotation track function
#' 
#' The function extracts the annotated features from the TxDb object and distributes them on different plotting layers to prevent overlaps. 
#' @param TxDb A TxDb object.
#' @param region A single GRanges object entry defining the genomic region to visualise.
#' @param min.unit An integer used to calculate the 'unit': width(region)/min.unit; the default is 100.
#' @param species [Optional] A character indicating either "Homo sapiens" or "Mus musculus". This is used only if the TxDb object is not provided.
#' @param version [Optional] A character specifying the annotation version (e.g., "27" or "M16"). This is used only if the TxDb object is not provided.
#' @return A three-levels list containing the information to plot the track.
#'
#' @import GenomicFeatures
#' @import data.table
#' @export
#' @examples
#' anno = createAnnoTrack(TxDb, region, min.unit=minimal.unit, species, version)


createAnnoTrack <- function(TxDb=NULL, region, min.unit=100, species=NULL, version=NULL){

   message("Check the annotation...")
   if( is.null(TxDb) ){
      if( is.null(species) | is.null(version) ) { stop("Either provide a TxDb object or species and version for the annotation!", call.=FALSE); return(1) }
      message("Annotation not found. Download GFF...")
      TxDb = getTxDb(species, version)
   }
   
   message("Extract features...")
   genes = genes(TxDb)
   genes = genes[!grepl("PAR", names(genes)),]
   
   exons = exonsBy(TxDb, by="gene")
   exons = exons[!grepl("PAR", names(exons)),]
   
   genes = genes[overlapsAny(genes, region),]
   
   genes = genes[order(width(genes), decreasing=TRUE),]
   
   exons = reduce(exons[names(exons)%in%names(genes)], min.gapwidth=min.unit/10)
   
   message("Apply layering...")
   genes$Y = as.numeric(NA)
   exons = lapply(exons,
                  function(x){
                     x$Y = as.numeric(NA);
                     x
                  })
   
   if( !identical(names(genes), names(exons)) )
      exons = exons[match(names(genes), names(exons))]
   
   overlaps = as.data.table(findOverlaps(genes, ignore.strand=TRUE))[!queryHits==subjectHits,]
   
   level = 0 
   genes[which(!seq_along(exons)%in%overlaps[, queryHits]),]$Y = level
   for( i in which(!seq_along(exons)%in%overlaps[, queryHits]) )
      exons[[i]]$Y = level 
   
   while( nrow(overlaps)>0 ){
      temp = copy(overlaps)
      temp = split(temp, temp[,queryHits])
      
      for( i in names(temp) ){
         if( i%in%names(temp) ){
            temp = temp[!names(temp)%in%as.character(temp[[i]][, subjectHits])]
         }
      }
      genes[as.numeric(names(temp)),]$Y = level
      for( i in as.numeric(names(temp)) )
         exons[[i]]$Y = level
      
      overlaps = overlaps[!queryHits%in%as.numeric(names(temp)),]
      level = level-0.25
   }
   
   message("Create track data...")
   myList = list()
   myList[["Genes"]] = data.table(as.data.frame(genes), key="gene_id")
   myList[["Exons"]] = rbindlist(lapply(exons, as.data.table), idcol="gene")
   
   arrows.pl = genes[width(genes)>min.unit,]
   arrows.pl = split(arrows.pl, names(arrows.pl))
   arrows.pl = lapply(arrows.pl, function(x) data.table(as.data.frame(tile(x, width=min.unit))))
   arrows.pl = rbindlist(arrows.pl, idcol="gene_id")
   setkey(arrows.pl, gene_id)
   myList[["Arrows"]] = arrows.pl[myList[["Genes"]], Y := i.Y]
   
   myList[["Genes"]] = myList[["Genes"]][!is.na(Y)]
   myList[["Exons"]] = myList[["Exons"]][!is.na(Y)]
   myList[["Arrows"]] = myList[["Arrows"]][!is.na(Y)]
   
   message("Done")
   return(myList)
}
