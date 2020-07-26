#' Create annotation track function
#' 
#' The function extracts the annotated features from the TxDb object and distributes them on different plotting layers to prevent overlaps. 
#' @param TxDb A TxDb object.
#' @param region A GRanges object with a single entry defining the genomic region to visualise.
#' @param min.unit An integer used to calculate the 'unit': width(region)/min.unit; the default is 100.
#' @param species [Optional] A character indicating either "Homo sapiens" or "Mus musculus". This is used only if the TxDb object is not provided.
#' @param genome [Optional] A character indicating the genome assembly (e.g. "hg19", "hg38" or "mm10"). This is used only if the TxDb object is not provided.
#' @param version [Optional] A character specifying the annotation version (e.g. "27" or "M16"). This is used only if the TxDb object is not provided.
#' @return A three-levels list containing the information to plot the track.
#'
#' @import GenomicFeatures
#' @import data.table
#' @export
#' @examples
#' TxDb = getTxDb("Homo sapiens", "hg19", "19", keep=FALSE)
#'
#' anno = createAnnoTrack(TxDb, region=GRanges("chr14", IRanges(50027696, 50277302)))

createAnnoTrack <- function(TxDb=NULL, region, min.unit=100, species=NULL, genome=NULL, version=NULL){

   message("Check the annotation...")
   if( is.null(TxDb) ){
      if( is.null(species) | is.null(genome) | is.null(version) ) { stop("Either provide a TxDb object or species, assembly and version for the annotation!", call.=FALSE); return(1) }
      message("Annotation not found. Download GFF...")
      TxDb = getTxDb(species, genome, version)
   }
   
   message("Extract features...")
   genes = suppressMessages(genes(TxDb))
   genes = genes[!grepl("PAR", names(genes)),]
   
   exons = exonsBy(TxDb, by="gene")
   exons = exons[!grepl("PAR", names(exons)),]
   
   genes = genes[overlapsAny(genes, region, ignore.strand=TRUE),]
   
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
   overlaps[, `:=`(queryIDs = genes[queryHits,]$gene_id, subjectIDs = genes[subjectHits,]$gene_id)]

   level = 0
   sel = which(!seq_along(exons)%in%overlaps[, queryHits])
   if( length(sel)>0 ){
      genes[sel,]$Y = level
      for( i in which(!seq_along(exons)%in%overlaps[, queryHits]) )
         exons[[i]]$Y = level 
   }
   while( nrow(overlaps)>0 ){
      temp = copy(overlaps)
      temp = split(temp, temp[, queryIDs])
      temp = temp[order(-sapply(temp, nrow))]

      for( i in names(temp) ){
         if( i%in%names(temp) ){
            temp = temp[!names(temp)%in%temp[[i]][, subjectIDs]]
         }
      }
      genes[names(temp),]$Y = level
      for( i in names(temp) )
         exons[[i]]$Y = level
      
      overlaps = overlaps[!queryIDs%in%names(temp),]
      level = level-0.25
   }
   
   message("Create track data...")
   myList = list()
   myList[["Genes"]] = data.table(as.data.frame(genes), key="gene_id")
   myList[["Exons"]] = rbindlist(lapply(exons, as.data.table), idcol="gene")
   
   if( !is.null(species) ){
      sp = paste(sapply(tstrsplit(species, " "), function(x) sapply(strsplit(x, ""), '[[', 1)), collapse="")
      pkg = paste("org", sp, "eg.db", sep=".")
      if( require(pkg, character.only=TRUE) ){
         dt = try( as.data.table(select(get(pkg), keys=myList[["Genes"]]$gene_id, keytype="ENTREZID", columns="SYMBOL")), silent = TRUE )
         if( class(dt)[1] == "try-error" ){
            dt = try( as.data.table(select(get(pkg), keys=gsub("\\.*.$", "", myList[["Genes"]]$gene_id), keytype="ENSEMBL", columns="SYMBOL")), silent = TRUE )
         }
         if( class(dt)[1] != "try-error" ){
            setnames(dt, c("ID", "SYMBOL"))
            myList[["Genes"]][, gene_name := dt[match(gene_id, ID), SYMBOL]]
         }
      }
   }

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
