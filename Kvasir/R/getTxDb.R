#' Create TxDb object function
#' 
#' The function checks if local versions of Gencode annotation GFF3 and the associated Txdb exist, otherwise it dowloads the required file and create the object. 
#' @param species A character indicating either "Homo sapiens" or "Mus musculus".
#' @param genome A character indicating either "Homo sapiens" or "Mus musculus".
#' @param version A character specifying the annotation version (e.g., "27" or "M16").
#' @param keep A logical specifying whether to keep the downloaded GFF3 file. Default is FALSE.
#' @return A TxDb object.
#'
#' @import GenomicFeatures
#' @import AnnotationDbi
#' @importFrom utils download.file
#' @export
#' @examples
#' TxDb = getTxDb("Homo sapiens", "hg38", "27", keep=FALSE)
#' TxDb = getTxDb("Mus musculus", "mm10", "M16", keep=FALSE)

getTxDb <- function(species, genome, version, keep=TRUE){
    if( !file.exists(paste0("gencode.", version, ".annotation.sqlite")) ){
        if( species %in% "Homo sapiens" ){
            if( !file.exists(paste0("gencode.", version, ".annotation.gff3.gz")) )
                download.file(paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                              paste0("gencode.", version, ".annotation.gff3.gz"), quiet = TRUE)
        }else{
            if( !file.exists(paste0("gencode.", version, ".annotation.gff3.gz")) )
                download.file(paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
                              paste0("gencode.", version, ".annotation.gff3.gz"), quiet = TRUE)
        }
        TxDb = makeTxDbFromGFF(file = paste0("gencode.", version, ".annotation.gff3.gz"),
            chrominfo=getChrSizes(genome),
            format = "gff3",
            dataSource = paste("Gencode version", version),
            organism = species)
        saveDb(TxDb, file=paste0("gencode.", version, ".annotation.sqlite"))
        if( file.exists(paste0("gencode.", version, ".annotation.gff3.gz")) & keep==FALSE )
            file.remove(paste0("gencode.", version, ".annotation.gff3.gz"))
    }
    return(loadDb(paste0("gencode.", version, ".annotation.sqlite")))
}
