#' Create TxDb object function
#' 
#' The function checks if local versions of Gencode annotation GFF3 and the associated Txdb exist, otherwise it dowloads the required file and create the object. 
#' @param species A character indicating either "Homo sapiens" or "Mus musculus".
#' @param version A character specifying the annotation version (e.g., "27" or "M16").
#' @return A TxDb object.
#'
#' @import GenomicFeatures
#' @import rtracklayer
#' @export
#' @examples
#' TxDb = getTxDb("Homo sapiens", "27")
#' TxDb = getTxDb("Mus musculus", "M16")

getTxDb <- function(species, version){
	if( !file.exists(paste0("gencode.v", version, ".annotation.sqlite")) ){
		if( species %in% "Homo sapiens" ){
			if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
				download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
					paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)
		}else{
			if( !file.exists(paste0("gencode.v", version, ".annotation.gff3.gz")) )
				download.file(paste0("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_", version, "/gencode.v", version, ".annotation.gff3.gz"),
					paste0("gencode.v", version, ".annotation.gff3.gz"), quiet = TRUE)
		}
		TxDb = makeTxDbFromGFF(file = paste0("gencode.v", version, ".annotation.gff3.gz"),
			format = "gff3",
			dataSource = paste("Gencode version", version),
			organism = species)
		saveDb(TxDb, file=paste0("gencode.v", version, ".annotation.sqlite"))
	}
	return(loadDb(paste0("gencode.v", version, ".annotation.sqlite")))
}
