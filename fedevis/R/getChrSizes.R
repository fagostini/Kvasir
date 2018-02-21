#' Create SeqInfo object function
#' 
#' The function retrieves the chromosome size information for the genome in input and create the Seqinfo object. 
#' @param genome A character indicating a genome reference (e.g, "hg19", "mm10").
#' @return A data.table object.
#'
#' @import GenomeInfoDb
#' @import data.table
#' @export
#' @examples
#' SeqInfo = getChrSizes("mm10")

getChrSizes <- function(genome){
	info = as.data.table(fetchExtendedChromInfoFromUCSC(genome))
	info = info[, list(UCSC_seqlevel, UCSC_seqlength, circular)]
	info = with(info, Seqinfo(seqnames=UCSC_seqlevel, seqlengths=UCSC_seqlength, isCircular=circular, genome=genome))
	return(info)
}
