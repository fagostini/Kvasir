#' Check List function
#' 
#' The function checks that the input argument is a list containing either strings or GRranges. 
#' @param x The object to test.
#' @return The List type.
#'
#' @export
#' @examples
#' List_1 = list("A", "B", "C") # This is correct
#' List_2 = list(2, 1, 4) # This will raise an error
#' checkList(List_1)
#' checkList(List_2)

checkList <- function(x){
  if(!is.list(x)){ stop("The object is not a list()", call.=FALSE); return(1) }
  tmp = unique(sapply(x, class))
  if(length(tmp)>1){ stop("Please provide a either a character or GRanges list", call.=FALSE); return(1) }
  return(tmp)
}