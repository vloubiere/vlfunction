#' reverseComplement
#'
#' Reverse complement DNA sequence character string
#'
#' @param DNA_char DNA character string
#' @examples 
#' vl_revComp("ATCG")
#' @return Reverse Complemented DNA sequence
#' @export

vl_revComp <- function(DNA_char, complement= T, reverse= T)
{
  if(length(DNA_char)>1)
    stop("length DNA_char should be 1!")
  .c <- strsplit(DNA_char, "")[[1]]
  if(!all(.c %in% c("A", "T", "C", "G")))
    stop("All characters should be one of A T C G")
  res <- .c
  if(complement)
  {
    res[.c=="A"] <- "T"
    res[.c=="T"] <- "A"
    res[.c=="C"] <- "G"
    res[.c=="G"] <- "C" 
  }
  if(reverse)
    res <- res[length(res):1]
  return(paste0(res, collapse = ""))
}