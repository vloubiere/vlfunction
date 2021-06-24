#' Digest sequence
#'
#' Simulate DNA digest from thermofisher "simple" enzymes
#'
#' @param seq character sequence to digest
#' @param enzyme Character vector containing the enzymes to use
#' @param keepsite If TRUE, then the enzymmatic consensu will be paste at cutting sites
#' @examples 
#' vl_digest(seq= "AAAAAAAAGGTACCTTTTTTTTTTTTGCGGCCGCAAAAAAAAAA",
#' enzyme= c("KpnI", "NotI"))
#' @return A vector containing the digested pieces
#' @export

vl_digest <- function(seq, enzyme, keepsite= F)
{
  seq <- toupper(seq)
  .c <- na.omit(vl_thermofisher_restriction_enzymes_table[enzyme, , on= "name1"])
  if(nrow(.c)!=length(enzyme))
    stop("Some enzyme(s) could not be found in table")
  if(!all(unlist(strsplit(.c$cutsite, "")) %in% c("A","T","C","G","^")))
    stop("Complex enzymes not handled yet")
  res <- seq
  if(keepsite)
    .c[, {res <<- gsub(consensus_F, paste0(consensus_F, "^", consensus_F), res)}, (.c)]
  else
    .c[, {res <<- gsub(consensus_F, cutsite, res)}, (.c)]
  res <- unlist(strsplit(res, "^", fixed= T))
  return(res)
}