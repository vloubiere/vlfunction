#' vOligoCalculator
#'
#' Calculates Oligo Tm and GC %
#'
#' @param seq Character vector containing DNA sequence
#' @examples 
#' vl_oligo_Tm("GCCGATTCTCGGGCAGTTCCTC")
#' @return pritns result
#' @export

vl_oligo_Tm <- function(seq)
{
  if(length(seq)>1)
    stop("seq should be length 1!")
  current <- data.table(c("A", "T", "C", "G"), N= 0, key= "V1")
  current[data.table(unlist(strsplit(seq, "")))[, .N, V1], N:= i.N, on= "V1"]

  Tm <- 2*(sum(current[c("A", "T"), N]))+4*(sum(current[c("G", "C"), N]))-7
  GC <- round(sum(current[c("G", "C"), N])/sum(current$N)*100, 1)
  return(list(Tm= Tm, 'GC%'= GC))
}