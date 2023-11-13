#' reverseComplement
#'
#' Reverse complement DNA sequence character string
#'
#' @param DNA_char DNA character string
#' 
#' @examples 
#' vl_revComp("ATCG")
#' 
#' @return Reverse Complemented DNA sequence
#' @export
vl_revComp <- function(DNA_char,
                       complement= T,
                       reverse= T)
{
  if(length(DNA_char)>1)
    stop("length DNA_char should be 1!")
  .c <- strsplit(DNA_char, "")[[1]]
  if(!all(.c %in% c("A", "T", "C", "G")))
    stop("All characters should be one of A T C G")
  res <- sapply(.c, function(x)
  {
    switch(x, 
           "A"= "T",
           "T"= "A",
           "C"= "G",
           "G"= "C") 
  })
  if(reverse)
    res <- rev(res)
  return(paste0(res, collapse = ""))
}

#' vOligoCalculator
#'
#' Calculates Oligo Tm and GC %
#'
#' @param seq Character vector containing DNA sequence
#' 
#' @examples 
#' vl_oligo_Tm("GCCGATTCTCGGGCAGTTCCTC")
#' 
#' @return Tm and GC% as a list
#' @export
vl_oligo_Tm <- function(seq)
{
  if(length(seq)>1)
    stop("seq should be length 1!")
  current <- data.table(c("A", "T", "C", "G"), N= 0, key= "V1")
  current[data.table(unlist(strsplit(seq, "")))[, .N, V1], N:= i.N, on= "V1"]
  
  Tm <- 2*(sum(current[c("A", "T"), N]))+4*(sum(current[c("G", "C"), N]))-7
  GC <- round(sum(current[c("G", "C"), N])/sum(current$N)*100, 1)
  return(list(Tm= Tm,
              'GC%'= GC))
}

#' Digest sequence
#'
#' Simulate DNA digest from thermofisher "simple" enzymes
#'
#' @param seq character sequence to digest
#' @param enzyme Character vector containing the enzymes to use
#' @param keepsite If TRUE, then the enzymmatic consensu will be paste at cutting sites
#' 
#' @examples 
#' vl_digest(seq= "AAAAAAAAGGTACCTTTTTTTTTTTTGCGGCCGCAAAAAAAAAA",
#' enzyme= c("KpnI", "NotI"))
#' 
#' @return A vector containing the digested pieces
#' @export
vl_digest <- function(seq,
                      enzyme,
                      keepsite= F)
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