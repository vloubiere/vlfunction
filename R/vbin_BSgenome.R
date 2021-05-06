#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param BSgenome BSgneome object to use.
#' @param bin_size bin size. default to 50bp
#' @examples 
#' vl_binBSgenome(BSgenome= BSgenome.Dmelanogaster.UCSC.dm3, bin_size= 50)
#' @return data.table containing bin coordinates
#' @export


vl_binBSgenome <- function(BSgenome,
                           bin_size= 50)
{
  if(!class(BSgenome)=="BSgenome")
    stop("genome should be a BSgenome object!")
  dat <- as.data.table(as.data.frame(GenomeInfoDb::seqinfo(BSgenome)), keep.rownames = "seqnames")
  bins <- dat[, .(start= seq(1, seqlengths, bin_size)), .(seqnames, seqlengths)]
  bins[, end:= start+bin_size-1]
  bins[end>seqlengths, end:= seqlengths]
  bins <- bins[end-start>0, !"seqlengths"]
  return(bins)
}
