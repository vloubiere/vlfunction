#' bw coverage
#'
#' Quantify bw file at a set of intervals. returns mean value/region
#'
#' @param bed Regions to plot. Can be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns"
#' @param bw Path to target bw file (character vector)
#' @export

vl_bw_coverage <- function(bed, 
                           bw)
{
  if(class(bed)[1]=="GRanges")
    bed <- data.table::as.data.table()
  if(!is.data.table(bed_DT))
    stop("!is.data.table(bed_DT)")
  if(!all(c("seqnames", "start", "end") %in% colnames(bed_DT)))
    stop("some bed keys are missing in bed colnames!")
  if(length(bw) != 1)
    stop("length(bw) != 1")
  
  # Format
  .b <- data.table::copy(bed)
  .b[, ID:= .I]
  
  # Import bw
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(bed), "score")
  var <- data.table::as.data.table(rtracklayer::import.bw(bw, selection= sel))
  keys <- c("seqnames", "start", "end")
  data.table::setkeyv(.b, keys)
  data.table::setkeyv(var, keys)
  
  # Compute counts
  ov <- data.table::foverlaps(var, .b, nomatch= NULL)
  ov[i.start<start, i.start:= start]
  ov[i.end>end, i.end:= end]
  ov[, width:= i.end-i.start+1]
  res <- data.table::data.table(score= merge(.b[, .(ID)], ov[, .(score= sum(width*score)/sum(width)), ID], by= "ID", all.x= T)$score)
  return(res)
}