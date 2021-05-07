#' Collapse ranges with a certain min dist
#'
#' Collapse ranges with a certain min dist between them
#'
#' @param bed GRanges or data.table file with "seqnames", "start", "end" (and optionally strand, set to * if absent)
#' @param mingap min gap ditance for merging. default is 1 (that is, touching coordinates)
#' @param stranded If set to true and strand column is provided, only merges coordinates with similar strand.
#' @examples 
#' ex <- GRanges(c("chr3L", "chr3L", "chr3L", "chr3L", "chr2R"), 
#' IRanges(c(1000, 2000, 2100, 2200, 2000), 
#' c(1500, 2099, 2199, 2299, 3000)),
#' strand= c("+","+","+","-","+"))
#' vl_collapse_DT_ranges(ex, mingap= 1, stranded= F)
#' vl_collapse_DT_ranges(ex, mingap= 1, stranded= T)
#' 
#' @return Collapse coor data.table
#' @export

vl_collapse_DT_ranges <- function(bed, 
                                  mingap= 1, 
                                  stranded= F)
{
  if(class(bed)[1]=="GRanges")
    bed <- data.table::as.data.table(bed)
  if(!data.table::is.data.table(bed) | !all(c("seqnames", "start", "end") %in% colnames(bed)))
    stop("bed must be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns")
  if(!"strand" %in% colnames(bed) | !stranded)
    bed[, strand:= "*"]
    
  DT <- copy(bed)
  setorderv(DT, c("seqnames", "start", "end", "strand"))
  i <- 0
  DT[, idx:= {
    i <<- i+1
    .idx <- c(i, sapply(.SD[-1, start]-.SD[-nrow(.SD), end], function(y) {
      if(y>mingap) 
        i <<- i+1
      return(i)
    }))
  }, .(seqnames, strand)]
  res <- DT[, .(min(start), max(end)), .(seqnames, idx, strand)]
  return(res)
}