#' Find closestBed idx
#'
#' For each line of a file, return the idx of the closest line in b (or non orverlapping line in a if b is null)
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself while excluding self-matching.
#' @param min_dist Min distance for closest feature 
#' @examples 
#' Example 1:
#' a <- b <- data.table(seqnames= c("chr4", "chr2R", "chr2R", "chr3L")
#' start= c(1e6, 1e6, 10e6, 1e6)
#' end= c(1.1e6, 1.1e6, 10.1e6, 1.1e6))
#' b <- a[-1]
#' b[, start:= start+5e5]
#' b[, end:= end+5e5]
#' closestBed(a, b)
#' closestBed(a, b, min_dist = 1e6)
#' closestBed(a, a)
#' closestBed(a, a, min_dist= 1)
#' closestBed(a)
#' 
#' Benchmarking:
#' require(rtracklayer)
#' tss <- as.data.table(import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf"))
#' tss <- tss[type=="gene", .(seqnames, start= ifelse(strand=="+", start+0, end-1)), gene_name]
#' tss[, end:= start+1]
#' This:
#' c1 <- closestBed(tss)
#' Is pretty similar to this:
#' c2 <- closestBed(tss, tss, 1)
#' Except for the cases where two genes have exact same start. Example check c1[67] and c2[67], and compare to tss[65:67]
#' closestBed(tss, tss)
#' closestBed(tss, min_dist = 100000)
#' 
#' @return For each line in a, returns the idx of its closest feature in b (vector)
#' @export

closestBed <- function(a, 
                       b= NULL, 
                       min_dist= 0)
{
  if(!is.data.table(a))
    a <- as.data.table(a)
  a[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
  a[, idx:= .I]
  if(is.null(b))
  {
    b <- a
    check_idx <- T
  }else
  {
    if(!is.data.table(b))
      b <- as.data.table(b)
    b[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
    b[, idx:= .I]
    check_idx <- F
  }
  # Compute
  res <- b[a, {
    dist <- abs(center-i.center)
    dist[dist<min_dist] <- NA
    if(check_idx)
      dist[idx==i.idx] <- NA
    if(all(is.na(dist)))
      as.integer(NA) else
        idx[which.min(dist)]
  }, .EACHI, on= "seqnames"]$V1
  return(res)
}


