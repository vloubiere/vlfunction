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
  dat <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  bins <- dat[, .(start= seq(1, end, bin_size)), .(seqnames, end, width)]
  bins[, end:= start+bin_size-1]
  bins[end>width, end:= width]
  bins <- bins[end-start>0, .(seqnames, start, end)]
  return(bins)
}

#' Generate random control regions
#'
#' Generate control regions with given width distribution
#'
#' @param BSgenome BSgneome object to use.
#' @param n Number of regions to sample
#' @param width Widths of regions to sample (length should either be 1 or equal n)
#' @param restrict_seqnames If specified, only the providsed seqnames will be used
#' @param n Number of regions to sample
#' 
#' @examples 
#' vl_control_regions_BSgenome(BSgenome = BSgenome.Dmelanogaster.UCSC.dm3, 
#' n= 5000, 
#' width= 1000, 
#' restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"))
#' 
#' @return data.table containing random control regions
#' @export


vl_control_regions_BSgenome <- function(BSgenome,
                                        n,
                                        width= 1, 
                                        restrict_seqnames= NULL)
{
  if(!class(BSgenome)=="BSgenome")
    stop("genome should be a BSgenome object!")
  if(width<1)
    stop("width cannot be smaller than 1")
  dat <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  if(is.null(restrict_seqnames))
    restrict_seqnames <- dat$seqnames
  
  dat <- dat[seqnames %in% restrict_seqnames]
  colnames(dat)[4] <- "chr_size"
  idx <- sample(nrow(dat), 
                replace = T, 
                size = n, 
                prob = dat$chr_size)
  rdm <- dat[idx, .(seqnames, chr_size)]
  rdm[, width:= width]
  rdm[, start:= sample(seq(width+1, chr_size-width-1), .N), .(chr_size, width)]
  rdm[, end:= start+width-1, .(chr_size, width)]
  
  return(rdm[, .(seqnames, start, end)])
}

#' Find closestBed idx
#'
#' For each line of a file, return the idx of the closest line in b (or non orverlapping line in a if b is null)
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself while excluding self-matching.
#' @param min_dist Min distance for closest feature 
#' @examples 
#' Example 1:
#' a <- b <- data.table(seqnames= c("chr4", "chr2R", "chr2R", "chr3L"),
#' start= c(1e6, 1e6, 10e6, 1e6),
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

vl_closestBed <- function(a, 
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
  res <- DT[, .(start= min(start), end= max(end)), .(seqnames, idx, strand)]
  return(res[, .(seqnames, start, end, strand)])
}