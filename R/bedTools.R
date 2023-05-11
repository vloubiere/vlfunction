#' Convert character to bed file
#'
#' @param coor Character vector of genomic coordinates, e.g chr3RHet:2281049-2281297:*
#' @examples
#' vl_toDTranges("chr3RHet:2281049-2281297:*")
#' @return DT ranges
#' @export
vl_toDTranges <- function(coor)
{
  if(!is.character(coor))
    stop("coor should be a character vector e.g chr3RHet:2281049-2281297:*")
  dat <- data.table::as.data.table(data.table::tstrsplit(coor, ":|-"))
  setnames(dat, c("seqnames", "start", "end", "strand")[seq(dat)])
  dat <- vl_importBed(dat)
  return(dat)
}

#' Import bed file
#'
#' Imports bed as data.table and check formats
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param extraCols Colnames for extra, non-canonical bed columns. "auto" tries to guess whether narrowPeak/broadPeak formats are used
#' @return DT ranges
#' @export
vl_importBed <- function(bed, ...) UseMethod("vl_importBed")

#' @describeIn vl_importBed for bed paths
#' @export
vl_importBed.character <- function(bed, 
                                   cols= c("seqnames", "start", "end", "name", "score", "strand"),
                                   extraCols= "auto")
{
  # Guess format if auto
  if(extraCols=="auto")
    extraCols <- if(grepl(".narrowPeak$", bed[1]))
      c("signalValue", "pValue", "qValue", "peak") else if(grepl(".broadPeak$", bed[1]))
        c("signalValue", "pValue", "qValue") else
          NULL
  # Columns
  if(!is.null(extraCols))
    cols <- c(cols, extraCols)
  # Fread
  bed <- rbindlist(lapply(bed, fread))
  # Name columns
  setnames(bed, cols[1:ncol(bed)])
  vl_importBed(bed)
}

#' @describeIn vl_importBed for GRanges
#' @export
vl_importBed.GRanges <- function(bed)
  vl_importBed(data.table::as.data.table(bed))

#' Default uses data.table format and does sanity checks
#' @export
vl_importBed.default <- function(bed)
{
  bed <- data.table::copy(bed)
  if("seqnames" %in% names(bed))
    bed[, seqnames:= as.character(seqnames)]
  if("start" %in% names(bed))
    bed[, start:= as.integer(start)]
  if("end" %in% names(bed))
    bed[, end:= as.integer(end)]
  if("name" %in% names(bed))
    bed[, name:= as.character(name)]
  if("score" %in% names(bed))
    bed[, score:= as.numeric(score)]
  if("strand" %in% names(bed))
  {
    bed[, strand:= as.character(strand)]
    bed[, strand:= ifelse(strand %in% c("+", "-"), strand, "*")]
  }
  
  if(any(bed[, start>end]))
    warning("bed file contains ranges with start>end -> malformed!")
   
  return(bed)
}

#' Find closestBed regions
#'
#' For each line of a file, returns the closest lines in b.
#' Similar to bedtools closest -a a.bed -b b.bed -D a
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself.
#' @param n Number of closest features to be reported. Default to 1
#' @param min_dist Min distance for closest feature 
#' @examples 
#' To all closest features
#' a <- data.table(seqnames= c("chr2R"), start= c(10000, 20000, 30000), end= c(10000, 20000, 30000), strand= c("+", "-", "*"))
#' b <- data.table(seqnames= c("chr2R"), start= c(11000, 21000, 31000), end= c(11000, 21000, 31000), strand= c("-", "*", "+"))
#' vl_closestBed(a,b)[]
#' To find closet yet non touching features
#' vl_closestBed(vl_STARR_DSCP_top_peaks, min_dist = 1)
#' To n closest features
#' vl_closestBed(vl_STARR_DSCP_top_peaks, min_dist = 1, n= 2)
#' 
#' @return Return "a" coor and closest "b" coordinates and distance in respect to A (When a is on the - strand, “upstream” means b has a higher (start,stop).)
#' @export
vl_closestBed <- function(a, 
                          b= NULL,
                          n= 1,
                          min_dist= 0)
{
  # Import
  a <- vl_importBed(a)
  if(is.null(b))
    b <- data.table::copy(a) else
      b <- vl_importBed(b)
  
  # Checks
  if(!"strand" %in% names(a))
    message("'a' does not contain strand column -> arbitrarily considered as +")else if("*" %in% a$strand)
      message("'a' strand column contains * -> arbitrarily considered as +!")

  # Closest
  idx <- b[a, {
    dist <- fcase(x.start>i.end, as.integer(x.start-i.end),
                  x.end<i.start, as.integer(i.start-x.end), 
                  default= 0L)
    sel <- between(dist, min_dist, unique(sort(dist[dist>=min_dist]))[n])
    .(.GRP, I= .I[sel], dist= dist[sel])
  }, .EACHI, on= "seqnames"]
  idx <- na.omit(idx)
  idx[I==0, I:= NA]
  
  # Return
  setnames(b, paste0(names(b), ".b"))
  res <- data.table(a[idx$GRP],
                    b[idx$I],
                    dist= idx$dist)
  if("strand" %in% names(res))
  {
    res[strand!="-" & start>end.b, dist:= -dist] # Meaning + or *
    res[strand=="-" & end<start.b, dist:= -dist]
  }else
    res[start>end.b, dist:= -dist] # In the case were no strand is provided for a, all considered as +
  
  return(res)
}

#' Resize bed
#'
#' Resize a bed file starting from specified origin
#'
#' @param bed Bed file to resize. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param center From where should the region be centered before extension. Either "center", "start" or "end"
#' @param upstream Upstream extension. default= 500L
#' @param downstream Downstream extension. default= 500L
#' @param ignore.strand Should the strand be considered when defininng start or end centering? Default= F
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accodtingly
#' @examples 
#' bed <- data.table(seqnames= "chr2L", start= 10000, end= 12000, strand= c("+", "-", "*"))
#' vl_resizeBed(bed, center= "start", upstream = 2000, downstream = 1000)[]
#' vl_resizeBed(bed, center= "end", upstream = 2000, downstream = 1000)[]
#' vl_resizeBed(bed, center= "center", upstream = 2000, downstream = 1000)[]
#' @return Resized DT ranges
#' @export
vl_resizeBed <- function(bed, 
                         center= "center",
                         upstream= 500,
                         downstream= 500, 
                         ignore.strand= F,
                         genome)
{
  # Checks
  if(!center %in% c("center", "start", "end"))
    stop("center should be one of center, start or end")
  bed <- vl_importBed(bed)
  if(!ignore.strand)
  {
    if(!"strand" %in% names(bed))
      message("'bed' does not contain strand column -> arbitrarily considered as +")else if("*" %in% bed$strand)
        message("'bed' strand column contains * -> arbitrarily considered as +!")
  }
  
  # define start
  if(center=="center") # Center does not depend on the strand
  {
    bed[, start:= round(rowMeans(.SD)), .SDcols= c("start", "end")]
  }else if(!ignore.strand && "strand" %in% names(bed)) # Meaning strand should and can be considered
  {
    if(center=="start")
      bed[strand=="-", start:= end]else if(center=="end")
        bed[strand!="-", start:= end]
  }else if(center=="end") # Note that if center=="start" and strand is not considered, start kept unchanged
    bed[, start:= end]

  # Ext
  if(!ignore.strand && "strand" %in% names(bed))
  {
    bed[, c("start", "end"):= .(start-ifelse(strand=="-", downstream, upstream),
                                start+ifelse(strand=="-", upstream, downstream))]
  }else
    bed[, c("start", "end"):= .(start-upstream, start+downstream)]

  # If genome is specified, resize accordingly
  if(!missing(genome))
  {
    chrSize <- GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))
    chrSize <- GenomicRanges::GRanges(chrSize)
    chrSize <- data.table::as.data.table(chrSize)
    bed[start<1, start:= 1]
    bed[end<1, end:= 1]
    bed[chrSize, start:= ifelse(start>i.end, i.end, start), on= "seqnames"]
    bed[chrSize, end:= ifelse(end>i.end, i.end, end), on= "seqnames"]
  }
  
  # return
  return(bed)
}

#' Collapse ranges with a certain min dist
#'
#' Collapse ranges with a certain min dist between them
#'
#' @param bed bed file to be collapsed. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames','start', 'end' columns. see ?vl_importBed()
#' @param mingap min gap ditance for merging. default is 1 (that is, touching coordinates)
#' @param return_idx_only If set to T, does not collapse regions but returns idx as an extra columns. default= F
#' @param ignore.strand Should strand be ignored? default= T
#' @examples 
#' gr <- GRanges(seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
#' ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
#' strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)))
#' bed <- vl_importBed(gr)
#' vl_collapseBed(bed)
#' vl_collapseBed(bed, return_idx_only = T)
#' vl_collapseBed(bed, ignore.strand = F)
#' vl_collapseBed(bed, ignore.strand = F, return_idx_only = T)
#' 
#' # Mingap
#' bed[7, c("start", "end"):= .(start+10, end+10)]
#' vl_collapseBed(bed)
#' vl_collapseBed(bed, mingap = 10)
#' vl_collapseBed(bed, return_idx_only = T)
#' 
#' @return Collapse coor data.table
#' @export
vl_collapseBed <- function(bed,
                           mingap= 1,
                           return_idx_only= F,
                           ignore.strand= T)
{
  # Hard copy of bed file
  DT <- vl_importBed(bed)
  DT[, init_ord:= .I]
  setorderv(DT, c("seqnames", "start", "end"))
  
  # Compute contig idx
  DT[, ord:= .I] 
  DT[, ext_end:= end+mingap] 
  idx <- DT$ord
  if(!ignore.strand && "strand" %in% names(DT))
    DT[DT, {idx[.GRP] <<- min(idx[.I])}, .EACHI, on= c("seqnames", "ext_end>=start", "ord<=ord", "strand")] else
      DT[DT, {idx[.GRP] <<- min(idx[.I])}, .EACHI, on= c("seqnames", "ext_end>=start", "ord<=ord")]
  DT[, idx:= data.table::rleid(idx)]

  # Collapse
  if(!return_idx_only)
  {
    if(!ignore.strand && "strand" %in% names(DT))
      res <- DT[, .(start= min(start), end= max(end)), .(seqnames, strand, idx)] else
        res <- DT[, .(start= min(start), end= max(end)), .(seqnames, idx)]
    setcolorder(res, c("seqnames", "start", "end"))
    setorderv(res, c("seqnames", "start", "end"))
  }else
    res <- DT[order(init_ord), idx]
  return(res)
}

#' Find closestBed regions
#'
#' Return regions "a" overlapping with regions "b"
#'
#' @param a Granges or data.table from which regions overlapping b have to be returned
#' @param b Set of regions of interest
#' @param ignore.strand Should strand be ignored? default= T
#' @examples 
#' a <- data.table(seqnames= "chr2L", start= c(1000,2000), end= c(2000, 3000))
#' b <- data.table(seqnames= "chr2L", start= 1500, end= 1600)
#' vl_intersectBed(a, b) 
#' 
#' @return Return 'a' ranges that overlap with any range in 'b'
#' @export
vl_intersectBed <- function(a, 
                            b, 
                            ignore.strand= T)
{
  # Import
  a <- vl_importBed(a)
  b <- vl_importBed(b)
  if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
    sel <- b[a, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start", "strand")]$N>0 else
      sel <- b[a, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>0
  a[sel]
}

#' Compute bed coverage
#'
#' For each bin, computes the number of overlapping reads from a bed file
#' @param a Ranges for which overlaps with b should be counted. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param b Ranges from which overlaps should be computed.
#' @param ignore.strand Should strand be ignored? default= T, meaning all overlapping elements in 'b'  will be considered
#' @examples
#' a <- data.table(seqnames= c("chr2R"), start= 10000, end= 20000, strand= c("+", "-"))
#' b <- data.table(seqnames= c("chr2R"), start= seq(10000, 19000, 1000), end= seq(11000, 20000, 1000), strand= rep(c("+", "-"), each= 10))
#' vl_covBed(a, b)
#' vl_covBed(a, b, ignore.strand = F)
#' 
#' @return For each range in 'a', reports the number of overlapping features in 'b'
#' @export
vl_covBed <- function(a,
                      b, 
                      ignore.strand= T)
{
  # Import reads
  a <- vl_importBed(a)
  b <- vl_importBed(b)
  # Count
  if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
    counts <- b[a, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start", "strand")]$N else
      counts <- b[a, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
  return(counts)
}
