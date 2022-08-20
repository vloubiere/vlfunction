#' Check if data.table file contains all bed canonical fields
#' If it does, make sure that class of bed columns is correct
#'
#' @param x Object to be tested
#' @return boolean
#' @export
vl_isDTranges <- function(x)
{
  if(is.data.table(x) && all(c("seqnames", "start", "end") %in% names(x)))
  {
    if(!is.factor(x$seqnames))
      x[, seqnames:= factor(seqnames, 
                            levels= sort(unique(as.character(seqnames))))]
    if(!is.integer(x$start))
      x[, start:= as.integer(start)]
    if(!is.integer(x$end))
      x[, end:= as.integer(end)]
    if("name" %in% names(x) && !is.factor(x$name))
      x[, name:= factor(name)]
    if("score" %in% names(x) && !is.numeric(x$score))
      x[, score:= as.numeric(score)]
    if("strand" %in% names(x) && !is.factor(x$strand))
      x[, strand:= factor(strand, levels = c("+", "-", "*"))]
    return(T)
  }else
    return(F)
}

#' Import bed file
#'
#' Imports bed as data.table and check formats
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @return Imported bed
#' @export
vl_importBed <- function(bed) UseMethod("vl_importBed")

#' @describeIn vl_importBed for bed paths
#' @export
vl_importBed.character <- function(bed)
{
  bed <- lapply(bed, function(x) fread(x, fill = T))
  bed <- data.table::rbindlist(bed)
  if(!vl_isDTranges(bed))
  {
    bedcols <- c("seqnames", "start", "end", "name", "score", "strand")
    if(ncol(bed)<6)
      bedcols <- bedcols[1:ncol(bed)]
    setnames(bed, names(bed)[seq(bedcols)], bedcols)
    vl_isDTranges(bed)
  }
  return(bed)
}

#' @describeIn vl_importBed for GRanges
#' @export
vl_importBed.GRanges <- function(bed)
{
  bed <- data.table::as.data.table(bed)
  vl_isDTranges(bed)
  return(bed)
}

#' Export bed file
#' @param bed Either a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @export
vl_exportBed <- function(bed, ...) UseMethod("vl_exportBed")

#' @describeIn vl_exportBed for GRanges
#' @export
vl_exportBed.GRanges <- function(bed, filename)
{
  fwrite(data.table::as.data.table(bed), 
         filename,
         sep= "\t", 
         col.names = F,
         quote= F,
         scipen = 20)
} 

#' @describeIn vl_exportBed for data.table
#' @export
vl_exportBed.data.table <- function(bed, filename)
{
  if(!vl_isDTranges(bed))
    stop("Could not find seqnames, start and end columns")
  cols <- c("seqnames", "start", "end", "name", "score", "strand")
  if(!"name" %in% names(bed))
    bed$name <- "."
  if(!"score" %in% names(bed))
    bed$score <- 0
  if(!"strand" %in% names(bed))
    bed$strand <- "."
  cols <- cols[cols %in% names(bed)]
  setcolorder(bed, cols)
  fwrite(bed, 
         filename,
         sep= "\t", 
         col.names = F,
         quote= F,
         scipen = 20)
} 

#' Find closestBed regions
#'
#' For each line of a file, returns the closest lines in b
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself.
#' @param n Number of closest features to be reported. Default to 1
#' @param min_dist Min distance for closest feature 
#' @examples 
#' a <- data.table(seqnames= "chr2L", start= sample(10000, 1000))
#' a[, end:= start+1000]
#' 
#' To all closest features
#' vl_closestBed(a, min_dist= 0)
#' 
#' To find closet yet non touching features
#' vl_closestBed(a, min_dist= 1)
#' 
#' @return Return "a" coor and closeet "b" coordinates together with distance
#' @export
vl_closestBed <- function(a, 
                          b= NULL,
                          n= 1,
                          min_dist= 0)
{
  if(!vl_isDTranges(a))
    a <- vl_importBed(a)
  if(is.null(b))
    b <- a else if(!vl_isDTranges(b))
      b <- vl_importBed(b)
  a <- data.table::copy(a)
  b <- data.table::copy(b)
  # Main function
  res <- b[a, {
    # Measure distance
    .c <- data.table(.SD, 
                     dist= fcase(x.start>i.end, x.start-i.end,
                                 x.end<i.start, x.start-i.end,
                                 default= 0L))# default means a & b overlap!
    # Order dist and apply dist & n cutoffs
    .c <- .c[abs(dist)>=min_dist]
    .c <- .c[data.table::first(order(abs(dist)), n)]
    setnames(.c, 
             c("start", "end"), 
             c("start.b", "end.b"))
    # Return a coordinates + closest b coor
    data.table(start= i.start, end= i.end, .c)
  }, .EACHI, on= "seqnames"]
  # Export
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
#' bed <- data.table(seqnames= "chr2L", start= 10000, end= 12000, strand= c("+", "-"))
#' vl_resizeBed(bed, center= "start", upstream = 2000, downstream = 1000)[]
#' @return Resized DT ranges
#' @export
vl_resizeBed <- function(bed, 
                         center= "center",
                         upstream= 500,
                         downstream= 500, 
                         ignore.strand= F,
                         genome)
{
  if(!ignore.strand && !"strand" %in% names(bed))
    message("ignore.strand=F but no strand column is detected -> ignored!")
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  regions <- data.table::copy(bed)
  
  if(!center %in% c("center", "start", "end"))
    stop("center should be one of center, start or end")
  
  # Check strand
  if(!ignore.strand & !("strand" %in% names(regions)))
    regions[, strand:= factor("*")]
    
  # define start
  if(center=="center")
  {
    regions[, start:= round(rowMeans(.SD)), .SDcols= c("start", "end")]
  }else if(center=="start")
  {
    regions[, start:= ifelse(strand=="-", end, start)]
  }else if(center=="end")
    regions[, start:= ifelse(strand=="-", start, end)]

  # Ext
  regions[, c("start", "end"):= .(start-ifelse(strand=="-", downstream, upstream),
                                  start+ifelse(strand=="-", upstream, downstream))]

  # If genome is specified, resize accordingly
  if(!missing(genome))
  {
    chrSize <- GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome)))
    chrSize <- data.table::as.data.table(chrSize)
    regions[chrSize, seqlength:= i.width, on= "seqnames"]
    regions[start<1, start:= 1]
    regions[end>seqlength, end:= seqlength]
  }
  
  # Final
  res <- data.table::copy(bed)
  res[, c("seqnames", "start", "end"):= regions[, .(seqnames, start, end)]]
  return(res)
}

#' Collapse ranges with a certain min dist
#'
#' Collapse ranges with a certain min dist between them
#'
#' @param bed bed file to be collapsed. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames','start', 'end' columns. see ?vl_importBed()
#' @param mingap min gap ditance for merging. default is 1 (that is, touching coordinates)
#' @param return_idx_only If set to T, does not collapse regions but returns idx as an extra columns. default= F
#' @examples 
#' bed <- data.table(seqnames= c("chr3L", "chr3L", "chr3L", "chr3L", "chr3L", "chr2R", "chr2R", "chr2R"),
#' start= c(1000, 2000, 2100, 2200, 2500, 2000, 2500, 5000),
#' end= c(1500, 2099, 2199, 2299, 3000, 3000, 3500, 6000))
#' vl_collapseBed(bed, mingap= 1)
#' @return Collapse coor data.table
#' @export
vl_collapseBed <- function(bed,
                           mingap= 1,
                           return_idx_only= F)
{
  # Hard copy of bed file
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  DT <- copy(bed)
  setkeyv(DT, c("seqnames", "start", "end"))
  
  # Compute contig idx
  DT[, idx:= cumsum(start>data.table::shift(end, fill= max(end))+mingap), seqnames]
  DT[, idx:= .GRP, .(seqnames, idx)]
  
  # Collapse
  if(!return_idx_only)
    DT <- DT[, .(start= min(start), end= max(end)), .(seqnames, idx)]

  return(DT)
}

#' Compute bed coverage
#'
#' For each bin, computes the number of overlapping reads from a bed file
#' @param bins bins for which enrichment has to be computed. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param bed bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @return numeric vector of the number of overlapping reads
#' @export
vl_covBed <- function(bins,
                      bed)
{
  # Import reads
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  if(!vl_isDTranges(bins))
    bins <- vl_importBed(bins)
  counts <- bed[bins, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
  
  return(counts)
}