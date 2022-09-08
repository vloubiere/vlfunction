#' Import bed file
#'
#' Imports bed as data.table and check formats
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param specialFormat "narrowPeak", "broadPeak"
#' @return Imported bed
#' @export
vl_importBed <- function(bed, ...) UseMethod("vl_importBed")

#' @describeIn vl_importBed for bed paths
#' @export
vl_importBed.character <- function(bed, 
                                   cols= c("seqnames", "start", "end", "name", "score", "strand"), 
                                   extraCols= NULL)
{
  # Fread
  bed <- rbindlist(lapply(bed, fread))
  # Name columns
  if(!is.null(extraCols))
  {
    cols <- if(extraCols=="narrowPeak") 
      c(cols, "signalValue", "pValue", "qValue", "peak") else if(extraCols=="broadPeak") 
        c(cols, "signalValue", "pValue", "qValue") else c(cols, extraCols)
  }
  setnames(bed, cols[1:ncol(bed)])
  if(!is.null(extraCols) && extraCols %in% c("narrowPeak", "broadPeak"))
  {
    bed[, signalValue:= as.numeric(signalValue)]
    bed[, pValue:= as.numeric(pValue)]
    bed[, qValue:= as.numeric(qValue)]
    if(extraCols=="narrowPeak")
      bed[, peak:= as.integer(peak)]
  }
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
    
  return(bed)
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
#' To all closest features
#' vl_closestBed(vl_STARR_DSCP_top_peaks)
#' To find closet yet non touching features
#' vl_closestBed(vl_STARR_DSCP_top_peaks, min_dist = 1)
#' To n closest features
#' vl_closestBed(vl_STARR_DSCP_top_peaks, min_dist = 1, n= 2)
#' 
#' @return Return "a" coor and closeet "b" coordinates together with distance
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
    if(all(c("strand", "strand.b") %in% names(res)))
    {
      res[strand=="+" & start>end.b, dist:= -dist]
      res[strand=="-" & end<start.b, dist:= -dist]
    }
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
  regions <- vl_importBed(bed)
  if(!center %in% c("center", "start", "end"))
    stop("center should be one of center, start or end")
  if(!ignore.strand & !("strand" %in% names(regions)))
    regions[, strand:= "*"]
  if(!ignore.strand && any(!regions$strand %in% c("+", "-")))
    message("ignore.strand=F but strand column absent or contains * ---")
    
  # define start
  if(center=="center")
    regions[, start:= round(rowMeans(.SD)), .SDcols= c("start", "end")] else if(center=="start")
      regions[, start:= ifelse(strand=="-", end, start)] else if(center=="end")
        regions[, start:= ifelse(strand=="-", start, end)]
  # Ext
  regions[, c("start", "end"):= .(start-ifelse(strand=="-", downstream, upstream),
                                  start+ifelse(strand=="-", upstream, downstream))]

  # If genome is specified, resize accordingly
  if(!missing(genome))
  {
    chrSize <- GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))
    chrSize <- GenomicRanges::GRanges(chrSize)
    chrSize <- data.table::as.data.table(chrSize)
    regions[start<1, start:= 1]
    regions[end<1, end:= 1]
    regions[chrSize, start:= ifelse(start>i.end, i.end, start), on= "seqnames"]
    regions[chrSize, end:= ifelse(end>i.end, i.end, end), on= "seqnames"]
  }
  
  # return
  return(regions)
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
  DT <- vl_importBed(bed)
  DT[, init_ord:= .I] 
  setorderv(DT, c("seqnames", "start", "end"))
  
  # Compute contig idx
  DT[, ord:= .I] 
  DT[, ext_end:= end+mingap] 
  idx <- DT$ord
  DT[DT, {idx[.GRP] <<- min(idx[.I])}, .EACHI, on= c("seqnames", "ext_end>=start", "ord<=ord")]
  DT[, idx:= data.table::rleid(idx)]

  # Collapse
  if(!return_idx_only)
    return(DT[, .(start= min(start), end= max(end)), .(seqnames, idx)]) else
      return(DT[order(init_ord), idx])
}

#' Find closestBed regions
#'
#' Return regions "a" overlapping with regions "b"
#'
#' @param a Granges or data.table from which regions overlapping b have to be returned
#' @param b Set of regions of interest
#' @examples 
#' a <- data.table(seqnames= "chr2L", start= c(1000,2000), end= c(2000, 3000))
#' b <- data.table(seqnames= "chr2L", start= 1500, end= 1600)
#' vl_intersectBed(a, b) 
#' @return Return "a" coor and closeet "b" coordinates together with distance
#' @export
vl_intersectBed <- function(a, b)
{
  # Import
  a <- vl_importBed(a)
  b <- vl_importBed(b)
  sel <- b[a, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>0
  a[sel]
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
  bed <- vl_importBed(bed)
  bins <- vl_importBed(bins)
  # Count
  counts <- bed[bins, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
  return(counts)
}
