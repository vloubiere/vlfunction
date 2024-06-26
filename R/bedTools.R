#' Convert character to bed file
#'
#' @param coor Character vector of genomic coordinates, e.g "chr3R:2281049-2281297:*"
#' @examples
#' vl_toDTranges("chr3R:2281049-2281297:*")
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
#' @examples
#' bed <- vl_importBed("/path/to/file.bed")
#' test <- data.table(seqnames= "chr3R",
#'                    start= 100000,
#'                    end= 110000)
#' bed <- vl_importBed(test)
#' 
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
  setnames(bed,
           cols[1:ncol(bed)])
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
    bed[, seqnames:= as.character(seqnames)] else
      stop("No seqnames column in bed file")
  if("start" %in% names(bed))
    bed[, start:= as.integer(start)] else
      stop("No start column in bed file")
  if("end" %in% names(bed))
    bed[, end:= as.integer(end)]else
      warning("No end column in bed file")
  if("name" %in% names(bed))
    bed[, name:= as.character(name)]
  if("score" %in% names(bed))
    bed[, score:= as.numeric(score)]
  if("strand" %in% names(bed))
  {
    bed[, strand:= as.character(strand)]
    bed[, strand:= ifelse(strand %in% c("+", "-"), strand, "*")]
  }
  
  if(any(bed[, start>end], na.rm = T))
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
#' @param min.dist Min distance for closest feature 
#' @examples 
#' a <- data.table(seqnames= "chr2R",
#'                 start= c(10000, 20000, 30000),
#'                 end= c(10000, 20000, 30000),
#'                 strand= c("+", "-", "*"))
#' b <- data.table(seqnames= "chr2R",
#'                 start= c(11000, 21000, 31000),
#'                 end= c(11000, 21000, 31000),
#'                 strand= c("-", "*", "+"))
#' 
#' # All closest features
#' vl_closestBed(a, b)[]
#' 
#' # To find closet yet non touching features
#' vl_closestBed(vl_STARR_DSCP_top_peaks,
#'               min.dist = 1)
#' # To n closest features
#' vl_closestBed(vl_STARR_DSCP_top_peaks,
#'               min.dist = 1,
#'               n= 2)
#' 
#' @return Return "a" coor and closest "b" coordinates and distance in respect to A (When a is on the - strand, “upstream” means b has a higher (start,stop).)
#' @export
vl_closestBed <- function(a, 
                          b= NULL,
                          n= 1,
                          min.dist= 0)
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
    sel <- between(dist, min.dist, unique(sort(dist[dist>=min.dist]))[n])
    .(.GRP, I= .I[sel], dist= dist[sel])
  }, .EACHI, on= "seqnames"]
  # idx <- na.omit(idx)
  idx[I==0, c("I", "dist"):= .(NA, NA)]
  
  # Return
  setnames(b, paste0(names(b), ".b"))
  res <- data.table(a[eval(idx$GRP, with= F)],
                    b[eval(idx$I, with= F)],
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
#' @param center From where should the region be centered before extension. Either "center", "start" or "end". Default= "center".
#' @param upstream Upstream extension. default= 500L
#' @param downstream Downstream extension. default= 500L
#' @param ignore.strand Should the strand be considered when defininng start or end centering? Default= F
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly
#' @examples 
#' bed <- data.table(seqnames= "chr2L",
#'                   start= 10000,
#'                   end= 12000,
#'                   strand= c("+", "-", "*"))
#' vl_resizeBed(bed, center= "start", upstream = 0, downstream = 0)[]
#' vl_resizeBed(bed, center= "end", upstream = 0, downstream = 0)[]
#' vl_resizeBed(bed, center= "center", upstream = 0, downstream = 0)[]
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
#' @param min.gap min gap ditance for merging. default is 1 (that is, touching coordinates)
#' @param return.idx.only If set to T, does not collapse regions but returns idx as an extra columns. default= F
#' @param ignore.strand Should strand be ignored? default= T
#' @examples 
#' bed <- data.table(seqnames= "chr3R",
#'                   start= c(10e3, 11e3, 12e3, 20e3, 21e3),
#'                   end= c(12e3, 13e3, 14e3, 22e3, 23e3),
#'                   strand= c("+", "+", "-", "+", "+"))
#'                   
#' vl_collapseBed(bed)
#' # Return mergin index only
#' vl_collapseBed(bed, return.idx.only = T)
#' # Only merge if strand is similar
#' vl_collapseBed(bed, ignore.strand = F)
#' # Allow a certain gap for merging
#' vl_collapseBed(bed, ming.gap = 1000)
#' vl_collapseBed(bed, ming.gap = 10000)
#' 
#' @return Collapse bed data.table
#' @export
vl_collapseBed <- function(bed,
                           min.gap= 1,
                           return.idx.only= F,
                           ignore.strand= T)
{
  # Hard copy of bed file
  DT <- vl_importBed(bed)
  DT <- DT[, names(DT) %in% c("seqnames", "start", "end", "strand"), with= F]
  DT[, init_ord:= .I]
  setorderv(DT, c("seqnames", "start", "end"))
  
  # Compute contig idx
  DT[, ord:= .I] 
  DT[, ext_end:= end+min.gap] 
  idx <- DT$ord
  if(!ignore.strand && "strand" %in% names(DT))
    DT[DT, {idx[.GRP] <<- min(idx[.I])}, .EACHI, on= c("seqnames", "ext_end>=start", "ord<=ord", "strand")] else
      DT[DT, {idx[.GRP] <<- min(idx[.I])}, .EACHI, on= c("seqnames", "ext_end>=start", "ord<=ord")]
  DT[, idx:= data.table::rleid(idx)]

  # Collapse
  if(!return.idx.only)
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
#' @param min.overlap A vector integer of length 1 or nrow(a) specifying the minimum overlap (bp) require to count overlaps. Default= 1L
#' @param invert If set to true, returns non-overlapping features. Default= F
#' @examples 
#' a <- data.table(seqnames= "chr2L",
#'                 start= c(1000,2000),
#'                 end= c(2000, 3000))
#' b <- data.table(seqnames= "chr2L",
#'                 start= 1500,
#'                 end= 1600)
#' vl_intersectBed(a, b)
#' 
#' @return Return 'a' ranges that overlap with any range in 'b'
#' @export
vl_intersectBed <- function(a, 
                            b,
                            ignore.strand= T,
                            min.overlap= 1L,
                            invert= F)
{
  # Import
  bed <- vl_importBed(a)
  a <- vl_importBed(a)
  a[, idx:= .I]
  a[, minOv:= min.overlap]
  b <- vl_importBed(b)
  # Prepare foverlaps
  if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
  {
    a <- a[, .(seqnames, start, end, strand, idx, minOv)]
    b <- b[, .(seqnames, start, end, strand)]
    setkeyv(a, c("seqnames", "strand", "start", "end"))
    setkeyv(b, c("seqnames", "strand", "start", "end"))
  }else
  {
    a <- a[, .(seqnames, start, end, idx, minOv)]
    b <- b[, .(seqnames, start, end)]
    setkeyv(a, c("seqnames", "start", "end"))
    setkeyv(b, c("seqnames", "start", "end"))
  }
  # Intersect
  inter <- foverlaps(b, a, nomatch= NULL)
  inter[, maxStart:= apply(.SD, 1, max), .SDcols= c("start", "i.start")]
  inter[, minEnd:= apply(.SD, 1, min), .SDcols= c("end", "i.end")]
  inter <- inter[, .(minEnd-maxStart+1>=minOv, idx)]
  sel <- unique(inter[(V1), idx])
  if(invert)
    sel <- setdiff(a$idx, sel)
  # Return
  return(bed[(sel)])
}

#' Compute bed coverage
#'
#' For each bin, computes the number of overlapping reads from a bed file
#' @param a Ranges for which overlaps with b should be counted. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param b Ranges from which overlaps should be computed.
#' @param ignore.strand Should strand be ignored? default= T, meaning all overlapping elements in 'b'  will be considered
#' @param min.overlap A vector integer of length 1 or nrow(a) specifying the minimum overlap (bp) require to count overlaps. Default= 1L
#' @examples
#' a <- data.table(seqnames= "chr2R",
#'                 start= 10000,
#'                 end= 20000,
#'                 strand= c("+", "-"))
#' b <- data.table(seqnames= "chr2R",
#'                 start= seq(10000, 19000, 1000),
#'                 end= seq(11000, 20000, 1000),
#'                 strand= rep(c("+", "-"), each= 10))
#' vl_covBed(a, b)
#' vl_covBed(a, b, ignore.strand = F)
#' 
#' @return For each range in 'a', reports the number of overlapping features in 'b'
#' @export
vl_covBed <- function(a,
                      b,
                      ignore.strand= T,
                      min.overlap= 1L)
{
  # Import
  a <- vl_importBed(a)
  a[, idx:= .I]
  a[, minOv:= min.overlap]
  b <- vl_importBed(b)
  # Prepare foverlaps
  if(!ignore.strand && "strand" %in% names(a) && "strand" %in% names(b))
  {
    a <- a[, .(seqnames, start, end, strand, idx, minOv)]
    b <- b[, .(seqnames, start, end, strand)]
    setkeyv(a, c("seqnames", "strand", "start", "end"))
    setkeyv(b, c("seqnames", "strand", "start", "end"))
  }else
  {
    a <- a[, .(seqnames, start, end, idx, minOv)]
    b <- b[, .(seqnames, start, end)]
    setkeyv(a, c("seqnames", "start", "end"))
    setkeyv(b, c("seqnames", "start", "end"))
  }
  # Intersect
  inter <- foverlaps(a, b)
  inter[, maxStart:= apply(.SD, 1, max), .SDcols= c("start", "i.start")]
  inter[, minEnd:= apply(.SD, 1, min), .SDcols= c("end", "i.end")]
  # Return
  return(inter[, sum(minEnd-maxStart+1>=minOv, na.rm= T), keyby= idx]$V1)
}