#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param genome BSgenome object to use. Example "dm3"
#' @param bins_width bins width default to 50bp
#' @param steps_width steps width separating each bin. default set to bins_width
#' @param restrict_seqnames If specified, bins are restricted to provided seqnames
#' @examples 
#' vl_binBSgenome(genome= "dm3", bins_width= 50)
#' @return data.table containing bin coordinates
#' @export
vl_binBSgenome <- function(genome,
                           bins_width= 50,
                           steps_width= bins_width,
                           restrict_seqnames= NULL)
{
  if(steps_width>bins_width)
    warning("steps_width>bins_width, meaning bins will not be contiguous")
  dat <- data.table::as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  if(!is.null(restrict_seqnames))
    dat <- dat[seqnames %in% restrict_seqnames]
  bins <- dat[, .(start= seq(1, end, steps_width)), .(seqnames, end, width)]
  bins[, end:= start+bins_width-1]
  bins[end>width, end:= width]
  bins <- bins[end-start>0, .(seqnames, start, end)]
  setkeyv(bins, c("seqnames", "start", "end"))
  return(bins)
}

#' Generate control regions
#'
#' Generate regions with similar dispersion and widths ditrib than provided bed file
#'
#' @param genome BSgneome object to use. ex: "dm3", "dm6"
#' @param bed Bed file used to produce similar control
#' @return data.table containing control regions
#' @export
vl_control_regions_BSgenome <- function(bed, genome)
{
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  # Format
  regions <- data.table::copy(bed)
  regions[, width:= end-start+1]
  dat <- data.table::as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  sample <- dat[regions, .(start= width, width= i.width), .EACHI, on= "seqnames"]
  # Random sampling
  sample[, start:= sample(start, .N), .(seqnames, start)]
  # Make sure stays within genome range
  sample[, start:= start-width+1]
  sample[start<1, start:= 1]
  sample[, end:= start+width-1]
  return(sample[,.(seqnames, start, end)])
}

#' random region 
#' 
#' Sample random regions from BSgenome
#'
#' @param n Number of regions to sample
#' @param width Widths of regions to sample (length should either be 1 or equal n)
#' @param restrict_seqnames If specified, only the provided seqnames will be used
#' @return data.table containing randomly sampled regions
#' @export
vl_random_regions_BSgenome <- function(genome,
                                       n,
                                       width= 1,
                                       restrict_seqnames= NULL)
{
  # Format
  dat <- data.table::as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  if(!is.null(restrict_seqnames))
    dat <- dat[seqnames %in% restrict_seqnames]
  # Random sampling
  dat <- dat[sample(x = seq(nrow(dat)), 
                    size = n, 
                    prob = dat$end-dat$start+1, replace= T)]
  dat$width <- width
  sample <- dat[, .(seqnames, start= end-width, width)]
  sample[, start:= sample(start, .N), .(seqnames, start)]
  # Make sure stays within genome range
  sample[, start:= start-width+1]
  sample[start<1, start:= 1]
  sample[, end:= start+width-1]
  return(sample[,.(seqnames, start, end)])
}
