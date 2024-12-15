#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param genome BSgenome object to use. Example "dm3"
#' @param bins.width bins width default to 50bp
#' @param steps.width steps width separating each bin. default set to bins.width
#' @param restrict.seqnames If specified, bins are restricted to provided seqnames
#' 
#' @examples 
#' vl_binBSgenome(genome= "dm3",
#'                bins.width= 50)
#'                
#' @return data.table containing bin coordinates
#' @export
vl_binBSgenome <- function(genome,
                           bins.width= 50L,
                           steps.width= bins.width,
                           restrict.seqnames= NULL)
{
  dat <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome, load.only = TRUE)))
  dat <- data.table::as.data.table(dat)
  # Restrict chromosomes
  if(!is.null(restrict.seqnames))
  {
    if(!any(dat$seqnames %in% restrict.seqnames))
      stop(paste("Provided seqnames do not exist in", genome)) else
        dat <- dat[as.character(seqnames) %chin% as.character(restrict.seqnames)]
  }
  dat <- dat[, .(seqnames, width)]
  # Compute bins start and end
  dat <- dat[, .(start= seq(1, max(width, 1), steps.width)), .(seqnames, width)]
  dat[, end:= start+bins.width-1]
  dat[end>width, end:= width]
  return(dat)
}

#' Generate control regions
#'
#' Generate regions with similar dispersion and widths ditrib than provided bed file
#'
#' @param bed Bed file used to produce similar control
#' @param genome BSgneome object to use. ex: "dm3", "dm6"
#' @param no.overlap If set to TRUE, avoids overlap between control sequences and the original bed file. Default= FALSE.
#' 
#' @examples 
#' vl_control_regions_BSgenome(vl_SUHW_top_peaks, "dm3")
#' 
#' @return data.table containing control regions
#' @export
vl_control_regions_BSgenome <- function(bed, genome, no.overlap= F)
{
  bed <- vl_importBed(bed)
  # Format
  regions <- data.table::copy(bed)
  regions[, width:= end-start+1]
  BS <- data.table::as.data.table(GenomicRanges::GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome, load.only = TRUE))))
  regions[BS, seqlength:= i.width, on= "seqnames"]
  # Random sampling
  regions[, start:= sample(seqlength-width, .N), .(seqlength, width)]
  regions[, end:= start+width-1]
  # Avoid overlaps between original bed file and control
  if(no.overlap)
  {
    regions[, cov:= vl_covBed(regions, bed)]
    i <- 1
    while(any(regions$cov))
    {
      set.seed(i)
      regions[cov>0, start:= sample(seqlength-width, .N), .(seqlength, width)]
      regions[cov>0, end:= start+width-1]
      regions[cov>0, cov:= vl_covBed(.SD, bed)]
      i <- i+1
    }
    regions$cov <- NULL
  }
  # Make sure stays within genome range
  return(regions)
}

#' random region 
#' 
#' Sample random regions from BSgenome
#'
#' @param n Number of regions to sample
#' @param width Widths of regions to sample (length should either be 1 or equal n)
#' @param restrict.seqnames If specified, only the provided seqnames will be used
#' 
#' @examples
#' vl_control_regions_BSgenome("dm3", 100, 1)
#' 
#' @return data.table containing randomly sampled regions
#' @export
vl_random_regions_BSgenome <- function(genome,
                                       n,
                                       width= 1,
                                       restrict.seqnames= NULL)
{
  # Format
  regions <- data.table::as.data.table(GenomicRanges::GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome, load.only = TRUE))))
  setnames(regions, "width", "seqlength")
  if(!is.null(restrict.seqnames))
    regions <- regions[seqnames %in% restrict.seqnames]
  regions <- regions[sample(nrow(regions), n, replace = T, prob = regions$width)]
  regions[, width:= width]
  # Random sampling
  regions[, start:= sample(seqlength-width, .N), .(seqlength, width)]
  regions[, end:= start+width-1]
  # Make sure stays within genome range
  return(regions)
}

#' Get genomic sequence
#' 
#' Returns the sequences of a DTranges
#'
#' @param bed Bed file for which regions have to be returned
#' @param genome BSgenome ID ("dm3", "dm6"...)
#' 
#' @examples
#' vl_getSequence(vl_SUHW_top_peaks, "dm3")
#' 
#' @return Character vector of sequences
#' @export
vl_getSequence <- function(bed,
                           genome)
{
  # Import ----
  bed <- vl_importBed(bed)
  if(!"strand" %in% names(bed))
  {
    message("'bed' strand set to unstranded (*)")
    bed[, strand:= "*"]
  }
  # Make sequence names ----
  names <- if("end" %in% names(bed))
    bed[, paste0(seqnames, ":", start, "-", end, ":", strand)] else
      bed[, paste0(seqnames, ":", start, ":", strand)]
  # Extract ----
  sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome, load.only = TRUE), 
                                names= bed$seqnames, 
                                start= bed$start, 
                                end= bed$end, 
                                strand= bed$strand, 
                                as.character= T)
  names(sequences) <- names
  # Return ----
  return(sequences)
}