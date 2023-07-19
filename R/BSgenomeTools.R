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
                           bins_width= 50L,
                           steps_width= bins_width,
                           restrict_seqnames= NULL)
{
  dat <- data.table::as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  # Restrict chromosomes
  if(!is.null(restrict_seqnames))
  {
    if(!any(dat$seqnames %in% restrict_seqnames))
      stop(paste("Provided seqnames do not exist in", genome)) else
        dat <- dat[as.character(seqnames) %chin% as.character(restrict_seqnames)]
  }
  dat <- dat[, .(seqnames, width)]
  # Compute bins start and end
  dat <- dat[, .(start= seq(1, max(c(width-bins_width, 1)), steps_width)), .(seqnames, width)]
  dat[, end:= start+bins_width-1]
  dat[end>width, end:= width]
  return(dat)
}

#' Generate control regions
#'
#' Generate regions with similar dispersion and widths ditrib than provided bed file
#'
#' @param genome BSgneome object to use. ex: "dm3", "dm6"
#' @param bed Bed file used to produce similar control
#' @return data.table containing control regions
#' @examples 
#' vl_control_regions_BSgenome(vl_SUHW_top_peaks, "dm3")
#' @export
vl_control_regions_BSgenome <- function(bed, genome)
{
  bed <- vl_importBed(bed)
  # Format
  regions <- data.table::copy(bed)
  regions[, width:= end-start+1]
  BS <- data.table::as.data.table(GenomicRanges::GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  regions[BS, seqlength:= i.width, on= "seqnames"]
  # Random sampling
  regions[, start:= sample(seqlength-width, .N), .(seqlength, width)]
  regions[, end:= start+width-1]
  # Make sure stays within genome range
  return(regions)
}

#' random region 
#' 
#' Sample random regions from BSgenome
#'
#' @param n Number of regions to sample
#' @param width Widths of regions to sample (length should either be 1 or equal n)
#' @param restrict_seqnames If specified, only the provided seqnames will be used
#' @return data.table containing randomly sampled regions
#' @examples
#' vl_control_regions_BSgenome("dm3", 100, 1)
#' @export
vl_random_regions_BSgenome <- function(genome,
                                       n,
                                       width= 1,
                                       restrict_seqnames= NULL)
{
  # Format
  regions <- data.table::as.data.table(GenomicRanges::GRanges(GenomeInfoDb::seqinfo(BSgenome::getBSgenome(genome))))
  setnames(regions, "width", "seqlength")
  if(!is.null(restrict_seqnames))
    regions <- regions[seqnames %in% restrict_seqnames]
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
#' @return Character vector of sequences
#' @examples
#' vl_getSequence(vl_SUHW_top_peaks, "dm3")
#' @export
vl_getSequence <- function(bed, genome)
{
  bed <- vl_importBed(bed)
  if(!"strand" %in% names(bed))
  {
    message("'bed' strand set to unstranded (*)")
    bed[, strand:= "*"]
  }
  sequences <- BSgenome::getSeq(BSgenome::getBSgenome(genome), 
                                names= bed$seqnames, 
                                start= bed$start, 
                                end= bed$end, 
                                strand= bed$strand, 
                                as.character= T)
  return(sequences)
}
