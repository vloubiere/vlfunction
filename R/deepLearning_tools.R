
#' Title
#'
#' @param h5 Path to an h5 file containing the contribution scores.
#' @param bed Bed file containing the regions for which contributions were computed. By default, the file is searched in the same folder as the h5 file.
#' @param selection An optional bed file containing a selected set of regions. If specified, only the regions in bed overlapping with the selection will be returned
#'
#' @return A contribution data.table containing, for each nucleotide, the corresponding base and its associated contribution score.
#' @export
#'
#' @examples
vl_importContrib <- function(h5,
                             bed= list.files(dirname(h5), ".bed$", full.names = TRUE),
                             selection)
{
  # Import contributions ----
  dat <- rhdf5::h5read(h5, "contrib_scores/class")
  
  # Import coordinates ----
  bed <- vl_importBed(bed)
  if(!"strand" %in% names(bed))
    bed[strand:= "*"]

  # Check if compatible ----
  if(nrow(bed) != dim(dat)[3])
    stop("Number of regions in bed file do not match the number of regions in the h5 contribution file.
         Please provide the path to the bed file containing the regions for which contributions were computed")
  
  # Restrict to bins of interest ----
  if(!missing(selection))
  {
    selRegions <- vl_importBed(selection)
    keep <- vl_covBed(bed, selRegions, ignore.strand = FALSE)>0
    bed <- bed[(keep)]
    dat <- dat[, , (keep), drop= FALSE]
  }
  # If no overlaps with selection, return empty result ----
  res <- if(!nrow(bed))
  {
    data.table(binIDX= NULL,
               seqnames= NULL,
               start= NULL,
               end= NULL,
               strand= NULL,
               score= NULL)
  }else
  {
    # Expand bins back to single nt ----
    bed[, binIDX:= .I]
    contrib <- if("strand" %in% names(bed))
    {
      bed[, {
        # Coordinates (depending on strand)
        coor <- if(strand=="-") seq(end, start) else seq(start, end)
        .(seqnames,
          start= coor,
          end= coor,
          strand= strand)
      }, binIDX]
    }else
      bed[, {
        # Coordinates (no strand)
        coor <- seq(start, end)
        .(seqnames,
          start= coor,
          end= coor)
      }, binIDX]

    # Add score ----
    contrib$score <- unlist(lapply(seq(dim(dat)[3]), function(i) rowSums(dat[,,i])))
    
    # Refine for the precise overlaps with selection ----
    if(!missing(selection))
      contrib <- vl_intersectBed(contrib, selRegions, ignore.strand = FALSE)
    
    # Clean ----
    contrib[, .(binIDX, seqnames, start, end, strand, score)]
  }
  
  # Return object ----
  return(res)
}

#' Plot contribution scores matrix
#'
#' @param bed A bed file containing a unique region for which contrib scores will be plotted.
#' @param h5 Path(s) to h5 files containing the contribution scores.
#' @param h5.bed Bed files containing the coordinates of the regions corresponding to provided h5 files.
#' @param genome The genome to be used.
#' @param agg.FUN In the case were several contribution scores would be found for a single nt, how should they be aggregated? Default= function(x) mean(x)
#' @param mot An optional bed file containing motifs to be added.
#' @param mot.name.column Name of the column containing the motif name.
#' @param xlab Default= "nt"
#' @param ylab Default= "Contribution"
#' @param xlim Default= sequence length
#' @param ylim Default= range(contrib)
#'
#' @return contrib plot
#' @export
vl_plot_contrib_logo <- function(bed,
                                 h5,
                                 h5.bed= sapply(dirname(h5), function(x) list.files(x, ".bed$", recursive = T, full.names = T)),
                                 genome,
                                 agg.FUN= function(x) mean(x),
                                 mot,
                                 mot.name.column= "motif_ID",
                                 xlab= "nt",
                                 ylab= "Contribution",
                                 xlim,
                                 ylim)
{
  # Import Bed
  bed <- vl_importBed(bed) # Very important!
  if(nrow(bed)>1)
    stop("Unique reigon should be specified")
  # Metadata ----
  meta <- data.table(h5= h5,
                     h5.bed= h5.bed)
  
  # Import contribution scores ----
  dat <- meta[, {
    vl_importContrib(h5,
                     bed = h5.bed,
                     selection= bed)
  }, .(h5, h5.bed)]
  
  # Aggregate if necessary ----
  if(uniqueN(dat[, .(seqnames, start)]) != nrow(dat))
  {
    message("Some nucleotides had >1 contribution score assigned to it, which will be aggregated using agg.FUN")
    dat <- dat[, .(score= agg.FUN(score)), .(seqnames, start, end)]
  }
  
  # Get sequence ----
  dat$base <- strsplit(vl_getSequence(bed, genome), "")[[1]]
  
  # Plotting vars ----
  dat[, xleft:= .I-1]
  if(missing(xlim))
    xlim <- c(0, nrow(dat))
  if(missing(ylim))
    ylim <- range(dat$score)

  # Plotting ----
  plot(NA,
       xlim= xlim,
       ylim= ylim,
       xlab= xlab,
       ylab= ylab,
       frame= FALSE)
  
  dat[, {
    vl_plotLetter(base,
                  xleft = xleft,
                  ytop= score,
                  width = 1,
                  height = score)
  }, (dat)]
  
  # Add motif boxes ----
  if(!missing(mot))
  {
    mot <- vl_importBed(mot)
    mot <- vl_intersectBed(mot, bed, ignore.strand= TRUE)
    if(nrow(mot))
      mot[, {
        xl <- start-bed$start
        xr <- end-bed$start
        rect(xleft = xl,
             ybottom = ylim[1],
             xright = xr,
             ytop = ylim[2])
        text((xl+xr)/2,
             ylim[2],
             get(mot.name.column)[1],
             pos= 3,
             xpd= T)
      }, (mot)]
  }
}


