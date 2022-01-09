#' Check if data.table file likely corresponds to a bed
#'
#' @param x Object to be testd
#' @return boolean
#' @export
vl_isDTranges <- function(x)
{
  all(c("seqnames", "start", "end") %in% names(x))
}

#-----------------------------------------------------------------------------------------------------------------------------------------#
# IMPORT
#-----------------------------------------------------------------------------------------------------------------------------------------#
#' Import bed file as data.table
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
  bed <- data.table::rbindlist(lapply(bed, function(x) fread(x, fill = T)))
  if(!vl_isDTranges(bed))
  {
    bedcols <- 1:ifelse(ncol(bed)>6, 6, ncol(bed))
    setnames(bed, 
             names(bed)[bedcols],
             c("seqnames", "start", "end", "name", "score", "strand")[bedcols])
  }
  return(bed)
}

#' @describeIn vl_importBed for GRanges
#' @export
vl_importBed.GRanges <- function(bed)
{
  return(data.table::as.data.table(bed))
}

#-----------------------------------------------------------------------------------------------------------------------------------------#
# EXPORT
#-----------------------------------------------------------------------------------------------------------------------------------------#
#' Export bed file
#'
#' @param bed Either a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @return Imported bed
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
    bed <- vl_importBed(bed)
  cols <- c("seqnames", "start", "end", "name", "score", "strand")
  if("name" %in% names(bed))
    bed$name <- "."
  if("score" %in% names(bed))
    bed$score <- 0
  if("strand" %in% names(bed))
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


#-----------------------------------------------------------------------------------------------------------------------------------------#
# GENERATE SETS OF REGIONS
#-----------------------------------------------------------------------------------------------------------------------------------------#
#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param BSgenome BSgenome object to use.
#' @param bins_width bins width default to 50bp
#' @param steps_width steps width separating each bin. default set to bins_width
#' @examples 
#' vl_binBSgenome(BSgenome= BSgenome.Dmelanogaster.UCSC.dm3, bins_width= 50)
#' @return data.table containing bin coordinates
#' @export

vl_binBSgenome <- function(BSgenome,
                           bins_width= 50,
                           steps_width= bins_width)
{
  if(!class(BSgenome)=="BSgenome")
    stop("genome should be a BSgenome object!")
  dat <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  bins <- dat[, .(start= seq(1, end, steps_width)), .(seqnames, end, width)]
  bins[, end:= start+bins_width-1]
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
#' @examples 
#' vl_control_regions_BSgenome(BSgenome = BSgenome.Dmelanogaster.UCSC.dm3, 
#' n= 5000, 
#' width= 1000, 
#' restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"))
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
  dat <- data.table::as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  if(is.null(restrict_seqnames))
    restrict_seqnames <- dat$seqnames
  
  dat <- dat[seqnames %in% restrict_seqnames]
  data.table::setnames(dat, "width", "chr_size")
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

#-----------------------------------------------------------------------------------------------------------------------------------------#
# BEDTOOLS EQUIVALENTS
#-----------------------------------------------------------------------------------------------------------------------------------------#
#' Find closestBed regions
#'
#' For each line of a file, returns the closest lines in b
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself.
#' @param k If set to "all" (default), return all features that match the min(|distance|). Else, returns the k first features
#' @param min_dist Min distance for closest feature 
#' @examples 
#' a <- data.table(chr= "chr2L", start= sample(10000, 1000))
#' a[, end:= start:1000]
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
                          k= "all",
                          min_dist= 0)
{
  if(k!="all")
    if(!is.numeric(k))
      stop("k should either be set to all or be a numeric value")
  if(!vl_isDTranges(a))
    a <- vl_importBed(a)
  if(is.null(b))
    b <- a else if(!vl_isDTranges(b))
      b <- vl_importBed(b)
  # Main function
  res <- b[a, {
    # Make res object
    .c <- data.table(start.a= i.start,
                     end.a= i.end,
                     start.b= start,
                     end.b= end)
    # Compute distance
    .c[, dist:= {
      if(end<i.start)
        end-i.start else if (start>i.end)
          start-i.end else
            0L
    }, .(start, end)]
    # Distance check AND order based on distance
    .c <- .c[abs(dist)>=min_dist][order(abs(dist))]
    # Filter k closest features if specified
    if(k=="all")
      .c <- .c[abs(dist)==abs(dist)[1]] else
        .c <- na.omit(.c[1:k])
  }, .EACHI, on= "seqnames"]
  # Export
  return(res)
}

#' Collapse ranges with a certain min dist
#'
#' Collapse ranges with a certain min dist between them
#'
#' @param bed bed file to be collapsed. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames','start', 'end' columns. see ?vl_importBed()
#' @param mingap min gap ditance for merging. default is 1 (that is, touching coordinates)
#' @param stranded If set to true and strand column is provided, only merges coordinates with similar strand.
#' @param return_idx_only If set to T, does not collapse regions but returns idx as an extra columns. default= F
#' @param compute_other_columns If set to T or to a vector of column names, extra columns will be compute by contig using aggregate.fun. 
#' @param aggregate.fun See compute_other_columns. Default to function(x) mean(x, na.rm= T)
#' @examples 
#' ex <- GRanges(c("chr3L", "chr3L", "chr3L", "chr3L", "chr2R"), 
#' IRanges(c(1000, 2000, 2100, 2200, 2000), 
#' c(1500, 2099, 2199, 2299, 3000))
#' vl_collapseBed(ex, mingap= 1)
#' @return Collapse coor data.table
#' @export

vl_collapseBed <- function(bed,
                           mingap= 1,
                           return_idx_only= F, 
                           compute_other_columns= NULL,
                           aggregate.fun= function(x) mean(x, na.rm= T))
{
  # Hard copy of bed file
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  DT <- copy(bed)
  setkeyv(DT, c("seqnames", "start"))
  
  # Checks
  if(is.null(compute_other_columns)) # Check if compute_other_columns corresponds to bed columns or bool
  {
    cols <- as.character()
    compute_other_columns <- F
  }else 
  {
    if(all(compute_other_columns %in% names(bed)))
      cols <- compute_other_columns else
        cols <- setdiff(names(bed), c("seqnames", "start", "end"))
      compute_other_columns <- T
  }
  
  # Compute contig idx
  DT[, idx:= cumsum(sapply(start-cummax(end)[c(1, seq(.N-1))], function(x) x>mingap)), seqnames]
  DT[, idx:= .GRP, .(seqnames, idx)] # Make idx unique
  if(!return_idx_only)
  {
    if(compute_other_columns) # compute other columns per group
      DT[, (cols):= lapply(.SD, aggregate.fun), idx, .SDcols= cols] else
        DT <- DT[, .(seqnames, start, end, idx)] # Otherwise only report collapsed intervals
    # Compute contigs start and end
    DT[, c("start", "end") := .(min(start), max(end)), .(seqnames, idx)]
    # Collapse
    DT <- unique(DT[, !"idx"])
  }

  return(DT)
}

#' Shuffle bed file
#'
#' Shuffle bed file in respect of covered regions, total reads and reads length
#'
#' @param bed bed file to shuffle. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @examples 
#' ChIP_bed <- "/test/H3K27me3_PH18_merge_uniq.bed"
#' res <- shuffleBed(ChIP_bed)
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export

vl_shuffleBed <- function(bed)
{
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  
  # Compute contigs 
  .coll <- vl_collapseBed(bed)
  # Compute effective chr sizes
  .coll[, chr_size:= sum(end-start+1), keyby= seqnames]
  # Compute sampling prob depending on the size
  .coll[, prob:= (end-start+1)/chr_size, seqnames]
  # Compute random reads that take into account covered chr 
  reads <- bed[, .N, .(seqnames, read_length= end-start)]
  .rdm <- reads[, {
    chr <- .coll[.BY] #seqnames key
    set.seed(.GRP)
    .s <- chr[sample(seq(.N), N, T, prob), # Sample regions depending on their size
              .(region_start= start, region_end= end-read_length)] # Resize regions depending on read length
    .s[, start:= {
      if(region_start==region_end) # sample rdm start positions
        rep(region_start, .N) else
        {
          set.seed(.GRP)
          sample(region_start:region_end, .N, replace = T)
        }
    },  .(region_start, region_end)]
    .s[, end:= start+read_length]
  }, .(seqnames, read_length, N)]
  # Clean
  .rdm <- .rdm[, .(seqnames, start, end)]
  setkeyv(.rdm, c("seqnames", "start"))
  
  return(.rdm)
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

#-----------------------------------------------------------------------------------------------------------------------------------------#
# Peak-calling related functions
#-----------------------------------------------------------------------------------------------------------------------------------------#
#' Compute ChIP enrichment
#'
#' Given a set of ChIP and Input counts, performs hypergeometric test to know if ChIP is enriched compared to Input
#'
#' @param regions Regions to analyse. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param ChIP_bed ChIP bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input_bed Input bed file.
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export

vl_enrichBed <- function(regions,
                         ChIP_bed,
                         Input_bed)
{
  # Hard copy regions
  if(!vl_isDTranges(regions))
    regions <- vl_importBed(regions)
  regions <- copy(regions)
  # ChIP coverage
  if(!vl_isDTranges(ChIP_bed))
    ChIP_bed <- vl_importBed(ChIP_bed)
  regions[, ChIP_counts:= vl_covBed(regions, ChIP_bed)]
  regions <- merge(regions,
                   unique(ChIP_bed[, .(ChIP_total_counts= .N), seqnames]))
  # Input coverage
  if(!vl_isDTranges(Input_bed))
    Input_bed <- vl_importBed(Input_bed)
  regions[, Input_counts:= vl_covBed(regions, Input_bed)]
  regions <- merge(regions,
                   unique(Input_bed[, .(Input_total_counts= .N), seqnames]))
  # Compute enrichment and pval
  regions[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, Input_counts, ChIP_total_counts, Input_total_counts)]
  regions[, padj:= p.adjust(pval, "fdr")]
  # Format narrowpeak file
  regions[pval==0, pval:= min(regions[pval>0, pval])]
  regions[padj==0, padj:= min(regions[padj>0, padj])]
  regions <- regions[, .(seqnames, start, end,
                         names= paste0("peak_", .I), 
                         score= round(OR/max(OR)*1000), 
                         strand= ".", 
                         signalValue= OR,
                         pValue= -log10(pval),
                         qValue= -log10(padj),
                         peak= -1)]

  return(regions)
}

#' Compute peak calling
#'
#' Compute peak calling using ChIP and Input
#'
#' @param ChIP_bed ChIP bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input_bed Input bed file.
#' @param average_function Function applied to ChIP read counts to define candidate peaks
#' @param gaussian_blur Applies gaussian blur to ChIP signal to identify peak candidates. Useful for noisy signal. default= FALSE
#' @param BSgenome A BSgenome object used for gw binning
#' @param bins_width bins width used for peak calling, Default to 100L (narrow Peaks). Use larger bins to call domains
#' @param steps_width Step between bins. default to round(bins_width/2)
#' @param bins_pval_cutoff The cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 0.05
#' @param bins_OR_cutoff Enrichment cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 1
#' @return significantly enriched peaks
#' @export

vl_peakCalling <- function(ChIP_bed,
                           Input_bed,
                           average_function= function(x) ceiling(mean(x)),
                           gaussian_blur= F,
                           BSgenome,
                           bins_width= 100,
                           steps_width= round(bins_width/2),
                           bins_pval_cutoff= 0.05,
                           bins_OR_cutoff= 1)
{
  #----------------------------#
  # Find potentially enriched bins gw and define peak candidates
  #----------------------------#
  # genome wide bins
  bins <- vl_binBSgenome(BSgenome, 
                         bins_width = bins_width, 
                         steps_width = steps_width)
  # Compute ChIP coverage and average signal
  if(!vl_isDTranges(ChIP_bed))
    ChIP_bed <- vl_importBed(ChIP_bed)
  bins[, ChIP_counts:= vl_covBed(bins, ChIP_bed)]
  bins <- merge(bins,
                unique(ChIP_bed[, .(ChIP_total_counts= .N), seqnames]))
  # Smooth signal?
  if(gaussian_blur)
    bins[, ChIP_counts:= round(vl_gaussian_blur(ChIP_counts))]
  # Average signal will be used to identify candidates
  bins[, average_counts:= average_function(ChIP_counts[ChIP_counts>0]), seqnames]
  bins[, average_total_counts := ChIP_total_counts, seqnames]
  # Remove bins that do not contain reads
  bins <- bins[ChIP_counts>0 & average_counts>0]
  # Enrichment
  bins[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, average_counts, ChIP_total_counts, average_total_counts)]
  # Only retain bins significantly enriched over background
  peaks <- bins[pval <= bins_pval_cutoff & OR >= bins_OR_cutoff]
  # Collapse touching bins into candidate peaks
  peaks <- vl_collapseBed(peaks)
  
  #----------------------------#
  # Test each candidate peaks and output significant ones
  #----------------------------#
  final <- vl_enrichBed(peaks, 
                        ChIP_bed, 
                        Input_bed)
  # Format and save significant peaks
  final <- final[qValue>(-log10(0.05))]
  
  return(final)
}
