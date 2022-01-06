#' Import bed file as data.table
#'
#' Imports bed as data.table and check formats
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @return Imported bed
#' @export
vl_importBed <- function(bed, ...) UseMethod("vl_importBed")

#' @describeIn vl_importBed for bed paths
#' @export
vl_importBed.character <- function(bed)
{
  bed <- data.table::rbindlist(lapply(bed, function(x) fread(x, fill = T)))
  if(!all(c("seqnames", "start", "end") %in% names(bed)))
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
  
#' @describeIn vl_importBed for data.table
#' @export
vl_importBed.data.table <- function(bed)
{
  if(!all(c("seqnames", "start", "end") %in% names(bed)))
    stop("If bed is a data.table, it should contain 'seqnames', 'start' and 'end' columns")
  return(bed)
} 

#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param BSgenome BSgenome object to use.
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

#' Find closestBed regions
#'
#' For each line of a file, returns the closest lines in b
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself.
#' @param k If set to "all" (default), return all features that match the min(|distance|). Else, returns the k first features
#' @param min_dist Min distance for closest feature 
#' @examples 
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
  a <- vl_importBed(a)
  if(is.null(b))
    b <- copy(a) else
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
#' c(1500, 2099, 2199, 2299, 3000)),
#' strand= c("+","+","+","-","+"))
#' vl_collapseBed(ex, mingap= 1, stranded= F)
#' vl_collapseBed(ex, mingap= 1, stranded= T)
#' @return Collapse coor data.table
#' @export

vl_collapseBed <- function(bed,
                           mingap= 1,
                           return_idx_only= F, 
                           compute_other_columns= NULL,
                           aggregate.fun= function(x) mean(x, na.rm= T))
{
  # Hard copy
  DT <- copy(vl_importBed(bed))
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
  bed <- vl_importBed(bed)
  
  # Compute contigs 
  .rdm <- vl_collapseBed(bed)
  # Compute effective chr sizes
  .rdm[, chr_size:= sum(end-start+1), keyby= seqnames]
  # Compute sampling prob depending on the size
  .rdm[, prob:= (end-start+1)/chr_size, seqnames]
  # Compute random reads that take into account covered chr 
  reads <- bed[, .N, .(seqnames, read_length= end-start)]
  .i <- reads[, {
    chr <- .rdm[.BY] #seqnames key
    set.seed(1)
    .s <- chr[sample(seq(.N), N, T, prob), # Sample regions depending on their size
              .(region_start= start, region_end= end-read_length)] # Resize regions depending on read length
    set.seed(1)
    .s[, start:= {
      if(region_start==region_end) # sample rdm start positions
        rep(region_start, .N) else
          sample(region_start:region_end, .N, replace = T)
    },  .(region_start, region_end)]
    .s[, end:= start+read_length]
  }, .(seqnames, read_length, N)]
  # Clean
  .i[, total_counts:= .N, seqnames]
  .i <- .i[, .(seqnames, start, end, total_counts)]
  setkeyv(.i, c("seqnames", "start"))
  
  return(.i)
}

#' Compute Enrichment
#'
#' Compute enrichment for a list of regions using fisher test to compare ChIP and Input
#'
#' @param bins bins for which enrichment has to be computed. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param ChIP_bed ChIP bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input_bed Input bed file. Default to NULL, meaning shuffle ChIP bed is used a background
#' @examples 
#' ChIP_bed <- "/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed"
#' res <- vl_enrichBed(ChIP_bed)
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export

vl_enrichBed <- function(bins,
                         ChIP_bed, 
                         Input_bed= NULL)
{
  # Count overlapping reads
  regions <- copy(vl_importBed(bins))
  # ChIP
  .c <- vl_importBed(ChIP_bed)
  .c[, total_counts:= .N, seqnames]
  # Input
  if(is.null(Input_bed))
    .i <- vl_shuffleBed(.c) else
    {
      .i <- vl_importBed(Input_bed)
      .i[, total_counts:= .N, seqnames]
    }
  
  # Count overlapping reads
  # ChIP
  regions$ChIP_counts <- .c[regions, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
  regions$ChIP_total_counts <- regions[, rep(.c[.BY, .N, on= "seqnames"], .N), seqnames]$V1
  # Input
  regions$Input_counts <- .i[regions, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
  regions$Input_total_counts <- regions[, rep(.i[.BY, .N, on= "seqnames"], .N), seqnames]$V1
  # Remove bins without reads
  regions <- na.omit(regions)
  
  # Compute enrichment over background
  regions[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, Input_counts, ChIP_total_counts, Input_total_counts)]
  regions[, padj:= p.adjust(pval, "fdr")]
  return(regions)
}

#' Compute peak calling
#'
#' Compute peak calling using ChIP and Input
#'
#' @param ChIP_bed ChIP bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input_bed Input bed file. Default to NULL, meaning shuffle ChIP bed is used a background
#' @param BSgenome A BSgenome object used for gw binning
#' @param binsize binsize used for peak calling, Default to 100 (narrow Peaks). Use larger bins to call domains
#' @examples 
#' @return significantly enriched peaks
#' @export

vl_peakCalling <- function(ChIP_bed,
                           Input_bed= NULL,
                           BSgenome,
                           binsize= 100)
{
  #----------------------------#
  # Import data and compute gw bins enrichment
  #----------------------------#
  # genome wide bins
  bins <- vl_binBSgenome(BSgenome, bin_size = binsize)
  # Import bed files
  .c <- vl_importBed(ChIP_bed)
  .c[, total_counts:= .N, seqnames]
  # Input
  if(is.null(Input_bed))
    .i <- vl_shuffleBed(.c) else
    {
      .i <- vl_importBed(Input_bed)
      .i[, total_counts:= .N, seqnames]
    }
  # Enrichment
  bins <- vl_enrichBed(ChIP_bed = .c, 
                       Input_bed = .i, 
                       bins = bins)
  
  #----------------------------#
  # Merge contiguous bins into peaks
  #----------------------------#
  # Only retain bins significantly enriched over background
  peaks <- bins[pval<1e-5]
  # Collapse touching bins into candidate peaks
  peaks <- vl_collapseBed(peaks)
  # Compute overall peaks enrihment
  peaks <- vl_enrichBed(ChIP_bed = .c, 
                        Input_bed = .i, 
                        bins = peaks)
  
  #----------------------------#
  # Clean and save narrowPeak object
  #----------------------------#
  final <- peaks[padj<0.05]
  final[pval==0, pval:= min(final[pval>0, pval])]
  final[padj==0, padj:= min(final[padj>0, padj])]
  final <- final[, .(seqnames, start, end,
                     names= paste0("peak_", .I), 
                     score= round(OR/max(OR)*1000), 
                     strand= ".", 
                     signalValue= OR,
                     pValue= -log10(pval),
                     qValue= -log10(padj),
                     peak= -1)]
  
  return(final)
}
