#' Bin BS genome
#'
#' Very fast using data.table
#'
#' @param BSgenome BSgneome object to use.
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
#' 
#' @examples 
#' vl_control_regions_BSgenome(BSgenome = BSgenome.Dmelanogaster.UCSC.dm3, 
#' n= 5000, 
#' width= 1000, 
#' restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"))
#' 
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
  dat <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  if(is.null(restrict_seqnames))
    restrict_seqnames <- dat$seqnames
  
  dat <- dat[seqnames %in% restrict_seqnames]
  colnames(dat)[4] <- "chr_size"
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

#' Find closestBed idx
#'
#' For each line of a file, return the idx of the closest line in b (or non orverlapping line in a if b is null)
#'
#' @param a Granges or data.table FOR which closest features have to be found.
#' @param b Granges or data.table FROM which closest features have to be found. If set to NULL (default), then a is matched to itself while excluding self-matching.
#' @param min_dist Min distance for closest feature 
#' @examples 
#' Example 1:
#' a <- b <- data.table(seqnames= c("chr4", "chr2R", "chr2R", "chr3L"),
#' start= c(1e6, 1e6, 10e6, 1e6),
#' end= c(1.1e6, 1.1e6, 10.1e6, 1.1e6))
#' b <- a[-1]
#' b[, start:= start+5e5]
#' b[, end:= end+5e5]
#' closestBed(a, b)
#' closestBed(a, b, min_dist = 1e6)
#' closestBed(a, a)
#' closestBed(a, a, min_dist= 1)
#' closestBed(a)
#' 
#' Benchmarking:
#' require(rtracklayer)
#' tss <- as.data.table(import("/groups/stark/vloubiere/genomes/ensembl/dm3/Drosophila_melanogaster.BDGP5.77.gtf"))
#' tss <- tss[type=="gene", .(seqnames, start= ifelse(strand=="+", start+0, end-1)), gene_name]
#' tss[, end:= start+1]
#' This:
#' c1 <- closestBed(tss)
#' Is pretty similar to this:
#' c2 <- closestBed(tss, tss, 1)
#' Except for the cases where two genes have exact same start. Example check c1[67] and c2[67], and compare to tss[65:67]
#' closestBed(tss, tss)
#' closestBed(tss, min_dist = 100000)
#' 
#' @return For each line in a, returns the idx of its closest feature in b (vector)
#' @export

vl_closestBed <- function(a, 
                          b= NULL, 
                          min_dist= 0)
{
  if(!is.data.table(a))
    a <- as.data.table(a)
  a[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
  a[, idx:= .I]
  if(is.null(b))
  {
    b <- a
    check_idx <- T
  }else
  {
    if(!is.data.table(b))
      b <- as.data.table(b)
    b[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
    b[, idx:= .I]
    check_idx <- F
  }
  # Compute
  res <- b[a, {
    dist <- abs(center-i.center)
    dist[dist<min_dist] <- NA
    if(check_idx)
      dist[idx==i.idx] <- NA
    if(all(is.na(dist)))
      as.integer(NA) else
        idx[which.min(dist)]
  }, .EACHI, on= "seqnames"]$V1
  return(res)
}

#' Collapse ranges with a certain min dist
#'
#' Collapse ranges with a certain min dist between them
#'
#' @param bed GRanges or data.table file with "seqnames", "start", "end" (and optionally strand, set to * if absent)
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
#' vl_collapse_DT_ranges(ex, mingap= 1, stranded= F)
#' vl_collapse_DT_ranges(ex, mingap= 1, stranded= T)
#' 
#' @return Collapse coor data.table
#' @export

vl_collapse_DT_ranges <- function(bed, 
                                  mingap= 1,
                                  return_idx_only= F, 
                                  compute_other_columns= NULL,
                                  aggregate.fun= function(x) mean(x, na.rm= T))
{
  if(class(bed)[1]=="GRanges")
    bed <- data.table::as.data.table(bed)
  if(!data.table::is.data.table(bed) | !all(c("seqnames", "start", "end") %in% names(bed)))
    stop("bed must be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns")
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
  
  # Hard copy
  DT <- copy(bed[, match(c("seqnames", "start", "end", cols), names(bed)), with= F])
  setkeyv(DT, c("seqnames", "start"))
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

#' Compute Enrichment
#'
#' Compute enrichment for a list of regions using fisher test to compare ChIP and Input
#'
#' @param bins Data.table containing bins (seqnames, start and end columns)
#' @param ChIP_bed Either a character vector poiting to bed file path OR a data.table containing ChIP reads
#' @param Input_bed Either a character vector poiting to bed file path OR a data.table containing Input reads. 
#' If set to NULL (default), use shuffled version of the ChIP reads as input.
#' @examples 
#' ChIP_bed <- "/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed"
#' res <- vl_bedEnrichment(ChIP_bed)
#' 
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export

vl_bedEnrichment <- function(bins,
                             ChIP_bed, 
                             Input_bed= NULL)
{
  # Function to Import bed files (if ChIP_bed and Input_bed are not already data.tables)
  bed_import <- function(x)
  {
    if(!is.data.table(bins))
      stop("bins should be a data.table containing seqnames, start and end columns")
    if(!all(c("seqnames", "start", "end") %in% names(bins)))
      stop("bins should contain seqnames, start and end columns")
    # Import function
    .c <- rbindlist(lapply(x, 
                           fread, 
                           sel= 1:3, 
                           col.names= c("seqnames", "start", "end"), 
                           key= c("seqnames", "start")))
    .c[, total_counts:= .N, seqnames]
    return(.c)
  }
  # ChIP
  if(is.character(ChIP_bed))
    .c <- bed_import(ChIP_bed) else if(is.data.table(ChIP_bed))
    {
      if(!all(c("seqnames", "start", "end") %in% names(ChIP_bed)))
        stop("Could not find seqnames, start, end in ChIP_bed data.table. Either provide such data.table or a .bed path")
      .c <- ChIP_bed
    }
  # Input
  if(is.character(Input_bed))
    Input_bed <- bed_import(Input_bed) else if (is.data.table(Input_bed))
    {
      if(!all(c("seqnames", "start", "end") %in% names(Input_bed)))
        stop("Could not find seqnames, start, end in Input_bed data.table. Either provide such data.table or a .bed path")
      .i <- Input_bed
    } else if (is.null(Input_bed)) # If Input not provided, use random shuffling of ChIP reads
    {
      # Compute contigs 
      .i <- vl_collapse_DT_ranges(bed = .c, 
                                  compute_other_columns = "total_counts", 
                                  aggregate.fun = function(x) max(x, na.rm=T))
      # resize depending on reads sizes (random starts will be sampled!)
      .i[, end:= end-min(.c[, end-start]), seqnames]
      # Compute effective chr sizes
      .i[, chr_size:= sum(end-start+1), seqnames]
      # Compute sample size fo each contig
      .i[, sample_size:= round(total_counts*((end-start+1)/chr_size))]
      # Sampling
      set.seed(1)
      .i <- .i[, .(start= sample(region_start:region_end, sample_size, replace = T)), 
               .(seqnames, region_start= start, region_end= end, sample_size, total_counts)]
      # Order and extend to match ChIP reads length
      setorderv(.i, c("seqnames", "start"))
      set.seed(1)
      .i[, end:= start+sample(.c[, end-start], .N)]
      # Clean
      .i <- .i[, .(seqnames, start, end, total_counts)]
    }

  # Count overlapping reads
  regions <- copy(bins)
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
#' @param ChIP_bed Either a character vector poiting to bed file path OR a data.table containing ChIP reads
#' @param Input_bed Either a character vector poiting to bed file path OR a data.table containing Input reads. 
#' @param BSgenome A BSgenome object used for gw binning
#' @param binsize binsize used for peak calling, Default to 100 (narrow Peaks). Use larger bins to call domains
#' @examples 
#' 
#' @return significantly enriched peaks
#' @export

vl_peakCalling <- function(ChIP_bed,
                           Input_bed= NULL,
                           BSgenome,
                           binsize= 100)
{
  if(!is.character(ChIP_bed))
    stop("ChIP_bed and Input_bed should both be character vectors pointing to bed files")
  if(!is.null(Input_bed))
    if(!is.character(Input_bed))
      stop("Input_bed should either be set to NULL or a bed file path")
  if(class(BSgenome)[1]!="BSgenome")
    stop("BSgenome should be a BSgenome object")
  
  # Small bins
  bins <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  bins <- bins[, .(start = seq(1, end, binsize)), .(seqnames, end, width)]
  bins[, end:= start+(binsize-1)]
  bins[end > width, end:= width]
  bins <- bins[end - start > 0, .(seqnames, start, end)]
  
  #----------------------------#
  # Import reads
  #----------------------------#
  bed_import <- function(x)
  {
    .c <- rbindlist(lapply(x, 
                           fread, 
                           sel= 1:3, 
                           col.names= c("seqnames", "start", "end"), 
                           key= c("seqnames", "start")))
    .c[, total_counts:= .N, seqnames]
    return(.c)
  }
  # ChIP
  .c <- bed_import(ChIP_bed)
  # Input
  if(is.character(Input_bed))
    .i <- bed_import(Input_bed) else # If Input not provided, use random shuffling of ChIP reads
    {
      # Compute contigs 
      .i <- vl_collapse_DT_ranges(bed = .c, 
                                  compute_other_columns = "total_counts", 
                                  aggregate.fun = function(x) max(x, na.rm=T))
      # Compute effective chr sizes
      .i[, chr_size:= sum(end-start+1), seqnames]
      # Compute prob vector in order to constrain shuffle reads into mappable regions AND match ChIP read length distribution
      prob_vector <- data.table(read_length= seq(max(.i[,end-start+1])))
      prob_vector <- merge(prob_vector,
                           .c[, .N, .(read_length= end-start+1)],
                           all.x= T)
      prob_vector[is.na(N), N:= 0]
      # Sampling
      .i[, sample_size:= round(total_counts*((end-start+1)/chr_size))]
      set.seed(1)
      .i <- .i[, .(start= sample(x = region_start:region_end, 
                                 size = sample_size,
                                 prob = rev(cumsum(prob_vector[read_length<=(region_end-region_start+1), N])),
                                 replace = T)), 
               .(seqnames, region_start= start, region_end= end, sample_size, total_counts)]
      
      # Sample
      .i[1, prob:= .(.(cumsum(prob_vector[read_length<=width, N]))), .(width= end-start+1)]
      .i[1, prob]
      # Extend table to get 1/row per read (depending on contig size)
      set.seed(1)
      .i <- .i[, .SD[sample(x= seq(.N), 
                            size= total_counts, 
                            prob = (end-start+1)/chr_size,
                            replace= T)], .(seqnames, total_counts)]
      setkeyv(.i, c("seqnames", "start"))
      # Sample read lengths from ChIP
      set.seed(1)
      .i[, read_length:= sample(.c[, end-start], .N)]
      # Compute shuffled reads start and end
      .i[, end:= end-read_length]
      set.seed(1)
      .i[, start:= sample(start:end, .N), .(seqnames, start, end)]
      .i[, sample_size:= round(total_counts*((end-start+1)/chr_size))]
      .i <- .i[rep(seq(.N), sample_size)]
      # resize depending on reads sizes (random starts will be sampled!)
      .i[, end:= end-min(.c[, end-start]), seqnames]
      
      
      # Sampling
      set.seed(1)
      .i <- .i[, .(start= sample(region_start:region_end, sample_size, replace = T)), 
               .(seqnames, region_start= start, region_end= end, sample_size, total_counts)]
      # Order and extend to match ChIP reads length
      setorderv(.i, c("seqnames", "start"))
      set.seed(1)
      .i[, end:= start+sample(.c[, end-start], .N)]
      # Clean
      .i <- .i[, .(seqnames, start, end, total_counts)]
    }
  
  #----------------------------#
  # Compute bins enrichment
  #----------------------------#
  bins <- vl_bedEnrichment(ChIP_bed = .c, 
                           Input_bed = .i, 
                           bins = bins)
  
  #----------------------------#
  # Merge contiguous bins into peaks
  #----------------------------#
  # Only retain bins significantly enriched over background
  peaks <- bins[pval<1e-5]
  # Collapse touching bins into candidate peaks
  i <- 0
  peaks[, idx:= {
    i <<- i+1
    .idx <- c(i, sapply(.SD[-1, start]-.SD[-nrow(.SD), end], function(y) {
      if(y>1) 
        i <<- i+1
      return(i)
    }))
  }, seqnames]
  peaks <- peaks[, .(start= start[1], end= end[.N]), .(seqnames, idx)]
  # Compute overall peaks enrihment
  peaks <- vl_bedEnrichment(ChIP_bed = .c, 
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