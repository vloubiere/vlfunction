#' Peak calling
#'
#' Compute peaks from bed files
#'
#' @param ChIP_bed Vector of CHIP bed files path
#' @param Input_bed Vector of INPUT bed files path. default= NULL uses local enrichment. Otherwise, 1 input file can be specified or as many files as for ChIP
#' @param binsize size of the bins used to call peaks
#' @param Nbins_test Number of INPUT bins used to average INPUT signal and perform Fishwer
#' @param BSgenome BSgenome object used for binning the data
#' @param collapse_touching_peaks If set to FALSE, return all bins with related padj and OR. Else, returns collapsed reads with max OR and -log10(padj) (Default).
#' @param min_N_replicates minimum number of replicates required to retain peak. Default= all replicates
#' @param cutoff_input_counts minimum number of input reads for the bin to be considered for testing. default= 5
#' @examples 
#' ChIP_bed <- c("../available_data_dm3/db/bed//GSE119708_ATAC_rep1_uniq.bed", "../available_data_dm3/db/bed//GSE119708_ATAC_rep2_uniq.bed")
#' peaks <- vl_peakCalling(ChIP_bed)
#' 
#' @return peaks
#' @export

vl_peakCalling <- function(ChIP_bed, 
                           Input_bed= NULL, 
                           binsize= 100,
                           Nbins_test= 100,
                           BSgenome= BSgenome.Dmelanogaster.UCSC.dm3,
                           collapse_touching_peaks= T,
                           min_N_replicates= length(ChIP_bed),
                           cutoff_input_counts= 5)
{
  if (!class(BSgenome) == "BSgenome") 
    stop("genome should be a BSgenome object!")
  if(!length(Input_bed) %in% c(0, 1, length(ChIP_bed)))
    stop("Input_bed should either be set to null (local enrichment) or have a length of 1 (single input) or match ChIP_bed files length!")
  if(is.null(Input_bed))
    Input_bed <- ChIP_bed
  if(length(ChIP_bed)>1 & length(Input_bed)==1)
    Input_bed <- rep(Input_bed, length(ChIP_bed))

  # Small bins
  bins <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  bins <- bins[, .(start = seq(1, end, binsize)), .(seqnames, end, width)]
  bins[, end:= start+(binsize-1)]
  bins[end > width, end:= width]
  bins <- bins[end - start > 0, .(seqnames, start, end, width)]
  bins[, bins_ID:= .I]
  
  #----------------------#
  # Count reads
  #----------------------#
  .q <- mclapply(c(ChIP_bed, Input_bed), function(x) 
  {
    .c <- fread(x)
    total_reads <- nrow(.c)
    .c <- .c[bins, .(bins_ID, .N), .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
    return(.c[, .(seqnames= V1, 
                  bins_ID, 
                  counts= N, 
                  total_counts= total_reads)])
  })
  names(.q) <- c(paste0("CHIP_", seq(ChIP_bed)), paste0("INPUT_", seq(Input_bed)))
  .q <- rbindlist(.q, idcol = T)
  
  #----------------------#
  # Compute Enrichment
  #----------------------#
  res <- dcast(.q, bins_ID+seqnames~.id, value.var = list("counts", "total_counts"))
  mav <- function(x,n){filter(x,rep(1/n,n), sides= 2)} # Rolling average function
  cols <- grep("^counts_INPUT_", colnames(res), value = T)
  # Cutoff minimum number input reads (all bins within +/- 5)
  res[, check_input:= { 
    check <- c(rep(0, 5), rowSums(.SD), rep(0, 5))
    rowSums(sapply(1:10, function(x) check[x:(x+.N-1)]))>=(10*cutoff_input_counts)
  }, .SDcols= patterns("^counts_INPUT")] 
  res[, (cols):= mclapply(.SD, function(x) ceiling(mav(x, Nbins_test+1))), seqnames, .SDcols= cols]
  res <- res[(check_input), !"check_input"] # Apply input cutoff check
  res <- melt(na.omit(res), 
              id.vars = "bins_ID", 
              measure.vars = patterns("^counts_CHIP", "^counts_INPUT", "^total_counts_CHIP", "^total_counts_INPUT"))
  res[, c("OR", "pval"):= {
    mat <- matrix(c(value1, value3-value1, value2, value4-value2), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, value1:value4]
  
  # Padj cutoff + keep only bins called in X replicates
  res[, padj:= p.adjust(pval, method = "fdr")]
  res <- res[padj<0.05, rep:= .N, bins_ID][rep==min_N_replicates]
  
  # Compute mean enrichments and merge
  peaks <- res[, .(OR= mean(OR), "-log10(padj)"= mean(-log10(padj))), bins_ID]
  peaks <- merge(bins, peaks)[, .(seqnames, start, end, OR, `-log10(padj)`)]
  peaks[, center:= round(rowMeans(.SD)), .SDcols= c("start", "end")]
  if(collapse_touching_peaks)
  {
    coll <- vl_collapse_DT_ranges(peaks)[, !"strand"]
    peaks <- peaks[coll, .(OR= max(OR), 
                           `-log10(padj)`= max(`-log10(padj)`), 
                           max_coor= center[which.max(OR)]), .EACHI, on= c("seqnames", "start<=end", "end>=start")]
    peaks <- cbind(coll, peaks[, .(OR, `-log10(padj)`, max_coor)])
  }
  return(peaks)
}



