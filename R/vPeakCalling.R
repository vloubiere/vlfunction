#' Peak calling
#'
#' Compute peaks from bed files
#'
#' @param ChIP_bed Vector of CHIP bed files path
#' @param Input_bed Vector of INPUT bed files path. default= NULL uses local enrichment. Otherwise, 1 input file can be specified or as many files as for ChIP
#' @param binsize size of the bins used to call peaks
#' @param BSgenome BSgenome object used for binning the data
#' @param collapse_touching_peaks If set to FALSE, return all bins with related padj and OR. Else, returns collapsed reads with max OR and -log10(padj) (Default).
#' @param collapse_mingap mingap for collapsing peaks. see ?vl_collapse_DT_ranges().
#' @param min_N_replicates minimum number of replicates required to retain peak. Default= all replicates
#' @examples 
#' ChIP_bed <- c("../available_data_dm3/db/bed//GSE119708_ATAC_rep1_uniq.bed", "../available_data_dm3/db/bed//GSE119708_ATAC_rep2_uniq.bed")
#' peaks <- vl_peakCalling(ChIP_bed)
#' 
#' @return peaks
#' @export

vl_peakCalling <- function(ChIP_bed, 
                           Input_bed= NULL, 
                           binsize= 100,
                           BSgenome= BSgenome.Dmelanogaster.UCSC.dm3,
                           collapse_touching_peaks= T,
                           collapse_mingap= 3*binsize+1,
                           min_N_replicates= length(ChIP_bed))
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
  bins <- bins[end - start > 0, .(seqnames, start, end)]
  
  #----------------------#
  # Count reads
  #----------------------#
  .q <- mclapply(c(ChIP_bed, Input_bed), function(x) 
  {
    .c <- fread(x)
    bins[, total_reads:= nrow(.c)]
    cbind(bins, counts= .c[bins, .N, .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]$N)
  })
  names(.q) <- c(paste0("CHIP_", seq(ChIP_bed)), paste0("INPUT_", seq(Input_bed)))
  .q <- rbindlist(.q, idcol = T)
  # Check that bins contain read for all INPUT replicates
  .q[, check:= all(counts[grep("^INPUT", .id)]>0), .(seqnames, start, end)]
  
  #----------------------#
  # Cast table
  #----------------------#
  res <- .q[(check), !"check"]
  res[, c(".id", "rep"):= tstrsplit(.id, "_")]
  res <- dcast(res, 
               seqnames+start+end+rep~.id, 
               value.var = c("counts", "total_reads"))
  res[, center:= round(rowMeans(.SD)), .SDcols= c("start", "end")]
  res <- na.omit(res)

  #----------------------#
  # Find significantly enriched bins
  #----------------------#
  # Compute average counts/bin/chromosome
  res[, input_average:= ceiling(mean(counts_INPUT)), .(seqnames, rep)]
  # Fisher test
  res[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(counts_CHIP, input_average, total_reads_CHIP, total_reads_INPUT)]
  res[, padj:= p.adjust(pval, method= "fdr")]
  # Check that N padj<0.05 is higher than min number of replicates required
  res[, check:= length(which(padj<0.05))>=min_N_replicates, .(seqnames, start, end)]

  #----------------------#
  # Collapse touching peaks
  #----------------------#
  # Filter peaks that pass thresholds
  coll <- vl_collapse_DT_ranges(res[(check)], mingap = collapse_mingap)
  # Compute Fisher on merged peaks
  setkeyv(coll, c("seqnames", "start", "end"))
  setkeyv(res, c("seqnames", "start", "end"))
  peaks <- foverlaps(res, 
                     coll, 
                     nomatch = NULL)
  peaks <- peaks[, .(counts_CHIP= sum(counts_CHIP),
                     input_average= sum(input_average)), 
                 .(seqnames, start, end, rep, total_reads_CHIP, total_reads_INPUT)]
  peaks[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(counts_CHIP, input_average, total_reads_CHIP, total_reads_INPUT)]
  peaks[, padj:= p.adjust(pval)]
  # Filter peaks that pass thresholds and export
  peaks[, check:= length(which(padj<0.05))>=min_N_replicates, .(seqnames, start, end)]
  
  return(peaks[(check), !"check"])
}



