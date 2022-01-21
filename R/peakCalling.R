#' Compute peak calling
#'
#' Compute peak calling using ChIP and Input
#'
#' @param ChIP ChIP bed file. Should be a vector of or bw file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input Input file (same format as ChIP).
#' @param read_length If ChIP and Input are bigwig files, this argument will be used to estimate total reads. default to 50
#' @param average_function Function applied to ChIP read counts to define candidate peaks
#' @param gaussian_blur Applies gaussian blur to ChIP signal to identify peak candidates. Useful for noisy signal. default= FALSE
#' @param BSgenome A BSgenome object used for gw binning
#' @param bins_width bins width used for peak calling, Default to 100L (narrow Peaks). Use larger bins to call domains
#' @param steps_width Step between bins. default to round(bins_width/2)
#' @param bins_pval_cutoff The cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 0.05
#' @param bins_OR_cutoff Enrichment cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 1
#' @return significantly enriched peaks
#' @export

vl_peakCalling <- function(ChIP,
                           Input,
                           read_length= 50,
                           average_function= function(x) ceiling(mean(x)),
                           gaussian_blur= F,
                           BSgenome,
                           bins_width= 100,
                           steps_width= round(bins_width/2),
                           bins_pval_cutoff= 0.05,
                           bins_OR_cutoff= 1)
{
  bw <- all(grepl(".bw$", ChIP))
  if(bw && !all(grepl(".bw$", Input)))
    stop("ChIP are bw files while Input are bed! Stick to one type only")
  
  #----------------------------#
  # Find potentially enriched bins gw and define peak candidates
  #----------------------------#
  # genome wide bins
  bins <- vl_binBSgenome(BSgenome, 
                         bins_width = bins_width, 
                         steps_width = steps_width)
  # Compute ChIP coverage and average signal
  if(bw)
  {
    bins[, ChIP_counts:= rowSums(do.call(cbind, lapply(ChIP, function(x) vl_bw_coverage(bins, x))), na.rm= T)]
    bins[vl_bw_totalReads(ChIP, read_length), ChIP_total_counts:= i.total_counts, on= "seqnames"]
    bins[, ChIP_counts:= round(ChIP_counts/sum(ChIP_counts)*ChIP_total_counts), seqnames]
  }else{
    if(!vl_isDTranges(ChIP))
      ChIP <- vl_importBed(ChIP)
    bins[, ChIP_counts:= vl_covBed(bins, ChIP)]
    bins <- merge(bins,
                  unique(ChIP[, .(ChIP_total_counts= .N), seqnames]))
  }
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
  if(bw)
  {
    final <- vl_enrichBed(peaks, 
                          ChIP, 
                          Input)
  }else{
    final <- vl_enrichBed(peaks, 
                          ChIP, 
                          Input)
  }
  
  # Format and save significant peaks
  final <- final[qValue>(-log10(0.05))]
  
  return(final)
}
