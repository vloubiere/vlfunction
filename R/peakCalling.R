#' bw coverage integer
#'
#' IN DEVELOPMENT!! -> Quantify bw file at a set of intervals and returns integer values
#'
#' Note that strand-specific overlap is not implemented! bw files do not contain strand info!?
#'
#' @param bed Regions to quantify. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param bw Path to target bw file (character vector)
#' @param read.length read length, used when integer_counts is set to TRUE
#' @param integer_counts If set to TRUE, returns integer counts using approximated total read counts based on read.length
#' @examples 
#' # Bins
#' bins <- vl_binBSgenome(genome= "dm3",
#'                        restrict.seqnames = "chr3R")
#' bins <- bins[start>=10.3e6 & end<=10.4e6]
#' 
#' # Coverage
#' cov1 <- vl_bw_coverage(bins,
#'                        bw= "../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw")
#' cov2 <- vl_bw_coverage_integer(bins,
#'                                bw= "../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw",
#'                                read.length = 50)
#'                                
#' # plot
#' plot(cov1, type= "l")
#' plot(cov2, type= "l")
#' @export
vl_bw_coverage_integer <- function(bed, 
                                   bw,
                                   read.length= 50)
{
  var <- vl_importBed::copy(bed)
  var$score <- vl_bw_coverage(var, 
                              bw, 
                              na.value = 0)
  
  # As integer
  total_counts <- sum(var[, score*(end-start+1)])/read.length
  total_width <- sum(var[, (end-start+1)])
  var[, score:= round(score*total_width/total_counts)]
  
  return(var$score)
}

#' bw total reads
#'
#' IN DEVELOPMENT!! -> Estimate total reads within bw file
#'
#' @param bw Path to target bw file(s) (character vector)
#' @param read.length length of the reads used to generate bw file
#' @export
vl_bw_totalReads <- function(bw,
                             read.length)
{
  # Import bw(s)
  var <- rbindlist(lapply(bw, function(x) data.table::as.data.table(rtracklayer::import.bw(x))))
  res <- var[, .(total_counts= round(sum(score*width, na.rm= T)/read.length)), seqnames]
  
  return(res)
}

#' Compute ChIP enrichment
#'
#' IN DEVELOPMENT!! -> compared to enrichBed, uses bw files as Input_bw
#'
#' @param regions Regions to analyse. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param ChIP_bw ChIP bw files paths. If several provided, cat
#' @param Input_bw Input bw files paths. If several provided, cat
#' @param read.length Used to estimate total reads. default to 50
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export
vl_bw_enrich <- function(regions,
                         ChIP_bw,
                         Input_bw,
                         read.length= 50)
{
  if(!all(grepl(".bw$", c(ChIP_bw, Input_bw))))
    stop("ChIP_bw and Input_bw should be bw files")
  
  # Hard copy regions
  regions <- vl_importBed(regions)
  
  # bw coverage
  count_wrap <- function(bw)
  {
    .c <- lapply(bw, function(x) 
    {
      vl_bw_coverage_integer(regions, 
                             x, 
                             read.length = read.length)
    })
    .c <- do.call(cbind, .c)
    .t <- vl_bw_totalReads(bw, 
                           read.length = read.length)
    setkeyv(.t, "seqnames")
    data.table(rowSums(.c, na.rm= T),
               .t[as.character(regions$seqnames), total_counts])
  }
  regions[, c("ChIP_counts", "ChIP_total_counts"):= count_wrap(ChIP_bw)]
  regions[, c("Input_counts", "Input_total_counts"):= count_wrap(Input_bw)]
  
  # Compute enrichment and pval
  check <- regions[, ChIP_counts>0 & Input_counts>0] # Only regions containing reads
  regions[(check), c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, Input_counts, ChIP_total_counts, Input_total_counts)]
  regions[, padj:= p.adjust(pval, "fdr")]
  # Format narrowpeak file
  regions[pval==0, pval:= min(regions[pval>0, pval])]
  regions[padj==0, padj:= min(regions[padj>0, padj])]
  regions <- regions[, .(seqnames, start, end,
                         name= paste0("peak_", .I), 
                         score= round(OR/max(OR, na.rm = T)*1000), 
                         strand= ".", 
                         signalValue= OR,
                         pValue= -log10(pval),
                         qValue= -log10(padj),
                         peak= -1)]
  
  return(regions)
}

#' Compute ChIP enrichment
#'
#' IN DEVELOPMENT!! -> Given a set of ChIP and Input counts, performs hypergeometric test to know if ChIP is enriched compared to Input
#'
#' @param regions Regions to analyse. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param ChIP.bed ChIP bed file. Should be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input.bed Input bed file.
#' @return original bins file with OR, pval and padj corresponding to fisher result
#' @export
vl_enrichBed <- function(regions,
                         ChIP.bed,
                         Input.bed)
{
  # Hard copy regions
  regions <- vl_importBed(regions)
  regions <- copy(regions)
  # ChIP coverage
  ChIP.bed <- vl_importBed(ChIP.bed)
  regions[, ChIP_counts:= vl_covBed(regions, ChIP.bed)]
  total <- unique(ChIP.bed[, .N, seqnames])
  regions[total, ChIP_total_counts:= i.N, on= "seqnames"]
  # Input coverage
  Input.bed <- vl_importBed(Input.bed)
  regions[, Input_counts:= vl_covBed(regions, Input.bed)]
  total <- unique(Input.bed[, .N, seqnames])
  regions[total, Input_total_counts:= i.N, on= "seqnames"]
  # Compute enrichment and pval
  check <- regions[, ChIP_counts>0 & Input_counts>0] # Only regions containing reads
  regions[(check), c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, Input_counts, ChIP_total_counts, Input_total_counts)]
  regions[, padj:= p.adjust(pval, "fdr")]
  # Format narrowpeak file
  regions[pval==0, pval:= min(regions[pval>0, pval])]
  regions[padj==0, padj:= min(regions[padj>0, padj])]
  regions <- regions[, .(seqnames, start, end,
                         name= paste0("peak_", .I), 
                         score= round(OR/max(OR, na.rm = T)*1000), 
                         strand= ".", 
                         signalValue= OR,
                         pValue= -log10(pval),
                         qValue= -log10(padj),
                         peak= -1)]
  
  return(regions)
}

#' Compute peak calling
#'
#' IN DEVELOPMENT!! -> Computes peak calling using ChIP and Input
#'
#' @param ChIP ChIP bed file. Should be a vector of or bw file paths, a GRange object or a data.table containing 'seqnames', 'start', 'end' columns. see ?vl_importBed()
#' @param Input Input file (same format as ChIP).
#' @param read.length If ChIP and Input are bigwig files, this argument will be used to estimate total reads. default to 50
#' @param fun.average Function applied to ChIP read counts to define candidate peaks
#' @param gaussian.blur Applies gaussian blur to ChIP signal to identify peak candidates. Useful for noisy signal. default= FALSE
#' @param BSgenome A BSgenome object used for gw binning
#' @param bins.width bins width used for peak calling, Default to 100L (narrow Peaks). Use larger bins to call domains
#' @param steps.width Step between bins. default to round(bins.width/2)
#' @param bins.pval.cutoff The cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 0.05
#' @param bins.OR.cutoff Enrichment cutoff that will be applied to find candidate bins and merge them into peaks. Less stringent cutoffs means broader regions. default to 1
#' @return significantly enriched peaks
#' @export
vl_peakCalling <- function(ChIP,
                           Input,
                           read.length= 50,
                           fun.average= function(x) ceiling(mean(x)),
                           gaussian.blur= F,
                           BSgenome,
                           bins.width= 100,
                           steps.width= round(bins.width/2),
                           bins.pval.cutoff= 0.05,
                           bins.OR.cutoff= 1)
{
  bw <- all(grepl(".bw$", ChIP))
  if(bw && !all(grepl(".bw$", Input)))
    stop("ChIP are bw files while Input are bed! Stick to one type only")
  
  # Find potentially enriched bins gw and define peak candidates ----
  # genome wide bins
  bins <- vl_binBSgenome(BSgenome, 
                         bins.width = bins.width, 
                         steps.width = steps.width)
  # Compute ChIP coverage and average signal
  if(bw)
  {
    bins[, ChIP_counts:= rowSums(do.call(cbind, lapply(ChIP, function(x) vl_bw_coverage(bins, x))), na.rm= T)]
    bins[vl_bw_totalReads(ChIP, read.length), ChIP_total_counts:= i.total_counts, on= "seqnames"]
    bins[, ChIP_counts:= round(ChIP_counts/sum(ChIP_counts)*ChIP_total_counts), seqnames]
  }else{
    ChIP <- vl_importBed(ChIP)
    bins[, ChIP_counts:= vl_covBed(bins, ChIP)]
    bins <- merge(bins,
                  unique(ChIP[, .(ChIP_total_counts= .N), seqnames]))
  }
  # Smooth signal?
  if(gaussian.blur)
    bins[, ChIP_counts:= round(vl_gaussian_blur(ChIP_counts))]
  # Average signal will be used to identify candidates
  bins[, average_counts:= fun.average(ChIP_counts[ChIP_counts>0]), seqnames]
  bins[, average_total_counts := ChIP_total_counts, seqnames]
  # Remove bins that do not contain reads
  bins <- bins[ChIP_counts>0 & average_counts>0]
  # Enrichment
  bins[, c("OR", "pval"):= {
    mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
    fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
  }, .(ChIP_counts, average_counts, ChIP_total_counts, average_total_counts)]
  # Only retain bins significantly enriched over background
  peaks <- bins[pval <= bins.pval.cutoff & OR >= bins.OR.cutoff]
  # Collapse touching bins into candidate peaks
  peaks <- vl_collapseBed(peaks)
  
  # Test each candidate peaks and output significant ones ----
  if(bw)
  {
    final <- vl_bw_enrich(peaks, 
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
