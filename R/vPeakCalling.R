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

ChIP_bed <- c("../available_data_dm3/db/bed//GSE119708_ATAC_rep1_uniq.bed", "../available_data_dm3/db/bed//GSE119708_ATAC_rep2_uniq.bed")
vl_peakCalling <- function(ChIP_bed, 
                           Input_bed= NULL, 
                           binsize= 100,
                           Nbins_test= 100,
                           BSgenome= BSgenome.Dmelanogaster.UCSC.dm3,
                           collapse_touching_peaks= T)
{
  if (!class(BSgenome) == "BSgenome") 
    stop("genome should be a BSgenome object!")
  if(!length(Input_bed) %in% c(0, 1, length(ChIP_bed)))
    stop("Input_bed should either be set to null (local enrichment) or have a length of 1 (single input) or match ChIP_bed files length!")
  
  # Small bins
  bins <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  bins <- bins[, .(start = seq(1, end, binsize)), .(seqnames, end, width)]
  bins[, end:= start+(binsize-1)]
  bins[end > width, end:= width]
  bins <- bins[end - start > 0, .(seqnames, start, end, width)]
  bins[, bins_ID:= .I]
  
  # Count ChIP reads
  .q <- mclapply(ChIP_bed, function(x) 
    {
    .c <- fread(x)
    total_reads <- nrow(.c)
    .c <- .c[bins, .(bins_ID, .N), 
             .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
    return(.c[, .(seqnames= V1, 
                  bins_ID, 
                  CHIP= N, 
                  CHIP_total= total_reads)])
  })
  .q <- rbindlist(.q, idcol = T)
  
  # Control reads
  if(is.null(Input_bed))
  {
    .i <- .q
    colnames(.q)[3:4] <- c("CTL", "CTL_total")
  }else if(length(Input_bed)==length(ChIP_bed))
  {
    .i <- mclapply(Input_bed, function(x) 
    {
      .c <- fread(x)
      total_reads <- nrow(.c)
      .c <- .c[bins, .(bins_ID, .N), 
               .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
      return(.c[, .(seqnames= V1, 
                    CTL= N, 
                    CTL_total= total_reads)])
    })
    .i <- rbindlist(.i, idcol = T)
  }else if(length(Input_bed)==1)
  {
    .c <- fread(Input_bed)
    total_reads <- nrow(.c)
    .c <- .c[bins, .(bins_ID, .N), 
             .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
    .c <- .c[, .(seqnames= V1, 
                 CTL= N, 
                 CTL_total= total_reads)]
    .i <- data.table(.id= rep(seq(ChIP_bed), each= nrow(.c)),
                     seqnames= rep(.c$seqnames, length(ChIP_bed)),
                     CTL= rep(.c$CTL, length(ChIP_bed)),
                     CTL_total= rep(.c$CTL_total, length(ChIP_bed)))
  }
  
  # Compute counts
  Nbins <- ceiling(round(Nbins_test/2))
  obj <- split(.q, by= c(".id", "seqnames"))
  obj <- mclapply(obj, function(x) 
  {
    # Compute input slidin window idx
    x[, input_min_line:= ifelse((.I-Nbins)<1, 1, .I-Nbins)]
    x[, input_max_line := ifelse((.I+Nbins)>.N, .N, .I+50)]
    x <- x[rep(seq(nrow(x)), input_max_line-input_min_line+1)]
    x[, input_lines:= input_min_line:input_max_line, .(input_min_line, input_max_line)]
    # Compute Input mean count in larger window
    input <- .i[.id==x$.id[1] & seqnames==x$seqnames[1]]
    x[, c("CTL", "CTL_total"):= input[input_lines, .(CTL, CTL_total)]]
    x <- x[, .(CTL= ceiling(mean(CTL))), .(.id, bins_ID, CHIP, CHIP_total, CTL_total)]
    # Compute pval and OR
    x[, c("OR", "pval"):= {
      mat <- matrix(c(CHIP, CHIP_total-CHIP, CTL, CTL_total-CTL), nrow= 2, byrow = T)
      fisher.test(mat, alternative = "greater")[c("estimate", "p.value")]
    }, .(CHIP, CHIP_total, CTL, CTL_total)]
    return(x)
  })
  obj <- rbindlist(obj)
  
  # Padj cutoff + keep only bins called in all replicates
  obj[, padj:= p.adjust(pval, method = "fdr")]
  obj <- obj[padj<0.05, rep:= .N, bins_ID][rep==length(ChIP_bed)]
  
  # Compute mean enrichments and merge
  peaks <- obj[, .(OR= mean(OR), "-log10(padj)"= mean(-log10(padj))), bins_ID]
  peaks <- merge(bins, peaks)[, .(seqnames, start, end, OR, `-log10(padj)`)]
  if(collapse_touching_peaks)
  {
    coll <- vl_collapse_DT_ranges(peaks)[, !"strand"]
    peaks <- peaks[coll, .(OR= max(OR), `-log10(padj)`= max(`-log10(padj)`)), .EACHI, on= c("seqnames", "start<=end", "end>=start")]
  }
  return(peaks)
}



