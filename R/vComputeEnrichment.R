#' Compyte enrichment
#'
#' Compute peaks from bed files
#'
#' @param ChIP_bed Vector of CHIP bed files path
#' @param Input_bed Vector of INPUT bed files path. default= NULL uses local enrichment. Otherwise, 1 input file can be specified or as many files as for ChIP
#' @param peaks DT containing peaks
#' @param ext_peaks size extension to compute enrichment at peak
#' @param BSgenome BSgenome object used for binning the data
#' @examples 
#' ChIP_bed <- c("../available_data_dm3/db/bed//GSE119708_ATAC_rep1_uniq.bed", "../available_data_dm3/db/bed//GSE119708_ATAC_rep2_uniq.bed")
#' 
#' @return enrichment DT
#' @export

vl_computeEnrichment <- function(ChIP_bed, 
                                Input_bed= NULL,
                                peaks,
                                ext_peaks= 10000,
                                BSgenome= BSgenome.Dmelanogaster.UCSC.dm3)
{
  if(!class(BSgenome) == "BSgenome")
    stop("genome should be a BSgenome object!")
  if(!length(Input_bed) %in% c(0, 1, length(ChIP_bed)))
    stop("Input_bed should either be set to null (local enrichment) or have a length of 1 (single input) or match ChIP_bed files length!")
  if(is.null(Input_bed))
    Input_bed <- ChIP_bed
  if(length(ChIP_bed)>1 & length(Input_bed)==1)
    Input_bed <- rep(Input_bed, length(ChIP_bed))
  if(nrow(unique(peaks)) != nrow(peaks))
    stop("nrow(unique(peaks)) != nrow(peaks), make sure all peaks file lines are unique")

  # Count reads ChIP
  .q <- mclapply(ChIP_bed, function(x)
  {
    .c <- fread(x)
    total_reads <- nrow(.c)
    .c <- .c[peaks, .N, .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
    return(cbind(peaks, .c[, .(counts= N, 
                               total_counts= total_reads)]))
  })
  names(.q) <- paste0("CHIP_", seq(ChIP_bed))
  .q <- rbindlist(.q, idcol = T)
  
  # Count reads CTL
  chr <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome)))
  epeaks <- peaks[, .(seqnames, start= start-ext_peaks, end= end+ext_peaks)]
  epeaks[chr, width:= i.width, on= "seqnames"]
  epeaks[start<1, start:= 1]
  epeaks[end>width, end:= width]
  ratio <- (peaks$end-peaks$start+1)/(epeaks$end-epeaks$start+1)
  .i <- mclapply(Input_bed, function(x)
  {
    .c <- fread(x)
    total_reads <- nrow(.c)
    .c <- .c[epeaks, .N, .EACHI, on= c("V1==seqnames", "V2<=end", "V3>=start")]
    return(cbind(peaks, .c[, .(counts= ceiling(N*ratio), 
                               total_counts= total_reads)]))
  })
  names(.i) <- paste0("INPUT_", seq(Input_bed))
  .i <- rbindlist(.i, idcol = T)
  
  # Compute Enrichment
  res <- rbind(.q, .i)
  res <- dcast(res, seqnames+start+end~.id, value.var = list("counts", "total_counts"))
  res <- melt(na.omit(res), 
              id.vars = c("seqnames", "start", "end"), 
              measure.vars = patterns("^counts_CHIP", "^counts_INPUT", "^total_counts_CHIP", "^total_counts_INPUT"))
  res[, enr:= log2((value1+1)/value3)-log2((value2+1)/value4)]
  res <- res[, .(log2_enr= log2(mean(2^enr))), .(seqnames, start, end)]
  res <- res[peaks, , on =c("seqnames", "start", "end")]
  
  return(res)
}

