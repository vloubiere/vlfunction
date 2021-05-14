#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Can be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns"
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param extend How much should the bed regions be extended (starting from the center). dEfault= c(-5000, 5000)
#' @param stranded Should the average track be stranded? If yes, - features are reversed :D
#' @param nbins Number of bins spanning the extended regions. Default= 501
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param plot Should the average track be ploted? default= T
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param ylim ylim for plotting. default= range(data)
#' @param col Colors to use. If specified, should match length(tracks)*length(unique(set_IDs)). default= NULL (redirects to rainbow)
#' @param axes Should the x axis be plotted? Default to T
#' @param legend Should the legend be plotted? default to T
#' @examples 
#' Unstranded 
#' STARR <- fread("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/peaks/DSCP_600bp_gw_cut_merged.peaks.txt", select = c(1,2))
#' STARR <- STARR[, .(seqnames= V1, start= V2, end= V2)]
#' obj1 <- vl_average_track(bed= STARR,
#'                          tracks= "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw", 
#'                          set_IDs = c(rep("High", 2000), rep("Low", 2781)),
#'                          #' extend = c(-5000, 5000),
#'                          stranded = F,
#'                          col= c("lightgrey", "red"), 
#'                          names = c("STARR-Seq"))
#'
#' Stranded 
#' prom <- readRDS("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/Analyses/Rdata/unique_proms_merged.df.RDS")
#' prom <- as.data.table(GRanges(unique(prom$Flybase.TSS)))
#' obj2 <- vl_average_track(bed= prom,
#' tracks= "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw",
#' extend = c(-5000, 5000),
#' stranded = T)
#' @return An object that can be used with the vl_average_track_plot_only() function.
#' @export

vl_average_track <- function(bed,
                             tracks,
                             set_IDs= NULL,
                             extend= c(-5000, 5000),
                             stranded= F,
                             nbins= 501, 
                             names= NULL,
                             plot= T,
                             xlab= "genomic distance",
                             ylab= "Enrichment",
                             ylim= NULL,
                             col= NULL, 
                             axes= T,
                             legend= T)
{
  if(class(bed)[1]=="GRanges")
    bed <- data.table::as.data.table(bed)
  if(!data.table::is.data.table(bed) | !all(c("seqnames", "start", "end") %in% colnames(bed)))
    stop("bed must be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns")
  if(is.null(set_IDs))
    set_IDs <- rep(1, nrow(bed))
  if(length(set_IDs)!=nrow(bed))
    stop("set_IDs length should match nrows(bed)!")
  if(stranded & !("strand" %in% colnames(bed)))
    stop("stranded= T but no 'strand' column in bed file")
  if(any(!file.exists(tracks)))
    stop("some track files could not be found! full paths prodvided?")
  if(any(!grepl(".bw$|.bedgraph$|.bedgraph.gz$|.bg$|.bg.gz$", tracks)))
    stop("some track files formats are not supported. Supported extensions: .bw, .bedgraph, .bedgraph.gz, .bg, .bg.gz")
  if(is.null(names))
    names <- gsub(".bw$", "", basename(tracks))
  if(length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length(", length(tracks), ")")
  if(is.null(col))
    col <- grDevices::rainbow(length(tracks)*length(unique(set_IDs)))
  if(length(col) != length(tracks)*length(unique(set_IDs)))
    stop("Colors vector length should be the same as bw files vector length * unique set_IDs (", 
         length(tracks)*length(unique(set_IDs)), ")")
  
  bins <- copy(bed)
  bins[, set_ID:= set_IDs]
  add <- seq(extend[1], extend[2], length.out = nbins+1)
  if(!stranded)
    bins[, strand:= "+"]
  bins[, center:= round(rowMeans(.SD)), .SDcols= c("start", "end")]
  bins[, region_ID:= .I]
  bins <- bins[, .(start= ceiling(center+add[-length(add)]),
                   end= floor(center+add[-1]), 
                   bin_ID= if(strand=="+") seq(nbins) else rev(seq(nbins))), 
                 .(set_ID, region_ID, seqnames, strand)]
  setkeyv(bins, c("seqnames", "start", "end"))
  
  #--------------------------#
  # Quantif tracks
  #--------------------------#
  sel <- GenomicRanges::GRanges(bins[, .(seqnames= seqnames[1], 
                                         start= min(start), 
                                         end= max(end)), .(region_ID)])
  sel <- rtracklayer::BigWigSelection(sel, "score")
  q <- parallel::mclapply(tracks, function(x)
  {
    if(grepl(".bw$", x))
      .c <- data.table::as.data.table(rtracklayer::import(x, selection= sel))
    if(grepl(".bedgraph$|.bedgraph.gz$|.bg$|.bg.gz$", x))
      .c <- data.table::fread(x, 
                              fill= T, 
                              col.names = c("seqnames", "start", "end", "score"))
    data.table::setkeyv(.c, c("seqnames", "start", "end"))
    res <- data.table::foverlaps(.c, bins, nomatch = 0)
    res <- res[, .(score= max(abs(score), na.rm= T)), .(set_ID, region_ID, bin_ID)]
    res <- res[, .(mean= mean(score, na.rm = T), se= sd(score, na.rm = T)/sqrt(length(score))), .(set_ID, bin_ID)]
  })
  names(q) <- tracks
  final <- data.table::rbindlist(q, idcol = "track")
  final[, Cc:= col[.GRP], .(track, set_ID)]
  final[, name:= names[.GRP], track]
  setorderv(final, "bin_ID")
  
  #--------------------------#
  # PLOT
  #--------------------------#
  if(plot)
  {
    if(is.null(ylim))
      ylim <- range(c(final[, mean-se], final[, mean+se]), na.rm = T)
    plot(NA, 
         xaxt= "n", 
         xlim = c(1, nbins), 
         ylim = ylim, 
         xlab= xlab,
         ylab= ylab)
    if(axes)
      axis(1, 
           at= seq(1, nbins, length.out = 3), 
           labels = seq(extend[1], extend[2], length.out = 3))
    
    final[, 
          {
            polygon(c(bin_ID, rev(bin_ID)), 
                    c(mean+se, rev(mean-se)), 
                    border= NA, 
                    col= adjustcolor(Cc, 0.5))
            lines(bin_ID, mean, 
                  col= Cc)
          }, .(track, set_ID, Cc)]
    
    if(legend)
    {
      .u <- unique(final[, .(name, set_ID, Cc)])
      if(length(unique(.u$set_ID))>1) # Only specify set_IDs if several sets are used!
        labels <- paste0(.u$name , " @ ", .u$set_ID) else 
          labels <- .u$name
      legend("topleft", 
             bty= "n", 
             fill = .u$Cc, 
             legend = labels)
    }
  }
  obj <- list(heatmap= final, 
              nbins= nbins,
              extend= extend)
  invisible(obj) 
}

#' bw Average tracks plot only
#'
#' Plots average tracks from a vl_average_track output object
#'
#' @param obj An object returned by the vl_average_track() function. The object can easily be modified to change colors and so on...
#' @param ylim ylim for plotting. default= range(data)
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param axes Should the x axis be plotted? Default to T
#' @param legend Should the legend be plotted? default to T
#' @examples 
#' Unstranded 
#' STARR <- fread("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/peaks/DSCP_600bp_gw_cut_merged.peaks.txt", select = c(1,2))
#' STARR <- STARR[, .(seqnames= V1, start= V2, end= V2)]
#' obj1 <- vl_average_track(bed= STARR,
#'                          set_IDs = c(rep("High", 2000), rep("Low", 2781)),
#'                          #' extend = c(-5000, 5000),
#'                          stranded = F,
#'                          tracks= "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw", 
#'                          col= c("lightgrey", "red"), 
#'                          names = c("STARR-Seq"), 
#'                          plot=F)
#' vl_average_track_plot_only(obj = obj1)
#'
#' Stranded 
#' prom <- readRDS("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/Analyses/Rdata/unique_proms_merged.df.RDS")
#' prom <- as.data.table(GRanges(unique(prom$Flybase.TSS)))
#' obj2 <- vl_average_track(bed= prom,
#'                          extend = c(-5000, 5000),
#'                          stranded = T,
#'                          tracks= "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw",
#'                          plot= F)
#' vl_average_track_plot_only(obj = obj2)
#' @export
vl_average_track_plot_only <- function(obj, 
                                       ylim= NULL,
                                       xlab= "genomic distance",
                                       ylab= "Enrichment",
                                       axes= T, 
                                       legend= T)
{
  if(is.null(ylim))
    ylim <- range(c(obj$heatmap[, mean-se], obj$heatmap[, mean+se]), na.rm = T)
  plot(NA, 
       xaxt= "n", 
       xlim = c(1, obj$nbins), 
       ylim = ylim, 
       xlab= xlab,
       ylab= ylab)
  if(axes)
    axis(1, 
         at= seq(1, obj$nbins, length.out = 3), 
         labels = seq(obj$extend[1], obj$extend[2], length.out = 3))
  
  obj$heatmap[, 
              {
                polygon(c(bin_ID, rev(bin_ID)), 
                        c(mean+se, rev(mean-se)), 
                        border= NA, 
                        col= adjustcolor(Cc, 0.5))
                lines(bin_ID, mean, 
                      col= Cc)
              }, .(track, set_ID, Cc)]

  
  if(legend)
  {
    .u <- unique(obj$heatmap[, .(name, set_ID, Cc)])
    if(length(unique(.u$set_ID))>1) # Only specify set_IDs if several sets are used!
      labels <- paste0(.u$name , " @ ", .u$set_ID) else 
        labels <- .u$name
      legend("topleft", 
             bty= "n", 
             fill = .u$Cc, 
             legend = labels)
  }
}
