#' bw Average tracks plot only
#'
#' Plots average tracks from a vl_average_bw_track output object
#'
#' @param obj An object returned by the vl_average_bw_track() function. The object can easily be modified to change colors and so on...
#' @param ylim ylim for plotting. default= range(data)
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param axes Should the x axis be plotted? Default to T
#' @param legend Should the legend be plotted? default to T
#' @examples 
#' Unstranded 
#' STARR <- fread("/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/peaks/DSCP_600bp_gw_cut_merged.peaks.txt", select = c(1,2))
#' STARR <- STARR[, .(seqnames= V1, start= V2, end= V2)]
#' obj1 <- vl_average_bw_track(bed= STARR,
#'                             set_IDs = c(rep("High", 2000), rep("Low", 2781)),
#'                             #' extend = c(-5000, 5000),
#'                             stranded = F,
#'                             tracks= "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_600bp_gw_cut_merged.bw", 
#'                             col= c("lightgrey", "red"), 
#'                             names = c("STARR-Seq"), 
#'                             plot=F)
#' vl_average_bw_track_plot_only(obj = obj1)
#'
#' Stranded 
#' prom <- readRDS("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/Analyses/Rdata/unique_proms_merged.df.RDS")
#' prom <- as.data.table(GRanges(unique(prom$Flybase.TSS)))
#' obj2 <- vl_average_bw_track(bed= prom,
#'                             extend = c(-5000, 5000),
#'                             stranded = T,
#'                             tracks= "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw",
#'                             plot= F)
#' vl_average_bw_track_plot_only(obj = obj2)
#' @export

vl_average_bw_track_plot_only <- function(obj, 
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
