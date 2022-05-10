#' bw coverage
#'
#' Quantify bw file at a set of intervals. returns mean value/region
#'
#' Note that strand-specific overlap is not implemented! bw files do not contain strand info!?
#'
#' @param bed Regions to quantify. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param bw Path to target bw file (character vector)
#' @param na_value Value to use for NAs. default to 0
#' @examples 
#' bins <- vl_binBSgenome("dm3", restrict_seqnames = "chr3R")
#' cov1 <- vl_bw_coverage(bins, "../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw")
#' 
#' # Compare to GRanges method
#' binned_average_function <- function(gr, file_bw_path){
#' file_bw <- import.bw(file_bw_path, as="RleList")
#' seqlevels(gr) <- names(file_bw)
#' bins1 <- binnedAverage(gr,file_bw, "average_score")
#' return(bins1$average_score)
#' }
#' cov2 <- binned_average_function(GRanges(bins), "../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw")
#' identical(cov1, cov2)
#' @export
vl_bw_coverage <- function(bed, 
                           bw,
                           na_value= 0)
{
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist! EXIT")
  
  # Format
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  .b <- data.table::copy(bed[, .(seqnames, start, end, .I)])

  # Import bw
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(na.omit(.b)), "score")
  var <- data.table::as.data.table(rtracklayer::import.bw(bw, selection= sel))
  
  # Overlap
  res <- var[.b, .(seqnames, x.start, x.end, score, i.start, i.end, I), on= c("seqnames", "start<=end", "end>=start")]
  
  # Clip bw ranges to bin ramges
  res[x.start<i.start, x.start:= i.start]
  res[x.end>i.end, x.end:= i.end]
  
  # Compute score
  res <- res[, sum(score*(x.end-x.start+1))/(i.end-i.start+1), .(seqnames, i.start, i.end, I)]$V1
  res[is.na(res)] <- na_value
  return(res)
}

#' Bins bed file and compute bw signal
#' 
#' Compared to vl_bw_coverage, this functions bins bed file before quantifying the signal. Useful to make heatmaps/average tracks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param set_IDs Set IDs specifying the groups as subsets of the bed file
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param stranded Should the average track be stranded?
#' @param nbins Number of bins spanning the extended regions. Default= 500
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @export
vl_bw_coverage_bins <- function(bed,
                                tracks,
                                set_IDs,
                                upstream, 
                                downstream,
                                stranded,
                                nbins, 
                                names)
{
  # Hard copy bed
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  # Checks
  if(!"strand" %in% names(bed))
    bed[, strand:= factor("*")]
  if(!identical(unique(tracks), tracks))
    stop("tracks should be unique")
  if(!identical(unique(names), names))
    stop("names should be unique")
  
  # Binning
  bins <- data.table(bed, set_IDs)
  bins[, region_ID:= .I]
  bins <- vl_resizeBed(bins, "center", upstream, downstream)
  bins <- bins[, {
    coor <- round(seq(start, end, length.out= nbins+1))
    coor <- data.table(start= coor[-length(coor)],
                       end= coor[-1])
    coor[-1, start:= start+1]
  }, .(seqnames,
       region_ID,
       set_IDs,
       strand)]
  bins[, bin.x:= rowid(region_ID)]
  if(stranded)
    bins[as.character(strand)=="-", bin.x:= nbins-bin.x+1]
  
  #--------------------------#
  # Quantif tracks
  #--------------------------#
  obj <- parallel::mclapply(tracks, function(x) {
    data.table(bins[, .(file= x, 
                        set_IDs, 
                        region_ID, 
                        bin.x)], 
               score= vl_bw_coverage(bins, x))
  }, mc.preschedule = T)
  names(obj) <- names
  obj <- rbindlist(obj, idcol= "name")
  return(obj)
}

#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param set_IDs Set IDs specifying the groups as subsets of the bed file
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param stranded Should the average track be stranded?
#' @param nbins Number of bins spanning the extended regions. Default= 500
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param plot Should the average track be ploted? default= T
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param ylim ylim for plotting. default= range(data)
#' @param col Color to be used for plotting. 
#' @param legend Should the legeng be plotted? default to T
#' @param legend.cex Legend cex. defaults to 1
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' vl_bw_average_track(bed, tracks, plot= T, upstream = 1000, downstream = 1000, set_IDs = sets)
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export
vl_bw_average_track <- function(bed,
                                tracks,
                                set_IDs= 1,
                                upstream= 5000,
                                downstream= 5000,
                                stranded= F,
                                nbins= 100, 
                                names= gsub(".bw$", "", basename(tracks)),
                                center_label= "Center",
                                plot= T,
                                xlab= "genomic distance",
                                ylab= "Enrichment",
                                ylim,
                                col= c("#E69F00","#68B1CB","#15A390","#96C954","#77AB7A","#4F6A6F","#D26429","#C57DA5","#999999"),
                                legend= T,
                                legend.cex= 1)
{
  obj <- vl_bw_coverage_bins(bed= bed,
                             tracks= tracks,
                             set_IDs= set_IDs,
                             upstream= upstream, 
                             downstream= downstream,
                             stranded= stranded,
                             nbins= nbins, 
                             names= names)
  obj[, col:= colorRampPalette(col)(.NGRP)[.GRP], keyby= .(name, set_IDs)]
  setattr(obj, "class", c("vl_bw_average_track", "data.table", "data.frame"))
  if(plot)
    plot.vl_bw_average_track(obj,
                             xlab= xlab,
                             xaxis= c(upstream, center_label, downstream),
                             ylab= ylab,
                             ylim= ylim,
                             legend= legend,
                             legend.cex= legend.cex)
  invisible(obj)
}

#' @describeIn vl_average_bw_track Method to plot average tracks
#' @export
plot.vl_bw_average_track <- function(obj,
                                     xlab= "genomic distance",
                                     xaxis= c("Upstream", "Center", "Downstream"),
                                     ylab= "Enrichment",
                                     ylim,
                                     legend= T,
                                     legend.cex= 1)
{
  pl <- obj[, .(mean= mean(score, na.rm= T), 
                se= sd(score, na.rm= T)/sqrt(.N)), .(name, col, set_IDs, bin.x)]
  plot(NA, 
       xlim= range(pl$bin.x),
       ylim= if(missing(ylim)) range(c(pl[, mean-se], pl[, mean+se])) else ylim,
       ylab= ylab,
       xlab= xlab,
       xaxt= "n")
  axis(1, 
       c(1, max(obj$bin.x)/2, max(obj$bin.x)),
       labels= xaxis)
  pl[, {
    polygon(c(bin.x, rev(bin.x)), 
            c(mean+se, rev(mean-se)),
            border= NA,
            col= adjustcolor(col[1], 0.5))
    lines(bin.x, mean, col= col[1])
  }, .(name, set_IDs, col)]
  # Legend
  if(legend)
  {
    leg <- unique(pl[, .(set_IDs, name, col)])
    if(length(unique(leg$set_IDs))>1 & length(unique(leg$name))>1)
      leg[, labels:= paste0(name, " @ ", set_IDs)] else if(length(unique(leg$set_IDs))>1)
        leg[, labels:= set_IDs] else if(length(unique(leg$name))>1)
          leg[, labels:= name]
    legend("topleft",
           legend= leg$labels,
           fill= leg$col,
           bty= "n",
           cex= legend.cex)
  }
}

#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param set_IDs Set IDs specifying the groups as subsets of the bed file
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param stranded Should the average track be stranded?
#' @param nbins Number of bins spanning the extended regions. Default= 100
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param plot Should the average track be ploted? default= T
#' @param col Vector of colors to be used for heatmap
#' @param order_cols Numeric vector giving the indexes of the tracks to be used for ordering
#' @param order_FUN Function used to aggregated per region and order heatmap. default= function(x) mean(x, na.rm= T)
#' @param max_FUN Function to be used for clipping. default= function(x) quantile(x, 0.995, na.rm= T). Using max -> no clipping
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' set_IDs <- c(rep("suhw", 100), rep("STARR", 1000))
#' test <- vl_bw_heatmap(bed, tracks, set_IDs= set_IDs, plot= T, upstream = 1000, downstream = 1000, order_FUN = mean, order_cols = 2)
#' plot(test)
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export
vl_bw_heatmap <- function(bed,
                          tracks,
                          set_IDs= 1,
                          upstream= 5000,
                          downstream= 5000,
                          stranded= F,
                          nbins= 100, 
                          names= gsub(".bw$", "", basename(tracks)),
                          plot= T,
                          col= c("blue", "yellow"),
                          order_cols= 1,
                          order_FUN= function(x) mean(x, na.rm= T),
                          max_FUN= function(x) quantile(x, 0.995, na.rm= T),
                          na_col= "lightgrey")
{
  obj <- vl_bw_coverage_bins(bed= bed,
                             tracks= tracks,
                             set_IDs= set_IDs,
                             upstream= upstream, 
                             downstream= downstream,
                             stranded= stranded,
                             nbins= nbins, 
                             names= names)
  setattr(obj, "class", c("vl_bw_heatmap", "data.table", "data.frame"))
  if(plot)
    plot.vl_bw_heatmap(obj,
                       col= col,
                       order_cols= order_cols,
                       max_FUN= max_FUN,
                       order_FUN= order_FUN,
                       na_col= na_col)
  invisible(obj)
}

#' @export
plot.vl_bw_heatmap <- function(obj, 
                               col= c("blue", "yellow"),
                               order_cols= 1,
                               order_FUN= function(x) mean(x, na.rm= T),
                               max_FUN= function(x) quantile(x, 0.995, na.rm= T),
                               na_col= "lightgrey")
{
  # Format 
  if(!is.factor(obj$name))
    obj[, name:= factor(name, unique(as.character(name)))]
  if(!is.factor(obj$set_IDs))
    obj[, set_IDs:= factor(set_IDs, unique(as.character(set_IDs)))]
  # Reorder region_ID if order cols specified
  if(is.numeric(order_cols))
  {
    order_cols <- levels(obj$name)[order_cols]
    ord <- dcast(obj[name %in% order_cols], 
                 region_ID~name, 
                 value.var = "score", 
                 fun.aggregate = order_FUN)
    setorderv(ord, order_cols, -1)
    ord[, order:= .I]
    obj[ord, region_ID:= i.order, on= "region_ID"]
  }
  # Clip outliers
  obj[, max:= max_FUN(score), name]
  obj[score>max, score:= max]
  # Dcast image
  dmat <- dcast(obj, set_IDs+region_ID~name+bin.x, value.var = "score")
  im <- mat <- as.matrix(dmat[, !c("set_IDs", "region_ID")])
  # Plotting parameters
  Cc <- circlize::colorRamp2(range(im, na.rm= T), col)
  im[!is.na(im)] <- Cc(im[!is.na(im)])
  im[is.na(im)] <- na_col
  Nbins <- max(obj$bin.x)
  track.names.x <- seq(1, ncol(mat)-Nbins, length.out= length(levels(obj$name)))+Nbins/2
  Sets.y <- c(0,cumsum(table(unique(obj[, .(set_IDs, region_ID)])$set_IDs)))
  if(length(Sets.y)>2)
  {
    Sets.lines.y <- nrow(mat)-Sets.y[-c(1, length(Sets.y))]
    Sets.names.y <- nrow(mat)-(Sets.y[-1]-diff(Sets.y)/2)
  }
  
  #--------------------------------#
  # PLOT
  #--------------------------------#
  plot.new()
  plot.window(xlim= c(0.5, ncol(mat)+0.5),
              ylim= c(0.5, nrow(mat)+0.5))
  rasterImage(xleft = 1, 
              ybottom = 1, 
              xright = ncol(im),
              ytop = nrow(im),
              im)
  abline(v= seq(0, ncol(im), Nbins), col= "white")
  rect(1,1,ncol(mat), nrow(mat), xpd= T)
  text(track.names.x,
       par("usr")[4], 
       levels(obj$name),
       xpd= T)
  if(length(Sets.y)>2)
  {
    segments(1, Sets.lines.y, ncol(mat), Sets.lines.y, xpd= T)
    text(0,
         Sets.names.y,
         levels(obj$set_IDs),
         pos= 2,
         xpd= T)
  }
}