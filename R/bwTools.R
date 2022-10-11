#' Merge several bigwig files into 1
#' 
#' Convenience function to merge directly into r
#'
#' @param bw A character vector containing the path of bw files to be merged
#' @param genome BSgenome. "dm3", "dm6", "mm10"
#' @param bins_width Bin size. default= 25L
#' @param output Output bw file path
#' @param scoreFUN A function to be applied to the score column before return
#' @export
vl_bw_merge <- function(tracks, genome, bins_width= 25L, output, scoreFUN= NULL)
{
  bw <- vl_binBSgenome(genome, bins_width = bins_width)
  bw[, score:= as.numeric(0)]
  # Compute value
  for(track in tracks)
  {
    bw[, score:= score+vl_bw_coverage(bw, track)]
    print(track)
  }
  # Apply function to score if specified
  if(!is.null(scoreFUN))
    bw[, score:= scoreFUN(score)]
  # Collapse score rle
  bw <- bw[, .(start= start[1], end= end[.N]), .(seqnames, data.table::rleid(score), score)]
  # Format GRanges
  .g <- GenomicRanges::GRanges(bw)
  BS <- getBSgenome(genome)
  GenomeInfoDb::seqlevels(.g) <- GenomeInfoDb::seqlevels(BS)
  GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BS)
  # save
  rtracklayer::export.bw(.g, 
                         con= output)
}

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
#' track <- "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw"
#' bins <- vl_binBSgenome("dm3", restrict_seqnames = "chr3R")
#' t1 <- Sys.time()
#' cov1 <- vl_bw_coverage(bins, track)
#' t1-Sys.time()
#' 
#' # Compare to GRanges method
#' binned_average_function <- function(gr, file_bw_path){
#' file_bw <- import.bw(file_bw_path, as="RleList")
#' seqlevels(gr) <- names(file_bw)
#' t1 <- Sys.time()
#' bins1 <- binnedAverage(gr,file_bw, "average_score")
#' return(bins1$average_score)
#' }
#' t1 <- Sys.time()
#' cov2 <- binned_average_function(GRanges(bins), track)
#' t1-Sys.time()
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
  .b <- vl_importBed(bed)[, .(seqnames, start, end, .I)]
  
  # Import bw
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(na.omit(.b)), "score")
  var <- data.table::as.data.table(rtracklayer::import.bw(bw, selection= sel))
  
  # Overlap
  res <- var[.b, .(seqnames, x.start, x.end, score, i.start, i.end, I), on= c("seqnames", "start<=end", "end>=start")]
  
  # Clip bw ranges to bin ramges
  res[x.start<i.start, x.start:= i.start]
  res[x.end>i.end, x.end:= i.end]
  
  # Compute score
  res <- res[, sum(score*(x.end-x.start+1))/width, .(I, width= i.end-i.start+1)]$V1
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
  regions <- vl_importBed(bed)
  # Checks
  if(stranded && !"strand" %in% names(bed))
  {
    message("stranded= TRUE but no strand was provided for 'bed' -> set to unstranded (*)")
    bed[, strand:= "*"]
  }
  if(!identical(unique(tracks), tracks))
    stop("tracks should be unique")
  if(!identical(unique(names), names))
    stop("names should be unique")
  
  # Binning
  bins <- data.table(regions, set_IDs)
  bins[, region_ID:= .I]
  bins[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
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
  bins[, bin.x:= seq(-upstream, downstream, length.out= nbins)[rowid(region_ID)]]
  if(stranded)
    bins[strand=="-", bin.x:= rev(bin.x), region_ID]
  
  #--------------------------#
  # Quantif tracks
  #--------------------------#
  obj <- parallel::mclapply(tracks, function(x) {
    data.table(bins[, .(set_IDs, region_ID, bin.x)], 
               score= vl_bw_coverage(bins, x))
  }, mc.preschedule = T)
  obj <- rbindlist(obj, idcol= "name")
  obj[, name:= names[name]]
  return(obj)
}

#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param names Track names to plot, further used for ordering. By default, bw basenames will be used (as factors).
#' @param set_IDs Set IDs specifying the groups as subsets of the bed file. Used for ordering.
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param stranded Should the average track be stranded?
#' @param nbins Number of bins spanning the extended regions. Default= 101L
#' @param plot Should the average track be ploted? default= T
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param ylim ylim for plotting. default= range(data)
#' @param col Color to be used for plotting. Dataset is ordered using keyby= .(names, set_IDs) before assigning colors
#' @param legend Should the legend be plotted? default to T
#' @param legend.cex Legend cex. defaults to 1
#' @param col.adj Opacity of polygons and lines. default= c(0.5,1)
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' vl_bw_average_track(bed, tracks, plot= T, upstream = 1000, downstream = 1000, set_IDs = sets)
#' @export
vl_bw_average_track <- function(bed,
                                names,
                                set_IDs,
                                tracks,
                                upstream= 5000,
                                downstream= 5000,
                                stranded= F,
                                nbins= 101L, 
                                center_label= "Center",
                                plot= T,
                                xlab= "genomic distance",
                                ylab= "Enrichment",
                                ylim,
                                col= c("#E69F00","#68B1CB","#15A390","#96C954","#77AB7A","#4F6A6F","#D26429","#C57DA5","#999999"),
                                legend= T,
                                legend.cex= 1,
                                col.adj= c(0.5, 1))
{
  # By default, preserve order bw tracks as specified in input
  if(missing(names))
  {
    names <- gsub(".bw$", "", basename(tracks))
    names <- factor(names, levels= unique(names))
  }
  if(missing(set_IDs))
    set_IDs <- 1
  obj <- vl_bw_coverage_bins(bed= bed,
                             tracks= tracks,
                             set_IDs= set_IDs,
                             upstream= upstream, 
                             downstream= downstream,
                             stranded= stranded,
                             nbins= nbins, 
                             names= names)
  obj[, col:= colorRampPalette(col)(.NGRP)[.GRP], keyby= .(name, set_IDs)]
  setattr(obj, 
          "class", 
          c("vl_bw_average_track", "data.table", "data.frame"))
  if(plot)
    plot.vl_bw_average_track(obj,
                             xlab= xlab,
                             xlab.at= c(-upstream, 0, downstream),
                             center_label= center_label,
                             ylab= ylab,
                             ylim= ylim,
                             legend= legend,
                             legend.cex= legend.cex,
                             col.adj= col.adj)
  invisible(obj)
}

#' @describeIn vl_bw_average_track Method to plot average tracks
#' @export
plot.vl_bw_average_track <- function(obj,
                                     xlab= "genomic distance",
                                     xlab.at= c(min(obj$bin.x), 0, max(obj$bin.x)),
                                     center_label= "Center",
                                     ylab= "Enrichment",
                                     ylim,
                                     legend= T,
                                     legend_pos= "topleft",
                                     legend.cex= 1,
                                     col.adj= c(0.5, 1))
{
  pl <- obj[, .(mean= mean(score, na.rm= T), 
                se= sd(score, na.rm= T)/sqrt(.N)), .(name, col, set_IDs, bin.x)]
  if(missing(ylim))
    ylim <- range(c(pl[, mean-se], pl[, mean+se]))
  plot(NA, 
       xlim= range(pl$bin.x),
       ylim= ylim,
       ylab= ylab,
       xlab= xlab,
       xaxt= "n")
  axis(1, 
       at = xlab.at,
       labels = c(xlab.at[1], center_label, xlab.at[3]))
  pl[, {
    polygon(c(bin.x, rev(bin.x)), 
            c(mean+se, rev(mean-se)),
            border= NA,
            col= adjustcolor(col, col.adj[1]))
    lines(bin.x, 
          mean, 
          col= adjustcolor(col, col.adj[2]))
  }, .(name, set_IDs, col)]
  # Legend
  if(legend)
  {
    leg <- unique(pl[, .(set_IDs, name, col)])
    if(length(unique(leg$set_IDs))>1 & length(unique(leg$name))>1)
      leg[, labels:= paste0(name, " @ ", set_IDs)] else if(length(unique(leg$set_IDs))>1)
        leg[, labels:= set_IDs] else if(length(unique(leg$name))>1)
          leg[, labels:= name]
    if("labels" %in% names(leg))
      legend(legend_pos,
             legend= leg$labels,
             fill= leg$col,
             bty= "n",
             cex= legend.cex)
  }
}

#' bw Average heatmap
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
#' tracks <- c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", 
#' "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' set_IDs <- c(rep("suhw", 100), rep("STARR", 1000))
#' test <- vl_bw_heatmap(bed, tracks, set_IDs= set_IDs, plot= T, upstream = 1000, downstream = 1000, order_FUN = mean, order_cols = 2)
#' plot(test, order_col= 2)
#' plot(test)
#' @export
vl_bw_heatmap <- function(bed,
                          tracks,
                          set_IDs= 1,
                          upstream= 5000,
                          downstream= 5000,
                          stranded= F,
                          nbins= 101L, 
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
  setattr(obj, 
          "class", 
          c("vl_bw_heatmap", "data.table", "data.frame"))
  if(plot)
    plot.vl_bw_heatmap(obj,
                       col= col,
                       order_cols= order_cols,
                       max_FUN= max_FUN,
                       order_FUN= order_FUN,
                       na_col= na_col)
  invisible(obj)
}

#' @describeIn vl_bw_heatmap Method to plot bw heatmaps
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