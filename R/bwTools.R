#' bw coverage
#'
#' Quantify bw file at a set of intervals. returns mean value/region
#'
#' Note that strand-specific overlap is not implemented! bw files do not contain strand info!?
#'
#' @param bed Regions to quantify. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param bw Path to target bw file (character vector)
#' @export

vl_bw_coverage <- function(bed, 
                           bw)
{
  if(length(bw) != 1)
    stop("length(bw) != 1")
  if(!file.exists(bw))
    stop("bw file does not exist! EXIT")
  
  # Format
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  .b <- data.table::copy(bed)
  .b[, .ID:= .I]
  
  # Import bw
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(.b), "score")
  var <- data.table::as.data.table(rtracklayer::import.bw(bw, selection= sel))
  keys <- c("seqnames", "start", "end")
  data.table::setkeyv(.b, keys)
  data.table::setkeyv(var, keys)
  
  # Compute counts
  ov <- data.table::foverlaps(var, .b)
  ov[i.start<start, i.start:= start]
  ov[i.end>end, i.end:= end]
  ov[, width:= i.end-i.start+1]
  res <- ov[, .(score= sum(score*width)/(end[1]-start[1]+1)), .ID]
  res <- res[.(seq(nrow(.b))), score, on= ".ID"]
  return(res)
}


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


#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
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
#' @param legend Should the legen be plotted? default to T
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
#'                             names = c("STARR-Seq"))
#'
#' Stranded 
#' prom <- readRDS("/groups/stark/hendy/Projects/Chromatin_remodeler_specificity/Analyses/Rdata/unique_proms_merged.df.RDS")
#' prom <- as.data.table(GRanges(unique(prom$Flybase.TSS)))
#' obj2 <- vl_average_bw_track(bed= prom,
#'                             extend = c(-5000, 5000),
#'                             stranded = T,
#'                             tracks= "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw")
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export

vl_average_bw_track <- function(bed,
                                set_IDs= NULL,
                                tracks,
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
  # Hard copy bed
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  bins <- data.table::copy(bed)
  if(is.null(set_IDs))
    set_IDs <- rep(1, nrow(bins))
  if(length(set_IDs)!=nrow(bins))
    stop("set_IDs length should match nrows(bed)!")
  if(stranded & !("strand" %in% colnames(bins)))
    stop("stranded= T but no 'strand' column in bed file")
  if(any(!file.exists(tracks)))
    stop("some bw files could not be found! full paths prodvided?")
  if(is.null(names))
    names <- gsub(".bw$", "", basename(tracks))
  if(length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length(", length(tracks), ")")
  if(is.null(col))
    col <- grDevices::rainbow(length(tracks)*length(unique(set_IDs)))
  if(length(col) != length(tracks)*length(unique(set_IDs)))
    stop("Colors vector length should be the same as bw files vector length * unique set_IDs (", 
         length(tracks)*length(unique(set_IDs)), ")")
  
  # Format bins
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
    .c <- data.table::as.data.table(rtracklayer::import.bw(x, selection= sel))
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
  
  #-------------------------#
  # Make object
  #-------------------------#
  obj <- list(heatmap= final, 
              nbins= nbins,
              extend= extend)
  
  #--------------------------#
  # PLOT
  #--------------------------#
  if(plot)
    vl_average_bw_track_plot_only(obj, 
                                  ylim = ylim,
                                  xlab= xlab,
                                  ylab= ylab, 
                                  axes= axes, 
                                  legend= legend)
  invisible(obj) 
}


#' bw heatmap plot only
#'
#' Plots heatmap from a vl_heatmap_bw_track output object
#'
#' @param obj An object returned by the vl_heatmap_bw_track() function. The object can easily be modified to change colors and so on...
#' @param center_label Label for the center of the heatmaps
#' @param col Color vector to be used for heatmaps. col= c("blue", "yellow", "white")
#' @export

vl_heatmap_bw_track_plot_only <- function(obj,
                                          center_label= "TSS",
                                          col= c("blue", "yellow", "white"))
{
  Ntracks <- length(unique(obj$track))
  
  # Compute plotting colors
  Cc <- colorRampPalette(col)(101)
  obj[, plot_Cc:= Cc[.GRP], keyby= col_idx]
  
  # Compute image
  im <- dcast(obj, 
              region_order~track+bin_ID, 
              value.var = "plot_Cc")
  mat <- as.matrix(im, 
                   1)
  par(xaxs= "i",
      yaxs= "i",
      mai= c(1.02, 
             max(strwidth(unique(obj$set_ID), "inches"))+0.5, 
             0.82, 
             0.42))
  
  # Plot heatmap
  plot.new()
  rasterImage(mat, 0,0,1,1, interpolate = F)
  
  # Plot lines
  .lv <- seq(0, 1, length.out = length(unique(obj$track))+1)
  abline(v= .lv, 
         lwd= 0.25)
  .lh <- obj[, 1-(min(region_order)-1)/max(obj$region_order), set_ID]$V1
  .lh <- c(.lh, 0)
  abline(h= .lh, 
         lwd= 0.25)
  
  # Plot axes
  axAt <- .lv[-length(.lv)]+diff(.lv)/2
  axis(1, 
       at= axAt, 
       labels = rep(center_label, Ntracks),
       lwd= 0, 
       lwd.ticks = 1, 
       line = 0)
  axis(1,
       at= .lv[-length(.lv)]+strwidth(obj$ext1[1], "user")/2, 
       labels= rep(obj$ext1[1], Ntracks),
       lwd= 0, 
       lwd.ticks = 0, 
       line = -0.5, 
       cex.axis= 0.8)
  axis(1,
       at= .lv[-1]-strwidth(obj$ext2[1], "user")/2, 
       labels= rep(obj$ext2[1], Ntracks),
       lwd= 0, 
       lwd.ticks = 0, 
       line = -0.5, 
       cex.axis= 0.8)
  axis(2, 
       at= .lh[-length(.lh)]+diff(.lh)/2, 
       labels = unique(sort(obj$set_ID)),
       las= 1,
       lwd= 0, 
       lwd.ticks= 0)
  axis(3, 
       at= axAt, 
       labels = obj[, unique(name), keyby= track]$V1,
       lwd= 0, 
       lwd.ticks = 0, 
       line = 0)
  
  # Plot legends
  rasterImage(matrix(colorRampPalette(col)(100), nrow= 1),
              xleft = axAt-1/Ntracks*0.3, 
              ybottom = grconvertY(1.35, "lines", "user"), 
              xright = axAt+1/Ntracks*0.3,
              ytop = grconvertY(1.65, "lines", "user"),
              xpd= T)
  rect(xleft = axAt-1/Ntracks*0.3, 
       ybottom = grconvertY(1.35, "lines", "user"), 
       xright = axAt+1/Ntracks*0.3,
       ytop = grconvertY(1.65, "lines", "user"),
       xpd= T, 
       lwd= 0.5)
  obj[, {
    text(axAt[.GRP]-1/Ntracks*0.3,
         grconvertY(1.5, "lines", "user"), 
         pos= 2,
         0,
         xpd= T, 
         cex= 0.8)
    text(axAt[.GRP]+1/Ntracks*0.3,
         grconvertY(1.5, "lines", "user"), 
         pos= 4,
         formatC(max, 
                 format = "g", 
                 digits = 2),
         xpd= T, 
         cex= 0.8)
  }, .(track, max)]
}


#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param set_IDs Vector specifying group of regions.
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param extend How much should the bed regions be extended (starting from the center). dEfault= c(-5000, 5000)
#' @param stranded Should the average track be stranded? If yes, - features are reversed :D
#' @param nbins Number of bins spanning the extended regions. Default= 501
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param plot Should the heatmap be ploted? default= T
#' @param center_label Label for the center of the heatmaps
#' @param max_FUN Function used to compute clipping max for each track. default= function(x) quantile(x, 0.995, na.rm =T)
#' @param orderTrackIdx yIdx of the tack to be used for ordering the heatmaps. Default=1 (1st left heatmap)
#' @param oder_FUN Function used for regions ordering. Default= function(x) mean(x, na.rm= T)
#' @param orderTrackIdx yIdx of the tack to be used for ordering the heatmaps. Default=1 (1st left heatmap)
#' @param col Color vector to be used for heatmaps. col= c("blue", "yellow", "white")
#' @return An object that can be used with the vl_heatmap_bw_track_plot_only() function.
#' @export

vl_heatmap_bw_track <- function(bed,
                                set_IDs= NULL,
                                tracks,
                                extend= c(-5000, 5000),
                                stranded= F,
                                nbins= 101, 
                                names= NULL,
                                plot= T,
                                center_label= "TSS",
                                max_FUN= function(x) quantile(x, 0.995, na.rm= T),
                                orderTrackIdx= 1,
                                order_FUN= function(x) mean(x, na.rm= T),
                                col= c("blue", "yellow", "white"))
{
  # Hard copy bed
  if(!vl_isDTranges(bed))
    bins <- copy(vl_importBed(bed))
  if(is.null(set_IDs))
    set_IDs <- rep(1, nrow(bins))
  if(length(set_IDs)!=nrow(bins))
    stop("set_IDs length should match nrows(bed)!")
  if(stranded & !("strand" %in% colnames(bins)))
    stop("stranded= T but no 'strand' column in bed file")
  if(any(!file.exists(tracks)))
    stop("some bw files could not be found! full paths prodvided?")
  if(is.null(names))
    names <- gsub(".bw$", "", basename(tracks))
  if(length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length(", length(tracks), ")")
  if(!is.numeric(orderTrackIdx) | length(orderTrackIdx)!=1)
    stop("orderTrackIdx should be numeric of length 1 specifying the idx of track to use for ordering")
  
  #--------------------------#
  # Compute bins
  #--------------------------#
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
    .c <- data.table::as.data.table(rtracklayer::import.bw(x, selection= sel))
    data.table::setkeyv(.c, c("seqnames", "start", "end"))
    res <- data.table::foverlaps(.c, bins, nomatch = 0)
    res <- res[, .(score= max(abs(score), na.rm= T)), .(set_ID, region_ID, bin_ID)]
  })
  names(q) <- tracks
  
  #--------------------------#
  # Final data.table
  #--------------------------#
  final <- data.table::rbindlist(q, idcol = "track")
  final[, name:= names[.GRP], track]
  final[,ext1:= extend[1]]
  final[,ext2:= extend[2]]
  # Make Set_IDs as factor to freeze ordering
  if(!is.factor(final$set_ID))
    final[, set_ID:= factor(set_ID, levels = rev(unique(set_ID)))]
  # Make tracks as factor to freeze ordering
  final[, track:= factor(track, levels = tracks)]
  # Compute max, norm values
  final[, max:= max_FUN(score), track]
  final[, norm:= score]
  final[score>max, norm:= 100]
  final[score<=max, norm:= score/max*100]
  # Compute colors idx
  final[, col_idx:= round(norm)]
  final[is.na(col_idx), col_idx:= 0]
  # Compute order
  final[, region_order:= -order_FUN(.SD[track==tracks[orderTrackIdx], score]), region_ID]
  final[, region_order:= .GRP, keyby= .(set_ID, region_order, region_ID)]
  
  #--------------------------#
  # PLOT
  #--------------------------#
  if(plot)
    vl_heatmap_bw_track_plot_only(obj= final,
                                  center_label= center_label,
                                  col= col)
  invisible(final) 
}


#' Merge bw files
#'
#' @param x fastq file path (read1)
#' @param output_folder Folder where to put processed files
#' @param output_prefix prefix for processed files
#' @examples 
#' vl_bw_merge(x= c("db/bw/cutnrun_reps/H3K27Ac_PH18_rep1.bw",
#' "db/bw/cutnrun_reps/H3K27Ac_PH18_rep2.bw"), 
#' output_folder = "db/bw/cutnrun_merge/",
#' output_prefix = "H3K27Ac_PH18_merge.bw")
#' @export

vl_bw_merge <- function(x, 
                        output_folder,
                        output_prefix)
{
  if(!is.character(x))
    stop("x should bw a character vector of bw paths")
  dat <- data.table(file= x)
  seqL <- list()
  dat <- dat[, {
    .c <- rtracklayer::import.bw(file)
    seqL[[.GRP]] <<- seqlengths(.c)
    as.data.table(.c)
  }, file]
  if(all(sapply(seqL, function(x) identical(x, seqL[[1]]))))
    seqL <- seqL[[1]] else
      stop("seqlenghts differ for some files!?!")
  res <- dat[, {
    .c <- matrix(sort(unique(c(start, end))), 
                 ncol = 2, 
                 byrow = T)
    colnames(.c) <- c("start", "end")
    as.data.table(.c)
  }, seqnames]
  res$score <- dat[res, sum(score), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
  res <- GRanges(res)
  seqlengths(res) <- seqL[match(seqlevels(res), names(seqL))]
  rtracklayer::export.bw(object = res, 
                         paste0(output_folder, "/", output_prefix, ".bw"))
}