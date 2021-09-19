#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Can be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns"
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
