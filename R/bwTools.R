#' Merge several bigwig files into 1
#' 
#' Convenience function to merge directly into r
#'
#' @param bw A character vector containing the path of bw files to be merged
#' @param genome BSgenome. "dm3", "dm6", "mm10"
#' @param bins.width Bin size. default= 25L
#' @param output Output bw file path
#' @param restrict.seqnames If specified, bins are restricted to provided seqnames. Default= NULL
#' @param scoreFUN A function to be applied to the score column before return
#' @export
vl_bw_merge <- function(tracks,
                        genome,
                        bins.width= 25L,
                        output,
                        restrict.seqnames= NULL,
                        scoreFUN= NULL)
{
  bw <- vl_binBSgenome(genome,
                       bins.width = bins.width,
                       restrict.seqnames = restrict.seqnames)
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
  BS <- BSgenome::getBSgenome(genome)
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
#' @param na.value Value to use for NAs. default to 0
#' 
#' @examples 
# track <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/cutnrun/H3K27Ac_PH18_merge.bw"
#' bins <- vl_binBSgenome(genome= "dm3",
#'                        restrict.seqnames = "chr3R")
#' cov <- vl_bw_coverage(bins, track)
#' 
#' # GRanges method
#' binned_average_function <- function(gr, file_bw_path){
#' file_bw <- import.bw(file_bw_path, as="RleList")
#' seqlevels(gr) <- names(file_bw)
#' t1 <- Sys.time()
#' bins1 <- binnedAverage(gr,file_bw, "average_score")
#' return(bins1$average_score)
#' }
#' cov <- binned_average_function(GRanges(bins), track)
#' 
#' @export
vl_bw_coverage <- function(bed,
                           bw,
                           na.value= 0)
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
  res[is.na(res)] <- na.value
  return(res)
}

#' Bins bed file and compute bw signal
#' 
#' Compared to vl_bw_coverage, this functions bins bed file before quantifying the signal. Useful to make heatmaps/average tracks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file. Default= 1L (no subsets).
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Factors use for ordering and naming of plots. Defaults to the basenames of .bw files (as factors).
#' @param anchor Anchor point to use. Must be one of 'center' or 'region' (meaning that the whole region will be used as anchor point, and extended by upstream/downstream values). Default= 'center'.
#' @param upstream Upstream  extension of bed regions (centered on center). Default= 5000L.
#' @param downstream Downstream  extension of bed regions (centered on center). Default= 5000L.
#' @param nbins A single integer value (anchor= 'center') or a vector of 3 integers (anchor= 'region') specifying the number of bins to be used. When anchor is set to 'regions', the 3 values will correspond to the number of bins before the region, the number of bins within the region and the number of bins after. Default= 101L (anchor= 'center') or c(20L, 61L, 21L) (anchor= 'region').
#' @param ignore.strand Should the strand be ignored? Default= FALSE.
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly
#' 
#' @examples
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' vl_bw_coverage_bins(bed, set.IDs = sets, tracks= tracks)[]
#' vl_bw_coverage_bins(bed, set.IDs = sets, tracks= tracks, anchor= "region")[]
#' 
#' @export
vl_bw_coverage_bins <- function(bed,
                                set.IDs= 1L,
                                tracks,
                                names= gsub(".bw$", "", basename(tracks)),
                                anchor= "center",
                                upstream= 5000L, 
                                downstream= 5000L,
                                nbins= if(anchor=="center") 101L else if(anchor=="region") c(20L, 61L, 20L),
                                ignore.strand= FALSE,
                                genome)
{
  # Hard copy bed
  bed <- vl_importBed(bed)
  cols <- intersect(c("seqnames", "start", "end", "strand"), names(bed))
  bed <- bed[, ..cols]
  bed[, setID:= set.IDs]
  bed[, regionID:= .I]
  
  # Checks
  if(!identical(unique(tracks), tracks))
    stop("tracks should be unique")
  if(any(!file.exists(tracks)))
    stop("Some provided bw file(s) do not exist")
  if(!identical(unique(names), names))
    stop("names should be unique")
  if(!is.factor(names))
    names <- factor(names, unique(names))
  
  # Resize
  if(anchor=="center")
  {
    bed <- vl_resizeBed(bed = bed,
                        center= "center",
                        upstream= upstream,
                        downstream= downstream,
                        ignore.strand = ignore.strand,
                        genome= genome)
  }else if(anchor=="region")
  {
    up <- vl_resizeBed(bed = bed,
                       center= "start",
                       upstream= upstream,
                       downstream= -1L,
                       ignore.strand = ignore.strand,
                       genome= genome)
    down <- vl_resizeBed(bed = bed,
                         center= "end",
                         upstream= -1L,
                         downstream= downstream,
                         ignore.strand = ignore.strand,
                         genome= genome)
  }else
    stop("anchor should be one of 'center' or 'region'")
  
  # Bin
  bins <- if(anchor=="center")
  {
    .c <- vl_binBed(bed, nbins = nbins)
    .c[, bin.x:= seq(-upstream, downstream, length.out= nbins)[rowid(regionID)]]
  }else if(anchor=="region")
  {
    b1 <- vl_binBed(up, nbins = nbins[1])
    b2 <- vl_binBed(bed, nbins = nbins[2])
    b2[, binIDX:= binIDX+max(b1$binIDX)]
    b3 <- vl_binBed(down, nbins = nbins[3])
    b3[, binIDX:= binIDX+max(b2$binIDX)]
    .c <- rbind(b1, b2, b3)
    setorderv(.c, c("regionID", "binIDX"))
    .c[, bin.x:= binIDX]
  }
  if(!ignore.strand & "strand" %in% names(bins))
    bins[strand=="-", bin.x:= rev(bin.x), regionID]
  
  # Quantif tracks ----
  obj <- parallel::mclapply(tracks, function(x) {
    data.table(bins[, .(setID, regionID, bin.x)], 
               score= vl_bw_coverage(bins, x))
  }, mc.preschedule = T)
  obj <- rbindlist(obj, idcol= "name")
  obj[, name:= names[name]]
  
  # Return
  return(obj)
}

#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file. Used for ordering.
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Factors used for ordering and naming of plots. By default, .bw basenames will be used (as factors).
#' @param anchor Anchor point to use. Must be one of 'center' or 'region' (meaning that the whole region will be used as anchor point, and extended by upstream/downstream values). Default= 'center'.
#' @param upstream Upstream  extension of bed regions (centered on center). Default= 5000L.
#' @param downstream Downstream  extension of bed regions (centered on center). Default= 5000L.
#' @param nbins A single integer value (anchor= 'center') or a vector of 3 integers (anchor= 'region') specifying the number of bins to be used. When anchor is set to 'regions', the 3 values will correspond to the number of bins upstream of the region, the number of bins within the region and the number of downstream bins. Default= 101L (anchor= 'center') or c(20L, 61L, 21L) (anchor= 'region').
#' @param ignore.strand Should the strand be ignored? Default= FALSE.
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly.
#' @param plot Should the average track be ploted? Default= TRUE.
#' @param col Color to be used for plotting (dataset is ordered by 1/ factorized names and 2/ set.IDs before assigning colors).
#' @param col.adj Opacity of polygons. default= 0.5.
#' @param xlab X label. Default= "Genomic distance".
#' @param ylab Y labels. Default= "Enrichment".
#' @param ylim ylim for plotting. Default= range(data).
#' @param xlab.at Position of the x labels. Default= c(-upstream, 0, downstream) when anchor= 'center' or c(1, nbins[1]+1, sum(nbins)/2, sum(nbins[1:2]+1), sum(nbins)) when anchor= 'region'.
#' @param xlab.labs Label(s) to write on the x axis. Default= c(-upstream, "Center, downstream) when anchor= 'center' or c(-upstream, "Start", "Region", "End", downstream) when anchor= "region".
#' @param legend Should the legend be plotted? default to TRUE.
#' @param legend.pos Legend position. Default to "topleft".
#' @param legend.cex Legend cex. Defaults= 1.
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' vl_bw_average_track(bed, tracks= tracks, plot= T, upstream = 1000, downstream = 1000, set.IDs = sets)
#' vl_bw_average_track(bed, tracks= tracks, plot= T, upstream = 1000, downstream = 1000, set.IDs = sets, anchor= "region")
#' @export
vl_bw_average_track <- function(bed,
                                set.IDs= 1L,
                                tracks,
                                names= gsub(".bw$", "", basename(tracks)),
                                anchor= "center",
                                upstream= 5000L,
                                downstream= 5000L,
                                nbins= if(anchor=="center") 101L else if(anchor=="region") c(20L, 61L, 20L), 
                                ignore.strand= FALSE,
                                genome,
                                plot= T,
                                col= c("#E69F00","#68B1CB","#15A390","#96C954","#77AB7A","#4F6A6F","#D26429","#C57DA5","#999999"),
                                col.adj= 0.5,
                                xlab= "Genomic distance",
                                ylab= "Enrichment",
                                ylim= NULL,
                                xlab.at= if(anchor=="center") c(min(obj$bin.x), 0, max(obj$bin.x)) else if(anchor=="region") c(1, nbins[1]+1, sum(nbins)/2, sum(nbins[1:2]+1), sum(nbins)),
                                xlab.labs= if(anchor=="center") c(xlab.at[1], "Center", xlab.at[3]) else if(anchor=="region") c(-upstream, "Start", "Region", "End", downstream),
                                legend= TRUE,
                                legend.pos= "topleft",
                                legend.cex= 1)
{
  # Compute signal ----
  obj <- vl_bw_coverage_bins(bed= bed,
                             set.IDs= set.IDs,
                             tracks= tracks,
                             names= names,
                             anchor= anchor,
                             upstream= upstream, 
                             downstream= downstream,
                             nbins= nbins, 
                             ignore.strand= ignore.strand,
                             genome= genome)
  
  # Add color and make object----
  obj[, col:= colorRampPalette(col)(.NGRP)[.GRP], keyby= .(name, setID)]
  setattr(obj, 
          "class", 
          c("vl_bw_average_track", "data.table", "data.frame"))
  
  # Plot ----
  if(plot)
    plot(obj= obj,
         xlab= xlab,
         ylab= ylab,
         ylim= ylim,
         xlab.at= xlab.at,
         xlab.labs= xlab.labs,
         legend= legend,
         legend.pos= legend.pos,
         legend.cex= legend.cex,
         col.adj= col.adj)
  invisible(obj)
}

#' @describeIn vl_bw_average_track Method to plot average tracks.
#' @export
plot.vl_bw_average_track <- function(obj,
                                     xlab= "Genomic distance",
                                     ylab= "Enrichment",
                                     ylim= NULL,
                                     xlab.at= c(min(obj$bin.x), 0, max(obj$bin.x)),
                                     xlab.labs= c(xlab.at[1], "Center", xlab.at[3]),
                                     legend= TRUE,
                                     legend.pos= "topleft",
                                     legend.cex= 1,
                                     col.adj= 0.5)
{
  # Compute mean and standard error ----
  pl <- obj[, .(mean= mean(score, na.rm= T), 
                se= sd(score, na.rm= T)/sqrt(.N)), .(name, col, setID, bin.x)]
  
  # Initiate plot ----
  if(is.null(ylim))
  {
    ylim <- range(c(pl[, mean-se], pl[, mean+se]))
    ylim[2] <- ylim[2]+diff(ylim)/5
  }
    
  plot(NA, 
       xlim= range(pl$bin.x),
       ylim= ylim,
       ylab= ylab,
       xlab= xlab,
       xaxt= "n",
       frame= F)
  axis(1, 
       at = xlab.at,
       labels = xlab.labs)
  
  # Plot traces ----
  pl[, {
    polygon(c(bin.x, rev(bin.x)), 
            c(mean+se, rev(mean-se)),
            border= NA,
            col= adjustcolor(col[1], col.adj))
    lines(bin.x, 
          mean, 
          col= col[1])
  }, .(name, setID, col)]
  
  # Legend ----
  if(legend)
  {
    leg <- unique(pl[, .(setID, name, col)])
    if(length(unique(leg$setID))>1 & length(unique(leg$name))>1)
      leg[, labels:= paste0(name, " @ ", setID)] else if(length(unique(leg$setID))>1)
        leg[, labels:= setID] else if(length(unique(leg$name))>1)
          leg[, labels:= name]
    if("labels" %in% names(leg))
      legend(legend.pos,
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
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file. Default= 1L (no subgroups).
#' @param tracks Vector of bw files to plot. Use full paths to avoid issues.
#' @param upstream Upstream  extension of bed regions (centered on center). Default= 5000L.
#' @param downstream Downstream  extension of bed regions (centered on center). Default= 5000L.
#' @param ignore.strand Should the strand information be ignored? Default= FALSE.
#' @param nbins Number of bins spanning the extended regions. Default= 101L.
#' @param names Factors used for ordering and naming of plots. By default, .bw basenames will be used (as factors).
#' @param plot Should the average track be ploted? Default= TRUE.
#' @param center.label Label center heatmap.
#' @param col Vector of colors to be used for heatmap.
#' @param order.col Index of the column(s) to be used for ordering. If set to FALSE, no ordering besides set.IDs. Default= 1L.
#' @param fun.order Function used to aggregated per region and order heatmap. Default= function(x) mean(x, na.rm= T)
#' @param max Allows to manually specify max values. Otherwise, max is computed using the fun.max function for each track.
#' @param fun.max Function to be used for clipping. Default= function(x) quantile(x, 0.995, na.rm= T). To avoid clipping, use: function(x) max(x, na.rm= T)
#' @param fun.max Function to be used for clipping. Default= function(x) quantile(x, 0.995, na.rm= T). To avoid clipping, use: function(x) max(x, na.rm= T)
#' @param cex.labels cex paramter for labels
#' 
#' @examples 
#' bed <- list(suhw= vl_SUHW_top_peaks,
#'             STARR= vl_STARR_DSCP_top_peaks)
#' bed <- rbindlist(bed,
#'                  idcol = T,
#'                  fill = T)
#' tracks <- c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", 
#'             "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' test <- vl_bw_heatmap(bed,
#'                       set.IDs= bed$.id,
#'                       tracks,
#'                       plot= T,
#'                       upstream = 1000,
#'                       downstream = 1000,
#'                       fun.order = mean,
#'                       order.col= 2)
#' plot(test)
#' 
#' @export
vl_bw_heatmap <- function(bed,
                          set.IDs= 1,
                          tracks,
                          upstream= 5000L,
                          downstream= 5000L,
                          ignore.strand= T,
                          nbins= 101L,
                          space= 20L,
                          names= gsub(".bw$", "", basename(tracks)),
                          plot= T,
                          center.label= "Center",
                          col= c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF"),
                          order.col= 1,
                          max= NULL,
                          fun.order= function(x) mean(x, na.rm= T),
                          fun.max= function(x) quantile(x, 0.99, na.rm= T),
                          na.col= "lightgrey",
                          cex.labels= 0.6)
{
  # Compute coverage
  hm <- vl_bw_coverage_bins(bed= bed,
                            set.IDs= set.IDs,
                            tracks= tracks,
                            upstream= upstream, 
                            downstream= downstream,
                            ignore.strand= ignore.strand,
                            nbins= nbins, 
                            names= names)
  # Format
  if(!is.factor(hm$name))
    hm[, name:= factor(name, unique(as.character(name)))]
  if(!is.factor(hm$setID))
    hm[, setID:= factor(setID, unique(as.character(setID)))]
  
  # Reorder regionID depending on order.col
  if(order.col)
  {
    ord <- dcast(hm[levels(hm$name)[order.col], on= "name"],
                 regionID+setID~name,
                 value.var = "score",
                 fun.aggregate = fun.order)
    setorderv(ord, names(ord)[-1], order = c(1, rep(-1, ncol(ord)-2)))
    hm[, regionID:= factor(regionID, ord$regionID)]
  }else
  {
    setorderv(hm, "setID")
    hm[, regionID:= factor(regionID, unique(regionID))]
  }
  
  # Compute max values
  if(is.null(max))
    hm[, max:= fun.max(score), name] else
      hm[, max:= max[.GRP], name]
  
  # SAVE
  obj <- ls()
  obj <- mget(c("hm", obj[obj!="hm"]))
  setattr(obj, "class", c("vl_bw_heatmap", "list"))
  if(plot)
    plot.vl_bw_heatmap(obj)
  
  invisible(obj)
}

#' @describeIn vl_bw_heatmap Method to plot bw heatmaps
#' @export
plot.vl_bw_heatmap <- function(obj)
{
  list2env(obj, environment())
  
  # Clip outliers
  clip <- copy(hm)
  clip[, clip:= cut(score, c(seq(0, max, length.out= 100), Inf), 1:100), max]
  clip[, clip:= as.numeric(clip)]
  clip[is.na(clip), clip:= 0]
  
  # Dcast image
  im <- dcast(clip,
              name+setID+regionID~bin.x,
              value.var = "clip")
  
  # Add white space
  white <- matrix(NA,
                  nrow= nrow(im),
                  ncol= space)
  white <- as.data.table(white)
  im <- cbind(im, white)
  set.IDs <- split(im, im$name)[[1]]$setID
  im <- split(im[, -c(1,2,3), with= F], im$name)
  im <- do.call(cbind, im)
  im <- as.matrix(im[, 1:(ncol(im)-space)])
  im <- im/100
  
  # Plot
  vl_heatmap(im,
             row.clusters= set.IDs,
             cluster.rows= F,
             cluster.cols= F,
             col= col,
             show.rownames= F,
             show.colnames= F,
             na.col= "white", 
             row.cluster.line.col= "white",
             row.clusters.pos= "left",
             box.lwd= 0.25,
             legend.title= "Score")
  
  # Add Title and center
  for(i in seq(levels(clip$name)))
  {
    # Title
    text(nbins/2+((i-1)*(nbins+space)),
         par("usr")[4]+diff(grconvertY(c(0, par("mgp")[1]+.5), "line", "user")),
         levels(clip$name)[i],
         xpd= NA,
         cex= par("cex.lab"))
    
    # Genomic distance axis
    at <- c(1, nbins)+((i-1)*(nbins+space))
    anchor <- upstream/(upstream+downstream)
    at <- c(at[1], at[1]+(at[2]-at[1])*anchor, at[2])
    axis(1,
         at,
         labels = c(upstream, center.label, downstream),
         gap.axis= 0)
  }
}