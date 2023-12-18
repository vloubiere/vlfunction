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
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param ignore.strand Should the strand be ignored?
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly
#' @param nbins Number of bins spanning the extended regions. Default= 500
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @export
vl_bw_coverage_bins <- function(bed,
                                set.IDs,
                                tracks,
                                upstream, 
                                downstream,
                                ignore.strand,
                                genome,
                                nbins, 
                                names)
{
  # Hard copy bed
  bed <- vl_importBed(bed)
  # Checks
  if(!identical(unique(tracks), tracks))
    stop("tracks should be unique")
  if(!identical(unique(names), names))
    stop("names should be unique")
  if(any(!file.exists(tracks)))
    stop("Some provided bw file(s) do not exist")
  
  # Binning
  cols <- c("seqnames", "start", "end")
  if(!ignore.strand & "strand" %in% names(bed))
    cols <- c(cols, "strand")
  bed <- bed[, cols, with= F]
  bed[, c("set.IDs", "region_ID"):= .(set.IDs, .I)]
  bed <- vl_resizeBed(bed,
                      "center",
                      upstream,
                      downstream,
                      ignore.strand = ignore.strand,
                      genome= genome)
  bins <- bed[, {
    coor <- round(seq(start, end, length.out= nbins+1))
    coor <- data.table(start= coor[-length(coor)],
                       end= coor[-1])
    coor[-1, start:= start+1]
  }, setdiff(names(bed), c("start", "end"))]
  bins[, bin.x:= seq(-upstream, downstream, length.out= nbins)[rowid(region_ID)]]
  if(!ignore.strand & "strand" %in% names(bins))
    bins[strand=="-", bin.x:= rev(bin.x), region_ID]
  
  # Quantif tracks ----
  obj <- parallel::mclapply(tracks, function(x) {
    data.table(bins[, .(set.IDs, region_ID, bin.x)], 
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
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file. Used for ordering.
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Track names to plot, further used for ordering. By default, bw basenames will be used (as factors).
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param ignore.strand Should the strand be ignored? Default= T
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly
#' @param nbins Number of bins spanning the extended regions. Default= 101L
#' @param plot Should the average track be ploted? default= T
#' @param xlab X label. default= "genomic distance"
#' @param ylab Y labels. default= "Enrichment"
#' @param ylim ylim for plotting. default= range(data)
#' @param col Color to be used for plotting. Dataset is ordered using keyby= .(names, set.IDs) before assigning colors
#' @param legend Should the legend be plotted? default to T
#' @param legend.pos Legend position. Default to "topleft"
#' @param legend.cex Legend cex. defaults to 1
#' @param col.adj Opacity of polygons. default= 0.5
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("../available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw", "../gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' vl_bw_average_track(bed, tracks= tracks, plot= T, upstream = 1000, downstream = 1000, set.IDs = sets)
#' @export
vl_bw_average_track <- function(bed,
                                set.IDs,
                                tracks,
                                names,
                                upstream= 5000,
                                downstream= 5000,
                                ignore.strand= T,
                                genome,
                                nbins= 101L, 
                                center.label= "Center",
                                plot= T,
                                xlab= "genomic distance",
                                ylab= "Enrichment",
                                ylim,
                                col= c("#E69F00","#68B1CB","#15A390","#96C954","#77AB7A","#4F6A6F","#D26429","#C57DA5","#999999"),
                                legend= T,
                                legend.pos= "topleft",
                                legend.cex= 1,
                                col.adj= 0.5)
{
  # By default, preserve order bw tracks as specified in input
  if(missing(names))
  {
    names <- gsub(".bw$", "", basename(tracks))
    names <- factor(names, levels= unique(names))
  }
  if(missing(set.IDs))
    set.IDs <- 1
  obj <- vl_bw_coverage_bins(bed= bed,
                             tracks= tracks,
                             set.IDs= set.IDs,
                             upstream= upstream, 
                             downstream= downstream,
                             ignore.strand= ignore.strand,
                             nbins= nbins, 
                             names= names,
                             genome= genome)
  obj[, col:= colorRampPalette(col)(.NGRP)[.GRP], keyby= .(name, set.IDs)]
  setattr(obj, 
          "class", 
          c("vl_bw_average_track", "data.table", "data.frame"))
  if(plot)
    plot.vl_bw_average_track(obj,
                             xlab= xlab,
                             xlab.at= c(-upstream, 0, downstream),
                             center.label= center.label,
                             ylab= ylab,
                             ylim= ylim,
                             legend= legend,
                             legend.pos= legend.pos,
                             legend.cex= legend.cex,
                             col.adj= col.adj)
  invisible(obj)
}

#' @describeIn vl_bw_average_track Method to plot average tracks
#' @export
plot.vl_bw_average_track <- function(obj,
                                     xlab= "genomic distance",
                                     xlab.at= c(min(obj$bin.x), 0, max(obj$bin.x)),
                                     center.label= "Center",
                                     ylab= "Enrichment",
                                     ylim,
                                     legend= T,
                                     legend.pos= "topleft",
                                     legend.cex= 1,
                                     col.adj= 0.5)
{
  pl <- obj[, .(mean= mean(score, na.rm= T), 
                se= sd(score, na.rm= T)/sqrt(.N)), .(name, col, set.IDs, bin.x)]
  if(missing(ylim))
    ylim <- range(c(pl[, mean-se], pl[, mean+se]))
  plot(NA, 
       xlim= range(pl$bin.x),
       ylim= ylim,
       ylab= ylab,
       xlab= xlab,
       xaxt= "n",
       frame= F)
  axis(1, 
       at = xlab.at,
       labels = c(xlab.at[1], center.label, xlab.at[3]))
  pl[, {
    polygon(c(bin.x, rev(bin.x)), 
            c(mean+se, rev(mean-se)),
            border= NA,
            col= adjustcolor(col[1], col.adj))
    lines(bin.x, 
          mean, 
          col= col[1])
  }, .(name, set.IDs, col)]
  # Legend
  if(legend)
  {
    leg <- unique(pl[, .(set.IDs, name, col)])
    if(length(unique(leg$set.IDs))>1 & length(unique(leg$name))>1)
      leg[, labels:= paste0(name, " @ ", set.IDs)] else if(length(unique(leg$set.IDs))>1)
        leg[, labels:= set.IDs] else if(length(unique(leg$name))>1)
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
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param upstream Upstream  extension of bed regions (centered on center)
#' @param downstream Downstream  extension of bed regions (centered on center)
#' @param ignore.strand Should the strande be ignored? Default= T
#' @param nbins Number of bins spanning the extended regions. Default= 100
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param plot Should the average track be ploted? default= T
#' @param venter_label Label center heatmap
#' @param col Vector of colors to be used for heatmap
#' @param order.col Index of the column(s) to be used for ordering. If set to FALSE, no ordering besides Set.IDs. Default= 1L
#' @param fun.order Function used to aggregated per region and order heatmap. default= function(x) mean(x, na.rm= T)
#' @param max Allows to manually specify max values. Otherwise, max is computed using the fun.max function for each track.
#' @param fun.max Function to be used for clipping. default= function(x) quantile(x, 0.995, na.rm= T). Using max -> no clipping
#' @param fun.max Function to be used for clipping. default= function(x) quantile(x, 0.995, na.rm= T). Using max -> no clipping
#' @param cex.labels cex paramter for labels
#' 
#' @examples 
#' bed <- list(suhw= vl_SUHW_top_peaks,
#'              STARR= vl_STARR_DSCP_top_peaks)
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
                          upstream= 5000,
                          downstream= 5000,
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
  if(!is.factor(hm$set.IDs))
    hm[, set.IDs:= factor(set.IDs, unique(as.character(set.IDs)))]
  
  # Reorder region_ID depending on order.col
  if(order.col)
  {
    ord <- dcast(hm[levels(hm$name)[order.col], on= "name"],
                 region_ID+set.IDs~name,
                 value.var = "score",
                 fun.aggregate = fun.order)
    setorderv(ord, names(ord)[-1], order = c(1, rep(-1, ncol(ord)-2)))
    hm[, region_ID:= factor(region_ID, ord$region_ID)]
  }else
  {
    setorderv(hm, "set.IDs")
    hm[, region_ID:= factor(region_ID, unique(region_ID))]
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
              name+set.IDs+region_ID~bin.x,
              value.var = "clip")
  # Add white space
  white <- matrix(NA,
                  nrow= nrow(im),
                  ncol= space)
  white <- as.data.table(white)
  im <- cbind(im, white)
  set.IDs <- split(im, im$name)[[1]]$set.IDs
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
    # browser()
    at <- c(1, nbins)+((i-1)*(nbins+space))
    anchor <- upstream/(upstream+downstream)
    at <- c(at[1], at[1]+(at[2]-at[1])*anchor, at[2])
    axis(1,
         at,
         labels = c(upstream, center.label, downstream),
         gap.axis= 0)
  }
}