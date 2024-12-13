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
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns.
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file. Default= 1L (unique subset).
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Track names. Defaults to the basenames of .bw files.
#' @param anchor Anchor point to use. Must be one of 'center' or 'region' (meaning that the whole region will be used as anchor point, and extended by upstream/downstream values). Default= 'center'.
#' @param upstream Upstream  extension of bed regions (centered on center). Default= 5000L.
#' @param downstream Downstream  extension of bed regions (centered on center). Default= 5000L.
#' @param nbins A single integer value (anchor= 'center') or a vector of 3 integers (anchor= 'region') specifying the number of bins to be used. When anchor is set to 'regions', the 3 values will correspond to the number of bins before the region, the number of bins within the region and the number of bins after. Default= 501L (anchor= 'center') or c(100L, 301L, 100L) (anchor= 'region').
#' @param ignore.strand Should the strand be ignored? Default= FALSE.
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly
#' 
#' @examples
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw",
#' "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
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
                                nbins= if(anchor=="center") 501L else if(anchor=="region") c(100L, 301L, 100L),
                                ignore.strand= FALSE,
                                genome)
{
  # Hard copy bed ----
  bed <- vl_importBed(bed)
  cols <- intersect(c("seqnames", "start", "end", "strand"), names(bed))
  bed <- bed[, ..cols]
  bed[, setID:= set.IDs]
  bed[, regionID:= .I]
  
  # Checks ----
  if(!identical(unique(tracks), tracks))
    stop("tracks should be unique")
  if(any(!file.exists(tracks)))
    stop("Some provided bw file(s) do not exist")
  if(!identical(unique(names), names))
    stop("names should be unique")
  
  # Binning ----
  bins <- if(anchor=="center")
  {
    ext <- vl_resizeBed(bed = bed,
                        center= "center",
                        upstream= upstream,
                        downstream= downstream,
                        ignore.strand = ignore.strand,
                        genome= genome)
    .c <- vl_binBed(ext,
                    nbins = nbins,
                    ignore.strand = ignore.strand)
    # Bin x position
    .c[, bin.x:= seq(-upstream, downstream, length.out= nbins)[rowid(regionID)]]
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
    # Binning
    b1 <- vl_binBed(up, nbins = nbins[1], ignore.strand = ignore.strand)
    b2 <- vl_binBed(bed, nbins = nbins[2], ignore.strand = ignore.strand)
    b3 <- vl_binBed(down, nbins = nbins[3], ignore.strand = ignore.strand)
    .c <- rbind(b1, b2, b3)
    # Bin x position
    setorderv(.c,
              c("regionID", "start"))
    .c[, bin.x:= seq(.N), regionID]
  }else
    stop("anchor should be one of 'center' or 'region'")
  
  # Correct for strand ----
  if(!ignore.strand & "strand" %in% names(bins))
    bins[strand=="-", bin.x:= rev(bin.x), regionID]
  
  # Quantif tracks ----
  obj <- parallel::mclapply(tracks, function(x) {
    data.table(bins[, .(setID, regionID, bin.x)], 
               score= vl_bw_coverage(bins, x))
  }, mc.preschedule = T)
  obj <- rbindlist(obj, idcol= "name")
  obj[, name:= names[name]]
  
  # Return ----
  return(obj)
}

#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param bed Regions to plot. Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param set.IDs Set IDs specifying the groups as subsets of the bed file.
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Tracks names. By default, .bw basenames will be used.
#' @param anchor Anchor point to use. Must be one of 'center' or 'region' (meaning that the whole region will be used as anchor point, and extended by upstream/downstream values). Default= 'center'.
#' @param upstream Upstream  extension of bed regions (centered on center). Default= 5000L.
#' @param downstream Downstream  extension of bed regions (centered on center). Default= 5000L.
#' @param nbins A single integer value (anchor= 'center', default= 501L) or a vector of 3 integers (anchor= 'region', default= c(100L, 301L, 100L)) specifying the number of bins to be used.
#' @param ignore.strand Should the strand be ignored? Default= FALSE.
#' @param genome BSgenome used to check limits of extended regions. Out of limit ranges will be resized accordingly.
#' @param plot Should the average track be ploted? Default= TRUE.
#' @param ylim ylim for plotting. Default= NULL, meaning range(data) will be used.
#' @param xlab X label. Default= "Genomic distance".
#' @param ylab Y labels. Default= "Enrichment".
#' @param col Color to be used for plotting. Note: dataset is first ordered by 1/ names and 2/ set.IDs before assigning colors.
#' @param col.adj Opacity of standard error polygons. default= 0.5.
#' @param legend Should the legend be plotted? default to TRUE.
#' @param legend.pos Legend position. Default to "topleft".
#' @param legend.cex Legend cex. Defaults= 1.
#' @param xlab.labs Label(s) to write on the x axis. Default= c(-upstream, "Center", downstream) when anchor= 'center' or c(-upstream, "Start", "Region", "End", downstream) when anchor= "region".
#' @param abline Should ablines be plotted to the center (anchor= 'center') or the start and end positions of anchored regions (anchor= 'region')? Default= TRUE.
#' @param abline.col Color of ablines. Default= "gold",
#' @param abline.lty Line type for ablines. Default= 3.
#' @param abline.lwd Line width of ablines. Default= 1.
#' @examples 
#' bed <- rbind(vl_SUHW_top_peaks, vl_STARR_DSCP_top_peaks, fill= T)
#' sets <- c(rep("suhw", 100), rep("STARR", 1000))
#' tracks <- c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE41354_SuHw_rep1_uniq.bw",
#'             "/groups/stark/vloubiere/projects/gw_STARRSeq_bernardo/db/bw/DSCP_200bp_gw.UMI_cut_merged.bw")
#' 
#' par(mfrow= c(2,2))
#' vl_bw_average_track(bed,
#'                     set.IDs = sets,
#'                     tracks= tracks,
#'                     upstream = 1000,
#'                     downstream = 1000)
#' pl <- vl_bw_average_track(bed,
#'                           col= c("blue", "red"),
#'                           set.IDs = sets,
#'                           tracks= tracks,
#'                           upstream = 1000,
#'                           downstream = 1000,
#'                           anchor= "region")
#' 
#' # Play with plotting parameters
#' plot(pl) # Plot with the parameters of the initial call
#' plot(pl,
#'      col= c("#E69F00","#68B1CB","#15A390","#96C954","#77AB7A","#4F6A6F","#D26429","#C57DA5","#999999"),
#'      ylim= c(0, 200),
#'      xlab.labs= c(-1000, "TSS", "Gene", "TTS", 1000))
#' 
#' @export
vl_bw_average_track <- function(bed,
                                set.IDs= 1L,
                                tracks,
                                names= gsub(".bw$", "", basename(tracks)),
                                anchor= "center",
                                upstream= 5000L,
                                downstream= 5000L,
                                nbins= if(anchor=="center") 501L else if(anchor=="region") c(100L, 301L, 100L), 
                                ignore.strand= FALSE,
                                genome,
                                plot= T,
                                ylim= NULL,
                                xlab= "Genomic distance",
                                ylab= "Enrichment",
                                col= c("#E69F00",
                                       "#68B1CB",
                                       "#15A390",
                                       "#96C954",
                                       "#77AB7A",
                                       "#4F6A6F",
                                       "#D26429",
                                       "#C57DA5",
                                       "#999999"),
                                col.adj= 0.5,
                                legend= TRUE,
                                legend.pos= "topleft",
                                legend.cex= 1,
                                xlab.labs= NULL,
                                abline= TRUE,
                                abline.col= "black",
                                abline.lty= 3,
                                abline.lwd= 1)
{
  # Compute signal ----
  quantif <- vl_bw_coverage_bins(bed= bed,
                                 set.IDs= set.IDs,
                                 tracks= tracks,
                                 names= names,
                                 anchor= anchor,
                                 upstream= upstream, 
                                 downstream= downstream,
                                 nbins= nbins, 
                                 ignore.strand= ignore.strand,
                                 genome= genome)
  
  # Make object ----
  if(!is.factor(quantif$name) && length(unique(quantif$name))>1)
  {
    message("'names' coerced to factors with levels:")
    quantif[, name:= factor(name, unique(name))]
    print(levels(quantif$name))
  }
  if(!is.factor(quantif$setID) && length(unique(quantif$setID))>1)
  {
    message("'setID' coerced to factors with levels:")
    quantif[, setID:= factor(setID, unique(setID))]
    print(levels(quantif$setID))
  }
  setorderv(quantif, c("name", "setID"))
  obj <- mget(ls(), envir = environment())
  setattr(obj, 
          "class", 
          "vl_bw_average_track")
  
  # Plot ----
  if(plot)
    plot(obj= obj,
         ylim= ylim,
         xlab= xlab,
         ylab= ylab,
         col= col,
         col.adj= col.adj,
         legend= legend,
         legend.pos= legend.pos,
         legend.cex= legend.cex,
         xlab.labs= xlab.labs,
         abline= abline,
         abline.col= abline.col,
         abline.lty= abline.lty,
         abline.lwd= abline.lwd)
  
  # Return object ----
  invisible(obj)
}

#' @describeIn vl_bw_average_track Method to plot average tracks.
#' @export
plot.vl_bw_average_track <- function(obj,
                                     ylim= obj$ylim,
                                     xlab= obj$xlab,
                                     ylab= obj$ylab,
                                     col= obj$col,
                                     col.adj= obj$col.adj,
                                     legend= obj$legend,
                                     legend.pos= obj$legend.pos,
                                     legend.cex= obj$legend.cex,
                                     xlab.labs= obj$xlab.labs,
                                     abline= obj$abline,
                                     abline.col= obj$abline.col,
                                     abline.lty= obj$abline.lty,
                                     abline.lwd= obj$abline.lwd)
{
  # Retrieve environment ----
  quantif <- obj$quantif
  anchor  <- obj$anchor 
  upstream <- obj$upstream
  downstream <- obj$downstream
  nbins <- obj$nbins
  
  # Check x labels and compute positions ----
  if(is.null(xlab.labs))
    xlab.labs <- if(anchor=="center") c(-upstream, "Center", upstream) else if(anchor=="region")
      c(-upstream, "Start", "Region", "End", downstream)
  if(anchor=="center" && length(xlab.labs)!=3)
    stop("When anchor is set to 'center', xlabl.labs should be of length 2")
  if(anchor=="region" && length(xlab.labs)!=5)
    stop("When anchor is set to 'region', xlabl.labs should be of length 5")
  xlab.at <- if(anchor=="center")
    c(-upstream, 0, downstream) else if(anchor=="region") # Start, middle and end
      c(0, nbins[1]/sum(nbins), .5, sum(nbins[1:2])/sum(nbins), 1)*(sum(nbins)-1)+1 # Upstream region, start, center and end of anchored region, downstream region
  
  # Compute mean and standard error ----
  pl <- quantif[, .(mean= mean(score, na.rm= T), 
                    se= sd(score, na.rm= T)/sqrt(.N)), .(name, setID, bin.x)]
  
  # Initiate plot ----
  if(is.null(ylim))
  {
    ylim <- range(c(pl[, mean-se], pl[, mean+se]), na.rm= T)
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
  pl[, col:= {
    .col <- colorRampPalette(col)(.NGRP)[.GRP]
    polygon(c(bin.x, rev(bin.x)), 
            c(mean+se, rev(mean-se)),
            border= NA,
            col= adjustcolor(.col, col.adj))
    lines(bin.x, 
          mean, 
          col= .col)
    .col
  }, .(name, setID)]
  
  # Plot ablines ----
  if(abline)
  {
    x <- if(length(xlab.at)==3)
    {
      xlab.at[2]
    }else if(length(xlab.at)==5)
    {
      xlab.at[c(2,4)]
    }
    print(x)
    segments(x0 = x,
             y0 = ylim[1],
             x1 = x,
             y1 = max(pl[, mean+se])+strheight("M"),
             lty= abline.lty,
             col= abline.col,
             lwd= abline.lwd)
  }
  
  # Legend ----
  if(legend)
  {
    # Unique values
    leg <- unique(pl[, .(setID, name, col)])
    # Make names track@setID
    if(length(unique(leg$setID))>1 & length(unique(leg$name))>1) # Several subsets and several tracks
    {
      leg[, labels:= paste0(name, " @ ", setID)]
      
    }else if(length(unique(leg$setID))>1) # Several subsets and one track
    {
      leg[, labels:= setID] 
    }else if(length(unique(leg$name))>1) # Several tracks and unique subset
      leg[, labels:= name]
    # Plot if at least two traces
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
                          col= c("#440154FF",
                                 "#482878FF",
                                 "#3E4A89FF",
                                 "#31688EFF",
                                 "#26828EFF",
                                 "#1F9E89FF",
                                 "#35B779FF",
                                 "#6DCD59FF",
                                 "#B4DE2CFF",
                                 "#FDE725FF"),
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
