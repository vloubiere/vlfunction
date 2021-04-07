#' bw screenshot
#'
#' This function plots a screenshot of bw files for a set of regions.
#'
#' @param bed Regions to plot. Can be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns"
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param max Max values to clip the data. If specified, must be the same length as bw vector. By default, uses max values/track.
#' @param col Ploting colors. If specified, must be the same length as bw vector. By default, uses "black"
#' @param n_genes Max number of genes to plot.
#' @param gband gband ratio. default is 0.1
#' @export

# bed <- data.table::data.table(seqnames= "chrX",
#                               start= 10e6,
#                               end= 10.1e6)
# tracks <- "desktop_genomic_files/GSM1257809_TBP_E2_4_518.bw"

vscreenshot <- function(bed,
                        tracks,
                        names= NULL,
                        max= NULL,
                        col= NULL,
                        n_genes= 4,
                        gband= 0.1)
{
  if(class(bed)[1]=="GRanges")
    bed <- data.table::as.data.table(bed)
  if(!data.table::is.data.table(bed) | !all(c("seqnames", "start", "end") %in% colnames(bed)))
    stop("bed must be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns")
  if(!is.null(names) & length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length")
  if(!is.null(col) & length(col) != length(tracks))
    stop("colors vector length should be the same as bw files vector length")
  if(any(!file.exists(tracks)))
    stop("some bw files could not be found! full paths prodvided?")

  #--------------------------#
  # Binning bed file
  #--------------------------#
  Nbins <- round(1000/nrow(bed))+1
  bins <- bed[, .(start= seq(start, end, length.out = Nbins)[-Nbins],
                  end= seq(start, end, length.out = Nbins)[-1],
                  x= seq(Nbins-1)+(Nbins+19)*(.GRP-1)), .(seqnames, region_ID= seq(nrow(bed)))]
  bins[, bin_ID:= .I] # Used to track unmatched bins
  data.table::setkeyv(bins, c("seqnames", "start", "end"))
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(bed), "score")

  #--------------------------#
  # Quantif tracks
  #--------------------------#
  q <- parallel::mclapply(tracks, function(x)
  {
    .c <- data.table::as.data.table(rtracklayer::import.bw(x, selection= sel))
    data.table::setkeyv(.c, c("seqnames", "start", "end"))
    res <- data.table::foverlaps(.c, bins, nomatch = 0)
    .strand <- res[which.max(abs(score)), sign(score)] # Compute strand
    add <- bins[!(bin_ID %in% res$bin_ID),
                .(seqnames,
                  start,
                  end,
                  region_ID,
                  x,
                  strand= .strand,
                  score= 0)] # Missing bins assigned 0
    res <- res[, .(strand= .strand,
                   score= max(abs(c(0, score)), na.rm= T)),
                 .(seqnames,
                   start,
                   end,
                   region_ID,
                   x)]
    rbind(res, add)
  })
  q <- data.table::rbindlist(q, idcol = "feature_ID")
  if(is.null(names))
    q[, name:= gsub(".bw$", "", basename(tracks))[.GRP], feature_ID] else
      q[, name:= names[.GRP], feature_ID]
  if(is.null(max))
    q[, max:= max(abs(c(0, score)), na.rm= T), feature_ID] else
      q[, max:= max[.GRP], feature_ID]
  if(is.null(col))
    q[, col:= "black", feature_ID] else
      q[, col:= col[.GRP], feature_ID]
  q <- q[, .(y= if(strand==1) rev(seq(102L)) else if (strand==-1) seq(102L), # Negative tracks will be plotted upside down
             value= as.character(ifelse(seq(0, 101)>score/max*100, "white", col))), (q)]
  q[, y:= as.integer(y+(.GRP-1L)*102L), feature_ID]

  #--------------------------#
  # Transcripts
  #--------------------------#
  TxDb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene
  .t <- data.table::as.data.table(GenomicFeatures::transcriptsByOverlaps(TxDb,
                                                                         GenomicRanges::GRanges(bed),
                                                                         columns= c("TXNAME", "GENEID")))
  if(nrow(.t)>0)
  {
    .t <- .t[, .(GENEID= as.character(unlist(GENEID))), seqnames:TXNAME]
    .t[, SYMBOL:= AnnotationDbi::mapIds(org.Dm.eg.db::org.Dm.eg.db,
                                        key= GENEID,
                                        column="SYMBOL",
                                        keytype="FLYBASE",
                                        multiVals="first")]
    .t[is.na(SYMBOL), SYMBOL:= GENEID]
    data.table::setkeyv(.t, c("seqnames", "start", "end"))
    ov <- data.table::foverlaps(.t, bins)
    ov[, width:= as.integer(ifelse(i.end>max(end), max(end), i.end)- # Correct width
                            ifelse(i.start<min(start), min(start), i.start)+1), region_ID]
    data.table::setorderv(ov, c("width", "TXNAME"), c(-1, 1))
    ov[, ord:= order(-width, TXNAME), x]
    ov[, ord:= max(ord), TXNAME]
    ov <- ov[ord<=n_genes] # Select given number of genes to plot
    # Add exons
    .e <- data.table::as.data.table(GenomicFeatures::exonsByOverlaps(TxDb,
                                                                     GenomicRanges::GRanges(bed),
                                                                     columns= "TXNAME"))
    .e <- .e[, .(TXNAME= as.character(unlist(TXNAME))), seqnames:strand]
    ov$exon <- .e[ov, .N>0, .EACHI, on= c("seqnames", "start<=end", "end>=start", "TXNAME")]$V1
    ov <- ov[, .(y= (1:10)+(ord-1)*10,
                 res= if(exon) c(rep(1, 7),0,0,0) else c(rep(1, 7),1,0,1)), (ov)]
    ov[strand=="+", value:= ifelse(res==0, "tomato", "white")]
    ov[strand=="-", value:= ifelse(res==0, "cornflowerblue", "white")]
  }

  #--------------------------#
  # PLOT tracks
  #--------------------------#
  empty <- data.table::CJ(x= seq(max(q$x))[!seq(max(q$x)) %in% q$x], y= unique(q$y))
  empty[, value:= "white"]
  res <- rbind(q, empty, fill= T)
  im <- as.matrix(data.table::dcast(res, y~x, value.var = "value"), 1)

  # image
  par(xaxs= "i", yaxs= "i")
  plot.new()
  rasterImage(im, 0, 0+gband, 1, 1, interpolate = F)

  # add labels
  adj <- (1-gband)/length(tracks)/15
  res[!is.na(feature_ID),
      {
        lab_y <- 1-mean(range(y))/nrow(im)
        min_y <- ifelse(strand==(-1), 1-min(y)/nrow(im)-adj, 1-max(y)/nrow(im)+adj)
        max_y <- ifelse(strand==1,    1-min(y)/nrow(im)-adj, 1-max(y)/nrow(im)+adj)
        y <- c(min_y, lab_y, max_y)
        text(x= 0,
             y= y/(1/(1-gband))+gband, # resize absolute y values to the band dedicated to track
             pos= 2,
             labels = c(0, name[1], round(max[1]*strand, 2)),
             cex= c(0.6, 1, 0.6),
             xpd= T)
      }, .(feature_ID, strand)]

  #--------------------------#
  # PLOT genes -> plotting it after allows constant tracks/gene band ratio
  #--------------------------#
  if(nrow(ov)>0)
  {
    empty <- data.table::CJ(x= seq(max(ov$x))[!seq(max(ov$x)) %in% ov$x], y= unique(ov$y))
    empty[, value:= "white"]
    gres <- rbind(ov, empty, fill= T)
    gim <- as.matrix(data.table::dcast(gres, y~x, value.var = "value"), 1)
    gim[1,] <- NA # Important to see bottom line of the last track
    rasterImage(gim, 0, 0, 1, gband, interpolate = F)

    # Add gene label
    text(0,
         gband/2,
         labels = "Genes",
         pos= 2,
         xpd= T)

    # add symbols
    ov[!is.na(SYMBOL),
       {
         if(diff(range(x))>=ncol(gim)/25)
            text(mean(range(x))/ncol(gim),
                 1-(min(y)+3)/nrow(gim)*gband-(1-gband),
                 SYMBOL[1],
                 cex= 0.8)
       }, .(TXNAME, region_ID)]
  }
}
