#' bw screenshot
#'
#' This function plots a screenshot of bw files for a set of regions.
#'
#' @param bed Regions to plot. Can be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns"
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param highlight_regions Granges or data.table object specifying regions to highlight
#' @param highlight_col Color to use for highlighting regions. default is "lightgrey"
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param max Max values to clip the data. If specified, must be the same length as bw vector. By default, uses max values/track.
#' @param col Ploting colors. If specified, must be the same length as bw vector. By default, uses "black"
#' @param n_genes Max number of genes to plot.
#' @param gband gband ratio. default is 0.1
#' @export

# source("/groups/stark/vloubiere/functions/my_screenshot_bw.R")
# # 
# bws <- list.files("/groups/stark/serebreni/ChIP_bws/", ".bw$", full.names = T)[-c(2,3)] # Remove 3rd track
# bws <- bws[c(1,2,15,16,4,3,6,5,8,7,10,9,12,11,14,13)]
# names <- c("BEAF-32", "DREF", "M1BP", "TAF1", "POLII", "POLII", "TBP", "TBP",
#            "TFIIA", "TFIIA", "TFIIB", "TFIIB", "TFIIF", "TFIIF", "TFIIH", "TFIIH")
# regions <- fread("/groups/stark/vloubiere/projects/MS_leo/db/pull_down_promoters/Pulldowns_CPs.txt")
# sel <- regions[Type=="DRE"][1:20]
# bed = sel[, .(seqnames= chr, start= pos-1000, end= pos+1000)]
# tracks = bws
# names = names
# col= colorRampPalette(c("cornflowerblue", "red"))(length(bws))
# n_genes = 9
# max = NULL
# gband = 0.15
# vscreenshot(bed= bed, tracks)
# file.show("/groups/stark/serebreni/ChIP_bws/screenshots_examples.pdf")

vl_screenshot <- function(bed,
                          tracks,
                          highlight_regions= NULL,
                          highlight_col= "lightgrey",
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
  if(any(!file.exists(tracks)))
    stop("some bw files could not be found! full paths prodvided?")
  if(is.null(names))
    names <- gsub(".bw$", "", basename(tracks))
  if(length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length(", length(tracks), ")")
  if(is.null(highlight_regions))
    highlight_regions <- data.table::as.data.table(GRanges())
  if(class(highlight_regions)[1]=="GRanges")
    highlight_regions <- data.table::as.data.table(highlight_regions)
  if(!data.table::is.data.table(highlight_regions) | !all(c("seqnames", "start", "end") %in% colnames(highlight_regions)))
    stop("highlight_regions must be a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns")
  if(is.null(col))
    col <- rep("black", length(tracks))
  if(length(col) != length(tracks))
    stop("Colors vector length should be the same as bw files vector length(", length(tracks), ")")
  if(!is.null(max) & length(max) != length(tracks))
    stop("max vector length should be the same as bw files vector length(", length(tracks), ")")

  #--------------------------#
  # Binning bed file
  #--------------------------#
  Nbins <- round(1000/nrow(bed))+1
  bins <- bed[,{
    binsize <- round(((end-start)+1)/Nbins)
    if(binsize==0)
      binsize <- 1
    .s <- seq(start, end, binsize)
    .e <- .s[-1]-1
    .s <- .s[-length(.s)]
    .(start= .s, end= .e, x= seq(.s))
  }, .(seqnames, region_ID= seq(nrow(bed)))]
  inter <- cumsum(bins[, max(x), region_ID]$V1+19) # Interpolate x values (19 empty between lines)
  inter <- inter-inter[1]
  bins[, x:= x+inter[.GRP], region_ID]
  bins[, bin_ID:= .I] # Used to track unmatched bins
  
  #--------------------------#
  # Highlight regions
  #--------------------------#
  data.table::setkeyv(bins, c("seqnames", "start", "end"))
  data.table::setkeyv(highlight_regions, c("seqnames", "start", "end"))
  bins$highlight_region <- highlight_regions[bins, .N>0, .EACHI, on= c("seqnames", "end>=start", "start<=end")]$V1

  #--------------------------#
  # Quantif tracks
  #--------------------------#
  sel <- rtracklayer::BigWigSelection(GenomicRanges::GRanges(bed), "score")
  q <- parallel::mclapply(tracks, function(x)
  {
    .c <- data.table::as.data.table(rtracklayer::import.bw(x, selection= sel))
    data.table::setkeyv(.c, c("seqnames", "start", "end"))
    res <- data.table::foverlaps(.c, bins, nomatch = 0)
    .strand <- res[which.max(abs(score)), sign(score)] # Compute strand
    # Compute result
    .c[bins, .(region_ID, 
               x, 
               score= max(abs(c(0, score)), na.rm= T), 
               strand= .strand, 
               highlight_region), 
       .EACHI, 
       on= c("seqnames", "start<=end", "end>=start")]
  })
  q <- data.table::rbindlist(q, idcol = "feature_ID")
  q[, name:= names[.GRP], feature_ID]
  q[, col:= col[.GRP], feature_ID]
  if(is.null(max))
    q[, max:= max(abs(c(0, score)), na.rm= T), feature_ID] else
      q[, max:= max[.GRP], feature_ID]
  q <- q[, .(y= if(strand== (-1)) seq(102L) else rev(seq(102L)), # Negative tracks will be plotted upside down
             value= as.character(ifelse(seq(0, 101)>score/max*100, ifelse(highlight_region, highlight_col, "white"), col))), (q)]
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
    .g <- data.table::foverlaps(.t, bins)
    .g[, width:= as.integer(ifelse(i.end>max(end), max(end), i.end)- # Correct width
                            ifelse(i.start<min(start), min(start), i.start)+1), region_ID]
    data.table::setorderv(.g, c("width", "TXNAME"), c(-1, 1))
    .g[TXNAME==TXNAME[1], ord:= 1]
    .g[, GRP:= .GRP, TXNAME]
    for(i in unique(.g$GRP)[-1])
    {
      .ov <- .g[GRP<i & x %in% .g[GRP==i, x]]
      .g[GRP==i, ord:= max(c(1, .ov$ord+1), na.rm= T)]
    }
    # unique(.g[, .(TXNAME, ord, SYMBOL)]) # Useful for debugging gene ordering
    .g <- .g[ord<=n_genes] # Select given number of genes to plot
    # Add exons
    .e <- data.table::as.data.table(GenomicFeatures::exonsByOverlaps(TxDb,
                                                                     GenomicRanges::GRanges(bed),
                                                                     columns= "TXNAME"))
    .e <- .e[, .(TXNAME= as.character(unlist(TXNAME))), seqnames:strand]
    .g$exon <- .e[.g, .N>0, .EACHI, on= c("seqnames", "start<=end", "end>=start", "TXNAME")]$V1
    .g <- .g[, .(y= (1:10)+(ord-1)*10,
                res= if(exon) c(rep(1, 7),0,0,0) else c(rep(1, 7),1,0,1)), (.g)]
    .g[strand=="+", value:= as.character(ifelse(res==0, "tomato", "white"))]
    .g[strand=="-", value:= as.character(ifelse(res==0, "cornflowerblue", "white"))]
  }

  #--------------------------#
  # PLOT tracks
  #--------------------------#
  empty <- data.table::CJ(x= seq(max(q$x))[!seq(max(q$x)) %in% q$x], y= unique(q$y))
  empty[, value:= "white"]
  res <- base::rbind(q, empty, fill= T)
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
  if(exists(".g"))
    if(nrow(.g)>0)
    {
      empty <- data.table::CJ(x= seq(max(q$x))[!seq(max(q$x)) %in% .g$x], y= unique(.g$y))
      empty[, value:= "white"]
      gres <- base::rbind(.g, empty, fill= T)
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
      .g[!is.na(SYMBOL) & !is.na(x),
         {
           if(diff(range(x))>=ncol(gim)/25)
             text(mean(range(x))/ncol(gim),
                  1-(min(y)+3)/nrow(gim)*gband-(1-gband),
                  SYMBOL[1],
                  cex= 0.8)
         }, .(TXNAME, region_ID)]
    }
}

