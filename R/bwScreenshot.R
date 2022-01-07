#' bw screenshot
#'
#' This function plots a screenshot of bw files for a set of regions.
#'
#' @param bed Regions to plot. Can be a vector of bed file paths, a GRange object or a data.table containing 'seqnames', 'start' and 'end' columns
#' @param tracks Vector of bw files to plot. Use full paths to avoid pbs.
#' @param genome "dm3", "dm6" or "mm10"
#' @param highlight_regions Path to a bed file OR Granges OR data.table object specifying regions to highlight
#' @param highlight_col Color to use for highlighting regions. default is "lightgrey"
#' @param names Track names to plot. If specified, must be the same length as bw vector. By default, bw basenames will be used.
#' @param max Max values to clip the data. If specified, must be the same length as bw vector. By default, uses max values/track.
#' @param col Ploting colors. If specified, must be the same length as bw vector. By default, uses "black"
#' @param n_genes Max number of genes to plot.
#' @param gband gband ratio. default is 0.1
#' @examples 
#' vl_screenshot(bed = GRanges("chr3R", IRanges(1+c(0, 5e6), 1e6+c(0, 5e6))),
#' genome= "dm3",
#' tracks = c("/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw", 
#' "/groups/stark/vloubiere/projects/pe_STARRSeq/db/peaks/H3K4me3_peaks.txt",
#' "/groups/stark/vloubiere/projects/available_data_dm3/db/bw/GSE81795_H3K4me3_rep1_uniq.bw"),
#' strand = c(1,-1,-1),
#' highlight_regions = "/groups/stark/vloubiere/projects/pe_STARRSeq/db/peaks/H3K4me3_peaks.txt", 
#' n_genes = 4)
#' 
#' @export

vl_screenshot <- function(bed,
                          tracks,
                          genome,
                          highlight_regions= NULL,
                          highlight_col= "lightgrey",
                          names= NULL,
                          max= NULL,
                          strand= NULL,
                          col= NULL,
                          n_genes= 4,
                          gband= 0.1)
{
  bed <- vl_importBed(bed)
  if(any(!is.character(tracks)))
    stop("tracks can only use file paths (either .bw, .bed or .txt)")
  if(any(!file.exists(tracks)))
    stop("some track files could not be found! full paths prodvided?")
  if(is.null(names))
    names <- gsub("(.*)[.].*", "\\1", basename(tracks))
  if(length(names) != length(tracks))
    stop("names vector length should be the same as bw files vector length(", length(tracks), ")")
  if(is.null(highlight_regions))
    highlight_regions <- data.table::as.data.table(GenomicRanges::GRanges())
  highlight_regions <- vl_importBed(highlight_regions)
  if(is.null(col))
    col <- rep("black", length(tracks))
  if(length(col) != length(tracks))
    stop("Colors vector length should be the same as bw files vector length(", length(tracks), ")")
  if(!is.null(max) & length(max) != length(tracks))
    stop("max vector length should be the same as bw files vector length(", length(tracks), ")")
  if(is.null(strand))
    strand <- rep(1, length(tracks))
  if(length(strand) != length(tracks))
    stop("strand vector length should be the same as bw files vector length(", length(tracks), ")")
  
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
  for(i in unique(bins$region_ID)[-1])
  {
    add <- bins[region_ID==(i-1), max(x)]+19*(i-1)# Interpolate x values (19 empty between lines)
    bins[region_ID==i, x:= x+add]
  }
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
  # q <- parallel::mclapply(tracks, function(x)
  q <- lapply(tracks, function(x)
  {
    # Import differs for bigiwig and bed/txt files (e.g to plot peaks)
    if(grepl(".bw|.bigWig|.bigwig", x))
      .c <- data.table::as.data.table(rtracklayer::import.bw(x, selection= sel)) else
      {
        .c <- vl_importBed(x) # bed files essentially imported as 0/1 scores
        .c[, c("strand", "score"):= .("*", 1)]
      }
    data.table::setkeyv(.c, c("seqnames", "start", "end"))
    res <- data.table::foverlaps(.c, bins, nomatch = 0)
    # Compute result
    final <- .c[bins, .(region_ID, 
                        x, 
                        score= max(abs(c(0, score)), na.rm= T),
                        highlight_region), 
                .EACHI, 
                on= c("seqnames", "start<=end", "end>=start")]
  })
  q <- data.table::rbindlist(q, idcol = "feature_ID")
  
  #--------------------------#
  # Compute vars used for plotting
  #--------------------------#
  # Add strand/type/names/colors
  q[, strand:= strand[feature_ID], feature_ID]
  q[, type:= ifelse(grepl(".bw|.bigWig|.bigwig", tracks)[feature_ID], 
                    yes = "bw", # ChIP or bed will be handled slighly differently
                    no = "bed"), feature_ID] 
  q[, name:= names[.GRP], feature_ID]
  q[, col:= col[.GRP], feature_ID]
  # Compute max for each track
  if(is.null(max))
    max <- q[, max(abs(c(0, score)), na.rm= T), feature_ID]$V1
  q[, max:= max[feature_ID], feature_ID]
  # Compute plotting variables (y axis)
  q[, y_Nbins:= ifelse(type=="bw", 100, 10)] # Number of bins used for plotting (-> track width)
  q[, y_start:= sum(q[.BY, y_Nbins+1, # Starting bin has to be unique for each feature (1:100, 101:200 [...])
                      .(feature_ID, y_Nbins), 
                      on= "feature_ID<feature_ID"]$V1, na.rm= T), feature_ID]
  # Replicate each line depending on y_Nbins (1 line/pixel)
  q <- q[rep(seq(nrow(q)), q$y_Nbins)]
  q[, y:= seq(y_start, y_start+y_Nbins-1), .(y_start, y_Nbins, x)] # y pos in the raster matrix
  q[, y_Npos:= round(score/max*y_Nbins)] # Number of y bins containing signal
  # Compute pixel coloring depending on overlaps
  q[, value:= "white"] # White background
  q[(highlight_region), value:= highlight_col] # highlighted regions coloring
  # Coloring of bins containing signal
  q[strand==1 & y_Nbins-(y-y_start+1)<=y_Npos, value := col]
  q[strand==(-1) & (y-y_start+1)<=y_Npos, value := col]
  # Add white band for bed features
  q[type=="bed" & !between(y, y_start+3, y_start+y_Nbins-4), value:= "white"]
  
  #--------------------------#
  # Transcripts
  #--------------------------#
  annotation <- switch(genome, 
                       "dm3"= list(TxDb= TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                                   org= org.Dm.eg.db::org.Dm.eg.db,
                                   Keytype= "FLYBASE"),
                       "dm6"= list(TxDb= TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                                   org= org.Dm.eg.db::org.Dm.eg.db,
                                   Keytype= "FLYBASE"),
                       "mm10"= list(TxDb= TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    org= org.Mm.eg.db::org.Mm.eg.db,
                                    Keytype= "ENTREZID"))
  .t <- data.table::as.data.table(GenomicFeatures::transcriptsByOverlaps(annotation$TxDb,
                                                                         GenomicRanges::GRanges(bed),
                                                                         columns= c("TXNAME", "GENEID")))
  if(nrow(.t)>0)
  {
    .t <- .t[, .(GENEID= as.character(unlist(GENEID))), seqnames:TXNAME]
    # Try to find symbols
    .symbols <- try(AnnotationDbi::mapIds(annotation$org,
                                          key= as.character(.t$GENEID),
                                          column= "SYMBOL",
                                          keytype= annotation$Keytype,
                                          multiVals= "first"), silent = T)
    # Fix error when no symbols can be ma[apped
    if(class(.symbols)=="try-error")
      .symbols <- .t$GENEID
    # Add symbols/FBgn
    .t[, SYMBOL:= .symbols]
    .t[is.na(SYMBOL), SYMBOL:= GENEID]
    data.table::setkeyv(.t, c("seqnames", "start", "end"))
    # Object
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
    .e <- data.table::as.data.table(GenomicFeatures::exonsByOverlaps(annotation$TxDb,
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
  if(identical(par("mar"), c(5.1, 4.1, 4.1, 2.1)))
    par(xaxs= "i", 
        yaxs= "i",
        mai= c(0.25,
               max(strwidth(unique(res$name), units = "inches"))+0.25,
               0.25,
               0.25))
  plot.new()
  rasterImage(im, 0, 0+gband, 1, 1, interpolate = F)
  
  # add labels
  adj <- 1*strheight("", cex = 0.6)
  res[!is.na(feature_ID),
      {
        lab_y <- 1-mean(range(y))/nrow(im)
        min_y <- max_y <- NA # Do not plot min and max if the feature is a "bed"
        if(type=="bw")
        {
          min_y <- ifelse(strand==(-1), 1-min(y)/nrow(im)-adj, 1-max(y)/nrow(im)+adj)
          max_y <- ifelse(strand==1,    1-min(y)/nrow(im)-adj, 1-max(y)/nrow(im)+adj)
        }
        y <- c(min_y, lab_y, max_y)
        text(x= 0,
             y= y/(1/(1-gband))+gband, # resize absolute y values to the band dedicated to track
             pos= 2,
             labels = c(0, name[1], round(max[1]*strand, 2)),
             cex= c(0.6, 1, 0.6),
             xpd= T)
      }, .(feature_ID, strand, type)]
  
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
           gband/2.5,
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