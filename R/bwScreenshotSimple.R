#' Simplified version. less options, faster
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks bw of bed (peaks)
#' @param space space between regions. Default is 30 px
#' @param widths Length 2 integer vector specifying the width of bw and bed tracks, respectively. Default= c(100L, 20L)
#' @param highlight_bed Regions to be highlighted
#' @param col Track colors
#' @param highlight_col Color used for highlighted region
#' @param max Max value for bw tracks 
#' @param names names for bw/bed files
#' @param genome Genome used to plot transcripts. Available: "dm3", "dm6", "mm10"
#'
#' @examples
#' regions <- data.table(seqnames= "chr3R",
#' start= c(2166691, 12049006),
#' end= c(2949651, 12932873),
#' name= c("ANTP-C", "BX-C"))
#' vl_isDTranges(regions)
#' 
#' bw <-c("../available_data_dm3/db/bw/GSE41440_H3K27ac_rep1_uniq.bw",
#' "../available_data_dm3/db/bw/GSE41440_H3K27me3_rep1_uniq.bw")
#' 
#' highlight_regions <- data.table(seqnames= "chr3R",
#' start= 12586652,
#' end= 12630183)
#' vl_isDTranges(highlight_regions)
#' 
#' peaks <- "db/peaks/ATAC_peaks.txt"
#' 
#' vl_screenshot(bed = regions, 
#' tracks = c(bw, peaks), 
#' highlight_bed = highlight_regions, 
#' col = c("black", "blue", "red"),
#' genome = "dm3")
#' @export
vl_screenshot <- function(bed,
                          tracks,
                          highlight_bed= NULL,
                          col= "black",
                          highlight_col= "lightgrey",
                          space= 30,
                          widths= c(100L, 20L),
                          names= NULL,
                          max= NULL,
                          genome)
{
  if(!vl_isDTranges(bed))
    bed <- vl_importBed(bed)
  if(!identical(bed, unique(bed)))
    stop("Non-unique regions in bed!")
  if(is.null(names))
    names <- gsub("(.*)[.].*", "\\1", basename(tracks))
  if(!identical(names, unique(names)))
    stop("Please provide unique names")
  if(!is.integer(widths) | any(widths<10))
    stop("widths should be integers >= 10")
  
  #----------------------------------#
  # Compute object
  #----------------------------------#
  # Binning
  bins <- bed[, {
    coor <- round(seq(start, 
                      end, 
                      length.out= round(1000/nrow(bed))+1))
    if(.GRP<.NGRP)  # Interpolate NAs between regions
      coor <- append(coor, rep(NA, space))
    res <- data.table(start= coor[-length(coor)], 
                      end= coor[-1])
    res[-1, start:= start+1]
  }, .(seqnames, 
       regionID= rleid(seqnames, start, end))]
  # Init obj
  obj <- data.table(file= tracks,
                    name= names,
                    type= ifelse(grepl(".bw", tracks), "bw", "bed"),
                    col= col)
  # Quantif signal
  obj <- obj[, {
    bins[, value:= as.numeric(switch(type,
                                     "bw"= vl_bw_coverage(bins, file),
                                     "bed"= vl_covBed(bins, file)>0))]
  }, (obj)]
  # background col
  obj[, bg:= "white"]
  if(!is.null(highlight_bed))
    obj[vl_covBed(obj, highlight_bed)>0, bg:= highlight_col]
  # Compute max
  if(is.null(max))
    obj[type=="bw", max:= max(value, na.rm= T), name] else if(length(max) == uniqueN(obj[type=="bw"], "name"))
      obj[type=="bw", max:= max[.GRP], name] else
        stop("length max should match the length of bw tracks (bed tracks handled automatically!)")
  obj[type=="bed", max:= as.numeric(1)]
  # Compute x,y pos and color
  obj[, x:= rowid(name)]
  obj <- obj[, {
    y <- seq(switch(type, 
                    "bw"= widths[1], 
                    "bed"= widths[2]))
    .(y, 
      Cc= as.character(ifelse(y>value/max*100, bg, col)))
  }, .(seqnames, start, end, regionID, x, type, name, max, bg, col)]
  # Borders
  obj[y==1 & Cc==bg & type=="bw" & !is.na(end), Cc:= col] # Full line at the bottom of bw tracks
  obj[(y<=5 | y>=widths[2]-5) & Cc!=bg & type=="bed", Cc:= bg] # Bg lines around bed tracks
  # Shift the different track bands in y
  yshift <- rev(cumsum(data.table::shift(rev(obj[, max(y), name]$V1), 1, fill = 0)))
  obj[, y:= y+yshift[.GRP], name]
  
  #----------------------------------#
  # PLOT
  #----------------------------------#
  im <- dcast(obj, -y~x, value.var = "Cc")
  im <- as.matrix(im, 1)
  plot.new()
  plot.window(xlim= c(1, ncol(im)),
              ylim= c(1, nrow(im)))
  rasterImage(im, 
              xleft = 1, 
              xright = ncol(im),
              ybottom = 1, 
              ytop = nrow(im))
  # Labels
  obj[, `:=`(lab.y= mean(y),
             lab.cex= switch(type,
                             "bw"= 1,
                             "bed"= (widths[2]-5)/strheight("M", "user")),
             max.y= max(y)-strheight(max, cex= 0.5),
             max.val= formatC(max, format= "g")), .(name, type, max)]
  obj[, {
    text(0,
         lab.y[1],
         name[1],
         pos= 2,
         xpd= T,
         cex= ifelse(lab.cex<1, lab.cex, 1)) #Adjust size of peak tracks
    if(type=="bw")
      text(0,
           max.y[1],
           max.val[1],
           pos= 2,
           xpd= T,
           cex= 0.5)
  }, .(name, type, lab.y, lab.cex, max.y, max.val)]
  
  #--------------------------#
  # Transcripts
  #--------------------------#
  if(!missing(genome))
  {
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
    # Extract overlapping transcripts
    transcripts <- as.data.table(GenomicFeatures::transcripts(annotation$TxDb, c("TXNAME", "GENEID")))
    transcripts <- transcripts[, .(GENEID= unlist(GENEID)), setdiff(names(transcripts), "GENEID")]
    if(nrow(transcripts)>0)
    {
      # Add symbols
      transcripts[, symbol:= GENEID]
      symbols <- try(AnnotationDbi::mapIds(annotation$org,
                                           key= unique(as.character(transcripts$GENEID)),
                                           column= "SYMBOL",
                                           keytype= annotation$Keytype,
                                           multiVals= "first"), silent = T)
      if(class(symbols)=="character")
        transcripts[GENEID %in% names(na.omit(symbols)), symbol:= symbols[data.table::chmatch(GENEID, names(symbols))]]
      # Add exons
      exons <- as.data.table(GenomicFeatures::exons(annotation$TxDb, 
                                                    filter= list("tx_name"= transcripts$TXNAME), 
                                                    columns= "TXNAME"))
      exons <- exons[, .(TXNAME= unlist(TXNAME)), setdiff(names(exons), "TXNAME")]
      transcripts <- rbind(transcripts[, type:= "transcript"],
                           exons[, type:= "exon"], fill= T)
      # Ovelrap with plotting regions
      overlap <- obj[!is.na(end), .(start= first(start), 
                                    end= last(end), 
                                    x0= first(x), 
                                    x1= last(x)), .(regionID, seqnames)]
      # Compute plotting positions
      pl <- transcripts[overlap, {
        .c <- data.table(x.start,
                         x.end,
                         i.start,
                         i.end,
                         x0,
                         x1,
                         type,
                         TXNAME,
                         symbol,
                         strand)
        # Clip borders
        .c[x.start<i.start, x.start:= i.start]
        .c[x.end>i.end, x.end:= i.end]
        # Compute x plotting positions
        .c[, seg.x0:= (x.start-i.start)/(i.end-i.start)*(x1-x0+1)+x0]
        .c[, seg.x1:= (x.end-i.start)/(i.end-i.start)*(x1-x0+1)+x0]
        # Adjust y position depending on overlaps
        .c <- .c[order(seg.x0-seg.x1, seg.x0)]
        .c[type=="transcript", idx:= .I]
        .c$y <- .c[.c, .N, .EACHI, on= c("idx<idx", "seg.x0<=seg.x1", "seg.x1>=seg.x0")]$N
        .c[, y:= y[type=="transcript"], TXNAME]
        .c[, y:= par("usr")[3]-strheight("M", cex= 1.8)*y]
      }, .EACHI, on= c("seqnames", "start<=end", "end>=start")]
      # Colors and lwd
      pl[, c("col", "lwd"):= 
           .(switch(as.character(strand), 
                    "-"= "cornflowerblue", 
                    "+"= "tomato", "grey"),
             switch(as.character(type), 
                    "transcript"= 1, 
                    "exon"= 3, 1)), .(strand, type)]
      # Plot
      segments(pl$seg.x0,
               pl$y,
               pl$seg.x1,
               pl$y,
               col= pl$col,
               xpd=T,
               lwd= pl$lwd,
               lend= 2)
      pl <- pl[seg.x1-seg.x0>50 & type=="transcript"]
      text(rowMeans(pl[, .(seg.x0, seg.x1)]),
           pl$y,
           pl$symbol,
           xpd=T,
           pos= 3,
           offset= 0.2,
           cex= 0.8)
    }
  }
}
