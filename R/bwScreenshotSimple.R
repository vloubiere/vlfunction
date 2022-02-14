#' Simplified version. less options, faster
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks bw of bed (peaks)
#' @param space space between regions. Default is 30 px
#' @param widths Length 2 integer vector specifying the width of bw and bed tracks, respectively. Default= c(100L, 20L)
#' @param names names for bw/bed files
#' @param gtf Optional gtf file(plot genes)
#'
#' @return
#' @export

vl_screenshot_simple <- function(bed,
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
  obj <- obj[, bins[], (obj)]
  # Quantif signal
  obj[!is.na(end), value:= {
    switch(type,
           "bw"= vl_bw_coverage(.SD, file),
           "bed"= vl_importBed(file)[.SD, ifelse(.N>0, 1, 0), 
                                     .EACHI, 
                                     on= c("seqnames", "start<=end", "end>=start")]$V1)
  }, .(file, type)]
  # background col
  obj[, bg:= "white"]
  if(!is.null(highlight_bed))
  {
    if(!vl_isDTranges(highlight_bed))
      highlight_bed <- vl_importBed(highlight_bed)
    obj[highlight_bed, bg:= highlight_col, on= c("seqnames", "start<=end", "end>=start")]
    
  }
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
  obj[y==1 & Cc==bg & type=="bw", Cc:= col] # Full line at the bottom of bw tracks
  obj[(y<=5 | y>=widths[2]-5) & Cc!=bg & type=="bed", Cc:= bg] # Bg lines around bed tracks
  # Shift the different track bands in y
  yshift <- rev(cumsum(data.table::shift(rev(obj[, max(y), name]$V1), 1, fill = 0)))
  obj[, y:= y+yshift[.GRP], name]
  
  #----------------------------------#
  # PLOT
  #----------------------------------#
  opar <- as.call(c(par, par()[c("mar", "xaxs", "yaxs")])) # Used to reinitialize plotting on exit
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
    transcripts <- transcripts[bed[transcripts, .N>0, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1]

    if(nrow(transcripts)>0)
    {
      # overlap with regions
      overlap <- obj[!is.na(end), .(start= first(start), 
                                    end= last(end), 
                                    x0= first(x), 
                                    x1= last(x)), .(regionID, seqnames)]
      setkeyv(overlap, c("seqnames", "start", "end"))
      #--------------------#
      # TRANSCRIPTS
      #--------------------#
      setkeyv(transcripts, c("seqnames", "start", "end"))
      transcripts <- foverlaps(transcripts, overlap, nomatch = NULL)
      # Compute x plotting positions
      transcripts[i.start<start, i.start:= start] # Cap window borders
      transcripts[i.end>end, i.end:= end]
      transcripts[, seg.x0:= (i.start-start)/(end-start)*(x1-x0+1)+x0] # compute gene segments
      transcripts[, seg.x1:= (i.end-start)/(end-start)*(x1-x0+1)+x0]
      # Adjust y position depending on overlaps
      transcripts <- transcripts[order(seg.x0-seg.x1, seg.x0)]
      transcripts[, idx:= .I]
      transcripts$y <- transcripts[transcripts, .N, .EACHI, on= c("idx<idx", "seg.x0<=seg.x1", "seg.x1>=seg.x0")]$N
      transcripts[, y:= par("usr")[3]-strheight("M", cex= 1.8)*y]
      # Color
      transcripts[, Cc:= switch(as.integer(strand), 
                                "+"= "tomato",
                                "+"= "cornflowerblue"), strand]
      # Plot
      transcripts[, segments(seg.x0,
                             y,
                             seg.x1,
                             y,
                             col= Cc,
                             xpd=T,
                             lwd= 3)]
      #--------------------#
      # SYMBOLS
      #--------------------#
      symbol <- transcripts[seg.x1-seg.x0>(par("usr")[2]-par("usr")[1])/20, # Restrict to long genes
                        .(x= rowMeans(.SD),
                          y= mean(y),
                          symbol= GENEID), 
                        .(TXNAME, GENEID), 
                        .SDcols= c("seg.x0", "seg.x1")]
      symbols <- try(AnnotationDbi::mapIds(annotation$org,
                                           key= as.character(symbol$GENEID),
                                           column= "SYMBOL",
                                           keytype= annotation$Keytype,
                                           multiVals= "first"), silent = T)
      if(class(symbols)=="character")
        symbol[, symbol:= symbols[data.table::chmatch(GENEID, names(symbols))]]
      # Plot
      text(symbol$x,
           symbol$y,
           symbol$symbol,
           xpd=T,
           pos= 3,
           offset= 0.2,
           cex= 0.8)
      #--------------------#
      # EXONS
      #--------------------#
      exons <- as.data.table(GenomicFeatures::exons(annotation$TxDb, 
                                                    filter= list("tx_name"= transcripts$TXNAME), 
                                                    columns= "TXNAME"))
      exons <- exons[, .(TXNAME= unlist(TXNAME)), setdiff(names(exons), "TXNAME")]
      setkeyv(exons, c("seqnames", "start", "end"))
      exons <- foverlaps(exons, overlap, nomatch = NULL)
      # Compute x/y plotting positions and color
      exons[i.start<start, i.start:= start] # Cap window borders
      exons[i.end>end, i.end:= end]
      exons[, seg.x0:= (i.start-start)/(end-start)*(x1-x0+1)+x0] # compute gene segments
      exons[, seg.x1:= (i.end-start)/(end-start)*(x1-x0+1)+x0]
      exons[transcripts, y:= i.y, on= c("TXNAME", "regionID")]
      exons[, Cc:= switch(as.integer(strand), 
                          "+"= "tomato",
                          "+"= "cornflowerblue"), strand]
      #plot
      exons[, segments(seg.x0,
                       y,
                       seg.x1,
                       y,
                       col= Cc,
                       xpd=T,
                       lwd= 3)]
    }
    on.exit(eval(opar), add=TRUE, after=FALSE)
  }
}
