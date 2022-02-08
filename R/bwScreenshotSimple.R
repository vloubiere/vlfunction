#' Simplified version. less options, faster
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks bw of bed (peaks)
#' @param space space between regions. Default is 30 px
#' @param names names for bw/bed files
#' @param gtf Optional gtf file(plot genes)
#'
#' @return
#' @export

vl_screenshot_simple <- function(bed,
                                 tracks,
                                 space= 30,
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
    .(start= coor[-length(coor)], 
      end= coor[-1])
  }, .(seqnames, 
       regionID= rleid(seqnames, start, end))]
  # Init obj
  obj <- data.table(file= tracks,
                    name= names,
                    type= ifelse(grepl(".bw", tracks), "bw", "bed"))
  # Quantif signal
  obj <- obj[, {
    signal <- bins[!is.na(end), 
                   value:= switch(type,
                                  "bw"= vl_bw_coverage(.SD, file),
                                  "bed"= vl_importBed(file)[.SD, ifelse(.N>0, 1, 0), 
                                                            .EACHI, 
                                                            on= c("seqnames", "start<=end", "end>=start")]$V1)]
  }, .(file, name, type)]
  # Compute max
  if(is.null(max))
    obj[type=="bw", max:= max(value, na.rm= T), name] else if(length(max) == uniqueN(obj[type=="bw"], "name"))
      obj[type=="bw", max:= max[.GRP], name] else
        stop("length max should match the length of bw tracks (bed tracks handled automatically!)")
  obj[type=="bed", max:= as.numeric(1)]
  # Compute x,y pos and color
  obj[, x:= rowid(name)]
  obj <- obj[, {
    y <- seq(switch(type, "bw"= 100, "bed"= 15))
    .(y, 
      Cc= as.character(ifelse(y>value/max*100, "white", "black")))
  }, .(seqnames, start, end, regionID, x, type, name, max)]
  obj[y==1 & Cc=="white" & type=="bw", Cc:= "black"] # Always keep black line at the bottom of bw tracks
  # Shift the different track bands in y
  yshift <- rev(cumsum(shift(rev(obj[, max(y), name]$V1), 1, fill = 0)))
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
  obj[, {
    text(0,
         mean(y),
         name[1],
         pos= 2,
         xpd= T,
         cex= 1)
    rect(0, 
         min(y)-strheight("M")/10,
         ncol(im)+1,
         min(y),
         border= NA,
         col= "white",
         lwd= 1)
    if(type=="bw")
      text(0,
           max(y)-strheight(max[1], cex= 0.5),
           formatC(max[1], format= "g"),
           pos= 2,
           xpd= T,
           cex= 0.5)
  }, .(name, type, max)]
  
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
    transcripts <- data.table::as.data.table(GenomicFeatures::transcriptsByOverlaps(annotation$TxDb,
                                                                                    GenomicRanges::GRanges(bed),
                                                                                    columns= c("TXNAME", "GENEID")))
    if(nrow(transcripts)>0)
    {
      transcripts <- transcripts[, .(GENEID= as.character(unlist(GENEID))), seqnames:TXNAME]
      # Try to find symbols
      transcripts[, SYMBOL:= {
        symbols <- try(AnnotationDbi::mapIds(annotation$org,
                                             key= as.character(GENEID),
                                             column= "SYMBOL",
                                             keytype= annotation$Keytype,
                                             multiVals= "first"), silent = T)
        # Fix error when no symbols can be ma[apped
        if(class(symbols)=="try-error")
          GENEID else
            symbols
      }]
      # Add exons
      exons <- as.data.table(GenomicFeatures::exons(annotation$TxDb, columns= "TXNAME", filter= list(tx_name= transcripts$TXNAME)))
      exons <- exons[, .(TXNAME= unlist(TXNAME)), setdiff(names(exons), "TXNAME")]
      transcripts <- rbindlist(list(transcript= transcripts,
                                    exon= exons), 
                               idcol = "type",
                               fill= T)
      # overlap with regions
      overlap <- obj[!is.na(end), .(start= first(start), 
                                    end= last(end), 
                                    x0= first(x), 
                                    x1= last(x)), .(regionID, seqnames)]
      setkeyv(overlap, c("seqnames", "start", "end"))
      setkeyv(transcripts, c("seqnames", "start", "end"))
      overlap <- foverlaps(transcripts, overlap, nomatch = NULL)
      # Compute x plotting positions
      overlap[i.start<start, i.start:= start] # Cap window borders
      overlap[i.end>end, i.end:= end]
      overlap[, seg.x0:= (i.start-start)/(end-start)*(x1-x0+1)+x0] # compute gene segments
      overlap[, seg.x1:= (i.end-start)/(end-start)*(x1-x0+1)+x0]
      overlap[, text.x:= rowMeans(.SD), .SDcols= c("seg.x0", "seg.x1")] # gene names pos
      overlap[seg.x1-seg.x0<(par("usr")[2]-par("usr")[1])/20, SYMBOL:= NA] # Do not write name small genes
      # Adjust y position depending on overlaps
      overlap <- overlap[order(seg.x0-seg.x1, seg.x0)]
      overlap[type=="transcript", idx:= .I]
      overlap$y <- overlap[overlap, .N, .EACHI, on= c("idx<idx", "seg.x0<=seg.x1", "seg.x1>=seg.x0")]$N
      overlap[, y:= y[type=="transcript"], TXNAME]
      overlap[, y:= par("usr")[3]-strheight("M", cex= 1.8)*y] 
      # Line width and col
      overlap[, lwd:= switch(type, 
                             "transcript"= 1,
                             "exon"= 3), type]
      overlap[, col:= switch(as.character(strand), 
                             "+"= "tomato",
                             "-"= "cornflowerblue"), strand]
      # Plot
      segments(overlap$seg.x0, 
               overlap$y,
               overlap$seg.x1, 
               overlap$y, 
               lwd= overlap$lwd,
               col= overlap$col,
               xpd= T,
               lend= 2)
      text(overlap$text.x,
           overlap$y,
           overlap$SYMBOL,
           xpd= T,
           pos= 3,
           offset= 0.2,
           cex= 0.8)
    }
    on.exit(eval(opar), add=TRUE, after=FALSE)
  }
}
  