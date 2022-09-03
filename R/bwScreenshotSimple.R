#' Simplified version. less options, faster
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks bw or bed (peaks) files.
#' @param space space between regions. Default is 30 px
#' @param widths Length 2 integer vector specifying the width of bw and bed tracks, respectively. Default= c(100L, 20L)
#' @param highlight_bed Regions to be highlighted
#' @param col Track colors
#' @param highlight_col Color used for highlighted region
#' @param bg_col Color used for background. default= "white"
#' @param density Type of track to plot. either "track" (default) or "density"
#' @param density_col Color palette used for density bw. default= viridis::viridis(10)
#' @param space= 30,
#' @param widths= c(100L, 20L),
#' @param min Min value for tracks. default to 0! (Setting to NA will compute min internally)
#' @param max Max value for tracks (will be internally reset to 1 for non-bw files)
#' @param names names for bw/bed files
#' @param genome Genome used to plot transcripts. Available: "dm3", "dm6", "mm10"
#' @param add Should only the tracks be added on top of existing plot? default= F
#'
#' @examples
#' regions <- data.table(seqnames= "chr3R",
#' start= c(2166691, 12049006),
#' end= c(2949651, 12932873),
#' name= c("ANTP-C", "BX-C"))
#' 
#' bw <-c("../available_data_dm3/db/bw/GSE41440_H3K27ac_rep1_uniq.bw",
#' "../available_data_dm3/db/bw/GSE41440_H3K27me3_rep1_uniq.bw")
#' 
#' highlight_regions <- data.table(seqnames= "chr3R",
#' start= 12586652,
#' end= 12630183)
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
                          bg_col= "white",
                          density= "track",
                          density_col= viridis::viridis(10),
                          space= 30,
                          widths= c(100L, 20L),
                          names= NULL,
                          min= 0,
                          max= as.numeric(NA),
                          genome,
                          add= F)
{
  bed <- vl_importBed(bed)
  if(!identical(bed, unique(bed)))
    stop("Non-unique regions in bed!")
  if(is.null(names))
  {
    names <- gsub("(.*)[.].*", "\\1", basename(tracks))
    names <- make.unique(names)
  }
  if(!identical(names, unique(names)))
    stop("Please provide unique names")
  if(!is.integer(widths) | any(widths<10))
    stop("widths should be integers >= 10")
  is_bw <- grepl(".bw", tracks)
  if(length(max)>1 && length(max)!=sum(is_bw))
    stop("max should either be length 1 or match the number of bw tracks")
  if(length(min)>1 && length(max)!=sum(is_bw))
    stop("min should either be length 1 or match the number of bw tracks")
  
  #----------------------------------#
  # Compute object
  #----------------------------------#
  # Binning and bg color
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
  bins[, bg:= bg_col]
  if(!is.null(highlight_bed))
    bins[vl_covBed(bins, highlight_bed)>0, bg:= highlight_col]
  # Init obj
  obj <- data.table(file= tracks,
                    name= names,
                    type= ifelse(is_bw, "bw", "bed"),
                    col= col,
                    min= min)
  obj[type=="bw", max:= max]
  obj[type!="bw", max:= 1]
  # Quantif signal
  obj <- obj[, {
    bins[, value:= as.numeric(switch(type,
                                     "bw"= vl_bw_coverage(bins, file),
                                     "bed"= vl_covBed(bins, file)>0))]
  }, (obj)]
  # Compute min/max and scale signal
  obj[type=="bw" & is.na(min), min:= min(value, na.rm= T), name]
  obj[type=="bw" & is.na(max), max:= max(value, na.rm= T), name]
  obj[type=="bed", min:= as.numeric(0)]
  obj[type=="bed", max:= as.numeric(1)]
  # Compute x,y pos and color
  obj[, x:= rowid(name)]
  bw_obj <- rbindlist(lapply(seq(widths[1]), function(x) obj[type=="bw"]), idcol = "y")
  bed_obj <- rbindlist(lapply(seq(widths[2]), function(x) obj[type=="bed"]), idcol = "y")
  obj <- rbind(bw_obj, bed_obj)
  obj[, Cc:= as.character(ifelse(y>(value-min)/(max-min)*100, bg, col))]
  if(density=="density")
    obj[type=="bw", Cc:= {
      idx <- round((value-min)/(max-min)*100)+1
      colorRampPalette(density_col)(101)[ifelse(idx>101, 101, idx)]}]
  obj[is.na(start), Cc:= bg]
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
  if(!add)
  {
    plot.new()
    plot.window(xlim= c(1, ncol(im)),
                ylim= c(1, nrow(im)))
  }
  rasterImage(im, 
              xleft = 1, 
              xright = ncol(im),
              ybottom = 1, 
              ytop = nrow(im))
  # Labels
  if(!add)
  {
    labs <- obj[, .(lab.y= mean(y),
                    lab.cex= switch(type,
                                    "bw"= 1,
                                    "bed"= (widths[2]-5)/strheight("M", "user")),
                    max.y= max(y)-strheight(max, cex= 0.5),
                    max.val= formatC(max, format= "g"),
                    min.y= min(y)+strheight(min, cex= 0.5),
                    min.val= formatC(min, format= "g")), .(name, type, max, min)]
    labs[, {
      text(0,
           lab.y[1],
           name[1],
           pos= 2,
           xpd= T,
           cex= ifelse(lab.cex<1, lab.cex, 1)) #Adjust size of peak tracks
      if(type=="bw")
        text(0,
             c(max.y, min.y),
             c(max.val, min.val),
             pos= 2,
             xpd= T,
             cex= 0.5)
    }, (labs)]
  }
  
  #----------------------------------#
  # Transcripts
  #----------------------------------#
  if(!missing(genome) & !add)
    vl_screenshot_transcripts(obj= obj, genome= genome)
}

#' @export
vl_screenshot_transcripts <- function(obj, genome)
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
  