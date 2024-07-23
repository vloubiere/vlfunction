#' Simplified version. less options, faster
#'
#' @param bed Either a vector of bed file paths, a GRanges object or a data.table containing 'seqnames', 'start', 'end' columns
#' @param tracks bw or bed (peaks) files.
#' @param space space between regions. Default is 30 px
#' @param widths Length 2 integer vector specifying the width of bw and bed tracks, respectively. Default= c(100L, 20L)
#' @param highlight.bed Regions to be highlighted
#' @param col Track colors
#' @param highlight.col Color used for highlighted region
#' @param bg.col Color used for background. default= "white"
#' @param density Type of track to plot. either "track" (default) or "density"
#' @param density.col Color palette used for density bw. default= viridis::viridis(10)
#' @param space= 30,
#' @param widths= c(100L, 20L),
#' @param min Min value for tracks. default to 0! (Setting to NA will compute min internally)
#' @param max Max value for tracks (will be internally reset to 1 for non-bw files)
#' @param names names for bw/bed files
#' @param genome Genome used to plot transcripts. Available: "dm3", "dm6", "mm10", "hg19
#' @param min.symbol Set the min size for symbols that will be plotted. Decreasing= more symbols. default= 1  
#' @param add Should only the tracks be added on top of existing plot? default= F
#'
#' @examples
#' bed <- data.table(seqnames= "chr3R",
#'                   start= c(2166691, 12049006),
#'                   end= c(2949651, 12932873),
#'                   name= c("ANTP-C", "BX-C"))
#' 
#' bw <-c("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/cutnrun/H3K27me3_PH18_merge.bw",
#'        "/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/cutnrun/H3K27Ac_PH18_merge.bw")
#' 
#' highlight_regions <- data.table(seqnames= "chr3R",
#'                                 start= 12586652,
#'                                 end= 12630183)
#' 
#' peaks <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak"
#' 
#' vl_screenshot(bed = bed, 
#'               tracks = c(bw, peaks), 
#'               highlight.bed = highlight_regions, 
#'               col = c("black", "blue", "red"),
#'               genome = "dm6")
#'
#' @export
vl_screenshot <- function(bed,
                          tracks,
                          highlight.bed= NULL,
                          col= "black",
                          highlight.col= "lightgrey",
                          bg.col= "white",
                          density= "track",
                          density.col= c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF"),
                          space= 30,
                          widths= c(100L, 20L),
                          names= NULL,
                          min= 0,
                          max= as.numeric(NA),
                          genome,
                          min.symbol= 1,
                          add= F)
{
  # Import bed file ----
  bed <- vl_importBed(bed)
  if(!identical(bed, unique(bed)))
    stop("Non-unique regions in bed!")
  if(is.null(names))
  {
    names <- gsub("(.*)[.].*", "\\1", basename(tracks))
    names <- make.unique(names)
  }
  if(!is.integer(widths) | any(widths<10))
    stop("widths should be integers >= 10")
  is_bw <- grepl(".bw", tracks)
  if(length(max)>1 && length(max)!=sum(is_bw))
    stop("max should either be length 1 or match the number of bw tracks")
  if(length(min)>1 && length(max)!=sum(is_bw))
    stop("min should either be length 1 or match the number of bw tracks")
  
  # Binning function ----
  binFUN <- function(start, end, length_of_vector) {
    integers <- start:end
    n_integers <- length(integers)
    reps <- length_of_vector %/% n_integers
    extra <- length_of_vector %% n_integers
    result_vector <- rep(integers, reps)
    if (extra > 0) {result_vector <- c(result_vector, integers[1:extra])}
    return(sort(result_vector))
  }
  
  # Binning and bg color ----
  bins <- bed[, {
    coor <- binFUN(start, end, round(1000/nrow(bed))+1)
    res <- data.table(start= coor[-length(coor)], 
                      end= coor[-1])
    res[.I>1 & end>start, start:= start+1]
    # Return
    res
  }, .(seqnames, 
       regionID= rleid(seqnames, start, end))]
  bins[, bg:= bg.col]
  if(!is.null(highlight.bed))
    bins[vl_covBed(bins, highlight.bed)>0, bg:= highlight.col]
  # Init obj ----
  obj <- data.table(ID= seq(tracks),
                    file= tracks,
                    name= names,
                    type= ifelse(is_bw, "bw", "bed"),
                    col= col,
                    min= min)
  obj[type=="bw", max:= max]
  obj[type!="bw", max:= 1]
  # Quantif signal ----
  obj <- obj[, {
    bins[, value:= as.numeric(switch(type,
                                     "bw"= vl_bw_coverage(bins, file),
                                     "bed"= vl_covBed(bins, file)>0))]
  }, (obj)]
  # Compute min/max and scale signal ----
  obj[type=="bw" & is.na(min), min:= min(value, na.rm= T), ID]
  obj[type=="bw" & is.na(max), max:= max(value, na.rm= T), ID]
  obj[type=="bed", min:= as.numeric(0)]
  obj[type=="bed", max:= as.numeric(1)]
  # Add white spaces (NA) between regions ----
  wsp <- data.table(seqnames= rep(as.character(NA), space),
                    start= as.numeric(NA),
                    end= as.numeric(NA),
                    value= as.numeric(NA))
  obj <- obj[, {
    if(regionID<max(obj$regionID))
      rbind(.SD, wsp) else
        .SD
  }, .(ID, file, name, type, col, min, max, regionID)]
  # Compute x,y pos and color ----
  obj[, x:= rowid(ID)]
  obj <- obj[, rbindlist(rep(list(.SD), 
                             switch(type, "bw"= widths[1], "bed"= widths[2])), 
                         idcol= "y"), .(file, type)]
  obj[, Cc:= as.character(ifelse(y>(value-min)/(max-min)*widths[1], bg, col))]
  if(density=="density")
    obj[type=="bw", Cc:= {
      idx <- round((value-min)/(max-min)*widths[1])+1
      colorRampPalette(density.col)(widths[1]+1)[ifelse(idx>widths[1]+1, widths[1]+1, idx)]}]
  obj[is.na(start), Cc:= bg]
  # Borders ----
  obj[y==1 & Cc==bg & type=="bw" & !is.na(end), Cc:= col] # Full line at the bottom of bw tracks
  obj[(y<=5 | y>=widths[2]-5) & Cc!=bg & type=="bed", Cc:= bg] # Bg lines around bed tracks
  # Shift the different track bands in y
  yshift <- rev(cumsum(data.table::shift(rev(obj[, max(y), ID]$V1), 1, fill = 0)))
  obj[, y:= y+yshift[.GRP], ID]
  
  # Plot ----
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
  # Labels ----
  if(!add)
  {
    labs <- obj[, .(lab.y= mean(y),
                    lab.cex= switch(type,
                                    "bw"= 1,
                                    "bed"= (widths[2]-5)/strheight("M", "user")),
                    max.y= max(y)-strheight(max, cex= 0.5),
                    max.val= formatC(max, format= "g"),
                    min.y= min(y)+strheight(min, cex= 0.5),
                    min.val= formatC(min, format= "g")), .(ID, name, type, max, min)]
    labs[, {
      axis(2,
           at = lab.y[1],
           labels = name[1],
           lwd = 0,
           cex.axis= par("cex.lab")*ifelse(lab.cex<1, lab.cex, 1)) #Adjust size of peak tracks
      if(type=="bw")
        axis(2,
             at = c(max.y, min.y),
             labels = c(max.val, min.val),
             lwd = 0,
             cex.axis= par("cex.axis")*0.5)
    }, (labs)]
  }
  
  # Transcripts ----
  if(!missing(genome) & !add)
    vl_screenshot_transcripts(obj= obj,
                              genome= genome,
                              min.symbol= min.symbol)
}

#' @export
vl_screenshot_transcripts <- function(obj,
                                      genome,
                                      min.symbol= 1)
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
                                    Keytype= "ENTREZID"),
                       "hg19"= list(TxDb= TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    org= org.Hs.eg.db::org.Hs.eg.db,
                                    Keytype= "ENTREZID"))
  # Extract overlapping transcripts ----
  transcripts <- GenomicFeatures::transcripts(annotation$TxDb, c("TXNAME", "GENEID"))
  transcripts <- as.data.table(transcripts)
  transcripts <- transcripts[, .(GENEID= unlist(GENEID)), setdiff(names(transcripts), "GENEID")]
  
  # Add symbols ----
  transcripts[, symbol:= GENEID]
  symbols <- try(AnnotationDbi::mapIds(annotation$org,
                                       key= unique(as.character(transcripts$GENEID)),
                                       column= "SYMBOL",
                                       keytype= annotation$Keytype,
                                       multiVals= "first"), silent = T)
  if(class(symbols)=="character")
    transcripts[GENEID %in% names(na.omit(symbols)), symbol:= symbols[data.table::chmatch(GENEID, names(symbols))]]
  # Add exons ----
  exons <- as.data.table(GenomicFeatures::exons(annotation$TxDb, 
                                                filter= list("tx_name"= transcripts$TXNAME), 
                                                columns= "TXNAME"))
  exons <- exons[, .(TXNAME= unlist(TXNAME)), setdiff(names(exons), "TXNAME")]
  transcripts <- rbind(transcripts[, type:= "transcript"],
                       exons[, type:= "exon"], fill= T)
  # Overlap with plotting regions ----
  overlap <- obj[!is.na(end), .(start= first(start), 
                                end= last(end), 
                                x0= first(x), 
                                x1= last(x)), .(regionID, seqnames)]
  # Compute plotting positions ----
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
    .c[, y:= par("usr")[3]-strheight("M", cex= 1.8)*(y+1)]
  }, .EACHI, on= c("seqnames", "start<=end", "end>=start")]
  # Colors and lwd ----
  pl[, c("col", "lwd"):= 
       .(switch(as.character(strand), 
                "-"= "cornflowerblue", 
                "+"= "tomato", "grey"),
         switch(as.character(type), 
                "transcript"= 1, 
                "exon"= 3, 1)), .(strand, type)]
  # Plot ----
  segments(pl$seg.x0,
           pl$y,
           pl$seg.x1,
           pl$y,
           col= pl$col,
           xpd=T,
           lwd= pl$lwd,
           lend= 2)
  pl <- pl[seg.x1-seg.x0>strwidth(symbol, cex= 0.8*min.symbol) & type=="transcript"]
  if(nrow(pl)>0)
    text(rowMeans(pl[, .(seg.x0, seg.x1)]),
         pl$y,
         pl$symbol,
         xpd=T,
         pos= 3,
         offset= 0.2,
         cex= 0.8)
}
  