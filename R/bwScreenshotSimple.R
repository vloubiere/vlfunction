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
                                 gtf= NULL)
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
  obj <- bed[, {
    coor <- round(seq(start, end, length.out= round(1000/nrow(bed))+1))
    .(start= coor[-length(coor)], 
      end= coor[-1]-c(rep(1, length(coor)-2),0))
  }, .(seqnames, 
       region_start= start, 
       region_end= end)]
  # quantif signal
  obj <- cbind(obj, 
               sapply(setNames(tracks, names), function(x) {
                 if(grepl(".bw$", x))
                   vl_bw_coverage(obj, x) else
                   {
                     x <- vl_importBed(x)
                     x[obj, ifelse(.N>0, 1, 0), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
                   }
               }))
  # Compute x pos
  obj[, x:= .I+((.GRP-1)*space), .(seqnames, region_start, region_end)]
  obj <- rbind(obj, 
               data.table(x= seq(max(obj$x))[!seq(max(obj$x)) %in% obj$x]), 
               fill= T)
  # Melt and compute tracks type
  obj <- melt(obj, measure.vars = names)
  obj[, type:= ifelse(grepl(".bw", tracks), "bw", "bed")[.GRP], variable]
  # Compute max
  if(is.null(max))
    obj[, max:= max(value, na.rm= T), variable] else
    {
      obj[type=="bw", max:= max[.GRP], variable]
      obj[type=="bed", max:= 1, variable]
    }
  # Compute y pos and col    
  obj <- obj[, {
    y <- seq(switch(type, "bw"= 100, "bed"= 15))
    .(y,
      Cc= as.character(ifelse(y>value/max*100, "white", "black")))
  }, .(seqnames, region_start, region_end, x, type, variable, max)]
  obj[, y:= y+100*(length(tracks)-.GRP), variable]
  
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
  text(0,
       obj[, mean(y), variable]$V1,
       unique(obj$variable),
       pos= 2,
       xpd= T,
       cex= 1)
  max_y <- grconvertY(obj[, max(y), variable]$V1, "user", "nfc")-grconvertY(0.5, "char", "nfc")
  text(0,
       grconvertY(max_y, "nfc", "user"),
       obj[, switch(type, 
                    "bw"= round(max, 1), 
                    bed= as.numeric(NA)), .(variable, type, max)]$V1,
       pos= 2,
       xpd= T,
       cex= 0.5)
  
  #----------------------------------#
  # Add genes
  #----------------------------------#
  if(!is.null(gtf))
  {
    if(!vl_isDTranges(gtf))
    {
      gtf <- rtracklayer::import(gtf)
      GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
      gtf <- as.data.table(gtf)
      gtf <- gtf[type=="gene"]
      gtf[, seqnames:= as.character(seqnames)]
      gtf[, start:= as.numeric(start)]
      gtf[, end:= as.numeric(end)]
      gtf[, strand:= as.character(strand)]
    }
    # Extract overlapping genes
    genes <- obj[, data.table(x_start= min(x), 
                              x_end= max(x),
                              gtf[seqnames==.BY[[1]] 
                                  & start<=region_end 
                                  & end>=region_start, .(start, 
                                                         end, 
                                                         strand, 
                                                         gene_symbol)]), 
                 .(seqnames, region_start, region_end)]
    # No arrow head for genes that finish outside 
    genes[, head:= T]
    genes[start<region_start, c("start", "head"):= .(region_start, F)]
    genes[end>region_end, c("end", "head"):= .(region_end, F)]
    # Compute arrows start and end
    genes[, x0:= (start-region_start)/(region_end-region_start)*(x_end-x_start)+x_start]
    genes[, x1:= (end-region_start)/(region_end-region_start)*(x_end-x_start)+x_start]
    genes <- genes[x1-x0>0]
    genes[, print:= x1-x0>(x_end-x_start)*0.25]
    # Compute Genes y plotting position
    ypos <- grconvertY(0, "npc", "nfc")
    ypos_g <- grconvertY(ypos-grconvertY(1, "char", "nfc"), "nfc", "user")
    ypos_gn <- grconvertY(ypos-grconvertY(2, "char", "nfc"), "nfc", "user")
    # plot genes and Symbols
    genes[, {
      suppressWarnings(
        arrows(ifelse(strand=="+", x0, x1),
               ypos_g,
               ifelse(strand=="+", x1, x0), 
               ypos_g, 
               xpd= T,
               col= ifelse(strand=="+", "tomato", "cornflowerblue"),
               lwd= 2, 
               lend= 2,
               length = ifelse(head, 0.025, 0))
      )
      if(print)
        text(x= (x0+x1)/2, y= ypos_gn, labels= gene_symbol, xpd= T)
    }, .(head, print)]
  }
}
