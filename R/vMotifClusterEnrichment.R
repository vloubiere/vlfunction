#' Counts motifs
#'
#' Counts motif occurences in a set of regions
#'
#' @param bed Either a GRanges object or a data.table that can be coerced to it
#' @param genome Genome to be used for coordinates ("dm6, "dm3")
#' @param resize Should the regions be resize according to extend arg? Default= T
#' @param extend Vector containing two positive integers (e.g c(500,500) indicating how the regions should be extend around their center
#' @param sel Either "Dmel" (convenient set of Dmel motifs), NULL (all motifs) or a list of IDs existing in vl_Dmel_motifs_DB$metadata$motif_name
#' 
#' @examples 
#' test <- cl_motif_counts(bed = GRanges("chr3R", IRanges(c(2e6, 3e6), c(2e6, 3e6)+1e3)),
#' genome= "dm6",
#' resize= T, 
#' extend = c(500, 500),
#' sel= vl_Dmel_motifs_DB$metadata$motif_name[1:10])
#' 
#' @return Network plot.
#' @export

vl_motif_cl_enrich <- function(obj, 
                               cl_column)
{
  DT <- copy(obj)
  if(any(DT[[cl_column]])==0)
    stop("cl of 0 is used internally and cannot be provided!")
  if(cl_column!="cl")
    names(DT)[names(DT)==cl_column] <- "cl"
  if(!all(c("motif", "motif_counts") %in% names(DT)))
    stop("Provided DT should contain c('motif', 'motif_counts') columns. see ?vl_motif_counts() output!")
  
  # Assign 0 to NAs
  DT[is.na(cl), cl:= 0]
  DT <- DT[, {
    .ccl <- data.table(V1= unique(cl))
    if(nrow(.ccl)>1)
      .ccl[, c("OR", "pval"):= {
        .f <- fisher.test(cl==V1, motif_counts>0)
        .f[c("estimate", "p.value")]
      }, V1]
  }, motif]
  DT[, padj:= p.adjust(pval, method = "fdr"), pval]
  DT[, log2OR:= log2(OR)]
  DT <- DT[V1!=0]
  names(DT)[2] <- "cl"
  return(DT)
}

vl_motif_cl_enrich_plot_only <- function(obj,
                                         padj_cutoff= 0.00001,
                                         N_top= Inf,
                                         auto_margin= T,
                                         cl_names)
{
  DT <- copy(obj)
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  sel <- DT[any(padj<=padj_cutoff & log2OR>0), motif]
  pl <- sel[motif %in% sel & padj<0.05 & log2OR>0]
  pl[, '-log10(pval)':= -log10(pval)]
  # Select top motif/cluster
  setorderv(pl, "-log10(pval)", order = -1)
  top_motifs <- pl[, motif[seq(.N)<=N_top], cl]$V1
  pl <- pl[motif %in% top_motifs]
  # Handle infinite log2OR
  if(any(is.infinite(pl$log2OR)))
    warning("There are suspcious infinite log2OR values found after padj cutoff")
  pl[log2OR==Inf, log2OR:= max(pl[is.finite(log2OR), log2OR])]
  pl[log2OR==(-Inf), log2OR:= min(pl[is.finite(log2OR), log2OR])]
  # Colors
  pval_lims <- range(pl$`-log10(pval)`)
  Cc <- colorRamp2(pval_lims, 
                   colors = c("blue", "red"))
  pl[, col:= Cc(`-log10(pval)`)]
  # Points scaling factor
  pl[, size:= (2^log2OR)]
  pl[, size:= par("cex")*log2(size+1)]
  size_lims <- range(pl$size)
  # Y coordinates
  setorderv(pl, 
            c("cl", "-log10(pval)", "log2OR", "motif"), 
            order = c(1, -1, -1, 1))
  pl[, y:= as.numeric(min(.I)), motif]
  pl[, y:= seq(1, 0, length.out = length(unique(pl$motif)))[.GRP], keyby= y]
  # X coordinates
  pl[, x:= seq(0, 1, length.out = length(unique(pl$cl)))[.GRP], keyby= cl]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  if(auto_margin)
    par(mai = c(max(strwidth(pl$cl, "inches"))+0.5,
                max(strwidth(pl$motif, "inches"))+0.5,
                0.5,
                strwidth("-log10(pval)", "inches")+0.25),
        xaxs= "i",
        yaxs= "i")
  
  # Lines
  plot.new()
  segments(unique(pl$x),
           par("usr")[3],
           unique(pl$x),
           par("usr")[4])
  segments(par("usr")[1],
           unique(pl$y),
           par("usr")[2],
           unique(pl$y))
  # Points
  points(pl$x, 
         pl$y, 
         cex= pl$size, 
         col= adjustcolor(pl$col, 0.8),
         pch= 16,
         xpd= T)
  axis(1,
       lty= 0,
       at = unique(pl$x), 
       labels = unique(pl$cl), 
       las= 2)
  axis(2,
       lty= 0,
       at = unique(pl$y), 
       labels = unique(pl$motif),
       las= 2)
  # Legend pval
  xleft <- 1-grconvertX(strwidth("-log10(pval)", "inches"), "inches", "ndc")
  xright <- xleft+grconvertX(1, "lines", "ndc")
  xleft <- grconvertX(xleft, "ndc", "npc")
  xright <- grconvertX(xright, "ndc", "npc")
  ybottom <- 0.7
  ytop <- 1-grconvertY(strheight("A", units = "inches")*2, "inches", "ndc")
  rasterImage(matrix(rev(Cc(seq(min(pval_lims), max(pval_lims), length.out = 101)))),
              xleft,
              ybottom,
              xright,
              ytop,
              xpd=T)
  ticks <- axisTicks(pval_lims, log=F)
  ymin.ticks <- ybottom+(min(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  ymax.ticks <- ybottom+(max(ticks)-pval_lims[1])/diff(pval_lims)*(ytop-ybottom)
  text(xright,
       seq(ymin.ticks, ymax.ticks, length.out = length(ticks)),
       labels = ticks,
       pos=4,
       xpd= T,
       cex= 0.6,
       offset= 0.25)
  text(xleft,
       1-0.5*strheight("A"),
       labels = "-log10(pval)",
       pos=4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  # Legend balloons
  scale <- axisTicks(size_lims, log= F)
  maxBalloonInch <- strheight("A", units = "inches", cex= max(scale))*0.75
  maxBalloonDiamX <- grconvertX(maxBalloonInch, "inches", "ndc")
  maxBalloonDiamY <- grconvertY(maxBalloonInch, "inches", "ndc")
  btop <- 0.6-maxBalloonDiamY/2
  bbot <- btop-maxBalloonDiamY*(length(scale)-1)
  xb <- grconvertX(grconvertX(xleft, "npc", "ndc")+(maxBalloonDiamX/2), "ndc", "npc")
  points(rep(xb, length(scale)),
         grconvertY(seq(bbot, btop, length.out = length(scale)), "ndc", "npc"), 
         xpd= T,
         col= "black",
         cex= scale,
         pch= 16)
  xb <- grconvertX(grconvertX(xb, "npc", "ndc")+(maxBalloonDiamX/2), "ndc", "npc")
  text(rep(xb, length(scale)),
       grconvertY(seq(bbot, btop, length.out = length(scale)), "ndc", "npc"), 
       labels= scale,
       pos= 4, 
       xpd= T, 
       offset= 0, 
       cex= 0.8)
  text(xleft,
       grconvertY(0.6+grconvertY(strheight("A", "inches", cex= 0.8), "inches", "ndc"), "ndc", "npc"),
       labels = "log2OR",
       pos= 4,
       xpd= T,
       cex= 0.8,
       offset= 0)
  
  
}
