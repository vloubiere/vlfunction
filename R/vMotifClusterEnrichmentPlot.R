#' Compute motif enrichment
#'
#' Compute motif enrichment for the cluster in cl_columns, using all the lines as background
#'
#' @param obj A fisher test matrix similar to ?vl_motif_cl_enrich() ouput.
#' @param padj_cutoff cutoff for ballons to be ploted
#' @param log2OR_cutoff cutoff for ballons to be ploted
#' @param N_top Select top enriched motifs/cluster
#' @param auto_margin computes optimal plotting margins
#' 
#' @examples 
#' 
#' @return Fisher test data.table.
#' @export


vl_motif_cl_enrich_plot_only <- function(obj,
                                         padj_cutoff= 0.00001,
                                         log2OR_cutoff= 0,
                                         N_top= Inf,
                                         auto_margin= T)
{
  DT <- copy(obj)
  if(!all(c("motif_name", "cl", "pval", "padj", "log2OR") %in% names(DT)))
    stop("Input should contain c('motif_name', 'cl', 'pval', 'padj', 'log2OR') columns. See ?vl_motif_cl_enrich() output")
  if(log2OR_cutoff<0)
    stop("log2OR_cutoff should be > 0")
  
  #----------------------------------#
  # Generate plot table
  #----------------------------------#
  # padj cutoff
  sel <- DT[, any(padj <= padj_cutoff & log2OR > log2OR_cutoff), motif_name][(V1), motif_name]
  if(length(sel)==0)
    stop("No enrichment found with provided paj cutoff!")
  # Plotting object
  pl <- DT[motif_name %in% sel & padj < 0.05 & log2OR > 0]
  pl[, '-log10(pval)':= -log10(pval)]
  # Select top motif/cluster
  setorderv(pl, "-log10(pval)", order = -1)
  top_motifs <- pl[, motif_name[seq(.N)<=N_top], cl]$V1
  pl <- pl[motif_name %in% top_motifs]
  # Handle infinite log2OR
  if(any(is.infinite(pl$log2OR)))
  {
    warning("There are suspcious infinite log2OR values found after padj cutoff")
    print(unique(pl[is.infinite(pl$log2OR), motif_name]))
  }
  pl[log2OR==Inf, log2OR:= max(pl[is.finite(log2OR), log2OR])]
  pl[log2OR==(-Inf), log2OR:= min(pl[is.finite(log2OR), log2OR])]
  # pval coloring
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
  pl[, y:= as.numeric(min(.I)), motif_name]
  pl[, y:= seq(1, 0, length.out = length(unique(pl$motif_name)))[.GRP], keyby= y]
  # X coordinates
  pl[, x:= seq(0, 1, length.out = length(unique(pl$cl)))[.GRP], keyby= cl]
  
  #-----------------------------#
  # PLOT
  #-----------------------------#
  if(auto_margin)
    par(mai = c(max(strwidth(pl$cl, "inches"))+0.5,
                max(strwidth(pl$motif_name, "inches"))+0.5,
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
       labels = unique(pl$motif_name),
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
  
  invisible(pl)
}