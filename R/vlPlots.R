#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars.only If set to TRU, then only plots */N.S. Default= FALSE
#' @param pos pos argument (label position). see ?text()
#' @param srt srt argument (rotation). see ?text()
#' @param ... Extra arguments passed to test function
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export
vl_plot_pval_text <- function(x, 
                              y, 
                              pval, 
                              stars.only= F,
                              pos= 3,
                              offset= 0,
                              ...)
{
  if(length(y)==1 & length(x)>1)
    y <- rep(y, length(x))
  star <- cut(pval, 
              breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
              labels = c("****", "***", "**", "*", "N.S"))
  star <- as.character(star)
  ns_val <- !is.na(star) & star=="N.S"
  value <- formatC(pval, format = "e", digits = 1)
  if(stars.only)
    value <- rep("", length(value))
  if(any(ns_val))
    text(x = x[ns_val], 
         y = y[ns_val], 
         labels= bquote(paste(.(value[ns_val]), ""^N.S)), 
         cex= 0.6,
         offset= offset-0.1,
         pos= pos,
         ...)
  if(!all(ns_val))
  text(x = x[!ns_val], 
       y = y[!ns_val], 
       labels= paste0(value[!ns_val], star[!ns_val]),
       offset= offset-0.2,
       pos= pos,
       ...)
}

#' Plots a heatkey
#'
#' @param breaks Legend breaks
#' @param col Color vectors with the same length as breaks
#' @param continuous Are the breaks supposed to represent continuous variable? If not, draws discrete categories. Default= T
#' @param left left pos
#' @param top top pos
#' @param height height
#' @param width width
#' @param main Title
#' @param ticks.cex cex expansion factor for ticks labels
#' @param main.cex cex expansion factor for legend title
#' @param main Title
#' @param border Border color. default= "black"
#'
#' @return Plots heatkey
#' @export
vl_heatkey <- function(breaks,
                       col,
                       continuous= T,
                       left= par("usr")[2],
                       top= par("usr")[4],
                       height= strheight("M")*6,
                       width= strwidth("M"),
                       ticks.cex= 0.6,
                       main.cex= 1,
                       main= NA,
                       border= "black")
{
  
  mat <- if(continuous)
  {
    Cc <- circlize::colorRamp2(breaks = breaks, 
                               colors = col)
    matrix(Cc(rev(seq(min(breaks, na.rm= T), 
                             max(breaks, na.rm= T), 
                             length.out= 100))),
                  ncol= 1)
  }else
    matrix(rev(col),
           ncol= 1)
  # Key
  rasterImage(mat, 
              xleft = left,
              ybottom = top-height, 
              xright = left+width, 
              ytop = top,
              interpolate = continuous,
              xpd= NA)
  rect(xleft = left,
       ybottom = top-height, 
       xright = left+width, 
       ytop = top,
       border= border,
       xpd= NA)
  # labels
  ticks <- if(continuous)
  {
    axisTicks(range(breaks, na.rm= T), 
              log= F, 
              nint = 4)
  }else
    breaks
  # labels position
  ypos <- if(continuous)
    (top-height)+height*((ticks-min(breaks, na.rm= T))/diff(range(breaks, na.rm=T))) else
      (top-height)+(height/length(breaks))*(seq(breaks)-0.5)
  # Plot labels
  text(left+width,
       ypos,
       ticks,
       cex= ticks.cex,
       pos= 4,
       xpd= NA)
  text(left,
       top+strheight("M"),
       main,
       cex= main.cex,
       pos= 4, 
       offset= 0,
       xpd= NA)
}

#' Plots balloons legend
#'
#' @param sizes Balloons sizes (cex)
#' @param labels Labels corresponding to each size
#' @param left left pos
#' @param top top pos
#' @param height height
#' @param main Title
#' @param ticks.cex cex expansion factor for ticks labels
#' @param main.cex cex expansion factor for legend title
#'
#' @return Plots heatkey
#' @export
vl_balloonskey <- function(sizes,
                           labels,
                           left= par("usr")[2], 
                           top= par("usr")[4],
                           ticks.cex= 0.6,
                           main.cex= 1,
                           main= NA)
{
  if(length(sizes) != length(labels))
    stop("sizes and labels should have the same legnth")
  width <- strwidth("M", cex= 0.55)*max(abs(sizes))
  center.top <- top-strheight("M", cex= 0.7)/2*max(abs(sizes))
  center.height <- strheight("M", cex= 0.7)*max(abs(sizes))*(length(sizes)-1)
  y <- seq(center.top-center.height,
           center.top,
           length.out=length(sizes))
  points(rep(left+width/2, length(sizes)), 
         y,
         cex= abs(sizes)+0.1,
         pch= ifelse(sizes>=0, 21, 22),
         xpd= T)
  text(rep(left+width+strwidth("M", cex= 0.7)/2, length(labels)),
       y,
       labels,
       pos= 4,
       offset= 0.25,
       xpd=T,
       cex= ticks.cex)
  text(left,
       top+strheight("M"),
       main, 
       pos= 4, 
       cex= main.cex,
       offset= 0,
       xpd= T)
}

#' Tilted x axis
#'
#' @param x x position
#' @param y Y position. Defaults to y axis position
#' @param labels labels to be plotted
#' @param srt Rotation. Default to 45
#' @param offset offset from the x pos. default= 0.25
#' @param pos Position. defaults to 2
#' @param xpd Default to NA
#' @param cex Default to par("cex.axis")
#' @param ... Extra arguments to be passed to text
#' @export
vl_tilt_xaxis <- function(x, 
                          y= rep(par("usr")[3], length(labels))-diff(grconvertY(c(0, par("mgp")[2]), "line", "user")), 
                          labels, 
                          srt= 45, 
                          offset= 0.25, 
                          pos= 2, 
                          xpd= NA, 
                          cex= par("cex.axis"), 
                          ...)
{
  text(x,
       y,
       labels= labels,
       srt= srt,
       offset= -offset,
       pos= pos,
       xpd= xpd,
       cex= cex)
}

#' Sets convenient par env
#'
#' @param bottom_strings Character vector. If specific, margin will be adjusted to afford longest string.
#' @param left_strings Character vector. If specific, margin will be adjusted to afford longest string.
#' @param top_strings Character vector. If specific, margin will be adjusted to afford longest string.
#' @param right_strings Character vector. If specific, margin will be adjusted to afford longest string.
#' @param mgp 
#' @param las 
#' @param tcl 
#' @param ... 
#'
#' @return Improves default par env
#' @export
vl_par <- function(bottom_strings,
                   left_strings,
                   top_strings,
                   right_strings,
                   mgp= c(1.5, 0.5, 0), 
                   las= 1, 
                   tcl= -0.2, 
                   ...)
{
  adj <- c(1.02, 0.82, 0.82, 0.42)
  lab_dist <- grconvertY(1, "lines", "inch")*(mgp[2]+1)
  if(!missing(bottom_strings))
    adj[1] <- max(strwidth(bottom_strings, units = "inches"))+lab_dist
  if(!missing(left_strings))
    adj[2] <- max(strwidth(left_strings, units = "inches"))+lab_dist
  if(!missing(top_strings))
    adj[3] <- strheight("M", units = "inches")*max(nchar(bottom_strings))+lab_dist
  if(!missing(right_strings))
    adj[4] <- max(strwidth(right_strings, units = "inches"))+lab_dist
  if(!identical(adj, c(1.02, 0.82, 0.82, 0.42)))
    par(mai= adj, mgp= mgp, las= las, tcl= tcl, ...) else
      par(mgp= mgp, las= las, tcl= tcl, ...)
}

#' Sets up my favorite par parameters. Optimized for 3x3 inch square.
#'
#' @param mai See ?par()
#' @param tcl 
#' @param mgp 
#' @param cex 
#' @param cex.lab 
#' @param cex.axis 
#' @param bty 
#' @param ... Extra arguments to be passed to ?par()
#'
#' @return Set up nice plotting parameters
#' @export
vl_par <- function(mai= c(.9, .9, .9, .9),
                   las= 1,
                   tcl= -0.1,
                   mgp= c(1.5, 0.35, 0),
                   cex= 1,
                   cex.lab= 9/12,
                   cex.axis= 7/12,
                   bty= "n",
                   ...)
{
  par(mai= mai,
      las= las,
      tcl= tcl,
      mgp= mgp,
      cex= cex,
      cex.lab= cex.lab,
      cex.axis= cex.axis,
      bty= bty,
      ...)
}

#' Plots Rsq or PCC coeff
#'
#' @param x 
#' @param value Rsquare or r value
#' @param type Type of value. one between "rsq" and "pcc"
#' @param adjusted Is the Rsq adjusted?
#' @param bty Box around the legend? default= F
#' @param digits Rounding digits. Default= 2 
#' @param ... 
#'
#' @return Add Rsquare value as a legend on aplot
#' @export
#'
#' @examples
#' plot(1,1)
#' vl_plot_R2(value= 0.01)
#' vl_plot_R2("topright", value= 0.01, adjusted= T)
vl_plot_coeff <- function(x= "topleft",
                          value,
                          type= "rsq",
                          adjusted= F,
                          bty= "n",
                          digits= 2,
                          ...)
{
  if(!type %in% c("pcc", "rsq"))
    stop("type has to be one of 'rsq', 'pcc'")
  if(type=="rsq")
  {
    if(adjusted)
    {
      legend(x, 
             legend= bquote(Adj.~R^2 == .(round(value, digits))),
             bty= bty,
             ...)
    }else
    {
      legend(x, 
             legend= bquote(R^2 == .(round(value, digits))),
             bty= bty,
             ...)
    }
  }
  if(type=="pcc")
    legend(x, 
           legend= bquote(italic(r) == .(round(value, digits))),
           bty= bty,
           ...)
  
}
