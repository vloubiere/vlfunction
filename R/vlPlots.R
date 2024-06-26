#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars Should stars be printed?
#' @param values Should p values be printed?
#' @param pos pos argument (label position). see ?text()
#' @param offset P values offest (default= 0)
#' @param cex cex expansion factor passed to text
#' @param ... Extra arguments passed to test function
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export
vl_plot_pval_text <- function(x, 
                              y, 
                              pval, 
                              stars= T,
                              values= F,
                              pos= 3,
                              offset= ifelse(values, -.2, -.35),
                              cex= .6,
                              ...)
{
  # Compute Stars ----
  star <- if(stars)
    cut(pval, 
        breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
        labels = c("****", "***", "**", "*", "N.S")) else
          rep("", length(pval))
  star <- as.character(star)
  
  # Format pval ----
  lab <- if(values)
    formatC(pval, digits = 1, format= "e") else
      rep("", length(pval))
  
  # Compute labels ----
  mapply(function(x, y, pos, offset, cex, p, l, s)
  {
    var <- if(values)
    {
      if(p<2.2e-308)
        bquote(italic(P) < "2.2e-308" * .(s)) else if(p>0.05)
          bquote(italic(P) == .(l)^"N.S") else
            bquote(italic(P) == .(l) * .(s))
    }else
    {
      if(p>0.05)
        bquote(.(l)^"N.S") else
          bquote(.(l) * .(s))
    }
    text(x,
         y,
         labels= var,
         offset= offset,
         pos= pos,
         cex= cex,
         xpd= NA,
         ...)
  }, x= x, y= y, pos= pos, offset= offset, cex= cex, p= pval, l= lab, s= star)
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
#' @param border Border color. default= "black". NA means no border
#' @param show.breaks Should breaks be shown?
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
                       border= "black",
                       show.breaks= T)
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
  if(!is.na(border))
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
  if(show.breaks)
  {
    text(left+width,
         ypos,
         ticks,
         cex= ticks.cex,
         pos= 4,
         xpd= NA)
    segments(left+width,
             ypos,
             left+width+strwidth("M")/5,
             ypos,
             xpd= T)
  }
  # Plot main
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
#' @examples
#' # example code
#' pdf("testbc1558fcb6ef.pdf", 3, 3)
#' plot(1, 1, main= "before")
#' vl_par()
#' plot(1, 1, main= "after")
#' dev.off()
#' 
#' file.show("testbc1558fcb6ef.pdf")
#' file.remove("testbc1558fcb6ef.pdf")
#' 
#'  path.R
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
                   lend= 2,
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
      lend= lend,
      ...)
}

#' Plots Rsq or PCC coeff
#'
#' @param x 
#' @param value Rsquare or r value. If numeric, will be rounded, while characters will be printed as is
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
  if(is.numeric(value))
    value <- round(value, digits)
  if(type=="rsq")
  {
    if(adjusted)
    {
      legend(x, 
             legend= bquote(Adj.~R^2 == .(value)),
             bty= bty,
             ...)
    }else
    {
      legend(x, 
             legend= bquote(R^2 == .(value)),
             bty= bty,
             ...)
    }
  }
  if(type=="pcc")
    legend(x, 
           legend= bquote(italic(r) == .(value)),
           bty= bty,
           ...)
  
}


vl_legend <- function(x= par("usr")[2],
                      y= par("usr")[4],
                      legend,
                      fill,
                      bty= "n",
                      xpd= T,
                      cex= 7/12,
                      border= NA,
                      ...)
{
  legend(x= x,
         y= y,
         legend= legend,
         fill= fill,
         bty= bty,
         xpd= xpd,
         cex= cex,
         border= border,
         ...)
}