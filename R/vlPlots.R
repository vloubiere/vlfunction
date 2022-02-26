#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars_only If set to TRU, then only plots */N.S. Default= FALSE
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export
vl_plot_pval_text <- function(x, 
                              y, 
                              pval, 
                              stars_only= F)
{
  if(length(y)==1 & length(x)>1)
    y <- rep(y, length(x))
  star <- cut(pval, 
              breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
              labels = c("****", "***", "**", "*", "N.S"))
  star <- as.character(star)
  ns_val <- !is.na(star) & star=="N.S"
  value <- formatC(pval, format = "e", digits = 1)
  if(stars_only)
    value <- rep("", length(value))
  if(any(ns_val))
    text(x = x[ns_val], 
         y = y[ns_val], 
         labels= bquote(paste(.(value[ns_val]), ""^N.S)), 
         cex= 0.6,
         offset= -0.1,
         pos= 3)
  if(!all(ns_val))
  text(x = x[!ns_val], 
       y = y[!ns_val], 
       labels= paste0(value[!ns_val], star[!ns_val]),
       offset= -0.2,
       pos= 3)
}

#' plot seqlogo rasterImage
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm pwm matrix
#' @param xleft left plot limit. Default= 0
#' @param ybottom bottom plot limit. Default= 0
#' @param xright right plot limit. Default= 1
#' @param ytop top plot limit. Default= 1
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' plot.new()
#' vl_seqlogo(pwm)
#' @export
vl_seqlogo <- function(pwm, 
                       xleft= 0, 
                       ybottom= 0, 
                       xright= 1, 
                       ytop= 1, 
                       add= T)
{
  if(!is.matrix(pwm))
  {
    stop("!is.matrix(pwm)")
  }
  require(fields)
  require(seqLogo)
  require(png)
  require(colorspace)
  
  tmp <- base::tempfile(fileext = ".png") 
  grDevices::png(tmp, type="cairo", width = 1000, height = 1000, units = "px")
  seqLogo::seqLogo(pwm, xaxis = F, yaxis = F)
  dev.off()
  im <- png::readPNG(tmp)
  res <- matrix(NA, nrow = nrow(im[,,1]), ncol = ncol(im[,,1]))
  res[im[,,1]==1 & im[,,2]==0 & im[,,3]==0] <- "firebrick1"
  res[im[,,1]==1 & im[,,2]>0.1 & im[,,2]<0.9 & im[,,3]==0] <- "goldenrod1"
  res[im[,,1]==0 & im[,,2]==1 & im[,,3]==0] <- "forestgreen"
  res[im[,,1]==0 & im[,,2]==0 & im[,,3]==1] <- "dodgerblue2"
  if(!add) plot.new()
  rasterImage(res[50:950,50:950], xleft= xleft, ybottom= ybottom, xright= xright, ytop= ytop, xpd= T)
}

#' figure label
#'
#' Plots a convenient label for paper figures
#'
#' @param text Text to plot.
#' @param region Can be "figure", "plot", "device"
#' @param cex Scaling factor
#' @param ... Extra args for text function
#' @export
vl_fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- grDevices::dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- graphics::grconvertX(c(0, ds[1]), from="in", to="user")
    y <- graphics::grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- graphics::par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- graphics::strwidth(text, cex=cex) * 60/100
  sh <- graphics::strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- graphics::par(xpd=NA)
  on.exit(graphics::par(old.par))
  graphics::text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

#' Heatkey plot
#'
#' @param breaks Legend breaks
#' @param col Color vectors with the same length as breaks
#' @param left left pos
#' @param top top pos
#' @param height height
#' @param width width
#' @param main Title
#' @param ticks.cex cex expansion factor for ticks labels
#' @param main.cex cex expansion factor for legend title
#'
#' @return Plots heatkey
#' @export
vl_heatkey <- function(breaks,
                       col,
                       left= par("usr")[2], 
                       top= par("usr")[4],
                       height= strheight("M")*6,
                       width= strwidth("M"),
                       ticks.cex= 0.6,
                       main.cex= 1,
                       main= NA)
{
  Cc <- circlize::colorRamp2(breaks = breaks, 
                             colors = col)
  mat <- matrix(Cc(rev(seq(min(breaks, na.rm= T), 
                           max(breaks, na.rm= T), 
                           length.out= 100))),
                ncol= 1)
  rasterImage(mat, 
              xleft = left,
              ybottom = top-height, 
              xright = left+width, 
              ytop = top,
              xpd= T)
  ticks <- axisTicks(range(breaks, na.rm= T), 
                     log= F, 
                     nint = 4)
  text(left+width,
       (top-height)+height*((ticks-min(breaks, na.rm= T))/diff(range(breaks, na.rm=T))),
       ticks,
       cex= ticks.cex,
       pos= 4,
       xpd= T)
  text(left,
       top+strheight("M"),
       main,
       cex= main.cex,
       pos= 4, 
       offset= 0,
       xpd= T)
}

#' Heatkey plot
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