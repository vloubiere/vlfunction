#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars_only If set to TRU, then only plots */N.S. Default= FALSE
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
                              stars_only= F,
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
  if(stars_only)
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

#' Plot help
#'
#' @return Simply prints an help chart showing how margins work
#' @export
vl_plot_help <- function()
{
  par(mar= c(5,4,4,2), las= 1, mgp= c(2, 0.5, 0), tcl= -0.2, oma= c(4.2,4.2,4.2,4.2))
  
  plot(0,0,xlab='xlab', ylab= "ylab", frame= F)
  
  # Plot area
  rect(grconvertX(0, "npc", "user"), grconvertY(0, "npc", "user"), grconvertX(1, "npc", "user"), grconvertY(1, "npc", "user"), border= "red", col= NA)
  text(0, par("usr")[4]-strheight("M", cex= 2), labels = "Plot", col= "red", cex= 2)
  arrows(0, 0.45, par("usr")[2], 0.45, length = 0.1, col= "red")
  arrows(0, 0.45, par("usr")[1], 0.45, length = 0.1, col= "red")
  text(0, 0.45, labels = "npc", col= "red", cex= 1, pos= 3, offset= 0.35)
  arrows(0.45, 0, 0.45, par("usr")[4], length = 0.1, col= "red")
  arrows(0.45, 0, 0.45, par("usr")[3], length = 0.1, col= "red")
  text(0.45, 0, labels = "npc", col= "red", cex= 1, pos= 4, srt= -90)
  
  # Margins area
  rect(grconvertX(0, "nfc", "user"), grconvertY(0, "nfc", "user"), grconvertX(1, "nfc", "user"), grconvertY(1, "nfc", "user"), border= "darkgreen", col= NA, xpd= NA)
  text(0, grconvertY(1, "nfc", "user")-strheight("M", cex= 2), labels = "Margins", col= "darkgreen", cex= 2, xpd= NA)
  text(0, grconvertY(1, "nfc", "user")-strheight("M", cex= 2), labels = "mar= c(5,4,4,2)", col= "darkgreen", xpd= NA, pos= 1, offset= 1)
  text(grconvertX(0, "nfc", "user"), par("usr")[3]-diff(grconvertY(c(0,1), "lines", "user")), "line 1", co= "darkgreen", xpd= NA, pos= 4)
  text(grconvertX(0, "nfc", "user"), par("usr")[3]-diff(grconvertY(c(0,2), "lines", "user")), "line 2", co= "darkgreen", xpd= NA, pos= 4)
  text(grconvertX(0, "nfc", "user"), par("usr")[3]-diff(grconvertY(c(0,3), "lines", "user")), "line 3", co= "darkgreen", xpd= NA, pos= 4)
  text(grconvertX(0, "nfc", "user"), par("usr")[3]-diff(grconvertY(c(0,4), "lines", "user")), "line 4", co= "darkgreen", xpd= NA, pos= 4)
  
  arrows(0, par("usr")[4]+strheight("M"), grconvertX(0, "nfc", "user"), par("usr")[4]+strheight("M"), length = 0.1, col= "darkgreen", xpd= NA)
  arrows(0, par("usr")[4]+strheight("M"), grconvertX(1, "nfc", "user"), par("usr")[4]+strheight("M"), length = 0.1, col= "darkgreen", xpd= NA)
  arrows(par("usr")[2]+strwidth("M"), 0, par("usr")[2]+strwidth("M"), grconvertY(0, "nfc", "user"), length = 0.1, col= "darkgreen", xpd= NA)
  arrows(par("usr")[2]+strwidth("M"), 0, par("usr")[2]+strwidth("M"), grconvertY(1, "nfc", "user"), length = 0.1, col= "darkgreen", xpd= NA)
  text(par("usr")[2], par("usr")[4]+strheight("M")*2, labels = "nfc", col= "darkgreen", cex= 1, pos= 2, xpd= NA, offset= 0)
  text(par("usr")[2]+strwidth("M")*2, par("usr")[4], labels = "nfc", col= "darkgreen", cex= 1, pos= 1, xpd= NA, srt= -90)
  
  # Outter margins
  xadj <- diff(grconvertX(c(0, .2), "lines", "user"))
  yadj <- diff(grconvertY(c(0, .2), "lines", "user"))
  rect(grconvertX(0, "ndc", "user")+xadj, grconvertY(0, "ndc", "user")+yadj, grconvertX(1, "ndc", "user")-xadj, grconvertY(1, "ndc", "user")-yadj, border= "blue", col= NA, xpd= NA, lwd= 2)
  text(0, grconvertY(1, "ndc", "user")-strheight("M", cex= 2)-yadj, labels = "Outer margin area", col= "blue", cex= 2, xpd= NA)
  text(0, grconvertY(1, "ndc", "user")-strheight("M", cex= 2)-yadj, labels = "oma= c(4,4,4,4)", col= "blue", xpd= NA, pos= 1, offset= 1)
  text(grconvertX(0, "ndc", "user"), grconvertY(0, "nfc", "user")-diff(grconvertY(c(0,1), "lines", "user")), "line 1", co= "blue", xpd= NA, pos= 4)
  text(grconvertX(0, "ndc", "user"), grconvertY(0, "nfc", "user")-diff(grconvertY(c(0,2), "lines", "user")), "line 2", co= "blue", xpd= NA, pos= 4)
  text(grconvertX(0, "ndc", "user"), grconvertY(0, "nfc", "user")-diff(grconvertY(c(0,3), "lines", "user")), "line 3", co= "blue", xpd= NA, pos= 4)
  
  arrows(0, grconvertY(1, "nfc", "user")+strheight("M"), grconvertX(0, "ndc", "user")+xadj, grconvertY(1, "nfc", "user")+strheight("M"), length = 0.1, col= "blue", xpd= NA)
  arrows(0, grconvertY(1, "nfc", "user")+strheight("M"), grconvertX(1, "ndc", "user")-xadj, grconvertY(1, "nfc", "user")+strwidth("M"), length = 0.1, col= "blue", xpd= NA)
  arrows(grconvertX(1, "nfc", "user")+strwidth("M"), 0, grconvertX(1, "nfc", "user")+strwidth("M"), grconvertY(0, "ndc", "user")+yadj, length = 0.1, col= "blue", xpd= NA)
  arrows(grconvertX(1, "nfc", "user")+strwidth("M"), 0, grconvertX(1, "nfc", "user")+strwidth("M"), grconvertY(1, "ndc", "user")-yadj, length = 0.1, col= "blue", xpd= NA)
  text(grconvertX(1, "nfc", "user"), grconvertY(1, "nfc", "user")+strheight("M")*2, labels = "ndc", col= "blue", cex= 1, pos= 2, xpd= NA, offset= 0)
  text(grconvertX(1, "nfc", "user")+strwidth("M")*2, grconvertY(1, "nfc", "user"), labels = "ndc", col= "blue", cex= 1, pos= 1, xpd= NA, srt= -90)
  
  text(grconvertX(0.5, "ndc", "user"), grconvertX(0, "nfc", "user")-diff(grconvertY(c(0,1), "lines", "user")), "par(mar= c(5,4,4,2), las= 1,\nmgp= c(2, 0.5, 0), tcl= -0.2,\noma= c(4,4,4,4))", xpd= NA, pos= 1, offset= 1.65)
}