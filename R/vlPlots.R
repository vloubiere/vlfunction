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

#' plot seqlogo letter
#'
#' See function vl_seqlogo
#'
#' @param letter "A", "T", "C" or "G"
#' @param xleft xleft position
#' @param ytop ytop position
vl_plotLetter <- function(letter, xleft, ytop, width, height)
{
  if(letter=="T")
  {
    x <- c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1
    y <- c(10, 10, 8.5, 8.5, 0, 0, 8.5, 8.5) * 0.1
    col <- "red"
  }else if(letter=="A")
  {
    x <- c(0, 4, 6, 2, 0, 4, 6, 10, 8, 4, 3.2, 6.8, 6.4, 3.6, 3.2) * 0.1
    y <- c(0, 10, 10, 0, 0, 10, 10, 0, 0, 10, 3, 3, 4, 4, 3) * 0.1
    col <- "forestgreen"
  }else if(letter=="C")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    col <- "dodgerblue2"
  }else if(letter=="G")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    h1 <- max(y.l1)
    r1 <- max(x.l1)
    h1 <- 0.4
    x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
    y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
    x <- c(rev(x), x.add)
    y <- c(rev(y), y.add)
    col <- "goldenrod1"
  }
  polygon(xleft+x*width, 
          ytop-(1-y)*height, 
          col= col,
          border= NA,
          xpd= T)
}

#' plot seqlogo rasterImage
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm List of pwm matrices
#' @param x x positions
#' @param y positions (centered)
#' @param pos eith 2 (left) of 4 (right)
#' @param cex.width width expansion factor applied before plotting motifs
#' @param cex.height height expansion factor applied before plotting motifs
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' plot.new()
#' vl_seqlogo(pwm)
#' @export
vl_seqlogo <- function(pwm, 
                       x,
                       y,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1,
                       add= T)
{
  if(!pos %in% c(2,4))
    stop("Unsupported pos value. Use either 2 (left) or 4 (right)")
  if(is.matrix(pwm))
    pwm <- list(pwm)

  # Make object and index
  obj <- data.table(pwm, x, y, cex.width, cex.height)
  obj[, idx:= .I]
  
  # Width only depends on cex
  obj[, width:= strwidth("M", cex= cex.width), cex.width]
  
  # For each base, compute xleft, xright, ytop, ybottom
  pl <- obj[, {
    # Import PWM and melt
    .c <- as.data.table(pwm[[1]], keep.rownames= "base")
    .c <- melt(.c, id.vars = "base")
    # xleft depends on the pos (2 or 4)
    if(pos==2)
    {
      setorderv(.c, "variable", -1)
      .c[, xleft:= x-(.GRP*width), variable]
    }else if(pos==4)
    {
      setorderv(.c, "variable")
      .c[, xleft:= x+((.GRP-1)*width), variable]
    }
    # Compute motif content per column and normalize importance
    .c[, content:= sum(value*log2(value/c(0.25, 0.25, 0.25, 0.25)), na.rm= T), variable]
    .c[, norm:= value*(content/max(content))]
    # Rank from lowest to biggest importance -> inscreasing ytop pos
    setorderv(.c, "value")
    .h <- strheight("M", cex= cex.height)
    .c[, c("height", "ytop"):= {
      heights <- cumsum(norm*.h)
      .(heights, (y-.h/2)+heights)
    }, variable]
  }, .(idx, y, width)]
  
  # Plot
  pl <- pl[height>0]
  pl[, vl_plotLetter(base[1], 
                     xleft[1], 
                     ytop[1], 
                     width[1], 
                     height[1]), .(base, xleft, ytop, width, height)]
  
  # Return object containing limits of each motif
  invisible(pl[, .(xleft= min(xleft), 
                   ybottom= min(ytop-height),
                   xright= max(xleft+width), 
                   ytop= max(ytop)), .(idx)])
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

#' Table plot
#'
#' Plots a data.table as text
#'
#' @param x data.table
#' @param colnames_col Color for col.names row
#' @param row_col Color for data rows
#' @return Plots heatkey
#'
#' @export
vl_plot_table <- function(x, colnames_col= "white", row_col= "white")
{
  if(!is.data.table(x))
    stop("x must be a data.table") else
      x[, lapply(.SD, as.character)]
  if(length(row_col)==1)
    row_col <- rep(row_col, nrow(x))
  if(length(row_col)!=nrow(x))
    stop("row_col should be lenght 1 or match the number of rows in x")
  if(length(colnames_col)!=1)
    stop("colnames_col should be lenght 1") 
  col <- c(colnames_col, row_col)
  x <- rbind(as.data.table(matrix(names(x), ncol= ncol(x))),
             x,
             use.names=FALSE)

  # Init plot to access str size
  plot.new()
  # Compute cols proportions and rows position
  grid.p <- apply(x, 2, function(x) max(strwidth(paste0(x, "MM"))))
  grid.p <- cumsum(grid.p/sum(grid.p))
  grid.x <- c(par("usr")[1], grid.p*diff(par("usr")[c(1,2)])+par("usr")[1])
  grid.y <- seq(par("usr")[1], par("usr")[2], length.out= nrow(x)+1)
  # Compute optimal cex for text
  cex.x <- max(diff(grid.x))/max(strwidth(paste0(unlist(x), "MM")))
  cex.y <- min(diff(grid.y))/max(strheight(unlist(x), cex= 2))
  cex.t <- min(c(cex.x, cex.y))
  # Adjust grid.x to text size
  ratio <- max(strwidth(paste0(unlist(x), "MM"), cex= cex.t))/max(diff(grid.x))
  half_size <- diff(par("usr")[c(1,2)])*ratio/2
  lims.x <- 0.5+half_size*c(-1,1)
  grid.x <- c(lims.x[1], grid.p*diff(lims.x)+lims.x[1])
  # Compute text position
  text.x <- rowMeans(matrix(c(grid.x[-length(grid.x)], grid.x[-1]), ncol = 2))
  text.y <- rowMeans(matrix(c(grid.y[-length(grid.y)], grid.y[-1]), ncol = 2))
  # Plot
  grid.y <- rev(grid.y)
  rect(min(grid.x),
       grid.y[-length(grid.y)],
       max(grid.x),
       grid.y[-1], 
       col= col,
       border= "white")
  segments(grid.x, par("usr")[1], grid.x, par("usr")[2], xpd= T)
  segments(min(grid.x), grid.y, max(grid.x), grid.y, xpd= T)
  text(matrix(text.x, nrow= nrow(x), ncol= ncol(x), byrow = T),
       matrix(text.y, nrow= nrow(x), ncol= ncol(x)),
       as.matrix(x)[nrow(x):1,],
       cex= cex.t,
       xpd= T)
}
