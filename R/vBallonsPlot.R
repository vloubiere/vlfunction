#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#' 
#' @param x Matrix specifyiung size of the balloons.
#' @param color_var Matrix specifyiung color of the balloons.
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param balloon_size_legend Title size legend
#' @param balloon_col_legend Title color legend
#' @param cex.balloons scaling factor for balloons size
#' @param auto_margins Use auto margins? Default= T
#' @examples
#' mat <- matrix(0:11, ncol = 3)
#' col <- matrix(6:-5, ncol = 3)
#' vl_balloons_plot(mat)
#' vl_balloons_plot(x= mat, color_var = col)
#' vl_balloons_plot(x= mat,
#' color_var = col,
#' balloon_size_legend = "test",
#' balloon_col_legend = "test2")
#' vl_balloons_plot(x= mat,
#' color_var = col,
#' balloon_size_legend = "test",
#' balloon_col_legend = "test2", 
#' cex.balloons = 8,
#' main= "Check!")
#' 
#' @return Balloon plot
#' @export
vl_balloons_plot <- function(x, ...)
  UseMethod("vl_balloons_plot")

#' @describeIn vl_balloons_plot matrix_method
#' @export
vl_balloons_plot.matrix <- function(x,
                                    color_var,
                                    x_breaks,
                                    color_breaks,
                                    col= c("cornflowerblue", "lightgrey", "tomato"),
                                    main= NA, 
                                    balloon_size_legend= NA,
                                    balloon_col_legend= NA,
                                    cex.balloons= 4,
                                    auto_margins= T)
{
  # Functions
  cex_scale <- function(x)
    (x-min(x_breaks, na.rm = T))/(max(x_breaks, na.rm = T)-min(x_breaks, na.rm = T))*cex.balloons+1
  compute_balloon_radius <- function(cex= cex.balloons)
    grconvertX(0.75/2, "chars", "in")*cex/2
  
  # Checks
  if(missing(x_breaks))
    x_breaks <- range(x, na.rm= T)
  if(missing(color_var))
  {
    color_var <- matrix(rep(1, nrow(x)*ncol(x)), 
                        nrow = nrow(x),
                        ncol = ncol(x))
    color_breaks <- c(0,1,2)
  }
  color_var <- color_var[nrow(color_var):1,]
  if(!identical(dim(x), dim(color_var)))
    stop("x and color_var matrices should have identical dimensions")
  if(missing(color_breaks))
    color_breaks <- seq(min(color_var, na.rm= T), 
                        max(color_var, na.rm= T), 
                        length.out= length(col))
  Cc <- circlize::colorRamp2(color_breaks, col)
  
  
  # Scale data
  scaled <- cex_scale(x)
  if(is.null(colnames(scaled)))
    colnames(scaled) <- seq(ncol(scaled))
  if(is.null(rownames(scaled)))
    rownames(scaled) <- seq(nrow(scaled))
  scaled <- scaled[nrow(scaled):1,]
  
  # Margins
  if(auto_margins)
  {
    bot <- 0.5+max(strwidth(colnames(x), "inches"))
    left <- 0.5+max(strwidth(rownames(x), "inches"))
    top <- 0.5+strheight(main, units = "inches")
    leg.width <- strwidth(c(balloon_size_legend, balloon_col_legend), "in")
    right <- 0.5+max(c(leg.width, 0.5))
    par(mai= c(bot, left, top, right),
        xaxs= "i",
        yaxs= "i")
  }
  # Init plot
  plot.new()
  # Compute xlim
  ball.radius <- compute_balloon_radius(cex= cex.balloons+1)
  plot.size.x <- grconvertX(1, "npc", "in")-grconvertX(0, "npc", "in")
  adj.x <- (ball.radius*1.3)/plot.size.x
  ext.x <- (ncol(scaled)-1)*adj.x
  xl <- c(1-ext.x, ncol(scaled)+ext.x)
  # Compute ylim
  plot.size.y <- grconvertY(1, "npc", "in")-grconvertY(0, "npc", "in")
  adj.y <- (ball.radius*1.3)/plot.size.y
  ext.y <- (nrow(scaled)-1)*adj.y
  yl <- c(1-ext.y, nrow(scaled)+ext.y)
  # Plot window
  plot.window(xlim = xl,
              ylim = yl)
  # Title
  title(main)
  # Grid
  segments(1, 
           seq(nrow(scaled)), 
           ncol(scaled),
           seq(nrow(scaled)))
  segments(seq(ncol(scaled)),
           1,
           seq(ncol(scaled)),
           nrow(scaled))
  # Points
  col_vec <- c(color_var)
  col_vec[!is.na(col_vec)] <- Cc(col_vec[!is.na(col_vec)])
  points(col(scaled),
         row(scaled),
         cex= c(scaled),
         pch= 21,
         bg= col_vec,
         col= "black")
  # Axes
  axis(1,
       at= seq(ncol(scaled)),
       labels = colnames(x),
       las= 2,
       lwd= 0,
       line= -0.5)
  axis(2,
       at= seq(nrow(scaled)),
       labels = rownames(scaled),
       las= 2,
       lwd= 0,
       line= -0.5)

  # Size legend
  x.leg <- grconvertX(1, "npc", "in")+grconvertX(0.5, "line", "in")+ball.radius
  x.leg.title <- grconvertX(x.leg-ball.radius, "in", "user")
  x.leg.text <- grconvertX(x.leg+ball.radius, "in", "user")
  x.leg <- grconvertX(x.leg, "in", "user")
  b.top <- grconvertY(nrow(scaled), "user", "in")-grconvertY(0.5, "line", "in")-ball.radius
  scale.values <- rev(axisTicks(x_breaks, log= F, nint = 4))
  scale.cex <- cex_scale(scale.values)
  scale.radius <- compute_balloon_radius(cex= scale.cex)
  scale.y <- sapply(seq(scale.cex[-1]), function(i) scale.radius[i]+scale.radius[i+1]+0.05)
  scale.y <- grconvertY(b.top-cumsum(c(0, scale.y)), "in", "user")
  points(rep(x.leg, length(scale.y)),
         scale.y, 
         cex= scale.cex, 
         xpd= T)
  text(x.leg.text, 
       scale.y, 
       scale.values,
       xpd= T, 
       pos= 4, 
       cex= 0.6,
       offset = 0.25)
  text(x.leg.title, 
       nrow(scaled), 
       balloon_size_legend,
       xpd= T, 
       pos= 4,
       offset = 0.5)
  
  # Color legend
  x.left <- grconvertX(x.leg, "user", "in")-ball.radius/2
  x.right <- x.left+grconvertX(1, "line", "in")
  y.leg <- grconvertY(min(scale.y), "user", "in")-grconvertY(2, "line", "in")
  y.top <- grconvertY(min(scale.y), "user", "in")-grconvertY(3, "line", "in")
  y.bot <- y.top-grconvertY(5, "line", "in")
  x.left <- grconvertX(x.left, "in", "user")
  y.bot <- grconvertY(y.bot, "in", "user")
  y.leg <- grconvertY(y.leg, "in", "user")
  x.right <- grconvertX(x.right, "in", "user")
  y.top <- grconvertY(y.top, "in", "user")
  col.im <- matrix(Cc(seq(max(color_breaks), min(color_breaks), length.out = 100)), ncol= 1)
  rasterImage(col.im,
              x.left, 
              y.bot, 
              x.right,
              y.top, 
              xpd= T)
  tick.lab <- axisTicks(range(color_breaks), log= F, nint = 4)
  tick.y <- y.bot+(y.top-y.bot)*(tick.lab-min(color_breaks))/(max(color_breaks)-min(color_breaks))
  text(rep(x.right, length(tick.y)), 
       tick.y,
       labels = tick.lab,
       cex= 0.6,
       pos= 4, 
       xpd= T, 
       offset = 0.25)
  text(x.leg.title,
       y.leg,
       labels = balloon_col_legend,
       pos= 4, 
       xpd= T, 
       offset = 0.5)
}
