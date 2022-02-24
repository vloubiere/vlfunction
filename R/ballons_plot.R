#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#' 
#' @param x Matrix specifyiung size of the balloons.
#' @param color_var Matrix specifyiung color of the balloons.
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param cex.balloons Expansion factor for balloons
#' @param main Title. Default= NA
#' @param balloon_size_legend Title size legend
#' @param balloon_col_legend Title color legend
#' @param auto_margins Use auto margins? Default= T
#' @examples
#' x <- matrix(-3:5, ncol= 3)
#' cols <- matrix(-3:5, ncol= 3)
#' vl_balloons_plot(x, cols, balloon_size_legend = "Size", balloon_col_legend = "Color")
#' 
#' @return Balloon plot
#' @export
#' 
vl_balloons_plot <- function(x, ...)
  UseMethod("vl_balloons_plot")


#' @describeIn vl_balloons_plot matrix_method
#' @export

vl_balloons_plot.matrix <- function(x,
                                    color_var,
                                    x_breaks,
                                    color_breaks,
                                    col= c("cornflowerblue", "white", "tomato"),
                                    cex.balloons= 1,
                                    main= NA, 
                                    balloon_size_legend= NA,
                                    balloon_col_legend= NA,
                                    auto_margins= T)
{
  # Checks
  if(missing(x_breaks))
    x_breaks <- axisTicks(range(x, na.rm= T), log= F)
  if(missing(color_breaks))
    color_breaks <- seq(min(color_var, na.rm= T),
                        max(color_var, na.rm= T),
                        length.out= length(col))
  
  # Revert matrices for plotting
  x <- x[nrow(x):1,,drop= F]
  color_var <- color_var[nrow(color_var):1,,drop= F]
  
  # Compute colors
  Cc <- circlize::colorRamp2(color_breaks, col)
  color_var[!is.na(color_var)] <- Cc(color_var[!is.na(color_var)])
  color_var[is.na(color_var)] <- "lightgrey"
      
  # Margins
  if(auto_margins)
  {
    bot <- 0.5+max(strwidth(colnames(x), "inches"))
    left <- 0.5+max(strwidth(rownames(x), "inches", cex = par("cex.axis")/par("cex")))
    top <- 0.5+strheight(main, units = "inches", cex = par("cex.axis")/par("cex"))
    leg.width <- strwidth(c(balloon_size_legend, balloon_col_legend), "in")
    right <- 0.5+max(c(leg.width, 0.5))
    par(mai= c(bot, left, top, right),
        xaxs= "i",
        yaxs= "i")
  }
  
  # Init plot
  plot.new()
  adj <- (par("cin")[1]*0.75/2*max(abs(x), na.rm=T)*cex.balloons*par("cex"))/par("pin")[1]
  adj.x <- adj*(ncol(x)-1)
  adj.y <- adj*(nrow(x)-1)
  plot.window(xlim = c(1-adj.x, ncol(x)+adj.x),
              ylim = c(1-adj.y, nrow(x)+adj.y))
  # Lines
  segments(1:ncol(x),
           1,
           1:ncol(x),
           nrow(x))
  segments(1,
           1:nrow(x),
           ncol(x),
           1:nrow(x))
  # Balloons
  points(rep(1:ncol(x), each= nrow(x)),
         rep(1:nrow(x), ncol(x)),
         bg= color_var,
         pch= ifelse(x>=0, 21, 22),
         cex= abs(x)*cex.balloons+0.1)
  # Legends
  axis(side= 1, 
       at= seq(ncol(x)),
       labels = colnames(x),
       lwd= NA,
       line= 0)
  axis(side= 2, 
       at= seq(nrow(x)),
       labels = rownames(x),
       lwd= NA,
       line= 0)
  vl_heatkey(color_breaks, 
             col, 
             top= nrow(x)-strheight("M"),
             main= balloon_col_legend)
  mtext(main)
  vl_balloonskey(sizes = x_breaks*cex.balloons,
                 labels = x_breaks,
                 top= nrow(x)-strheight("M")*10, 
                 main = balloon_size_legend)
}