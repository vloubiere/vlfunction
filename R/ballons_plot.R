#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#' 
#' @param x Matrix specifyiung size of the balloons.
#' @param color_var Matrix specifyiung color of the balloons.
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param cex.balloons Expansion factor for balloons. default= 1.
#' @param main Title. Default= NA
#' @param balloon_size_legend Title size legend
#' @param balloon_col_legend Title color legend
#' @examples
#' x <- matrix(-3:5, ncol= 3)
#' cols <- matrix(-3:5, ncol= 3)
#' vl_balloons_plot(x, cols, balloon_size_legend = "Size", balloon_col_legend = "Color")
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
                                    col= c("cornflowerblue", "white", "tomato"),
                                    cex.balloons= 1,
                                    main= NA, 
                                    balloon_size_legend= NA,
                                    balloon_col_legend= NA,
                                    auto_margins= T)
{
  # Checks
  if(missing(x_breaks))
    x_breaks <- axisTicks(range(x, na.rm= T), log= F, nint = 4)
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

  # Init plot
  plot.new()
  plot.window(xlim = c(1, ncol(x)),
              ylim = c(1, nrow(x)),
              xaxs= "r")
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
         cex= abs(x)*cex.balloons,
         xpd= T)
  # Axes
  axis(side= 1, 
       at= seq(ncol(x)),
       labels = colnames(x),
       lwd= NA)
  axis(side= 2, 
       at= seq(nrow(x)),
       labels = rownames(x),
       lwd= NA)
  # Legends
  left <- par("usr")[2]+(grconvertX(adj, "in", "user")-grconvertX(0, "in", "user"))
  vl_heatkey(color_breaks, 
             left= left,
             col, 
             top= nrow(x)-strheight("M"),
             main= balloon_col_legend)
  title(main= main)
  vl_balloonskey(sizes = x_breaks*cex.balloons,
                 labels = x_breaks,
                 left= left,
                 top= nrow(x)-strheight("M")*10, 
                 main = balloon_size_legend)
}