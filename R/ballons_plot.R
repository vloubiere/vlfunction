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
#' @param legend_left_pos If specified, overrides default legend left user coordinates
#' @examples
#' x <- matrix(-3:5, ncol= 3)
#' cols <- matrix(-3:5, ncol= 3)
#' par(mar= c(3,3,2,5), las= 1)
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
                                    col= c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF"),
                                    cex.balloons= 1,
                                    main= NA, 
                                    balloon_size_legend= NA,
                                    balloon_col_legend= NA,
                                    legend_left_pos,
                                    gap.axis= 0)
{
  # Checks
  if(missing(x_breaks))
  {
    x_breaks <- range(x, na.rm= T)
    if(length(unique(x_breaks))==1)
      x_breaks <- x_breaks+c(-0.5, 0.5)
    x_breaks <- axisTicks(x_breaks, log= F, nint = 4)
  }
  if(missing(color_breaks))
  {
    color_breaks <- range(color_var, na.rm= T)
    if(length(unique(color_breaks))==1)
      color_breaks <- color_breaks+c(-0.5, 0.5)
    color_breaks <- seq(min(color_breaks, na.rm= T),
                        max(color_breaks, na.rm= T),
                        length.out= length(col))
  }
  
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
              xaxs= "i",
              yaxs= "i")
  # Lines
  segments(1:ncol(x),
           1,
           1:ncol(x),
           nrow(x),
           xpd=T)
  segments(1,
           1:nrow(x),
           ncol(x),
           1:nrow(x),
           xpd=T)
  # Balloons
  points(rep(1:ncol(x), each= nrow(x)),
         rep(1:nrow(x), ncol(x)),
         col= color_var,
         pch= ifelse(x>=0, 19, 15),
         cex= abs(x)*cex.balloons,
         xpd= T)
  # Axes
  vl_tilt_xaxis(seq(ncol(x)),
       labels = colnames(x))
  axis(side= 2, 
       at= seq(nrow(x)),
       labels = rownames(x),
       gap.axis= gap.axis,
       lwd= NA)
  # Title
  title(main= main)
  # Legends
  if(missing(legend_left_pos))
  {
    legend_left_pos <- strwidth("M")*0.75/2*max(abs(x), na.rm= T)*cex.balloons # max balloon radius
    legend_left_pos <- grconvertX(1, "npc", "user")+legend_left_pos # add to right plot limit
  }
  vl_heatkey(color_breaks, 
             left= legend_left_pos,
             col, 
             top= nrow(x)-strheight("M"),
             main= balloon_col_legend)
  vl_balloonskey(sizes = x_breaks*cex.balloons,
                 labels = x_breaks,
                 left= legend_left_pos,
                 top= nrow(x)-strheight("M")*10, 
                 main = balloon_size_legend)
}