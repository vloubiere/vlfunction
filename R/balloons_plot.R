#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#' 
#' @param x Matrix specifyiung size of the balloons.
#' @param color.var Matrix specifyiung color of the balloons.
#' @param color.breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param cex.balloons Expansion factor for balloons. default= 1.
#' @param main Title. Default= NA
#' @param xaxt see ?par()
#' @param yaxt see ?par()
#' @param balloon.size.legend Title size legend
#' @param balloon.col.legend Title color legend
#' @param legend.left.pos If specified, overrides default legend left user coordinates
#' @examples
#' x <- color.var <- matrix(-3:5, ncol= 3)
#' vl_balloons_plot(x= x,
#'                  color.var= color.var)
#' 
#' @return Balloon plot
#' @export
vl_balloons_plot <- function(x, ...)
  UseMethod("vl_balloons_plot")


#' @describeIn vl_balloons_plot matrix_method
#' @export
vl_balloons_plot.matrix <- function(x,
                                    color.var,
                                    x.breaks,
                                    color.breaks,
                                    col= c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF"),
                                    cex.balloons= 1,
                                    main,
                                    xaxt= "s",
                                    yaxt= "s",
                                    balloon.size.legend= "Size",
                                    balloon.col.legend= "Color",
                                    legend.left.pos,
                                    gap.axis= 0)
{
  # Checks
  if(missing(x.breaks))
  {
    x.breaks <- range(x, na.rm= T)
    if(length(unique(x.breaks))==1)
      x.breaks <- x.breaks+c(-0.5, 0.5)
    x.breaks <- axisTicks(x.breaks, log= F, nint = 4)
  }
  if(missing(color.breaks))
  {
    color.breaks <- range(color.var, na.rm= T)
    if(length(unique(color.breaks))==1)
      color.breaks <- color.breaks+c(-0.5, 0.5)
    color.breaks <- seq(min(color.breaks, na.rm= T),
                        max(color.breaks, na.rm= T),
                        length.out= length(col))
  }
  
  # Revert matrices for plotting
  x <- x[nrow(x):1,,drop= F]
  color.var <- color.var[nrow(color.var):1,,drop= F]
  
  # Compute colors
  Cc <- circlize::colorRamp2(color.breaks, col)
  color.var[!is.na(color.var)] <- Cc(color.var[!is.na(color.var)])
  color.var[is.na(color.var)] <- "lightgrey"

  # Init plot
  plot.new()
  plot.window(xlim = c(1, ncol(x)),
              ylim = c(1, nrow(x)),
              xaxs= "i",
              yaxs= "i",
              xaxt= xaxt,
              yaxt= yaxt)
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
         col= color.var,
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
  if(!missing(main))
    title(main= main)
  # Legends
  if(missing(legend.left.pos))
  {
    legend.left.pos <- strwidth("M")*0.75/2*max(abs(x), na.rm= T)*cex.balloons # max balloon radius
    legend.left.pos <- grconvertX(1, "npc", "user")+legend.left.pos # add to right plot limit
  }
  vl_heatkey(color.breaks, 
             left= legend.left.pos,
             col, 
             top= nrow(x)-strheight("M"),
             main= balloon.col.legend)
  vl_balloonskey(sizes = x.breaks*cex.balloons,
                 labels = x.breaks,
                 left= legend.left.pos,
                 top= nrow(x)-strheight("M")*10, 
                 main = balloon.size.legend)
}