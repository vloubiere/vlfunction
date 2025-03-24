require(data.table)


#' Scatterplot with density lines on top and on the right
#'
#' @param x x variable
#' @param y y variable
#' @param label Labels to split the data
#' @param xlim x axis lims
#' @param ylim y axis lims
#' @param col Colors to be used (must match number of levels in label)
#' @param xlab x axis labels
#' @param ylab y axis labels
#' @param dens.lw densities width in lines. Default= 1L.
#' @param plot.legend Should the legend be plotted?
#' @param legend.x x position for the legend. Default= par("usr")[2].
#' @param legend.y y position for the legend. Default= par("usr")[4]+strheight("M")*5*dens.lw.
#' @param legend.cex Expansion factor for the legend. Default= .7.
#' @param useRaster Default= T
#' @param ... Extra arguments to be passed to plot().
#'
#' @return A scatterplot with density lines on the side
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- data.table(x= rnorm(1000, 1, 0.5),
#' y= rnorm(1000, 1, 0.5),
#' label= sample(rep(c("A", "B"), each= 500), 1000))
#' par(mai= c(0.9, 0.9, 0.9, 0.9),
#' tcl= -0.1,
#' mgp= c(2, .5, 0),
#' las= 1)
#' vl_densScatterplot(dat$x, dat$y, dat$label)
#' 
vl_densScatterplot <- function(x,
                               y,
                               label,
                               xlim= range(x, na.rm = T),
                               ylim= range(y, na.rm = T),
                               col= rainbow(length(unique(label))),
                               xlab= deparse1(substitute(x)),
                               ylab= deparse1(substitute(y)),
                               dens.lw= 1,
                               plot.legend= TRUE,
                               legend.x= par("usr")[2],
                               legend.y= par("usr")[4]+strheight("M")*5*dens.lw,
                               legend.cex= .7,
                               frame= FALSE,
                               useRaster= FALSE,
                               ...)
{
  # Checks
  if(!is.factor(label))
  {
    label <- factor(label, unique(label))
    message("label converted to factor")
  }
  if(length(col)==length(levels(label)))
    col <- col[label]
  
  # Make data.table
  dat <- data.table(x= x,
                    y= y,
                    label= label,
                    Cc= col)
  # Scatterplot
  dat[, {
    if(useRaster)
    {
      vl_rasterScatterplot(x= x,
                           y= y,
                           col= Cc,
                           xlab= xlab,
                           ylab= ylab,
                           xlim= xlim,
                           ylim= ylim,
                           frame= frame,
                           ...)
    }else
    {
      plot(x,
           y,
           col= Cc,
           xlab= xlab,
           ylab= ylab,
           xlim= xlim,
           ylim= ylim,
           frame= frame,
           ...)
    }
    
    if(length(x)>1 && !is.na(dens.lw)) {
      # Densities
      xlim <- par("usr")[1:2]
      ylim <- par("usr")[3:4]
      densities <- .SD[, .(dx= .(density(x, from= xlim[1], to= xlim[2], na.rm= TRUE)),
                           dy= .(density(y, from= ylim[1], to= ylim[2], na.rm= TRUE))), .(label, Cc)]
      
      densities[, {
        # x
        .c <- dx[[1]]
        .c$y <- .c$y/max(.c$y)
        lines(.c$x,
              ylim[2]+(.c$y*strheight("M")*5*dens.lw),
              col= Cc[1],
              xpd= T)
        # y
        .c <- dy[[1]]
        .c$y <- .c$y/max(.c$y)
        lines(xlim[2]+(.c$y*strwidth("M")*5*dens.lw),
              .c$x,
              col= Cc[1],
              xpd= T)
        .SD
      }, .(Cc, label)]
    }
  }, ]
  
  # Legend
  if(plot.legend)
    legend(par("usr")[2],
           par("usr")[4]+strheight("M")*5*dens.lw,
           fill= unique(col[order(label)]),
           legend = levels(label),
           xpd= T,
           bty= "n",
           border= NA)
}