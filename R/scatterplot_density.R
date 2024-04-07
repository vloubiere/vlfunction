require(data.table)


#' Title
#'
#' @param x x variable
#' @param y y variable
#' @param label Labels to split the data
#' @param xlim x axis lims
#' @param ylim y axis lims
#' @param col Colors to be used (must match number of levels in label)
#' @param xlab x axis labels
#' @param ylab y axis labels
#' @param dens.lw Width of the density lines to be drawn on the top and right of the scatterplot
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
                               frame= F,
                               ...)
{
  if(!is.factor(label))
    label <- factor(label)
  if(length(col) != length(levels(label)))
    stop("Number of colors and label levels should be equal")
  dat <- data.table(x= x,
                    y= y,
                    label= label,
                    Cc= col[label])
  
  dat[, {
    plot(x,
         y,
         col= Cc,
         xlab= xlab,
         ylab= ylab,
         frame= frame,
         ...)
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    densities <- .SD[, .(dx= .(density(x, from= xlim[1], to= xlim[2])),
                         dy= .(density(y, from= ylim[1], to= ylim[2]))), .(label, Cc)]
    
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
  }, ]
  legend(par("usr")[2],
         par("usr")[4]+strheight("M")*5*dens.lw,
         fill= col,
         legend = levels(label),
         xpd= T,
         bty= "n")
}