#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param dat Data.table or matrix containing 1 column/condition
#' @param col Colors to use for bars. Wrapped with colorRampPalette
#' @param alpha.f Alpha value for transition colors
#' @param space Space between bars. See ?barplot()
#' @param width Bars width. See ?barplot()
#' @param ... Extra parameters to be passed to barplot
#' @examples 
#' test <- data.table(c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
#'                    c("B", "A", "A", "B", "C", "B", "B", "C", "C", "A", "A"),
#'                    c("A", "A", "B", "B", "B", "B", "B", "B", "B", "C", "A"))
#' vl_alluvial_plot(test)
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export
vl_alluvial_plot <- function(dat,
                             col= rainbow(length(unique(unlist(dat)))),
                             show.bar.labels= TRUE,
                             labels.cex= 1,
                             show.legend= FALSE,
                             alpha.f= 0.5,
                             space= 1,
                             width= 0.5,
                             ...)
{
  # Make object ----
  DT <- melt(dat, 
             measure.vars = names(dat))

  # Bars ----
  setorderv(DT, c("variable", "value"))
  DT[, variable:= factor(variable, unique(variable))]
  DT[, value:= factor(value, unique(value))]
  mat <- dcast(DT, value~variable, fun.aggregate = length)
  mat <- as.matrix(mat, 1)
  mat <- mat[nrow(mat):1, , drop= FALSE]
  pl <- barplot(mat,
                col= col,
                space= space,
                width= width,
                ...)
  
  # Connections ----
  for(i in 1:(ncol(dat)-1))
  {
    pols <- dat[, i:(i+1), with= F]
    setnames(pols, c("V1", "V2"))
    setorderv(pols, c("V1", "V2"))
    pols[, c("y1", "y2"):= as.list(range(.I)), .(V1, V2)]
    setorderv(pols, c("V2", "V1"))
    pols[, c("y3", "y4"):= as.list(rev(range(.I))), .(V2, V1)]
    pols[, col:= adjustcolor(col[.GRP], alpha.f = alpha.f), V1]
    browser()
    pols[, polygon(c(pl[i]+width/2, pl[i]+width/2, pl[i+1]-width/2, pl[i+1]-width/2), 
                   c(y1[1]-1, y2[1], y3[1], y4[1]-1),
                   col= col[1]), .(y1, y2, y3, y4, col)]
  }
  
  # Print names
  if(show.bar.labels) {
    for(i in seq(ncol(mat))) {
      .c <- na.omit(mat[mat[,i]>0, i])
      text(pl[i],
           cumsum(.c)-.c/2,
           names(.c),
           cex= labels.cex)
    }
  }

  # Print legend ----
  if(show.legend) {
    legend(par("usr")[2],
           par("usr")[4],
           fill= col,
           legend = rownames(mat),
           bty= "n",
           xpd= T)
  }
  
  # Return ----
  invisible(pl)
}