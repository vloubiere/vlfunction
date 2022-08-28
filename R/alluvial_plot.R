#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param dat Data.table or matrix containing 1 column/condition
#' @param col Colors to use for bars. Wrapped with colorRampPalette
#' @param alpha.f Alpha value for transition colors
#' @param ... Extra parameters to be passed to barplot
#' @examples 
#' test <- data.table(c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
#' c("A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"),
#' c("A", "A", "B", "B", "B", "B", "B", "B", "B", "C", "C"))
#' vl_alluvial_plot(test)
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export
vl_alluvial_plot <- function(dat,
                             col= vl_palette_many_categ(length(unique(unlist(dat)))),
                             alpha.f= 0.5,
                             ...)
{
  # Make object
  DT <- melt(dat, 
             measure.vars = names(dat))
  # Bars
  setorderv(DT, c("variable", "value"))
  mat <- as.matrix(dcast(DT, value~variable, fun.aggregate = length), 1)
  pl <- barplot(mat,
                col= col,
                width = 0.5, 
                space= 1)
  for(i in 1:(ncol(dat)-1))
  {
    pols <- dat[, i:(i+1), with= F]
    setnames(pols, c("V1", "V2"))
    setorderv(pols, c("V1", "V2"))
    pols[, c("y1", "y2"):= as.list(range(.I)), .(V1, V2)]
    setorderv(pols, c("V2", "V1"))
    pols[, c("y3", "y4"):= as.list(rev(range(.I))), .(V2, V1)]
    pols[, col:= adjustcolor(col[.GRP], alpha.f = alpha.f), V2]
    pols[, polygon(c(i, i, i+0.5, i+0.5), 
                   c(y1-1, y2, y3, y4-1),
                   col= col[1]), .(y1, y2, y3, y4, col)]
  }
  legend(par("usr")[2],
         par("usr")[4],
         fill= col,
         legend = rownames(mat),
         bty= "n",
         xpd= T)
  
  invisible(pl)
}