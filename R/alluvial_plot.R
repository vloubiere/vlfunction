#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param dat Data.table containing the variables to plot
#' @param col Colors to use for bars. Wrapped with colorRampPalette
#' @param ... Extra parameters to be passed to barplot
#' @examples 
#' test <- data.table(c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
#' c("A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"),
#' c("A", "A", "B", "B", "B", "B", "B", "B", "B", "C", "C"))
#' vl_alluvial_plot(test,
#' class_levels = c("C", "B", "A"))
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export

vl_alluvial_plot <- function(dat,
                             col= vl_palette_categ2(length(unique(unlist(dat)))),
                             ...)
{
  # Make object
  res <- melt(dat, 
              measure.vars = names(dat),  
              value.factor = T)
  # Bars
  setorderv(res, c("variable", "value"))
  mat <- as.matrix(dcast(res, value~variable, fun.aggregate = length), 1)
  if(is.null(col))
    col <- vl_palette_categ2(nrow(mat))
  par(mai= c(1.02,
             0.82,
             0.82,
             max(strwidth(names(dat), "inch"))+1))
  barplot(mat, 
          col= col,
          width = 0.5, 
          space= 1)
  cols <- data.table(class= rownames(mat), 
                     Cc= col)
  for(i in 1:(ncol(dat)-1))
  {
    pols <- dat[, i:(i+1), with= F]
    setnames(pols, c("V1", "V2"))
    setkeyv(pols, c("V1", "V2"))
    pols[, c("y0", "y3"):= .(.I-1, .I)]
    pols[, c("y0", "y3"):= .(data.table::first(y0), 
                             data.table::last(y3)), .(V1, V2)]
    setkeyv(pols, c("V2", "V1"))
    pols[, c("y1", "y2"):= .(.I-1, .I)]
    pols[, c("y1", "y2"):= .(data.table::first(y1), 
                             data.table::last(y2)), .(V2, V1)]
    pols[cols, Cc:= adjustcolor(i.Cc, 0.5), on= "V1==class"]
    pols[, polygon(c(i, i+0.5, i+0.5, i), 
                   c(y0, y1, y2, y3),
                   col= Cc[1]), .(y0, y1, y2, y3, Cc)]
  }
  legend(par("usr")[2],
         par("usr")[4],
         fill= cols$Cc,
         legend = cols$class,
         bty= "n",
         xpd= T)
}