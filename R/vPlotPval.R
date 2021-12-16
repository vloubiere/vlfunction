#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars_only If set to TRU, then only plots */N.S. Default= FALSE
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export

vl_plot_pval_text <- function(x, 
                              y, 
                              pval, 
                              stars_only= F)
{
  if(length(y)==1 & length(x)>1)
    y <- rep(y, length(x))
  star <- cut(pval, 
              breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
              labels = c("****", "***", "**", "*", "N.S"))
  star <- as.character(star)
  ns_val <- star=="N.S"
  value <- formatC(pval, format = "e", digits = 1)
  if(stars_only)
    value <- rep("", length(value))
  if(any(ns_val))
    text(x = x[ns_val], 
         y = y[ns_val], 
         labels= bquote(paste(.(value[ns_val]), ""^N.S)), 
         cex= 0.6,
         offset= -0.1,
         pos= 3)
  if(!all(ns_val))
  text(x = x[!ns_val], 
       y = y[!ns_val], 
       labels= paste0(value[!ns_val], star[!ns_val]),
       offset= -0.2,
       pos= 3)
}

#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x list corresponding to the different groups
#' @param compare A list of length 2 vectors specifying the groups to compare (can be either integers or group names)
#' @param adj Ratio of the boxplots' range that will be used to space pval segments
#' @param outline Should oultiers be ploted?
#' @param las boxplot las
#' @param xlab boxplot xlab
#' @param ylim boxplot ylim
#' @param staplewex boxplot staplewex
#' @param whisklty boxplot whisklty
#' @param boxwex boxplot boxwex
#' @param ... Extra parameters to be passed to boxplot
#' @examples 
#' y <- split(InsectSprays$count, InsectSprays$spray)
#' vl_boxplot_pval(y, 
#' compare = list(c(1,2), c(1,3), c(1,5), c(4,5), c(5, 6)))
#' vl_boxplot_pval(y, 
#' compare = list(c("A","B"), c("A","C"), c("A","E"), c("D","E"), c("E", "F")))
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export
vl_boxplot_pval <- function(x,
                            compare,
                            adj= 10,
                            outline= F, 
                            las= 2,
                            xlab= NA,
                            ylim= c(min(box$stats, na.rm = T), max(obj$y)+inter),
                            staplewex = NA, 
                            whisklty = 1, 
                            boxwex = 0.5, ...)
{
  # Checks
  if(!is.list(x))
    stop("x should be a (named) list of values")
  if(is.numeric(unlist(compare)))
  {
    if(!all(between(unlist(compare), 1, length(x))))
      stop("compare contains integer either<1 or >length(x) that could not be matched to any sublist of x")
  }else if(!all(unlist(compare) %in% names(x)))
    stop("Values in compare do not all match names(x)")
  
  # Compute boxplot stats
  box <- boxplot(x, 
                 plot= F)
  
  # Make object with necessary checks
  obj <- data.table(name= box$names,
                    max= box$stats[5,])
  
  # If compare was specified using group names, find corresponding idx in x and vice versa
  obj <- rbindlist(lapply(compare, function(i) 
  {
    if(!is.numeric(i))
      i <- match(i, obj$name)
    i <- sort(i)
    .c <- obj[i[1]:i[2], .(x= i, cdition= name[c(1, .N)], max= max(max))]
    .c$var <- x[i]
    return(.c)
  }), idcol = T)
  obj[, idx:= seq(.N), .id]
  obj <- dcast(obj, .id+max~idx, value.var = list("cdition", "var", "x"))
  
  # Wilcoxon
  obj[, pval:= wilcox.test(unlist(var_1), unlist(var_2))$p.value, .id]
  
  # Compute Y pos
  setorderv(obj, "max")
  inter <- (max(box$stats)-min(box$stats))/adj
  obj[1, y:= max+inter]
  for(i in 2:nrow(obj))
  {
    overlap <- obj$y[obj[i, x_1]<=obj$x_2 & obj[i, x_2]>=obj$x_1]
    obj[i, y:= max(c(overlap, max), na.rm = T)+inter]
  }

  # Plot
  boxplot(x, 
          outline= outline, 
          las= las,
          xlab= xlab,
          ylim= ylim,
          staplewex = staplewex, 
          whisklty = whisklty, 
          boxwex = boxwex, ...)
  obj[, {
    segments(x_1[1], y, x_2[1], y)
    vl_plot_pval_text(mean(c(x_1, x_2)), 
                      y,
                      pval,
                      stars_only= T)
  }, .id]
  invisible(obj)
}
