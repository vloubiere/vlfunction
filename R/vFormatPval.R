#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars_only If set to TRU, then only plots */N.S. Default= FALSE
#' @param ... extra arguments for text() function (used for plotting)
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' vl_plot_pval(1:6, 0, pval, cex= 0.6)
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export

vl_plot_pval <- function(x, y, pval, stars_only= F, ...)
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
  text(x = x[ns_val], 
       y = y[ns_val], 
       labels= parse(text= paste0(value[ns_val], "^N.S")), 
       ...)
  text(x = x[!ns_val], 
       y = y[!ns_val], 
       labels= paste0(value[!ns_val], star[!ns_val]), 
       ...)
}
