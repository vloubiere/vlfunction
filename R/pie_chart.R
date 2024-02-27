
#' Title
#'
#' @param x Variable to be summarized
#' @param col colors to be used
#' @param ... extra arguments passed to pie(x, ...)
#'
#' @return Plots a pie chart
#' @export
#'
#' @examples
#' vl_pie(rep(1:3, 1:3))
vl_pie <- function(x, col= adjustcolor(rainbow(length(unique(x))), .3), ...)
{
  .c <- table(x)
  names(.c) <- paste0(names(.c), " (", formatC(.c, big.mark = ","), ")")
  pie(.c, col= col, ...)
}