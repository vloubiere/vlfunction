#' box plot
#'
#' This function is just a wrapper around default boxplot function with a nicer layout
#'
#' @param x Matrix, data.atble or data.frame containing data
#' @param las.x las x axis labels
#' @param las.y las y axis labels
#' @param xlab x axis label
#' @param ylab y axis label
#' @examples
#' test <- data.table(test= c(1:10, 50, 50, 90:100),
#' test1= c(1:10, 50, 50, 90:100)+10)
#' vl_boxplot(test, xlab = "test", ylab = "test")
#' @export

vl_boxplot <- function(x, las.x= 1, las.y= 1, xlab= NA, ylab= NA, ...)
{
  
  par(mgp= c(3,0.5,1))
  b <- boxplot(x,
               lwd.axis= 0,
               xaxt= "n", 
               yaxt= "n", 
               staplewex= 0,
               whisklty= 1, 
               boxwex= 0.5, ...)
  axis(side = 1, 
       at = seq_along(b$names), 
       labels = b$names,
       las= las.x,
       tick = FALSE)
  mtext(xlab, 1, 2)
  axis(2, 
       line= 0, 
       lwd= 0,
       lwd.ticks= 1,
       tck= -0.015,
       las= las.y)
  mtext(ylab, 2, 2)
  par(mgp= c(3, 1, 0))
}