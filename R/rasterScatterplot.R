#' Title
#'
#' @param x 
#' @param y default= NULL
#' @param type default= "p"
#' @param frame Default= F
#' @param size Size of the png image. Default= 2000
#' @param res Nominal resolution of png image. default= 600
#' @param ... Extra plotting parameters
#'
#' @return Raster scatterplot
#' @export
#'
#' @examples
#' vl_rasterScatterplot(1:3)
vl_rasterScatterplot <- function(x,
                                 y= NULL,
                                 type= "p",
                                 frame= F,
                                 res= 600L,
                                 size= 2000L,
                                 xlab= ifelse(is.null(y), "Index", deparse1(substitute(x))),
                                 ylab= ifelse(is.null(y), deparse1(substitute(x)), deparse1(substitute(y))),
                                 add= F,
                                 ...)
{
  # Initialize plot ----
  if(!add)
    plot(x= x,
         y= y, 
         frame= frame,
         type= "n",
         xlab= xlab,
         ylab= ylab,
         ...)
  
  # Extract plot area in both user and physical coordinates ----
  coords <- par("usr")
  gx <- grconvertX(coords[1:2], "user", "inches")
  gy <- grconvertY(coords[3:4], "user", "inches")
  width <- diff(gx)
  height <- diff(gy)
  ratio <- round(c(width, height)/max(c(width, height))*size)
  
  # Save as png ----
  tmp <- tempfile(fileext = "png")
  png(tmp,
      width = ratio[1],
      height = ratio[2],
      units = "px",
      res = res,
      type="cairo",
      bg = "transparent")
  par(mar = c(0,0,0,0))
  plot.new()
  plot.window(coords[1:2],
              coords[3:4],
              xaxs = "i",
              yaxs = "i")
  points(x= x, y = y, type = type, ...)
  lines(x= x, y = y, type = type, ...)
  dev.off()
  
  # Plot png ----
  panel <- png::readPNG(tmp)
  rasterImage(panel,
              coords[1],
              coords[3],
              coords[2],
              coords[4])
}