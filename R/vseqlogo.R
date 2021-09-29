#' plot seqlogo rasterImage
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm pwm matrix
#' @param xleft left plot limit. Default= 0
#' @param ybottom bottom plot limit. Default= 0
#' @param xright right plot limit. Default= 1
#' @param ytop top plot limit. Default= 1
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' @export


# pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)

vl_seqlogo <- function(pwm, 
                       xleft= 0, 
                       ybottom= 0, 
                       xright= 1, 
                       ytop= 1, 
                       add= T)
{
  if(!is.matrix(pwm))
  {
    stop("!is.matrix(pwm)")
  }
  require(fields)
  require(seqLogo)
  require(png)
  require(colorspace)
  
  tmp <- base::tempfile(fileext = ".png") 
  grDevices::png(tmp, type="cairo", width = 1000, height = 1000, units = "px")
  seqLogo::seqLogo(pwm, xaxis = F, yaxis = F)
  dev.off()
  im <- png::readPNG(tmp)
  res <- matrix(NA, nrow = nrow(im[,,1]), ncol = ncol(im[,,1]))
  res[im[,,1]==1 & im[,,2]==0 & im[,,3]==0] <- "firebrick1"
  res[im[,,1]==1 & im[,,2]>0.1 & im[,,2]<0.9 & im[,,3]==0] <- "goldenrod1"
  res[im[,,1]==0 & im[,,2]==1 & im[,,3]==0] <- "forestgreen"
  res[im[,,1]==0 & im[,,2]==0 & im[,,3]==1] <- "dodgerblue2"
  if(!add) plot.new()
  rasterImage(res[50:950,50:950], xleft= xleft, ybottom= ybottom, xright= xright, ytop= ytop, xpd= T)
}