#' applies a gasussian blur
#'
#' radius= +/- 10 bins
#'
#' @param signal A numeric vector of the signal to be blurred
#' @return Smoothed signal
#' @export
vl_gaussian_blur <- function(signal) {
  # Kernel (Pascal triangle)
  kernel <- c(1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1)
  # Smooth signal
  var <- c(rep(0, 10), signal, rep(0, 10))
  smooth <- sapply(seq(signal), function(i) sum(var[i:(i+20)]*kernel, na.rm= T)/sum(kernel, na.rm= T))
  return(smooth)
}

#' Title
#'
#' @param x number of colours ro return
#' @return vector of colors of length x
#' @export
vl_palette_categ1 <- function(x){
  colorRampPalette(c("orchid1", "darkorchid1", "purple", "darkorchid4", 
                     "olivedrab1", "limegreen", "olivedrab3", "olivedrab4", 
                     "lightsteelblue1", "cornflowerblue", "blue", "navy", 
                     "pink", "pink3", "indianred2", "red2", 
                     "black", "gold", "goldenrod", "goldenrod4", "sienna4"))(x)
}

#' Title
#'
#' @param x number of colours ro return
#' @return vector of colors of length x
#' @export
vl_palette_categ2 <- function(x){
  colorRampPalette(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                     "#0072B2", "#D55E00", "#CC79A7", "#999999"))(x)
}
