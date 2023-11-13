#' applies a gasussian blur
#'
#' radius= +/- 10 bins
#'
#' @param x A numeric vector of the signal to be blurred
#' @return Smoothed signal
#' @export
vl_gaussian_blur <- function(x) {
  # Kernel (Pascal triangle)
  kernel <- c(1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1)
  # Smooth signal
  var <- c(rep(0, 10), x, rep(0, 10))
  smooth <- sapply(seq(x), function(i) sum(var[i:(i+20)]*kernel, na.rm= T)/sum(kernel, na.rm= T))
  return(smooth)
}

