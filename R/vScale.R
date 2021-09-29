#' Scaling
#'
#' Scale the values of a vector between 0 and 1
#'
#' @param x vector
#' @return Scaled vector.
#' @export

vl_scale01 <- function(x){(x-min(x))/(max(x)-min(x))}