#' model equation
#'
#' Returns the equation from lm model. The code was found at https://stats.stackexchange.com/questions/63600/how-to-translate-the-results-from-lm-to-an-equation
#'
#' @param model Model
#' @param ... Arguments to 'format', such as 'digits' and 'trim'.
#' @return Equation
#' @export

vl_model_equation <- function(model, ...) {
  format_args <- list(...)
  
  model_coeff <- model$coefficients
  format_args$x <- abs(model$coefficients)
  model_coeff_sign <- sign(model_coeff)
  model_coeff_prefix <- dplyr::case_when(model_coeff_sign == -1 ~ " - ",
                                         model_coeff_sign == 1 ~ " + ",
                                         model_coeff_sign == 0 ~ " + ")
  model_eqn <- paste(strsplit(as.character(model$call$formula), "~")[[2]], # 'y'
                     "=",
                     paste(dplyr::if_else(model_coeff[1]<0, "- ", ""),
                           do.call(formatC, format_args)[1],
                           paste(model_coeff_prefix[-1],
                                 do.call(formatC, format_args)[-1],
                                 " * ",
                                 names(model_coeff[-1]),
                                 sep = "", collapse = ""),
                           sep = ""))
  return(model_eqn)
}

#' applies a gasussian blur
#'
#' radius= +/- 10 bins
#'
#' @param signal A numeric vector of the signal to be blurred
#' @param ... Arguments to 'format', such as 'digits' and 'trim'.
#' @return Equation
#' @export

vl_gaussian_blur <- function(signal) {
  # Kernel (Pascal triangle)
  kernel <- c(1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1)
  # Smooth signal
  var <- c(rep(0, 10), signal, rep(0, 10))
  smooth <- sapply(seq(signal), function(i) sum(var[i:(i+20)]*kernel, na.rm= T)/sum(kernel, na.rm= T))
  return(smooth)
}
