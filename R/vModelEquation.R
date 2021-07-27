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