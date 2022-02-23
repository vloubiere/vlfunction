#' model equation
#'
#' Returns the equation from lm model. The code was found at https://stats.stackexchange.com/questions/63600/how-to-translate-the-results-from-lm-to-an-equation
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

#' Compute RMSE and R2 from observed and predicted values
#'
#' Comes from https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
#'
#' @param observed vector of observed values
#' @param predicted vector of predicted values
#'
#' @return RMSE and R2
#' @export 
vl_model_eval <- function(observed, predicted) 
{
  if(length(predicted)!=length(observed))
    print("Observed and predicted should have the same length")
  SSE <- sum((predicted - observed)^2)
  SST <- sum((observed - mean(observed))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/length(observed))
  
  # Model performance metrics
  data.table(RMSE = RMSE,
             Rsquare = R_square)
}

#' Compute explained variance
#' 
#' From Bernardo, computes the % explained variance for each predictor
#'
#' @param model An object containing the results returned by a model fitting function (e.g., lm or glm).
#' @return % of explained variance
#' @export
vl_model_expVar <- function(model)
{
  af <- stats::anova(model)
  af$PctExp <- af$"Sum Sq"/sum(af$"Sum Sq")*100
  return(as.data.table(af, keep.rownames = "vars")[, .(vars, PctExp)])
}
