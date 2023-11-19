#' LASSO regression using TF motifs
#'
#' @param response Variable to predict
#' @param counts Matrix of motif counts
#'
#' @return Return a LASSO model
#' 
#' @examples
#' test <- vl_trainLASSO(ATACSeq_FC, counts)
#' 
#' @export
vl_trainLASSO <- function(response,
                          counts)
{
  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -.1)
  lasso_reg <- cv.glmnet(counts,
                         response,
                         alpha = 1,
                         lambda = lambdas,
                         standardize = TRUE,
                         nfolds = 5)
  # Best  lambda
  lambda_best <- lasso_reg$lambda.min
  # Modelling
  model <- glmnet(counts,
                  response,
                  alpha = 1,
                  lambda = lambda_best,
                  standardize = TRUE)
  # Add Rsq
  predict_test <- predict(model, 
                          s = lambda_best, 
                          newx = counts)
  model$rsq <- vl_model_eval(observed = response,
                             predicted = predict_test)$Rsquare
  return(model)
}