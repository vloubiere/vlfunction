#' ROC AUC
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param plot Should the ROC be plotted? Default= FALSE.
#' @param xlab xlab. Default= "FALSE positive rate".
#' @param ylab ylab. Default= "TRUE positive rate".
#' @param type The type of plot, when add= FALSE. Default= "l".
#' @param add When plot is set to TRUE, should only the line be added to an existing plot?
#' @param ... Extra arguments to be passed to plot (when add= FALSE) or lines (when add= TRUE).
#'
#' @return
#' @export
#'
#' @examples
vl_ROC_AUC <- function(predicted,
                       label,
                       plot= FALSE,
                       xlab= "False Positive Rate",
                       ylab= "True Positive Rate",
                       type= "l",
                       add= FALSE, ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))
  
  # Make data table ----
  dat <- data.table(label= label,
                    predicted= predicted)
  
  # Order ----
  setorderv(dat, "predicted", -1)
  
  # Compute AUC ----
  dat[, FPR:= cumsum(!label)/sum(!label)]
  dat[, TPR:= cumsum(label)/sum(label)]
  dat[, ROC_AUC:= c(0, sapply(seq(.N)[-1], function(i) (FPR[i] - FPR[i-1]) * (TPR[i] + TPR[i-1]) / 2))]
  
  # Plot ----
  if(plot)
  {
    if(add)
    {
      lines(dat$FPR,
            dat$TPR,
            ...)
    }else
    {
      lines(dat$FPR,
            dat$TPR,
            xlab= xlab,
            ylab= ylab,
            type= type,
            ...)
    }
  }
  
  # AUC function
  # AUC <- pROC::roc(AUC$label, AUC$score)$auc
  
  return(round(sum(dat$ROC_AUC), 2))
}

#' PR AUC
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param plot Should the ROC be plotted? Default= FALSE.
#' @param xlab xlab. Default= "FALSE positive rate".
#' @param ylab ylab. Default= "TRUE positive rate".
#' @param type The type of plot, when add= FALSE. Default= "l".
#' @param add When plot is set to TRUE, should only the line be added to an existing plot?
#' @param ... Extra arguments to be passed to plot (when add= FALSE) or lines (when add= TRUE).
#'
#' @return
#' @export
#'
#' @examples
vl_PR_AUC <- function(predicted,
                      label,
                      plot= FALSE,
                      xlab= "True Positive Rate",
                      ylab= "Positive Predictive Value",
                      type= "l",
                      add= FALSE, ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))
  
  # Make data table ----
  dat <- data.table(label= label,
                    pred= predicted)
  
  # Order ----
  setorderv(dat, "pred", -1)
  
  # Compute AUC ----
  dat[, TPR:= cumsum(label)/sum(label)]
  dat[, PPV:= cumsum(label)/seq(.N)]
  dat[, PR_AUC:= c(0, sapply(seq(.N)[-1], function(i) (TPR[i] - TPR[i-1]) * (PPV[i] + PPV[i-1]) / 2))]
  
  # Plot ----
  if(plot)
  {
    if(add)
    {
      lines(dat$TPR,
            dat$PPV,
            ...)
    }else
    {
      plot(dat$TPR,
           dat$PPV,
           xlab= xlab,
           ylab= ylab,
           type= type,
           ...)
    }
  }
  
  # AUC function
  # AUC <- pROC::roc(AUC$label, AUC$score)$auc
  
  return(round(sum(dat$PR_AUC), 2))
}


#' Positive Predicted Value curve
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#' @param Nleft Number of enhancers left before cutoff. Default= 100.
#' @param plot Should the PPV be plotted?
#' @param xlim x limits for plotting. Default= NULL.
#' @param ylim y limits for plotting. Default= NULL.
#' @param xlab x label. Default= "Prediction score".
#' @param ylab y label. Default= "Positive pred. value (%)".
#' @param lty.1 Line type before Nleft. Default= 1.
#' @param lty.2 Line type after Nleft. Default= 3.
#' @param col Color.
#' @param cex.max.PPV cex for the max PPV text label.
#' @param add When plot is set to TRUE, should only the lines/label be added to an existing plot?
#' @param ... Extra arguments to be passed to lines
#'
#' @return
#' @export
#'
#' @examples
vl_PPV <- function(predicted,
                   label,
                   Nleft= 100,
                   plot= FALSE,
                   xlim= NULL,
                   ylim= NULL,
                   xlab= "Prediction score",
                   ylab= "Positive pred. value (%)",
                   lty.1= 1,
                   lty.2= 3,
                   col= "black",
                   cex.max.PPV= .7,
                   add= FALSE,
                   ...)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))
  
  # Create a data table with observed and predicted values
  dat <- data.table(label = label,
                    pred = predicted)
  
  # Sort the data table by the predicted values in descending order
  setorderv(dat, cols = "pred", order = -1)
  
  # Calculate the cumulative percentage of actually positive values
  dat[, TP:= cumsum(label)] # TRUE positive
  dat[, FP:= cumsum(!label)] # FALSE positive
  dat[, PPV:= TP/(TP+FP)*100]
  setorderv(dat,
            cols = "pred")
  
  # Fit a smooth spline and find peaks
  x <- dat[1:(.N-Nleft), pred]
  y <- dat[1:(.N-Nleft), PPV]
  spline_fit <- smooth.spline(x, y)
  spline_derivative <- predict(spline_fit, deriv = 1)
  
  # Find peaks
  peaks <- which(diff(sign(spline_derivative$y)) == -2) + 1
  peak_x_values <- spline_derivative$x[peaks]
  peak_y_values <- predict(spline_fit, x = peak_x_values)$y
  peak_x_values <- c(peak_x_values, last(x)) # value at .N-Nleft
  peak_y_values <- c(peak_y_values, last(y)) # value at .N-Nleft 
  
  # Select ideal cutoff
  sel <- min(which(peak_y_values>max(0.95*peak_y_values)))
  x_cutoff <- peak_x_values[sel]
  y_cutoff <- peak_y_values[sel]
  # Plot PPV and cutoffs
  if(plot)
  {
    # Initiate plot
    if(!add)
    {
      
      plot(dat$pred,
           dat$PPV+diff(range(dat$PPV, na.rm= T))/15,
           type = "n",
           xlab = "Prediction score",
           ylab = "Positive pred. value (%)",
           xlim= xlim,
           ylim= ylim,
           ...)
    }
    
    # Plot PPV lines
    dat[,{
      # PPV
      lines(pred[1:(.N-Nleft)],
            PPV[1:(.N-Nleft)],
            lty= lty.1,
            col= col,
            ...)
      lines(pred[(.N-Nleft):.N],
            PPV[(.N-Nleft):.N],
            lty = lty.2,
            col= col,
            ...)
      
      # Plot the smooth spline
      # lines(spline_fit, col = "blue", lwd = 2)
      
      # Plot top PPV point
      points(x_cutoff,
             y_cutoff,
             col = adjustcolor(col, .7),
             pch = 19)
      text(x_cutoff,
           y_cutoff,
           paste0(round(y_cutoff, 1), "%"),
           pos= 3,
           offset= .5,
           col= col,
           cex= cex.max.PPV)
      
      # Add segments
      # segments(x_cutoff,
      #          0,
      #          x_cutoff,
      #          y_cutoff,
      #          lty= "33")
      # segments(0,
      #          y_cutoff,
      #          x_cutoff,
      #          y_cutoff,
      #          lty= "33")
      # text(x_cutoff,
      #      y_cutoff/2,
      #      round(x_cutoff, 2),
      #      pos= 4,
      #      offset= 1)
      # text(x_cutoff/2,
      #      y_cutoff,
      #      paste0(round(y_cutoff, 1), "%"),
      #      pos= 3,
      #      offset= 1)
    }]
  }
  
  # Return cutoffs
  return(list(min_PPV= dat$PPV[1],
              predict_cutoff= x_cutoff,
              PPV_at_cutoff= y_cutoff))
}

#' Compute Matthew's PCC
#'
#' @param predicted Predicted values from the model (ranging from 0 to 1).
#' @param label A vector of logical labels (or that can be coerced to logical).
#'
#' @return
#' @export
#'
#' @examples
vl_mPCC <- function(predicted,
                   label)
{
  if(!is.logical(label))
    label <- as.logical(label)
  if(sum(label)==0)
    warning(paste0(length(label), "/", length(label), " labels are set to FALSE"))
  
  # Make table
  conf_matrix <- table(pred= factor(predicted>0.5, c(FALSE, TRUE)),
                       obs= factor(as.logical(label), c(FALSE, TRUE)))
  TP <- as.numeric(conf_matrix[2, 2])
  TN <- as.numeric(conf_matrix[1, 1])
  FP <- as.numeric(conf_matrix[2, 1])
  FN <- as.numeric(conf_matrix[1, 2])
  # Compute
  mcPCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  # Return
  return(mcPCC)
}