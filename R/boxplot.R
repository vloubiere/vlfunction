#' violin plot
#'
#' This function plots a violin plot from either matrix,
#' data.table, data.frame or using formula.
#'
#' @param x Matrix, data.atble or data.frame containing data
#' @param x list data
#' @param compute_pval A list containing sublists of pairwise x indexes to compute pvalues from. example-> list(c(1,2), c(1,3))
#' @param outline show outliers? default=F
#' @param xlab 
#' @param ylab 
#' @param xlim 
#' @param ylim 
#' @param names group labels (x axis). By default, names= names(x). if set to false, no labels are ploted
#' @param boxcol box color
#' @param boxwex box width
#' @param violin If set to FALSE, the violin is not computed/shown
#' @param violcol violin color
#' @param violwex violin width
#' @param at x position for plotting
#' @param trim Should violin be trimmed?
#' @param add If TRUE, plot on the top of actual plot
#' @param ... Extra arguments set to plot function
#' @examples
#' # Create test matrix
#' set.seed(1234)
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' 
#' vl_boxplot(test, violin= T, outline= T, compute_pval= list(c(1,2),c(1,3),c(3,5)))
#' 
#' @export
vl_boxplot <- function(x, ...) UseMethod("vl_boxplot")

#' @describeIn vl_boxplot method for matrices
#' @export
vl_boxplot.matrix <- function(x, ...)
{
  x <- lapply(seq(ncol(x)), function(i) x[,i])
  vl_boxplot.default(x, ...)
}

#' @describeIn vl_boxplot method for formula
#' @export
vl_boxplot.formula <- function(formula, data = NULL, ...)
{
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"), "response")
  vl_boxplot.default(split(mf[[response]], mf[-response]), ...)
}

#' @describeIn vl_boxplot method for data.table
#' @export
vl_boxplot.data.table <- function(x, ...)
{
  x <- as.list(x)
  vl_boxplot.default(x, ...)
}

#' @describeIn vl_boxplot default method
#' @export
vl_boxplot.default <- function(x,
                               compute_pval,
                               outline= F,
                               xlab= NA,
                               ylab= NA,
                               xlim,
                               ylim,
                               names= NULL,
                               boxcol= "white",
                               boxwex= 0.15,
                               violin= F,
                               violcol= "white",
                               violwex= 0.35,
                               at,
                               trim= T,
                               add= F,
                               ...)
{
  # Format data list
  if(!is.list(x))
    x <- list(x)
  if(is.character(names))
    names(x) <- names
  if(is.null(names(x)))
    names(x) <- seq(x)
  # Format pval list
  if(!missing(compute_pval))
  {
    if(max(unlist(compute_pval))>length(x) | !is.numeric(unlist(compute_pval)))
      stop("All elements in compute_pval should be numeric and be smaller than length(x)")
    if(!is.list(compute_pval))
      compute_pval <- list(compute_pval)
    if(!all(lengths(compute_pval)==2))
      stop("compute_pval should be a list of vectors of length two containing pairwise x indexes to be compared")
  }
  # Check
  if(missing(at))
    at <- seq(x)
  
  # Format as data.table
  obj <- data.table(variable= names(x),
                    value= lapply(x, na.omit),
                    at= at,
                    boxCc= boxcol,
                    violCc= violcol)
  obj[, N_obs:= lengths(value)]
  
  #--------------------------#
  # Compute results
  #--------------------------#
  # Box
  obj[, c("Q1", "Q2", "Q3"):= lapply(c(0.25, 0.5, 0.75), function(q) sapply(value, function(x) quantile(x, q)))]
  obj[, c("min", "max"):= .(Q1-1.5*sapply(value, IQR),
                                   Q3+1.5*sapply(value, IQR))]
  # Outliers
  if(outline)
    obj[, outliers:= mapply(function(min, max, value) value[value>max | value<min], min, max, value, SIMPLIFY= F)]
  
  # Violins
  if(violin)
  {
    if(trim) # trim violin plot to data ranges
      obj[N_obs>2, dens:= lapply(value, function(x) density(x, from= min(x, na.rm= T), to= max(x, na.rm= T)))] else
        obj[N_obs>2, dens:= lapply(value, function(x) density(x))]
    obj[N_obs>2, x:= lapply(dens, function(x) x$y/max(x$y)*violwex)]
    obj[N_obs>2, y:= lapply(dens, function(x) x$x)]
  }
  
  # Compute min/max ploted values for each var
  cols <- intersect(c("min", "max", "outliers", "y"), names(obj))
  obj[N_obs>0, y_min:= min(unlist(.SD), na.rm= T), .SDcols= cols]
  obj[N_obs>0, y_max:= max(unlist(.SD), na.rm= T), .SDcols= cols]
  adj <- diff(c(min(obj$y_min, na.rm = T), # Adjust factor used for plotting (4% of range)
                max(obj$y_max, na.rm = T)))*0.04
  
  # Compute pvalues
  if(!missing(compute_pval))
  {
    pval <- as.data.table(do.call(rbind, compute_pval))
    # Remove combinations for which there are not enough values
    enough_obs <- pval[, .(check= all(lengths(obj[c(V1,V2), value])>0)), .(V1, V2)]$check
    if(any(!enough_obs))
      print("Some pval groups did not contain enough finite obs and were discarded")
    pval <- pval[(enough_obs)]
    if(nrow(pval)>0)
    {
      # Wilcoxon
      pval[, pval:= wilcox.test(obj[V1, unlist(value)],
                                obj[V2, unlist(value)])$p.value, .(V1, V2)]
      # x values correspond to obj "at" position
      pval[, c("x0", "x1"):= as.list(range(obj[c(V1, V2), at])), .(V1, V2)]
      # y position is the max of crossing boxes + adjustment
      pval$y <- obj[pval, max(y_max+adj*1.5, na.rm = T), .EACHI, on= c("at<=x1", "at>=x0")]$V1
      # Avoid overlapping segments by requiring some space between consecutive y values
      if(nrow(pval)>1)
      {
        pval <- pval[order(x0, x1)]
        pval[, overlap:= cumsum(x0-cummax(x1)[c(1, seq(.N-1))]>0)]# Make groups of overlapping segments
        pval <- pval[order(y, x1-x0, x0)]
        pval[, y:= {
          for(i in 2:(.N)) 
            if(y[i]<y[i-1]+adj*1.5)
              y[i] <- y[i-1]+adj*1.5
          y
        }, overlap]
      }
    }
  }
  #------------------####
  # PLOT
  #------------------####
  if(missing(xlim))
    xlim <- range(obj$at)+c(-0.5,0.5)
  if(missing(ylim))
  {
    ylim <- range(obj[, .(y_min-adj, y_max+adj)], na.rm= T) # adjust fact 4% range (see higher)
    if(!missing(compute_pval) && nrow(pval)>0 && max(pval[,y+adj])>ylim[2])
      ylim[2] <- max(pval[,y+adj])
  }
    
  # plot
  if(!add)
  {
    plot(NA,
         xlim = xlim,
         ylim = ylim,
         xlab= xlab,
         ylab= ylab, 
         xaxt= "n",
         ...)
    if(!isFALSE(names))
      axis(1,
           at= unique(obj$at), 
           labels= obj$variable)
  }

  # violins
  if(violin)
  {
    obj[N_obs>2, 
        mapply(function(x, y, at, violCc)
        {
          x <- c(at+x, at-rev(x))
          y <- c(y, rev(y))
          polygon(x, y, col = violCc)
        }, x, y, at, violCc)]
  }
  
  # boxplots
  obj[, segments(at, min, at, max)]
  obj[, rect(at-boxwex, Q1, at+boxwex, Q3, col= boxCc)]
  obj[, segments(at-boxwex, Q2, at+boxwex, Q2, lwd= 2*par("lwd"), lend= 2)]
  
  # outliers
  if(outline)
    obj[, mapply(function(at, min, max, outliers) {
      x <- jitter(rep(at, length(outliers)))
      points(x, 
             outliers, 
             pch= 16, 
             col= adjustcolor("lightgrey", 0.6))
    }, at, min, max, outliers)]
  
  # pval
  if(!missing(compute_pval) && nrow(pval)>0)
  {
    segments(pval$x0,
             pval$y,
             pval$x1,
             pval$y)
    vl_plot_pval_text(rowMeans(pval[, .(x0, x1)]),
                      pval$y,
                      pval$pval, 
                      stars_only = T)
  }
  
  # Returns object
  if(!missing(compute_pval) && nrow(pval)>0)
    obj <- list(obj= obj, pval= pval)
  invisible(obj)
}
