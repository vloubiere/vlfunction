#' violin plot
#'
#' This function plots a violin plot from either matrix,
#' data.table, data.frame or using formula.
#'
#' @param x Matrix, data.atble or data.frame containing data
#' @param formula Formula to apply to the data
#' @param outline plot outlier values? default= F
#' @param xlab x axis label
#' @param ylab y axis label
#' @param plot_labels Should box labels be plotted? Default= T
#' @param xlim x axis limits
#' @param ylim y axis limits
#' @export

# # Example matrix
# set.seed(123)
# nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
# nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
# mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
#             rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
#             rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
#                   matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
#                   matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3)))
# mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
# rownames(mat) = paste0("row", seq_len(nr))
# colnames(mat) = paste0("column", seq_len(nc))

vl_vioplot <- function(x,
                     formula,
                     outline= F,
                     xlab= "",
                     ylab= "",
                     plot_labels= T,
                     xlim= c(0.5, max(obj$at)+0.5),
                     ylim= if(outline) range(unlist(obj$value)) else range(c(obj$Q1, obj$Q5)),
                     col= "white", 
                     at= NULL)
{
  # Format data to a bindable DT list
  if(!missing(formula))
  {
    if(length(formula) != 3L)
      stop("'formula' incorrect")
    mf <- stats::model.frame.default(formula, x)
    response <- attr(attr(mf, "terms"), "response")
    x <- split(mf[[response]], mf[-response])
  }
  if(is.matrix(x) | is.data.frame(x) | is.numeric(x))
  {
    dat <- data.table::as.data.table(x)
    dat <- data.table::melt.data.table(dat, measure.vars = colnames(dat))
  }
  if(is.list(x))
    dat <- data.table::rbindlist(lapply(x, function(x) data.table::data.table(value= x)), idcol = "variable")
  nvar <- length(unique(dat$variable))
  if(length(col)==1)
    col <- rep(col, nvar) else if(length(col) != nvar)
      stop(paste0("col vector should be length 1 or match variables length (",  nvar, ")"))
  if(is.null(at))
    at <- seq(nvar)
  if(length(at) != nvar)
    stop(paste0("at vector should be length 1 or match variables length (",  nvar, ")"))

  # Compute boxplots
  obj <- dat[, .(value= list(na.omit(value))), variable]
  obj[, Cc:= col[.GRP], variable]
  obj[, at:= at[.GRP], variable]
  obj[, c("Q1", "Q2", "Q3", "Q4", "Q5", "outliers", "too_few_values", "few_values"):=
        {
          .c <- unlist(value)
          .b <- boxplot(.c, plot= F)
          list(.b$stats[1,1],
               .b$stats[2,1],
               .b$stats[3,1],
               .b$stats[4,1],
               .b$stats[5,1], 
               list(if(length(.c)>2) .b$out else NA),
               ifelse(length(.c)>2, F, T),
               list(if(length(.c)>2) NA else .c))
        }, variable]
  
  # Compute violins
  obj[, c("x", "y"):= 
        {
          if(!too_few_values)
          {
            .c <- density(unlist(value))
            keep <- .c$x>=Q1 & .c$x<=Q5
            .x <- (.c$y/max(.c$y)*0.35)[keep]
            .x <- c(.x[1], .x, .x[length(.x)])
            .x <- c(at-.x, at+rev(.x))
            .y <- .c$x[keep]
            .y <- c(Q1, .y, Q5)
            .y <- c(.y, rev(.y))
          }else
          {
            .x <- .y <- NA
          }
          .(list(.x),
            list(.y))
        }, .(variable, at, Q1, Q5, too_few_values)]

  #------------------####
  # PLOT
  #------------------####
  # plot
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       xlab= xlab,
       ylab= ylab, 
       xaxt= "n")
  axis(1, 
       at= unique(obj$at), 
       labels= if(plot_labels) unique(obj$variable) else rep(NA, length(obj$at)))
  
  # violins
  obj[, polygon(unlist(x), unlist(y), col= Cc[1]), .(Cc, variable)]
  
  # boxplots
  obj[!(too_few_values), rect(at-0.1, Q2[1], at+0.1, Q4[1], col= "white"), .(variable, at, Q1, Q5)]
  obj[!(too_few_values), segments(at-0.1, Q3[1], at+0.1, Q3[1], lwd= 2), .(variable, at, Q3)]
  
  # outliers
  if(outline)
    obj[, points(jitter(rep(at, length(unlist(outliers)))),
                 unlist(outliers),
                 cex= 0.5,
                 col= adjustcolor("lightgrey", 0.8),
                 pch= 19), .(variable, at)]
  
  # Too few values treated differently
  obj[, {
    .c <- na.omit(unlist(few_values[[1]]))
    if(length(.c)>0)
    {
      if(length(.c)>1) 
      {
        segments(at, # sd
                 mean(.c)-sd(.c),
                 at,
                 mean(.c)+sd(.c),
                 cex= 0.8,
                 pch= 19,
                 col= Cc[1])
        points(at, # mean
               mean(.c),
               cex= 0.8, 
               col= Cc[1])
      }
      points(jitter(rep(at, length(.c))),
             .c,
             cex= 0.5,
             col= adjustcolor("lightgrey", 0.8),
             pch= 19)
    }
  }, .(variable, at, Cc)]
  
  # Return object
  invisible(obj)
}

