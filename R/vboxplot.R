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

vboxplot <- function(x,
                     formula,
                     outline= F,
                     xlab= "",
                     ylab= "",
                     xlim= c(1, length(unique(obj$variable))),
                     ylim= if(outline) range(unlist(obj$value)) else range(c(obj$Q1, obj$Q5)))
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
  if(is.matrix(x) | is.data.frame(x))
    dat <- data.table::melt.data.table(data.table::as.data.table(x), measure.vars = colnames(x))
  if(is.list(x))
    dat <- data.table::rbindlist(lapply(x, function(x) data.table::data.table(value= x)), idcol = "variable")

  # Compute violins
  obj <- dat[, .(value= list(value)), variable]
  obj[, c("Q1", "Q2", "Q3", "Q4", "Q5", "outliers"):=
        {
          .b <- boxplot(unlist(value), plot= F)
          append(as.list(.b$stats[,1]), as.list(.b$out))
        }, variable]
  obj <- obj[, c("x", "y"):=
               {
                 .c <- density(unlist(value))
                 keep <- .c$x>=Q1 & .c$x<=Q5
                 .x <- (.c$y/max(.c$y)*0.35)[keep]
                 .x <- c(.x[1], .x, .x[length(.x)])
                 .x <- c(.GRP-.x, .GRP+rev(.x))
                 .y <- .c$x[keep]
                 .y <- c(Q1, .y, Q5)
                 .y <- c(.y, rev(.y))
                 .(list(.x),
                   list(.y))
               }, .(variable, Q1, Q5)]

  # Plot
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       xlab= xlab,
       ylab= ylab)

  obj[, polygon(unlist(x), unlist(y)), variable]
  obj[, rect(.GRP-0.1, Q2[1], .GRP+0.1, Q4[1]), .(variable, Q1, Q5)]
  obj[, segments(.GRP-0.1, Q3[1], .GRP+0.1, Q3[1], lwd= 2), .(variable, Q3)]
  if(outline)
    obj[, points(jitter(.GRP),
                 unlist(outliers),
                 cex= 0.5,
                 col= adjustcolor("lightgrey", 0.6),
                 pch= 19), variable]

  invisible(obj)
}

