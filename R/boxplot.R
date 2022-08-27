#' @title Boxplot
#' @description Just a wrapper around boxplot that makes it nicer and allows to add wilcox pvals
#' @param x list of variables to be plotted
#' @param compute_pval list of vectors of length two containing pairwise x indexes to be compared
#' @param pval_offset offset for pval plotting. Defaults to 0.04 (fraction of ylim)
#' @examples
#' set.seed(1234)
#' vl_boxplot(formula= len~supp*dose, data=ToothGrowth, notch=TRUE,
#' col= c("gold","darkgreen"),
#' main="Tooth Growth", xlab="Suppliment and Dose",
#' compute_pval= list(c(1,2), c(5,6), c(1,6)))
#' @export
vl_boxplot <- function(x, ...) UseMethod("vl_boxplot")

#' @describeIn vl_boxplot default method
#' @export
vl_boxplot.default <- function(x, ..., 
                               compute_pval= NULL, pval_offset= 0.08,
                               range = 1.5, width = NULL, varwidth = FALSE,
                               notch = FALSE, outline = FALSE, names, plot = TRUE,
                               border = par("fg"), col = NULL, log = "",
                               pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
                               horizontal = FALSE, add = FALSE, at = NULL,
                               frame= FALSE, whisklty = 1)
  {
    args <- list(x, ...)
    namedargs <-
      if(!is.null(attributes(args)$names)) attributes(args)$names != ""
    else rep_len(FALSE, length(args))
    ## pars <- c(args[namedargs], pars)
    groups <- if(is.list(x)) x else args[!namedargs]
    if(0L == (n <- length(groups)))
      stop("invalid first argument")
    if(length(class(groups)))
      groups <- unclass(groups)
    if(!missing(names))
      attr(groups, "names") <- names
    else {
      if(is.null(attr(groups, "names")))
        attr(groups, "names") <- 1L:n
      names <- attr(groups, "names")
    }
    cls <- sapply(groups, function(x) class(x)[1L])
    cl <- if(all(cls == cls[1L])) cls[1L] else NULL
    #<<<<<<<<<<<<<<<<<<<< vl chunck 1 copy groups for pval
    groups.copy <- groups
    #<<<<<<<<<<<<<<<<<<<< vl chunck 1
    for(i in 1L:n)
      groups[i] <- list(boxplot.stats(unclass(groups[[i]]), range)) # do.conf=notch)
    stats <- matrix(0, nrow = 5L, ncol = n)
    conf  <- matrix(0, nrow = 2L, ncol = n)
    ng <- out <- group <- numeric(0L)
    ct <- 1
    for(i in groups) {
      stats[,ct] <- i$stats
      conf [,ct] <- i$conf
      ng <- c(ng, i$n)
      if((lo <- length(i$out))) {
        out	  <- c(out,i$out)
        group <- c(group, rep.int(ct, lo))
      }
      ct <- ct+1
    }
    if(length(cl) && cl != "numeric") oldClass(stats) <- cl
    z <- list(stats = stats, n = ng, conf = conf, out = out, group = group,
              names = names)
    #<<<<<<<<<<<<<<<<<<<< vl chunck 2 pval object
    if(!is.null(compute_pval)) 
    {
      if(!is.list(compute_pval) | !all(lengths(compute_pval)==2))
        stop("compute_pval list of vectors of length two containing pairwise x indexes to be compared")
      # Make pairs object
      dat <- data.table::data.table(dat= groups.copy,
                                    x= if(is.null(at)) seq(groups.copy) else at) # Retrieve at positions
      if(outline)
      {
        dat[, max:= sapply(dat, max, na.rm= T)]
        ymin <- min(unlist(dat))
      }else{
        dat[, max:= z$stats[5,]]
        ymin <- min(z$stats[1,])
      }
      dat[, max:= as.numeric(max)] # Needed to not lose y precision later
      pval <- cbind(dat[sapply(compute_pval, min), !"max"], # Order comparison pairs! (x0<=x1)
                    dat[sapply(compute_pval, max), !"max"])
      data.table::setnames(pval, c("dat0", "x0", "dat1", "x1"))
      # Compute wilcox pval
      pval[, wilcox:= mapply(function(x, y) wilcox.test(unlist(x), unlist(y))$p.value, x= dat0, y= dat1)]
      # Compute x pos text
      pval[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
      # Compute y pos
      adj <- pval_offset*(max(dat$max)-ymin)
      pval[, y:= max(dat[.BY, max, on= c("x>=x0", "x<=x1")]), .(x0, x1)] # Max depending on overlapping boxes
      data.table::setorderv(pval, c("y", "x0", "x1"))
      pval[, idx:= .I] # Index to check previous bars
      for(i in seq(nrow(pval)))
      {
        .c <- sort(pval[pval[i], y, on= c("x0<=x1", "x1>=x0", "idx<=idx")]) # y pos overlapping lines
        .c <- .c[diff(c(.c, Inf))>2*adj  & .c>=pval[i,y]] # y values with enough space
        pval[i, y:= data.table::first(.c)+adj] # Smaller y value + adj
      }
      # Horizontal plot
      pval[, y0:= y] # Convenient if horizontal= T
      pval[, y1:= y] # Convenient if horizontal= T
      pos <- ifelse(horizontal, 4, 3)
      offset <- ifelse(horizontal, 0.3, 0)
      z <- c(z,
             list(ylim= c(ymin, max(c(dat$max, pval$y0+adj)))), #ylim
             list(pval= pval))
    }
    #>>>>>>>>>>>>>>>>>>>> vl chunck 2
    if(plot) {
      if(is.null(pars$boxfill) && is.null(args$boxfill)) pars$boxfill <- col
      #<<<<<<<<<<<<<<<<<<<< vl chunck 3 extract args and adjust ylim before ploting
      args <- c(list(z, notch = notch, width = width, varwidth = varwidth,
                     log = log, border = border, pars = pars,
                     outline = outline, horizontal = horizontal, add = add,
                     at = at), args[namedargs])
      if(!hasArg(ylim) && "ylim" %in% names(z))
        args <- c(args, list(ylim= z$ylim))
      do.call("bxp", args)
      #>>>>>>>>>>>>>>>>>>>> vl chunck 3
      #<<<<<<<<<<<<<<<<<<<< vl chunck 4 plot pval
      if("pval" %in% names(z))
      {
        # Return if horizontal
        if(horizontal)
          setnames(pval, 
                   c("x", "y", "x0", "x1", "y0", "y1"), 
                   c("y", "x", "y0", "y1", "x0", "x1"))
        pval[, {
          segments(x0, y0, x1, y1)
          vl_plot_pval_text(x,
                            y,
                            wilcox,
                            stars_only = T,
                            pos= pos,
                            offset= offset)
        }]
      }
      #>>>>>>>>>>>>>>>>>>>> vl chunck 4
      invisible(z)
    }
    else z
  }
#' @describeIn vl_boxplot method for matrices
#' @export
vl_boxplot.matrix <- function(x, use.cols = TRUE, ...)
{
  ## Purpose: Boxplot for each column or row [use.cols= TRUE / FALSE] of a matrix
  ## -------------------------------------------------------------------------
  ## Arguments: x: a numeric matrix; use.cols: logical, columns (T) or rows (F)
  ## <FIXME split.matrix>
  groups <- if(use.cols) {
    split(c(x), rep.int(1L:ncol(x), rep.int(nrow(x), ncol(x))))
  } else split(c(x), seq(nrow(x)))
  ## Make use of col/row names if present
  if (length(nam <- dimnames(x)[[1+use.cols]])) names(groups) <- nam
  invisible(vl_boxplot(groups, ...))
}

#' @describeIn vl_boxplot method for formula
#' @export
vl_boxplot.formula <- function(formula, data = NULL, ..., subset, na.action = NULL)
{
  if(missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$na.action <- na.action # force use of default for this method
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"), "response")
  vl_boxplot(split(mf[[response]], mf[-response]), ...)
}