#' @title Boxplot
#' @description Just a wrapper around boxplot that makes it nicer and allows to add wilcox pvals
#' @param x list of variables to be plotted
#' @param compute_pval list of vectors of length two containing pairwise x indexes to be compared
#' @param pval_offset offset for pval plotting. Defaults to 0.04 (fraction of ylim)
#' @param tilt.names Should names be tilted (ignored if horizontal= TRUE)
#' @param ... Extra parameters for boxplot()
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
vl_boxplot.default <-
  function(x, ..., 
           compute_pval= NULL, pval_offset= 0.08, tilt.names= F,
           range = 1.5, width = NULL, varwidth = FALSE,
           notch = FALSE, outline = FALSE, names, plot = TRUE,
           border = par("fg"), col = NULL, log = "",
           pars = list(boxwex = 0.4, staplewex = NA, outwex = NA),
           horizontal = FALSE, add = FALSE, at = NULL,
           frame= F, whisklty = 1)
  {
    # Call
    pl <- substitute(boxplot(x, ..., range = range, width = width, varwidth = varwidth,
                             notch = notch, outline = outline, names= names, plot = FALSE,
                             border = border, col = col, log = log,
                             pars = pars,
                             horizontal = horizontal, add = add, at = at,
                             frame= F, whisklty = 1))
    
    # Boxplot stats
    box <- eval(pl, parent.frame())
    
    # Compute pval
    if(!is.null(compute_pval))
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
      pval <- vl_compute_bxp_pval(groups= groups,
                                  box= box,
                                  compute_pval= compute_pval, 
                                  pval_offset= pval_offset,
                                  outline= outline,
                                  at= at,
                                  horizontal= horizontal)
      pl$ylim <- pval$ylim
    }
    
    # Plot boxplot
    pl$plot <- plot
    if(tilt.names && !horizontal)
      pl$names <- NA
    eval(pl, parent.frame())
    # Plot pval
    vl_plot_bxp_pval(pval = pval$pval, 
                     horizontal = pval$horizontal, 
                     pos = pval$pos,
                     offset = pval$offset)
    # Plot tilted names
    if(tilt.names && !horizontal)
      text(if(is.null(at)) seq(box$names) else at,
           rep(par("usr")[3], length(box$names))-par("lheight"),
           box$names,
           srt= 45,
           offset= -0.35,
           pos= 2,
           xpd= T,
           cex= par("cex")*par("cex.axis"))
    
    # Return
    invisible(c(box, pval))
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

#' @export
vl_compute_bxp_pval <- function(groups, box, compute_pval, pval_offset, outline, at, horizontal)
{
  if(!is.list(compute_pval) | !all(lengths(compute_pval)==2))
    stop("compute_pval list of vectors of length two containing pairwise x indexes to be compared")
  # Make pairs object
  dat <- data.table::data.table(dat= groups,
                                x= if(is.null(at)) seq(groups) else at) # Retrieve at positionsbrowser()
  if(outline)
  {
    dat[, max:= sapply(dat, max, na.rm= T)]
    ymin <- min(unlist(dat))
  }else{
    dat[, max:= box$stats[5,]]
    ymin <- min(box$stats[1,])
  }
  dat[, max:= as.numeric(max)] # Needed to not lose y precision later
  # Order comparison pairs! (x0<=x1)
  compute_pval <- lapply(compute_pval, function(i) i[order(dat[i, x])])
  pval <- cbind(dat[sapply(compute_pval, `[`, 1), !"max"], 
                dat[sapply(compute_pval, `[`, 2), !"max"])
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
  return(list(pval= pval,
              ylim= c(ymin, max(c(dat$max, pval$y0+adj))),
              horizontal= horizontal,
              pos= pos,
              offset= offset))
}

#' @export
vl_plot_bxp_pval <- function(pval, 
                             horizontal,
                             pos,
                             offset)
{
  # Rotate if horizontal
  if(horizontal)
  {
    setnames(pval, 
             c("x", "y", "x0", "x1", "y0", "y1"), 
             c("y", "x", "y0", "y1", "x0", "x1"))
  }
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