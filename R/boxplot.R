#' @title Boxplot
#' @description Just a wrapper around boxplot that makes it nicer and allows to add wilcox pvals
#' @param x list of variables to be plotted
#' @param compute_pval list of vectors of length two containing pairwise x indexes to be compared
#' @param pval_cex cex value for pval plotting. Default= .8.
#' @param names Names to plot under boxplot. If function specified, applied to names before plotting
#' @param tilt.names Should names be tilted (ignored if horizontal= TRUE)
#' @param srt rotation angle for titled names
#' @param violin Should violins be added?
#' @param viocol Violin colors. Default to transparent
#' @param viowex Expansion factor for violins
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
           compute_pval= NULL, pval_cex= .8, tilt.names= F, srt= 45,
           range = 1.5, width = NULL, varwidth = FALSE,
           notch = FALSE, outline = FALSE, names, plot = TRUE,
           border = par("fg"), col = NULL, log = "",
           pars = list(boxwex = ifelse(violin, .2, .4), staplewex = NA, outwex = NA),
           horizontal = FALSE, add = FALSE, at = NULL,
           frame= F, whisklty = ifelse(violin, 2, 1), ylim= NULL, xaxt= "s",
           violin= FALSE, viocol = "#FFFFFF00", viowex= 0.4)
  {
    # Boxplot stats
    if(!missing(names) && is.function(names))
    {
      box <- boxplot(x, ..., plot = F)
      box$names <- names(box$names)
    }else
      box <- boxplot(x, ..., names= names, plot = F)
    
    # Compute groups
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
    attr(groups, "names") <- box$names
    
    # Plot boxplot
    if(plot)
    {
      boxplot(x, ..., range = range, width = width, varwidth = varwidth,
              notch = notch, outline = outline,
              names= if(tilt.names && !horizontal) NA else box$names,
              plot = plot, border = if(violin) NA else border, col = if(violin) NA else col, log = log,
              pars = pars, horizontal = horizontal, add = add, at = at,
              frame= frame, whisklty = if(violin) 0 else whisklty, ylim= ylim, xaxt= xaxt)
      if(violin)
      {
        # Violins
        xpos <- if(is.null(at)) seq(groups) else at
        viocols <- rep(viocol, length.out= length(groups))
        lapply(seq(groups), function(i) 
        {
          .d <- density(groups[[i]],
                        from= box$stats[1,i],
                        to= box$stats[5,i], 
                        na.rm= T)
          x <- .d$y
          x <- x/max(x)*viowex/2
          x <- xpos[i]-c(x, rev(-x))
          y <- .d$x
          y <- c(y, rev(y))
          if(horizontal)
            polygon(y, x, col= viocols[i]) else
              polygon(x, y, col= viocols[i])
        })
        # Add boxes
        boxplot(x, ..., range = range, width = width, varwidth = varwidth,
                notch = notch, outline = outline,
                names= if(tilt.names && !horizontal) NA else box$names,
                plot = plot, border = border, col = col, log = log,
                pars = pars, horizontal = horizontal, add = T, at = at,
                frame= frame, whisklty = whisklty, ylim= ylim, xaxt= "n", yaxt= "n")
      }
      # Plot pval
      if(!is.null(compute_pval))
      {
        pval <- vl_compute_bxp_pval(groups= groups,
                                    box= box,
                                    compute_pval= compute_pval,
                                    outline= outline,
                                    at= at,
                                    horizontal= horizontal)
        if(nrow(pval))
          vl_plot_bxp_pval(pval = pval,
                           horizontal = horizontal,
                           pval_cex= pval_cex)
      }
      # Plot tilted names
      if(tilt.names && !horizontal && xaxt!="n")
        vl_tilt_xaxis(if(is.null(at)) seq(box$names) else at,
                      grconvertY(grconvertY(0, "npc", "inch")-grconvertY(par("mgp")[2], "line", "inch"), "inch", "user"),
                      box$names,
                      srt= srt,
                      offset= 0.25,
                      pos= 2,
                      xpd= NA,
                      cex= par("cex.axis"))
    }
    
    # Return
    if(is.null(compute_pval))
      invisible(box) else
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
vl_compute_bxp_pval <- function(groups, box, compute_pval, outline, at, horizontal)
{
  if(!is.list(compute_pval) | !all(lengths(compute_pval)==2))
    stop("compute_pval list of vectors of length two containing pairwise x indexes to be compared")
  if(any(unlist(compute_pval)>length(groups)))
    stop("Some indexes provided in compute_pval are bigger than the number of groups in current boxplot")
  # Make pairs object
  dat <- data.table::data.table(dat= groups,
                                x= if(is.null(at)) seq(groups) else at) # Retrieve at positions
  dat[, dat:= lapply(dat, na.omit)]
  if(outline)
    dat[, max:= sapply(dat, max, na.rm= T)] else
      dat[, max:= box$stats[5,]]
  # Order comparison pairs! (x0<=x1)
  compute_pval <- lapply(compute_pval, function(i) i[order(dat[i, x])])
  pval <- cbind(dat[sapply(compute_pval, `[`, 1), !"max"], 
                dat[sapply(compute_pval, `[`, 2), !"max"])
  data.table::setnames(pval, c("dat0", "x0", "dat1", "x1"))
  # Compute wilcox pval
  pval[, wilcox:= mapply(function(x, y) wilcox.test(unlist(x), unlist(y))$p.value, x= dat0, y= dat1)]
  # Compute x pos text
  pval[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
  # Compute y pos (in inches in case axis is logged)
  adj <- if(horizontal)
    strwidth("M", units = "inch") else
      strheight("M", units = "inch")
  # Overlapping boxes
  pval[, y:= max(dat[.BY, max, on= c("x>=x0", "x<=x1")]), .(x0, x1)] 
  if(horizontal)
    pval[, y:= grconvertX(y, "user", "inch"), .(x0, x1)] else
      pval[, y:= grconvertY(y, "user", "inch"), .(x0, x1)]
  # Adjust based on overlaps
  data.table::setorderv(pval, c("y", "x0", "x1"))
  pval[, idx:= .I] # Index to check previous bars
  for(i in seq(nrow(pval)))
  {
    .c <- sort(pval[pval[i], y, on= c("x0<=x1", "x1>=x0", "idx<=idx")]) # y pos overlapping lines
    .c <- .c[diff(c(.c, Inf))>2*adj  & .c>=pval[i,y]] # y values with enough space
    pval[i, y:= data.table::first(.c)+adj] # Smaller y value + adj
  }
  # Return
  return(pval)
}

#' @export
vl_plot_bxp_pval <- function(pval, 
                             horizontal,
                             pval_cex)
{
  # Convert to users coordinates
  if(horizontal)
  {
    pval[, y0:= grconvertX(y, "inch", "user")]# Convenient if horizontal= T
    pval[, y1:= y0] # Convenient if horizontal= T
    pval[, y:= grconvertX(y+strheight("M", "inch")*0.45, "inch", "user")]
    setnames(pval, # Rotate if horizontal
             c("x", "y", "x0", "x1", "y0", "y1"), 
             c("y", "x", "y0", "y1", "x0", "x1"))
  }else
  {
    pval[, y0:= grconvertY(y, "inch", "user")]# Convenient if horizontal= T
    pval[, y1:= y0] # Convenient if horizontal= T
    pval[, y:= grconvertY(y+strwidth("M", "inch")*0.45, "inch", "user")]
  }
  pval[, {
    segments(x0, y0, x1, y1, xpd= NA)
    .c <- cut(wilcox,
              c(-Inf, 0.00001, 0.001, 0.01, 0.05, Inf),
              c("****", "***", "**", "*", "N.S"),
              include.lowest= T)
    text(x, 
         y, 
         .c,
         cex= pval_cex*ifelse(.c=="N.S", 0.5, 1),
         srt= ifelse(horizontal, -90, 0),
         offset= 0,
         xpd= NA)
  }]
}