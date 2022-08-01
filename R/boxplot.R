#' @title Boxplot
#' @description Just a wrapper around boxplot that makes it nicer and allows to add violins and wilcox pvals
#' 
#' @param x list of variables to be plotted
#' @param compute_pval list of vectors of length two containing pairwise x indexes to be compared
#' @param outline Should outliers be plotted?
#' @param xlim xlim
#' @param ylim ylim
#' @param xlab x label
#' @param ylab y label
#' @param names box names (x axis)
#' @param boxcol Box colors (recycled)
#' @param boxwex Boxes expansion factor. default to 0.4, 0.25 if violin= T
#' @param box.lty Line type used for boxplots, default= 1
#' @param staplewex Staples expansion factor. default= NA
#' @param violin Should violins be plotted?
#' @param violcol Violin colors (recycled) 
#' @param violwex violins expansion factor. default to 0.4
#' @param wilcox.alternative When compute_pval is specified, alternative of the wilcox.test. default= "two.sided"
#' @param pval_offset offset for pval plotting. Defaults to 0.04 (fraction of ylim)
#' @param horizontal Whould the plot be made horizontal?
#' @param main main title
#' @param tilt.names Should boxes names be tilted (works only if horizontal= T)
#' @param axes Should the box around the plot be drawn? Default= F
#' @param xaxt Should the x axis be plotted or not ("n") ? default= "o"
#' @param yaxt Should the y axis be plotted or not ("n") ? default= "o"
#' @param las style of axis label. see ?par()
#' @param plot Should the boxlplot be ploted? default to TRUE
#' @param add Should the boxlplot be added to existing plot? default to FALSE
#' @param ... Extra parameters passed to boxplot, such as las, lwd... 
#' @examples
#' # Create test matrix
#' set.seed(1234)
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' vl_boxplot(test, compute_pval= list(c(1,2), c(1,4), c(1,5), c(5,6), c(1,10), c(9,10)), outline= F, violcol= c("red", "yellow"), boxcol= c("green", "purple"))
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
                               xlim, 
                               ylim,
                               main= NA,
                               xlab= NA,
                               ylab= NA,
                               names= NULL,
                               boxcol= "white",
                               boxwex= ifelse(violin, 0.25, 0.4),
                               box.lty= 1,
                               staplewex= NA,
                               violin= F,
                               violcol= "white",
                               violwex= 0.4,
                               wilcox.alternative= "two.sided",
                               pval_offset= 0.04,
                               horizontal= F,
                               tilt.names= F,
                               axes= F,
                               xaxt= "o",
                               yaxt= "o",
                               las= par("las"),
                               plot= T,
                               add= F,
                               ...)
{
  if(!missing(compute_pval) && !is.list(compute_pval))
    stop("compute_pval should be a list of vectors of length two containing pairwise x indexes to be compared")

  # Format data list
  if(!is.list(x))
    x <- list(x)
  if(is.character(names))
    names(x) <- names
  if(is.null(names(x)))
    names(x) <- seq(x)
  
  # Compute box stats and lims
  box <- sapply(x, 
                boxplot.stats,
                do.out= outline | violin)
  xrange <- c(0.5, length(x)+0.5)
  yrange <- range(box[c("out", "stats"),], na.rm= T)
  
  # Compute violins
  if(violin)
  {
    violcols <- c(matrix(violcol, ncol= 1, nrow= length(x)))
    viols <- lapply(seq(x), function(i) {
      var <- x[[i]]
      if(length(var)>2)
      {
        dens <- density(var,
                        from= min(var, na.rm= T),
                        to= max(var, na.rm= T), 
                        na.rm= T)
        xp <- dens$y/max(dens$y)*violwex
        xp <- c(i+xp, i-rev(xp))
        yp <- c(dens$x, rev(dens$x))
        list(x= xp, 
             y= yp,
             col= violcols[i])
      }else
        list()
    })
    viols <- rbindlist(viols[lengths(viols)>0],
                       idcol = T)
    viols <- viols[, .(x= .(x),
                       y= .(y)), .(.id, col)]
    # Adjust max
    if(min(unlist(viols$y))<yrange[1])
      yrange[1] <- min(unlist(viols$y))
    if(max(unlist(viols$y))>yrange[2])
      yrange[2] <- min(unlist(viols$y))
  }
  
  # Compute pvals
  if(!missing(compute_pval))
  {
    if(missing(ylim)) # plotting adjust
      adj <- diff(yrange)*pval_offset else
        adj <- diff(ylim)*pval_offset
    pval <- matrix(sapply(compute_pval, sort), ncol= 2, byrow = T)
    pval <- as.data.table(pval)
    setnames(pval, c("x0", "x1"))
    pval[, c("var1", "var2"):= .(x[x0], x[x1])]
    pval <- pval[lengths(var1)>0 & lengths(var2)>0]
    # Compute pvals
    if(nrow(pval)>0)
    {
      pval[, pval:= wilcox.test(unlist(var1),
                                unlist(var2),
                                alternative= wilcox.alternative)$p.value, .(x0, x1)]
      pval[, max:= max(unlist(box[c("out", "stats"), x0:x1]), na.rm= T), .(x0, x1)]
      setorderv(pval, c("x0", "x1", "max"))
      # Compute contig idx and adjust y
      pval[, idx:= cumsum(x0>data.table::shift(x1, fill= max(x1)))]
      setorderv(pval, c("idx", "max"))
      pval[, y:= max+cumsum(c(0, diff(max))<adj)*adj, idx]
      pval[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
      # Adjust max
      if(max(pval$y+adj)>yrange[2])
        yrange[2] <- max(pval$y+adj)
    }
  }

  #------------------------#
  # Plot
  #------------------------#
  if(missing(xlim))
    xlim <- xrange
  if(missing(ylim))
    ylim <- yrange
  if(plot)
  {
    if(!add)
    {
      plot.new()
      plot.window(xlim= if(horizontal) ylim else xlim,
                  ylim= if(horizontal) xlim else ylim)
    }
    title(main= main, 
          xlab= xlab, 
          ylab= ylab)
    if(violin && nrow(viols)>0)
      viols[, {
        polygon(if(horizontal) unlist(y) else unlist(x), 
                if(horizontal) unlist(x) else unlist(y), 
                col= col[1])
      }, .(.id, col)]
    boxplot(x,
            pch= NA,
            outline= F,
            add= T,
            staplewex= staplewex, 
            lty= box.lty, 
            boxwex= boxwex,
            col= boxcol,
            names= NA,
            horizontal= horizontal,
            axes= axes,
            xaxt= "n",
            yaxt= "n",
            ...)
    if(xaxt!="n")
      if(horizontal)
        axis(1, 
             lwd= ifelse(axes, 0, 1),
             lwd.ticks= 1,
             las= las) else
               axis(1,
                    at= seq(x),
                    labels= if(tilt.names) rep(NA, length(x)) else names(x),
                    lwd= ifelse(axes, 0, 1),
                    lwd.ticks= 1,
                    las= las)
    if(yaxt!="n")
      if(horizontal)
        axis(2,
             at= seq(x),
             labels= names(x),
             lwd= ifelse(axes, 0, 1),
             lwd.ticks= 1,
             las= las) else
               axis(2, 
                    lwd= ifelse(axes, 0, 1),
                    lwd.ticks= 1,
                    las= las)
    if(tilt.names && xaxt!="n" && !horizontal)
      text(seq(x),
           par("usr")[3]-(grconvertY(par("mgp")[2], "line", "user")-grconvertY(0, "line", "user")),
           names(x),
           srt= 45,
           offset= -0.35,
           pos= 2,
           xpd= T,
           cex= par("cex.axis"))
    if(outline)
    {
      points(if(horizontal) unlist(box["out",]) else jitter(rep(seq(x), lengths(box["out",]))),
             if(horizontal) jitter(rep(seq(x), lengths(box["out",]))) else unlist(box["out",]),
             pch= 16,
             col= adjustcolor("grey", 0.5))
    }
    
    if(!missing(compute_pval) && nrow(pval)>0)
      pval[, {
        segments(if(horizontal) y else x0, 
                 if(horizontal) x0 else y, 
                 if(horizontal) y else x1,
                 if(horizontal) x1 else y)
        vl_plot_pval_text(if(horizontal) y else x, 
                          if(horizontal) x else y, 
                          pval, 
                          stars_only = T,
                          pos= ifelse(horizontal, 4, 3),
                          srt= ifelse(horizontal, 270, 0))
      }]
  }
  
  obj <- list(stats= as.data.table(box["stats",]),
              xlim= if(horizontal) ylim else xlim,
              ylim= if(horizontal) xlim else ylim)
  if(violin && nrow(viols)>0)
    obj <- c(obj, list(violins= viols))
  if(!missing(compute_pval) && nrow(pval)>0)
    obj <- c(obj, list(pval= pval))
  invisible(obj)
}
