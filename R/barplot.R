#' Title
#'
#' @param height A vector of values describing the bars which make up the plot. 
#' @param sd sd values to be plotted as whiskers around bars
#' @param xlim 
#' @param ylim 
#' @param individual.var List of individual variables to be plotted as points.
#' @param ind.pch pch value for individual variables
#' @param ind.col color value for individual variables
#' @param ind.jitter jitter value for individual variables
#' @param compute.diff List of pairwise bar pairs to be compared (V2/V11 ratio)
#' @param diff.cex cex value for pairwise comparisons
#' @param xpd xpd value for pairwaise comparisons
#' @param horiz Not supported atm
#' @param ... Extra arguments to be passed to barplot
#'
#' @return A simple barplot
#' @export
#'
#' @examples
#' vl_barplot(1:3, rep(.2, 3))
vl_barplot <- function(height,
                       sd= NULL,
                       xlim= NULL,
                       ylim= NULL,
                       individual.var= NULL,
                       ind.pch= 16,
                       ind.col= adjustcolor("lightgrey", .7),
                       ind.jitter= .2,
                       compute.diff= NULL,
                       diff.cex= .8,
                       xpd= NA,
                       horiz= F,
                       ...)
{
  if(!is.vector(height))
    stop("Only height vectors are supported for now")
  if(horiz)
    stop("horiz not supported yet ;)")
  
  # Create data table to keep track of all values ----
  dat <- data.table(height= height,
                    sd.min= height-sd,
                    sd.max= height+sd,
                    individual.var= individual.var)
  dat[, min:= apply(.SD, 1, function(x) min(unlist(x), na.rm= T))]
  dat[, max:= apply(.SD, 1, function(x) max(unlist(x), na.rm= T))]
  
  # Compute ylim that takes into account outliers ----
  if(is.null(ylim))
  {
    ylim <- range(c(dat$min, dat$max), na.rm = T)
    if(ylim[1]>0)
      ylim[1] <- 0
    if(ylim[2]<0)
      ylim[2] <- 0
  }
  
  # Initiate barplot ----
  bar <- barplot(height, xlim= xlim, ylim= ylim, ...)
  dat[, x:= bar]
  
  # Add sd arrows if specified ----
  if(!is.null(sd))
  {
    if(length(sd)!=length(height))
      stop("sd should be a vector or a of the same length as height")
    arrows(bar,
           height,
           bar,
           height+sd,
           angle = 90,
           length = .05,
           xpd= xpd)
    arrows(bar,
           height,
           bar,
           height-sd,
           xpd= T,
           angle = 90,
           length = .05,
           xpd= xpd)
  }
  
  # Add individual measurements if specified ----
  if(!is.null(individual.var))
  {
    individual.var <- as.list(individual.var)
    if(length(individual.var)!=length(height))
      stop("individual.var should be a vector or a list of the same length as height")
    x <- rep(bar, lengths(individual.var))
    y <- unlist(individual.var)
    points(jitter(x, amount = ind.jitter),
           y,
           pch= ind.pch,
           col= ind.col,
           xpd= xpd)
  }
  
  # Add mean diff if specified ----
  if(!is.null(compute.diff))
  {
    # Pairs to compare
    comp <- do.call(rbind, compute.diff)
    comp <- data.table(V1= comp[,1],
                       V2= comp[,2],
                       var1= dat$height[comp[,1]],
                       var2= dat$height[comp[,2]],
                       x0= bar[comp[,1]],
                       x1= bar[comp[,2]])
    comp[, x:= rowMeans(.SD), .SDcols= c("x0", "x1")]
    # Compute FC and max/ values (y position)
    comp[, var:= paste0("x", formatC(var2/var1, digits = 2))]
    comp[, max:= max(dat$max[V1:V2]), .(V1, V2)]
    setorderv(comp, "max")
    comp[, y:= max(dat[.BY, max, on= c("x>=x0", "x<=x1")]), .(x0, x1)] 
    # Adjust based on overlaps
    data.table::setorderv(comp, c("y", "x0", "x1"))
    comp[, idx:= .I] # Index to check previous bars
    adj <- strheight("M")*1.2
    for(i in seq(nrow(comp)))
    {
      .c <- sort(comp[comp[i], y, on= c("x0<=x1", "x1>=x0", "idx<=idx")]) # y pos overlapping lines
      .c <- .c[diff(c(.c, Inf))>2*adj  & .c>=comp[i,y]] # y values with enough space
      comp[i, y:= data.table::first(.c)+adj] # Smaller y value + adj
    }
    comp[, {
      segments(x0, y, x1, y, xpd= xpd)
      text(x, 
           y, 
           var,
           cex= diff.cex,
           offset= 0.1,
           xpd= xpd,
           pos= 3)
    }]
  }
  invisible(bar)
}