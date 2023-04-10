#' Title
#'
#' @param dat A named list of variable of length 2 or 3
#' @param col Circles colors. Default to c("cornflowerblue", "tomato", "gold") and alpha= 0.3
#' @param border circle borders
#'
#' @return Plots a venn diagram for 2/3 groups
#' @export
#'
#' @examples
#' par(mfrow= c(1,2))
#' dat <- list(A= c("A", "B", "C", "D", "E", "F", "G", "H"),
#' B= c("I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "A"),
#' C= c("A", "B", "C", "K", "L", "M", "U", "V", "W", "X", "Y", "Z"))
#' vl_venn(dat)
#' dat <- list(A= c("A", "B", "C", "D", "E", "F", "G", "H"),
#' B= c("A", "B", "C", "D", "K", "L", "M", "U", "V", "W", "X", "Y", "Z"))
#' vl_venn(dat)
vl_venn <- function(dat, 
                    col= adjustcolor(c("cornflowerblue", "tomato", "gold"), 0.3), 
                    border= NA)
{
  stop("Function not working yet with more than two elements :(")
  
  if(!is.list(dat))
    stop("dat should be a named list")
  if(is.null(names(dat)))
    stop("dat should be a named list")
  if(!length(dat) %in% c(2,3))
    stop("dat should be a named list of length 2 or 3")
  if(length(col)!=length(dat))
    col <- rep(col, length(dat))
  if(length(border)!=length(dat))
    border <- rep(border, length(dat))
  
  
  # Functions
  circle_intersection <- function(x1, y1, r1, x2, y2, r2){
    rr1 <- r1 * r1
    rr2 <- r2 * r2
    d <- sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
    
    if (d > r2 + r1) # Circles do not overlap
    {
      return(0)
    } else if (d <= abs(r1 - r2) && r1 >= r2){ # Circle2 is completely inside circle1  
      return(pi*rr2)
    } else if (d <= abs(r1 - r2) && r1 < r2){ # Circle1 is completely inside circle2
      return(pi*rr1)
    } else { # Circles partially overlap
      phi <- (acos((rr1 + (d * d) - rr2) / (2 * r1 * d))) * 2
      theta <- (acos((rr2 + (d * d) - rr1) / (2 * r2 * d))) * 2
      area2 <- 0.5 * theta * rr2 - 0.5 * rr2 * sin(theta)
      area1 <- 0.5 * phi * rr1 - 0.5 * rr1 * sin(phi)
      return(area1 + area2)
    }
  }
  circle_coor <- function(x, y, r)
  {
    px <- cos(seq(0, 2*pi, length.out= 359))
    px <- px*r+x
    py <- sin(seq(0, 2*pi, length.out= 359))
    py <- py*r+y
    return(list(x= px, y= py))
  }
  intersect <- function(o1, o2, top= 2)
  {
    dat <- CJ(1:359, 1:359)
    dat <- cbind(data.table(x1= o1$x, y1= o1$y, ord1= 1:359)[dat[[1]]],
                 data.table(x2= o2$x, y2= o2$y)[dat[[2]]])
    dat[, dist:= sqrt((x1-x2)^2+(y1-y2)^2)]
    setorderv(dat, "dist")
    dat <- dat[seq(top)]
    setorderv(dat, "ord1")
    return(dat[, .(x= x1, y= y1, ord1)])
  }
  
  # Prepare data
  dat <- lapply(dat, unique)
  
  # Compute cirles
  r1 <- .25
  x1 <- .25
  y1 <- .25
  c1 <- circle_coor(x1, y1, r1)
  a1 <- pi*r1*r1
  
  r2 <- sqrt(length(dat[[2]])/length(dat[[1]])*r1*r1)
  x2 <- seq(x1+r2-r1, x1+r2*2, length.out = 100)
  y2 <- .25
  a2 <- length(dat[[2]])/length(dat[[1]])*a1
  inter <- sapply(x2, function(x) circle_intersection(x1, y1, r1, x, y1, r2))
  x2 <- x2[which.min(abs((a1*sum(dat[[1]] %in% dat[[2]])/length(dat[[1]]))-inter))]
  c2 <- circle_coor(x2, y2, r2)
  
  if(length(dat)==3)
  {
    r3 <- sqrt(length(dat[[3]])/length(dat[[1]])*r1*r1)
    y3 <- seq(y1+r2-r1, y1+r2*2, length.out = 100)
    a3 <- length(dat[[3]])/length(dat[[1]])*a1
    inter <- sapply(y3, function(x) circle_intersection(x1, y1, r1, x1, x, r3))
    y3 <- y3[which.min(abs((a1*sum(dat[[1]] %in% dat[[3]])/length(dat[[1]]))-inter))]
    coors <- circle_coor(x1, y1, y3-y1)
    inter <- mapply(function(x,y) {
      circle_intersection(x2, y2, r2, x, y, r3)
    }, x= coors$x[1:180], y= coors$y[1:180])
    x3 <- coors$x[which.min(abs((a2*sum(dat[[2]] %in% dat[[3]])/length(dat[[2]]))-inter))]
    y3 <- coors$y[which.min(abs((a2*sum(dat[[2]] %in% dat[[3]])/length(dat[[2]]))-inter))]
    c3 <- circle_coor(x3, y3, r3)
  }
  
  # Adjust the size
  if(length(dat)==3)
    x.adj <- max(c(c1$x, c2$x, c3$x)) else
      x.adj <- max(c(c1$x, c2$x))
  c1$x <- c1$x/x.adj
  c2$x <- c2$x/x.adj
  if(length(dat)==3)
    c3$x <- c3$x/x.adj
  if(length(dat)==3)
    y.min <- min(c(c1$y, c2$y, c3$y)) else
      y.min <- min(c(c1$y, c2$y))
  c1$y <- c1$y-y.min
  c2$y <- c2$y-y.min
  if(length(dat)==3)
    c3$y <- c3$y-y.min
  if(length(dat)==3)
    y.adj <- max(c(c1$y, c2$y, c3$y)) else
      y.adj <- max(c(c1$y, c2$y))
  c1$y <- c1$y/y.adj
  c2$y <- c2$y/y.adj
  if(length(dat)==3)
    c3$y <- c3$y/y.adj
  
  # Plot
  plot.new()
  polygon(c1$x, c1$y, col = col[1], border= border[1])
  polygon(c2$x, c2$y, col = col[2], border= border[2])
  if(length(dat)==3)
    polygon(c3$x, c3$y, col = col[3], border= border[3])
  
  # Add labels
  if(length(dat)==3)
  {
    int12 <- intersect(c1, c2)
    seg12 <- int12[, .(x= seq(x[1], x[2], length.out = 359), 
                       y= seq(y[1], y[2], length.out = 359))]
    cr12 <- intersect(seg12, c3, top = 1)
    text(mean(c(int12[2,x], cr12[1,x])), 
         mean(c(int12[2,y], cr12[1,y])),
         sum(dat[[1]] %in% dat[[2]] & !dat[[1]] %in% dat[[3]]))
    int13 <- intersect(c1, c3)
    seg13 <- int13[, .(x= seq(x[1], x[2], length.out = 359), 
                       y= seq(y[1], y[2], length.out = 359))]
    cr13 <- intersect(seg13, c2, top = 1)
    text(mean(c(int13[1,x], cr13[1,x])), 
         mean(c(int13[1,y], cr13[1,y])),
         sum(dat[[1]] %in% dat[[3]] & !dat[[1]] %in% dat[[2]]))
    int23 <- intersect(c2, c3)
    seg23 <- int23[, .(x= seq(x[1], x[2], length.out = 359), 
                       y= seq(y[1], y[2], length.out = 359))]
    cr23 <- intersect(seg23, c1, top = 1)
    text(mean(c(int23[1,x], cr23[1,x])), 
         mean(c(int23[1,y], cr23[1,y])),
         sum(dat[[2]] %in% dat[[3]] & !dat[[2]] %in% dat[[1]]))
    int123 <- intersect(seg12, seg23, 1)
    text(cr12$x, 
         mean(c(cr12$y, int12[1, y])),
         sum(dat[[1]] %in% dat[[2]] & dat[[1]] %in% dat[[3]]))
    text(c1$x[205],
         c1$y[205],
         pos= 4,
         paste0(names(dat)[1], "\n", sum(!dat[[1]] %in% unlist(dat[c(2, 3)]))))
    text(c2$x[335],
         c2$y[335],
         pos= 2,
         paste0(names(dat)[2], "\n", sum(!dat[[2]] %in% unlist(dat[c(1, 3)]))))
    text(c3$x[90],
         c3$y[90],
         pos= 1,
         paste0(names(dat)[3], "\n", sum(!dat[[3]] %in% unlist(dat[c(1, 2)]))))
  }else
  {
    text(mean(c(c1$x[180], c2$x[180])), 
         c2$y[180],
         paste0(names(dat)[1], "\n", sum(!dat[[1]] %in% dat[[2]])))
    text(mean(c(c1$x[1], c2$x[1])), 
         c1$y[1],
         paste0(names(dat)[2], "\n", sum(!dat[[2]] %in% dat[[1]])))
    text(mean(c(c2$x[180], c1$x[1])), 
         c1$y[1],
         sum(dat[[1]] %in% dat[[2]]))
  }
  
  # Return intersection
  stat <- rbindlist(lapply(dat, as.data.table), idcol = T)
  stat <- stat[, .(name= paste0(.id, collapse = "+")), V1]
  invisible(split(stat$V1, stat$name))
}
