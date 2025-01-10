#' Add motifs to enrichment plots
#'
#' @param DT DT object output from vl_motif_enrich or vl_motif_cl_enrich
#' @param cex.width expansion factor for motif widths
#' @param cex.height expansion factor for motif heights
#' 
#' @examples 
#' For vl_enr object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl)
#' 
#' For vl_enr_cl object
#' pl <- plot(vl_enr)
#' vl_add_motif(pl$DT)
#' 
#' @export
vl_add_motifs <- function(DT,
                          cex.width= 1,
                          cex.height= 1,
                          lwd= 0.1)
{
  # Extract PWMs
  DT <- unique(DT[, .(variable, name, y)])
  mats <- vl_Dmel_motifs_DB_full[DT, pwms_perc, on= "motif_ID==variable"]
  mats <- lapply(mats, TFBSTools::as.matrix)
  # Compute plotting coordinates
  ax.width <- max(strwidth(DT$name, cex= par("cex.axis")))+diff(grconvertX(c(0, par("mgp")[2]+0.5), "lines", "user"))
  coor <- vl_seqlogo(pwm = mats, 
                     x = par("usr")[1]-ax.width, 
                     y = DT$y, 
                     cex.width = cex.width,
                     cex.height = cex.height,
                     min_content= 0.05)
  # Plot
  coor[, segments(xleft, 
                  ybottom, 
                  xright, 
                  ybottom, 
                  xpd= T, 
                  col= "grey60",
                  lwd= lwd*par("lwd"))]
  invisible(coor)
}

#' plot seqlogo
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm List of pwm matrices
#' @param x x positions
#' @param y positions (centered)
#' @param pos eith 2 (left) of 4 (right)
#' @param cex.width width expansion factor applied before plotting motifs
#' @param cex.height height expansion factor applied before plotting motifs
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' 
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' rownames(pwm) <- c("A", "C", "G", "T")
#' plot.new()
#' vl_seqlogo(pwm)
#' 
#' @export
vl_seqlogo <- function(pwm, 
                       x,
                       y,
                       pos= 2,
                       cex.width= 1,
                       cex.height= 1,
                       add= T,
                       min_content= 0)
{
  if(!pos %in% c(2,4))
    stop("Unsupported pos value. Use either 2 (left) or 4 (right)")
  if(is.matrix(pwm))
    pwm <- list(pwm)
  
  # Make object and index
  obj <- data.table(pwm, x, y, cex.width, cex.height)
  obj[, idx:= .I]
  
  # Width only depends on cex
  obj[, width:= strwidth("M", cex= cex.width), cex.width]
  
  # For each base, compute xleft, xright, ytop, ybottom
  pl <- obj[, {
    # Import PWM and melt
    .c <- as.data.table(pwm[[1]], keep.rownames= "base")
    .c <- melt(.c, id.vars = "base")
    # Compute motif content per column and normalize importance
    .c[, content:= sum(value*log2(value/c(0.25, 0.25, 0.25, 0.25)), na.rm= T), variable]
    .c[, norm:= value*(content/max(content))]
    # Remove flanks with little content (5% max)
    setorderv(.c, "variable")
    .c <- .c[min(which(norm>min_content)):max(which(norm>min_content))]
    # xleft depends on the pos (2 or 4)
    if(pos==4)
    {
      # Already correclty sorted earlier
      .c[, xleft:= x+((.GRP-1)*width), variable]
    }else if(pos==2)
    {
      setorderv(.c, "variable", -1)
      .c[, xleft:= x-(.GRP*width), variable]
    }
    # Rank from lowest to biggest importance -> inscreasing ytop pos
    setorderv(.c, "norm")
    .h <- strheight("M", cex= cex.height)
    .c[, c("height", "ytop"):= {
      heights <- norm*.h
      .(heights, (y-.h/2)+cumsum(heights))
    }, variable]
  }, .(idx, y, width)]
  
  # Plot
  pl[, vl_plotLetter(base[1], 
                     xleft[1], 
                     ytop[1], 
                     width[1], 
                     height[1]), .(base, xleft, ytop, width, height)]
  
  # Return object containing limits of each motif
  invisible(pl[, .(xleft= min(xleft), 
                   ybottom= min(ytop-height),
                   xright= max(xleft+width), 
                   ytop= max(ytop)), .(idx)])
}

#' plot seqlogo letter
#'
#' See function vl_seqlogo
#'
#' @param letter "A", "T", "C" or "G"
#' @param xleft xleft position
#' @param ytop ytop position
#'
#' @export
vl_plotLetter <- function(letter, xleft, ytop, width, height)
{
  letter <- toupper(letter)
  if(letter=="T")
  {
    x <- c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1
    y <- c(10, 10, 8.5, 8.5, 0, 0, 8.5, 8.5) * 0.1
    col <- "red"
  }else if(letter=="A")
  {
    x <- c(0, 4, 6, 2, 0, 4, 6, 10, 8, 4, 3.2, 6.8, 6.4, 3.6, 3.2) * 0.1
    y <- c(0, 10, 10, 0, 0, 10, 10, 0, 0, 10, 3, 3, 4, 4, 3) * 0.1
    col <- "forestgreen"
  }else if(letter=="C")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    col <- "dodgerblue2"
  }else if(letter=="G")
  {
    angle1 <- seq(0.3 + pi / 2, pi, length = 100)
    angle2 <- seq(pi, 1.5 * pi, length = 100)
    x.l1 <- 0.5 + 0.5 * sin(angle1)
    y.l1 <- 0.5 + 0.5 * cos(angle1)
    x.l2 <- 0.5 + 0.5 * sin(angle2)
    y.l2 <- 0.5 + 0.5 * cos(angle2)
    x.l <- c(x.l1, x.l2)
    y.l <- c(y.l1, y.l2)
    x <- c(x.l, rev(x.l))
    y <- c(y.l, 1 - rev(y.l))
    x.i1 <- 0.5 + 0.35 * sin(angle1)
    y.i1 <- 0.5 + 0.35 * cos(angle1)
    x.i1 <- x.i1[y.i1 <= max(y.l1)]
    y.i1 <- y.i1[y.i1 <= max(y.l1)]
    y.i1[1] <- max(y.l1)
    x.i2 <- 0.5 + 0.35 * sin(angle2)
    y.i2 <- 0.5 + 0.35 * cos(angle2)
    x.i <- c(x.i1, x.i2)
    y.i <- c(y.i1, y.i2)
    x1 <- c(x.i, rev(x.i))
    y1 <- c(y.i, 1 - rev(y.i))
    x <- c(x, rev(x1))
    y <- c(y, rev(y1))
    h1 <- max(y.l1)
    r1 <- max(x.l1)
    h1 <- 0.4
    x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
    y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
    x <- c(rev(x), x.add)
    y <- c(rev(y), y.add)
    col <- "goldenrod1"
  }else
    stop("DNA letter not recognized. Other than ATCG?")
  polygon(x = xleft+x*width, 
          ytop-(1-y)*height, 
          col= col,
          border= NA,
          xpd= NA)
}