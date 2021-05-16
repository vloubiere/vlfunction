#' upset Plot
#'
#' Plots an upset plot convenient to visualize intersections
#'
#' @param dat_list List containing the names to intersect (see examples)
#' @param ylab Ylab main barplot
#' @examples 
#' test <- list(A= 1:1000, B= 1:1500, C= 1000:1750)
#' vl_upset_plot(test)
#' @export

vl_upset_plot <- function(dat_list, ylab= "Intersection size")
{
  if(is.null(names(dat_list)))
    stop("The list should be named!")
  if(any(grepl("\\|", names(dat_list))))
    stop("list names should not contain any '|' cause they are use internally")
  
  # Initiate plotting area
  init <- par(no.readonly=TRUE)
  plot.new()
  par(mar = c(length(dat_list)*2+1,
              grconvertX(max(strwidth(names(dat_list))), to= "lines"),
              1,
              1))
  
  # Format dat
  dat <- rbindlist(lapply(dat_list, as.data.table), idcol = T)
  colnames(dat)[2] <- "intersect"
  dat <- dat[, .(.id= paste(.id, collapse = "|")), intersect]
  dat <- dat[, .(N= .N), .id]
  setorderv(dat, "N", -1)
  dat$x <- barplot(dat$N, plot = F)
  dat <- dat[, .(all_IDs= unlist(tstrsplit(.id, "\\|"))), dat]
  sets <- dat[, sum(N), all_IDs]
  setorderv(sets, "V1", 1)
  sets$y <- barplot(-sets$V1,
                    plot = F)
  dat[sets, y:= i.y, on= "all_IDs"]
  dat <- dat[, .(y= .(y)), .(.id, N, x)]
  
  # Plot main barplot
  barplot(dat$N,
          las= 1,
          border= NA,
          col= "grey20",
          ylab= ylab)
  text(dat$x,
       dat$N,
       dat$N,
       cex=0.5,
       xpd=T,
       pos=3,
       offset= 0.25)
  dat[, x:= grconvertX(x, to= "ndc")]
  inset <- c(grconvertX(0.5, from= "lines", to= "ndc"),
             grconvertX(0, to= "ndc")-max(strwidth(names(dat_list), "figure"))-grconvertX(1, "lines", "ndc"),
             0,
             grconvertY(0, "user", to= "ndc")-grconvertY(0.75, "lines", "ndc"))
  
  # Plot set barplot
  par(fig = inset,
      mar= c(1,1,0,0),
      xaxs= "i",
      yaxs= "i",
      new = T)
  barplot(-sets$V1,
          horiz = T,
          axes= F,
          xpd= T,
          border= NA,
          col= "grey20")
  ticks <- axisTicks(c(0, max(sets$V1)), log= F)
  axis(3,
       at= -ticks,
       labels= NA,
       line = 0.5,
       tcl= -0.1)
  mtext(text = range(ticks),
        at= -range(ticks),
        side= 3,
        line= 0.75,
        cex= 0.6,
        xpd= T)
  mtext("Set size",
        line = 1.5,
        cex= 0.6)
  axis(4,
       at= sets$y,
       label= sets$all_IDs,
       las= 1,
       lwd= NA)
  dat[, y:= .(list(grconvertY(unlist(y), to= "ndc"))), .id]
  
  # Plot Intersection scheme
  par(fig = c(0,1,0,1))
  dat[, x:= grconvertX(x, "ndc", "user")]
  dat[, y:= .(list(grconvertY(unlist(y), "ndc", "user"))), .id]
  
  frame <- CJ(unique(dat$x), unique(unlist(dat$y)))
  points(frame$V1,
         frame$V2,
         col= "grey90",
         pch= 19,
         cex= 2)
  dat[, segments(x[1],
                 min(unlist(y)),
                 x[1],
                 max(unlist(y)),
                 xpd= T,
                 lwd = 2,
                 col= "grey20"), .id]
  dat[, points(rep(x, lengths(y)),
               unlist(y),
               col= "grey20",
               pch= 19,
               cex= 2), .id]
  
  on.exit(par(init),add=TRUE,after=FALSE)
}
