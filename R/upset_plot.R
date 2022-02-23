#' upset Plot
#'
#' Plots an upset plot convenient to visualize intersections
#'
#' @param dat_list List containing the names to intersect (see examples)
#' @param intersection_cutoff min cutoff for selections to be shown
#' @examples 
#' test <- list(A= 1:1000, B= 1:1500, C= 1000:1750)
#' vl_upset_plot(test)
#' @export

vl_upset_plot <- function(dat_list, 
                          intersection_cutoff= 0)
{
  if(is.null(names(dat_list)))
    stop("The list should be named!")
  if(!identical(unique(names(dat_list)), names(dat_list)))
    stop("dat_list names should be unique!")
  
  # Format
  dat <- rbindlist(lapply(dat_list, as.data.table), idcol = T)
  setnames(dat, "V1", "var")
  dat[, .id:= factor(.id,
                     levels= names(dat_list)[order(lengths(dat_list))])]
  # Intersections
  inter <- dcast(dat, 
                 var~.id, 
                 value.var = "var",
                 fun.aggregate = function(x) ifelse(any(!is.na(x)), 1, 0))
  inter <- inter[, .N, setdiff(names(inter), "var")][N>intersection_cutoff]
  setorderv(inter, "N", -1)
  # sets
  set <- dat[, .N, keyby= .id]
  
  #------------------#
  # PLOT
  #------------------#
  par(mai= c(grconvertY(nrow(set)+1, "lines", "inch"),
             0.5+max(strwidth(set$.id, "inch"))+strwidth("M", "inch")*10,
             0.42, 
             0.42))
  # Barplot intersections
  inter[, x:= barplot(N, 
                      las= 1,
                      ylab= "Intersection size",
                      border= NA,
                      col= "grey20")]
  inter[, text(x, N, N, cex= 0.8, pos= 3, xpd= T)]
  # Barplot sets
  right <- par("usr")[1]-strwidth("MM", "user")-max(strwidth(set$.id, "user"))
  plot.width <- right-grconvertX(par("fig")[1]+0.015, "ndc", "user")
  xaxs.max <- max(axisTicks(range(set$N), log= F))
  xaxs.w <- xaxs.max/max(set$N)*plot.width
  set[, top:= grconvertY(grconvertY(.I+0.5, "lines", "ndc")+par("fig")[3], "ndc", "user")]
  set[, height:= grconvertY(1, "lines", "user")-grconvertY(0, "lines", "user")]
  set[, width:= N/max(N)*plot.width]
  set[, {
    rect(right-width, 
         top-height, 
         right, 
         top, 
         xpd= T, 
         border= "white", 
         col= "grey20")
    lines(c(right[1]-xaxs.w, right[1]-xaxs.w, right[1], right[1]),
          par("usr")[3]+c(strheight("M")/10, 0, 0, strheight("M")/10),
             xpd= T)
    text(c(right[1], right[1]-xaxs.w/2, right[1]-xaxs.w),
         par("usr")[3],
         pos= 3,
         labels = c(0, "Set size", xaxs.max),
         xpd= T,
         offset= 0.25,
         cex= c(0.5, 0.8, 0.5))
    text(par("usr")[1],
         top-height/2,
         pos= 2,
         .id,
         xpd= T)
  }]
  # Grid 
  grid <- melt(inter[, !"N"], "x")[set[, .(y= top-height/2), .id], on= "variable==.id"]
  grid[, col:= ifelse(value==0,
                      "lightgrey",
                      "grey20")]
  grid[, points(x, 
           y, 
           xpd= T,
           col= col,
           pch= 19, 
           cex= 2)]
  grid[value==1, segments(x[1], 
                          min(y), 
                          x[1], 
                          max(y), 
                          xpd= T, 
                          col= "grey20", 
                          lwd= 2), x]
}
