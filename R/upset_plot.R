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
  # Intersections
  inter <- dcast(dat, 
                 var~.id, 
                 value.var = "var",
                 fun.aggregate = function(x) ifelse(any(!is.na(x)), 1, 0))
  inter <- inter[, .N, setdiff(names(inter), "var")]
  inter <- inter[N>=intersection_cutoff]
  setorderv(inter, "N", -1)
  # sets
  set <- dat[, .N, keyby= .id]
  setorderv(set, "N", -1)
  set[, I:= .I]
  
  #------------------#
  # PLOT
  #------------------#
  # Barplot intersections
  inter[, x:= barplot(N, 
                      las= 1,
                      ylab= "Intersection size",
                      border= NA,
                      col= "grey20",
                      xaxt= "n")]
  browser()
  inter[, text(x, N, N, cex= 0.8, pos= 3, xpd= TRUE)]
  # Add grid
  grid <- melt(inter[, !"N"], id.vars = "x")
  grid[set, y:= par("usr")[3]-strheight("M", cex= 2)*i.I, on= "variable==.id"]
  grid[, col:= ifelse(value==0, "lightgrey", "grey20")]
  grid[, points(x, 
                y, 
                xpd= TRUE, 
                cex= 2, 
                pch= 19, 
                col= col)]
  grid[value==1, segments(x[1], 
                          min(y), 
                          x[1], 
                          max(y), 
                          col= "grey20", 
                          xpd= TRUE,
                          lwd= 2), x]
  # Sets names
  set[grid, y:= i.y, on= ".id==variable", mult= "first"]
  set[, text(par("usr")[1], y, .id, xpd= TRUE, pos= 2)]
  # Sets barplot
  tot_width <- strwidth(paste0(max(set$N), " Set size  0"))
  set[, top:= y+strheight("M")*0.75]
  set[, height:= strheight("M")*1.5]
  set[, right:= par("usr")[1]-max(strwidth(paste0("M", .id)))]
  set[, width:= N/max(N)*tot_width]
  set[, rect(right-width, 
             top-height, 
             right, 
             top, 
             border= NA, 
             col= "grey20",
             xpd= TRUE)]
  # Sets barplot axis
  s.max <- max(axisTicks(c(0, max(set$N)), log= FALSE))
  x <- c(set[1,right], set[1,right]-(s.max/max(set$N)*tot_width))
  y <- min(set[, y-height])
  y <- c(y-strheight("M")*0.1, y)
  lines(rep(x, each= 2),
        c(y, rev(y)),
        xpd= TRUE)
  text(mean(x),
       min(set[, y-height]),
       offset= 0.8,
       "Set size",
       pos= 1,
       xpd= TRUE)
  text(x[c(1, 2)],
       y[c(1, 1)],
       c(0, s.max),
       pos= 1,
       xpd= TRUE,
       cex= 0.5,
       offset= 0.2)
  
  invisible(list(inter= inter, grid= grid, set= set))
}
