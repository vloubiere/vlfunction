#' upset Plot
#'
#' Plots an upset plot convenient to visualize intersections
#'
#' @param dat.list List containing the names to intersect (see examples)
#' @param col Color used for intersection bars
#' @param col.set Color used for set size bars
#' @param cex.grid Cex for grid points
#' @param grid.hex Height of grid lines
#' @param cex.inter Cex size intersection labels
#' @param show.empty Should empty instersections be shown? default= F
#' @param lwd line width
#' @param xaxt Should x axis be drawn? default= "n
#' @param intersection.cutoff min cutoff for selections to be shown
#' @param border Default to NA
#' @param ylab ylab
#' @param las las
#' @examples 
#' test <- list(A= 1:1000, B= 1:1500, C= 1000:1750)
#' vl_upset_plot(test)
#' @export
vl_upset_plot <- function(dat.list,
                          col= "grey20",
                          col.set= "grey20",
                          cex.grid= 2,
                          grid.hex= 2,
                          cex.inter= .8,
                          show.empty= F,
                          lwd= 1,
                          xaxt= "n",
                          intersection.cutoff= 0,
                          border= NA,
                          ylab= "Intersection size",
                          las= par("las"))
{
  if(is.null(names(dat.list)))
    stop("The list should be named!")
  if(!identical(unique(names(dat.list)), names(dat.list)))
    stop("dat.list names should be unique!")
  dat.list <- lapply(dat.list, unique)
  
  # Format ----
  dat <- rbindlist(lapply(dat.list, as.data.table), idcol = T)
  dat[, .id:= factor(.id, names(dat.list))]
  setnames(dat, "V1", "var")
  # Intersections ----
  inter <- dcast(dat, 
                 var~.id, 
                 value.var = "var",
                 fun.aggregate = function(x) ifelse(any(!is.na(x)), 1, 0),
                 drop = !(show.empty))
  inter <- inter[, .N, setdiff(names(inter), "var")]
  inter <- inter[N>=intersection.cutoff]
  setorderv(inter, "N", -1)
  # sets
  set <- as.data.table(table(dat$.id))
  setnames(set, c(".id", "N"))
  setorderv(set, "N", -1)
  set[, I:= .I-1]
  
  # PLOT ----
  ## Barplot intersections ----
  inter[, x:= barplot(N,
                      las= las,
                      ylab= ylab,
                      border= border,
                      col= col,
                      xaxt= xaxt)]
  inter[, text(mean(x), N, N[1], cex= cex.inter, pos= 3, xpd= TRUE), N]
  ## Add grid ----
  grid <- melt(inter[, !"N"], id.vars = "x")
  grid[, y:= par("usr")[3]-strheight("M", cex= grid.hex*1.25)]
  grid[set, y:= y-strheight("M", cex= grid.hex*1.25)*i.I, on= "variable==.id"]
  grid[, col:= ifelse(value==0, "lightgrey", "grey20")]
  grid[, points(x, 
                y, 
                xpd= TRUE, 
                cex= cex.grid, 
                pch= 19, 
                col= col)]
  grid[value==1, segments(x[1], 
                          min(y), 
                          x[1], 
                          max(y), 
                          col= "grey20", 
                          xpd= TRUE,
                          lwd= lwd*par("lwd")), x]
  # Sets names ----
  set[grid, y:= i.y, on= ".id==variable", mult= "first"]
  prev <- par("cex.axis")
  par(cex.axis= par("cex.lab"))
  axis(2,
       set$y,
       set$.id,
       xpd= NA,
       lwd= 0,
       line = 0,
       gap.axis = 0)
  par(cex.axis= prev)
  # Sets barplot ----
  tot_width <- strwidth(paste0(max(set$N), "Set size"), cex= par("cex.lab"))
  set[, top:= y+strheight("M", cex= grid.hex)/2]
  set[, height:= strheight("M", cex= grid.hex)]
  set[, right:= par("usr")[1]-diff(grconvertX(c(0, par("mgp")[2]), "lines", "user"))]
  set[, right:= right-max(strwidth(paste0("M", .id), cex= par("cex.lab")))]
  set[, width:= N/max(N)*tot_width]
  set[, rect(right-width, 
             top-height, 
             right, 
             top, 
             border= NA, 
             col= col.set,
             xpd= TRUE)]
  # Sets barplot axis ----
  s.max <- max(axisTicks(c(0, max(set$N)), log= FALSE))
  x <- c(set[1,right], set[1,right]-(s.max/max(set$N)*tot_width))
  y <- min(set[, y-height], na.rm = T)
  y <- c(y-strheight("M")*0.1, y)
  lines(rep(x, each= 2),
        c(y, rev(y)),
        xpd= TRUE)
  text(mean(x),
       min(set[, y-height], na.rm = T),
       offset= par("mgp")[1],
       "Set size",
       pos= 1,
       cex= par("cex.lab"),
       xpd= TRUE)
  text(x[c(1, 2)],
       y[c(1, 1)],
       c(0, s.max),
       pos= 1,
       xpd= TRUE,
       cex= par("cex.axis"),
       offset= par("mgp")[2])
  
  invisible(list(inter= inter,
                 grid= grid,
                 set= set))
}
