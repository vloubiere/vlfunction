#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param dat Data.table containing the variables to plot
#' @param class_levels Levels used for ordering. e.g c("Up", "Unaffected", "Down").
#' @param col Colors to use for bars. Wrapped with colorRampPalette
#' @examples 
#' test <- data.table(c("A", "A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
#' c("A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"),
#' c("A", "A", "B", "B", "B", "B", "B", "B", "B", "C", "C"))
#' vl_alluvial_plot(test,
#' class_levels = c("C", "B", "A"))
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export



vl_alluvial_plot <- function(dat,
                             class_levels= NULL,
                             col= c("cornflowerblue", "lightgrey", "tomato"))
{
  if(is.null(class_levels))
    class_levels <- as.character(unique(unlist(dat)))
  DT <- dat[, lapply(.SD, factor, levels= class_levels)]
  
  # Make object
  res <- melt(DT, measure.vars = names(DT))
  res[, value:= factor(value, class_levels)]
  # Compute rectangles
  setorderv(res, c("variable", "value"))
  res <- res[, .N, .(variable, value)]
  res[, cumsum:= cumsum(N), variable]
  res[, Cc1:= colorRampPalette(col)(length(class_levels))[.GRP], value]
  res[, top:= cumsum/sum(N), variable]
  res[, bottom:= c(0, top[-(.N)]), variable]
  res[, left:=  seq(0, 1, length.out= length(DT)*2)[(.GRP-1)*2+1], variable]
  res[, right:= seq(0, 1, length.out= length(DT)*2)[.GRP*2], variable]
  # Compute connections
  res[, to_variable:= unique(res$variable)[c(2:length(unique(res$variable)), 1)][.GRP], variable]
  res <- res[, .(to_value= unique(res$value)), (res)]
  res[, count:= length(which(DT[[as.character(variable)]]==value & # all transition counts
                             DT[[as.character(to_variable)]]==to_value)), .(variable, to_variable, value, to_value)]
  # Compute polygons
  res[, Cc2:= colorRampPalette(col)(length(class_levels))[.GRP], to_value]
  res[, Cc:= colorRampPalette(c(Cc1, Cc2))(3)[2], .(Cc1, Cc2)]
  res[, y1:= cumsum(count)/sum(count), variable]
  res[, y4:= c(0, y1[-(.N)]), variable]
  setorderv(res, "to_value")
  res[, y2:= cumsum(count)/sum(count), to_variable]
  res[, y3:= c(0, y2[-(.N)]), to_variable]
  res[, left1:=  seq(0, 1, length.out= length(DT)*2)[.GRP*2], variable]
  res[, right1:= seq(0, 1, length.out= length(DT)*2)[.GRP*2+1], variable]
  
  # PLOT
  plot.new()
  res[, rect(left[1], 
             bottom[1], 
             right[1], 
             top[1], 
             col= Cc1[1], 
             border= NA), .(top, bottom, left, right, Cc1)]
  res[, polygon(c(left1[1], 
                  right1[1], 
                  right1[1], 
                  left1[1]), 
                c(y1[1], 
                  y2[1], 
                  y3[1], 
                  y4[1]),
                col= adjustcolor(Cc[1], 0.5),
                border= NA), .(left1, right1, y1, y2, y3, y4, Cc)]
  res[, text(mean(c(left, right)), 
             1, 
             variable, 
             pos= 3, 
             xpd=T), .(variable, left, right)]
  res[, text(mean(c(left, right)), 
             mean(c(top, bottom)), 
             paste0(value, "\n(", N,")"), 
             xpd=T), .(value, N, left, right, top, bottom)]
  
  invisible(res)
}