#' bw heatmap plot only
#'
#' Plots heatmap from a vl_heatmap_bw_track output object
#'
#' @param obj An object returned by the vl_heatmap_bw_track() function. The object can easily be modified to change colors and so on...
#' @param center_label Label for the center of the heatmaps
#' @param col Color vector to be used for heatmaps. col= c("blue", "yellow", "white")
#' @export

vl_heatmap_bw_track_plot_only <- function(obj,
                                          center_label= "TSS",
                                          col= c("blue", "yellow", "white"))
{
  Ntracks <- length(unique(obj$track))
  
  # Compute plotting colors
  Cc <- colorRampPalette(col)(101)
  obj[, plot_Cc:= Cc[.GRP], keyby= col_idx]
  
  # Compute image
  im <- dcast(obj, 
              region_order~track+bin_ID, 
              value.var = "plot_Cc")
  mat <- as.matrix(im, 
                   1)
  par(xaxs= "i",
      yaxs= "i",
      mai= c(1.02, 
             max(strwidth(unique(obj$set_ID), "inches"))+0.5, 
             0.82, 
             0.42))
  
  # Plot heatmap
  plot.new()
  rasterImage(mat, 0,0,1,1, interpolate = F)
  
  # Plot lines
  .lv <- seq(0, 1, length.out = length(unique(obj$track))+1)
  abline(v= .lv, 
         lwd= 0.25)
  .lh <- obj[, 1-(min(region_order)-1)/max(obj$region_order), set_ID]$V1
  .lh <- c(.lh, 0)
  abline(h= .lh, 
         lwd= 0.25)
  
  # Plot axes
  axAt <- .lv[-length(.lv)]+diff(.lv)/2
  axis(1, 
       at= axAt, 
       labels = rep(center_label, Ntracks),
       lwd= 0, 
       lwd.ticks = 1, 
       line = 0)
  axis(1,
       at= .lv[-length(.lv)]+strwidth(obj$ext1[1], "user")/2, 
       labels= rep(obj$ext1[1], Ntracks),
       lwd= 0, 
       lwd.ticks = 0, 
       line = -0.5, 
       cex.axis= 0.8)
  axis(1,
       at= .lv[-1]-strwidth(obj$ext2[1], "user")/2, 
       labels= rep(obj$ext2[1], Ntracks),
       lwd= 0, 
       lwd.ticks = 0, 
       line = -0.5, 
       cex.axis= 0.8)
  axis(2, 
       at= .lh[-length(.lh)]+diff(.lh)/2, 
       labels = unique(sort(obj$set_ID)),
       las= 1,
       lwd= 0, 
       lwd.ticks= 0)
  axis(3, 
       at= axAt, 
       labels = obj[, unique(name), keyby= track]$V1,
       lwd= 0, 
       lwd.ticks = 0, 
       line = 0)
  
  # Plot legends
  rasterImage(matrix(colorRampPalette(col)(100), nrow= 1),
              xleft = axAt-1/Ntracks*0.3, 
              ybottom = grconvertY(1.35, "lines", "user"), 
              xright = axAt+1/Ntracks*0.3,
              ytop = grconvertY(1.65, "lines", "user"),
              xpd= T)
  rect(xleft = axAt-1/Ntracks*0.3, 
       ybottom = grconvertY(1.35, "lines", "user"), 
       xright = axAt+1/Ntracks*0.3,
       ytop = grconvertY(1.65, "lines", "user"),
       xpd= T, 
       lwd= 0.5)
  obj[, {
    text(axAt[.GRP]-1/Ntracks*0.3,
         grconvertY(1.5, "lines", "user"), 
         pos= 2,
         0,
         xpd= T, 
         cex= 0.8)
    text(axAt[.GRP]+1/Ntracks*0.3,
         grconvertY(1.5, "lines", "user"), 
         pos= 4,
         formatC(max, 
                 format = "g", 
                 digits = 2),
         xpd= T, 
         cex= 0.8)
  }, .(track, max)]
}