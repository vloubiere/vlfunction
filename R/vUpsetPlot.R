#' upset Plot
#'
#' Plots an upset plot convenient to visualize intersections
#'
#' @param dat_list List containing the names to intersect (see examples)
#' @param ylab Ylab main barplot
#' @param intersection_cutoff min cutoff for selections to be shown
#' @examples 
#' test <- list(A= 1:1000, B= 1:1500, C= 1000:1750)
#' vl_upset_plot(test)
#' @export

vl_upset_plot <- function(dat_list, 
                          ylab= "Intersection size",
                          intersection_cutoff= 0)
{
  if(is.null(names(dat_list)))
    stop("The list should be named!")
  if(any(grepl("\\|", names(dat_list))))
    stop("list names should not contain any '|' cause they are use internally")
  if(max(nchar(names(dat_list)))<3)
    names(dat_list) <- paste0("  ", names(dat_list), "  ")
  # if(par("mfrow"))
    
  #------------------------#
  # Main object
  #------------------------#  
  dat <- rbindlist(lapply(dat_list, as.data.table), idcol = T)
  colnames(dat)[2] <- "intersect"
  dat <- dat[, .(.id= paste(.id, collapse = "|")), intersect]
  dat <- dat[, .(N= .N), .id]
  dat <- dat[N>=intersection_cutoff]
  setorderv(dat, "N", -1)
  
  #------------------------#
  # Initiate plotting area
  #------------------------#
  N_cditions <- length(unique(na.omit(unlist(tstrsplit(dat$.id, "\\|")))))
  plot.new()
  par(mar = c(grconvertY(N_cditions*1.5+2, "chars", to= "lines"),
              6+grconvertX(max(strwidth(names(dat_list), "inches")), "inches", to= "lines"),
              2,
              1))
  
  #-------------------------#
  # Main barplot
  #-------------------------#
  width <- 0.9/nrow(dat)
  space <- 0.1/nrow(dat)
  dat[, left:= grconvertX(space*.I+width*(.I-1), "npc", "user")]
  dat[, right:= grconvertX(space*(.I-1)+width*.I, "npc", "user")]
  dat[, top:= grconvertY(N/max(N), "npc", "user")]
  dat[, bottom:= grconvertY(0, "npc", "user")]
  dat[, x:= rowMeans(.SD), .SDcols= c("left", "right")]
  dat[, rect(left[1],
             bottom[1],
             right[1],
             top[1], 
             border = NA, 
             col= "grey20",
             xpd= T), (dat)]
  # Print N on top
  dat[, text(x[1], 
             top[1],
             N[1],
             xpd= T,
             cex= 0.6,
             pos=3), (dat)]
  
  # Y axis
  ticks <- axisTicks(c(0, max(dat$N)), log= F)
  at <- grconvertY(ticks/max(dat$N), "npc", "user")
  axis(2,
       at= at,
       labels= NA,
       line = 0.5,
       tcl= -0.1,
       xpd= T)
  mtext(text = ticks,
        at= at,
        side= 2,
        line= 1,
        cex= 0.8,
        xpd= T,
        las= 1)
  title(ylab= ylab)
  
  #-------------------------#
  # Intersections
  #-------------------------#
  sets <- dat[, .(all_IDs= unlist(tstrsplit(.id, "\\|"))), dat]
  sets <- sets[, .(N= sum(N)), all_IDs]
  setorderv(sets, "N", 1)
  sets[, y:= grconvertY(1+.I*1.5, "chars", "user")]
  sets[, all_x:= .(.(dat$x))]
  tab <- dat[, .(check_IDs= unlist(tstrsplit(.id, "\\|")), x), .id]
  sets[, inter_x:= .(.(unlist(all_x) %in% tab[check_IDs %in% all_IDs, x])), all_IDs]
  # points
  sets[, {
    cx <- unlist(all_x)
    cy <- rep(y[1], length(cx))
    cxi <- unlist(inter_x)
    points(cx,
           cy, 
           pch=19, 
           cex=2,
           col= ifelse(cxi, "grey20", "grey80"),
           xpd= T)
  }, .(all_IDs)]
  # Segments
  seg <- data.table(x= unlist(sets$all_x), 
                    y= rep(sets$y, lengths(sets$all_x)))[unlist(sets$inter_x)]
  seg[, segments(x[1], 
                 min(y), 
                 x[1], 
                 max(y), 
                 xpd= T, 
                 lwd=2,
                 col= "grey20"), x]
  
  #-------------------------#
  # Sets barplot
  #-------------------------#
  sets[, text(grconvertX(0, "npc", "user"), 
              y[1], 
              labels = all_IDs[1], 
              pos= 2, 
              xpd= T), .(all_IDs, y)]
  plot_left <- grconvertX(0.02, from = "nfc", "user")
  plot_right <- grconvertX(0, from = "npc", "user")-max(strwidth(paste0("   ", names(dat_list))))
  sets[, right:= plot_right]
  sets[, left:= right-(N/max(N)*(plot_right-plot_left))]
  sets[, top:= y+strheight("A", "user")]
  sets[, bottom:= y-strheight("A", "user")]
  # Barplot
  sets[, rect(left[1],
              bottom[1],
              right[1],
              top[1], 
              border = NA, 
              col= "grey20",
              xpd= T), right:bottom]
  # Axis
  ticks <- axisTicks(c(0, max(sets$N)), log= F)
  ticks <- ticks[c(1, length(ticks))]
  at <- c(plot_right,
          plot_right-(ticks[2]/max(sets$N)*(plot_right-plot_left)))
  segments(at[1], 
           grconvertY(0, "npc", "user"), 
           at[2], 
           grconvertY(0, "npc", "user"), 
           xpd=T)
  segments(at, 
           grconvertY(0, "npc", "user"), 
           at, 
           grconvertY(grconvertY(0.1, "chars", "ndc"), "npc", "user"), 
           xpd=T)
  text(x= at, 
       grconvertY(0, "npc", "user"),
       ticks, 
       cex= 0.5, 
       xpd= T, 
       pos= 3)
  text(x= mean(at), 
       grconvertY(0, "npc", "user"),
       labels = "Size", 
       cex= 0.8, 
       xpd= T, 
       pos= 3)
  
  # on.exit(par(init), add=TRUE, after=FALSE)
  invisible(sets)
}
