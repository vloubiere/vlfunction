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


#' Plot pval
#'
#' Plots formated pval data as text
#'
#' @param x X position for plotting
#' @param y Y position for plotting
#' @param pval Pval to be plotted
#' @param stars_only If set to TRU, then only plots */N.S. Default= FALSE
#' @examples 
#' pval <- c(1e-10, 1e-5, 1e-3, 1e-2, 1e-1, 1)
#' plot(NA, xlim= c(0,7), ylim= c(-1,1))
#' vl_plot_pval(1:6, 0.5, pval)
#' @export

vl_plot_pval_text <- function(x, 
                              y, 
                              pval, 
                              stars_only= F)
{
  if(length(y)==1 & length(x)>1)
    y <- rep(y, length(x))
  star <- cut(pval, 
              breaks = c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), 
              labels = c("****", "***", "**", "*", "N.S"))
  star <- as.character(star)
  ns_val <- !is.na(star) & star=="N.S"
  value <- formatC(pval, format = "e", digits = 1)
  if(stars_only)
    value <- rep("", length(value))
  if(any(ns_val))
    text(x = x[ns_val], 
         y = y[ns_val], 
         labels= bquote(paste(.(value[ns_val]), ""^N.S)), 
         cex= 0.6,
         offset= -0.1,
         pos= 3)
  if(!all(ns_val))
  text(x = x[!ns_val], 
       y = y[!ns_val], 
       labels= paste0(value[!ns_val], star[!ns_val]),
       offset= -0.2,
       pos= 3)
}

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
  dat <- dat[, .(N= .N, vars= paste0(intersect, collapse = ",")), .id]
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
             xpd= T), .(left, bottom, right, top)]
  # Print N on top
  dat[, text(x[1], 
             top[1],
             N[1],
             xpd= T,
             cex= 0.6,
             pos=3), .(x, top, N)]
  
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
  sets[, y:= grconvertY(1+.I*1.5+grconvertY(0, "nfc", "chars"), "chars", "user")]
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
  invisible(list(dat=dat,
                 sets= sets))
}


#' plot seqlogo rasterImage
#'
#' plot seqlogo from pwm matrix
#'
#' @param pwm pwm matrix
#' @param xleft left plot limit. Default= 0
#' @param ybottom bottom plot limit. Default= 0
#' @param xright right plot limit. Default= 1
#' @param ytop top plot limit. Default= 1
#' @param add Should the pwm be plot on the top of opened device? Default= T
#' @examples
#' pwm <- matrix(c(0.33,0.21,0.33,0.13,0.1,0.42,0.38,0.1,0.26,0.26,0.27,0.21,0,0.03,0.19,0.78,0.1,0.05,0.1,0.75,0.24,0.05,0.18,0.53,0.8,0.04,0,0.16,0.13,0.16,0.02,0.69,0.04,0.05,0.7,0.21,0.24,0.09,0.57,0.1,0.02,0.8,0.15,0.03,0.22,0.28,0.31,0.19,0.35,0.26,0.26,0.13,0.19,0.33,0.26,0.22), nrow= 4)
#' plot.new()
#' vl_seqlogo(pwm)
#' @export

vl_seqlogo <- function(pwm, 
                       xleft= 0, 
                       ybottom= 0, 
                       xright= 1, 
                       ytop= 1, 
                       add= T)
{
  if(!is.matrix(pwm))
  {
    stop("!is.matrix(pwm)")
  }
  require(fields)
  require(seqLogo)
  require(png)
  require(colorspace)
  
  tmp <- base::tempfile(fileext = ".png") 
  grDevices::png(tmp, type="cairo", width = 1000, height = 1000, units = "px")
  seqLogo::seqLogo(pwm, xaxis = F, yaxis = F)
  dev.off()
  im <- png::readPNG(tmp)
  res <- matrix(NA, nrow = nrow(im[,,1]), ncol = ncol(im[,,1]))
  res[im[,,1]==1 & im[,,2]==0 & im[,,3]==0] <- "firebrick1"
  res[im[,,1]==1 & im[,,2]>0.1 & im[,,2]<0.9 & im[,,3]==0] <- "goldenrod1"
  res[im[,,1]==0 & im[,,2]==1 & im[,,3]==0] <- "forestgreen"
  res[im[,,1]==0 & im[,,2]==0 & im[,,3]==1] <- "dodgerblue2"
  if(!add) plot.new()
  rasterImage(res[50:950,50:950], xleft= xleft, ybottom= ybottom, xright= xright, ytop= ytop, xpd= T)
}

#' figure label
#'
#' Plots a convenient label for paper figures
#'
#' @param text Text to plot.
#' @param region Can be "figure", "plot", "device"
#' @param cex Scaling factor
#' @param ... Extra args for text function
#' @export

vl_fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- grDevices::dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- graphics::grconvertX(c(0, ds[1]), from="in", to="user")
    y <- graphics::grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- graphics::par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- graphics::strwidth(text, cex=cex) * 60/100
  sh <- graphics::strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- graphics::par(xpd=NA)
  on.exit(graphics::par(old.par))
  graphics::text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


#' Ballons plot
#'
#' Starting from two matrices (dot size and colors, respectively), generates a balloon plot
#' 
#' @param x Matrix specifyiung size of the balloons.
#' @param color_var Matrix specifyiung color of the balloons.
#' @param color_breaks Color breaks used for coloring
#' @param col Vector of colors used for coloring
#' @param main Title. Default= NA
#' @param balloon_size_legend Title size legend
#' @param balloon_col_legend Title color legend
#' @param cex.balloons scaling factor for balloons size
#' @param auto_margins Use auto margins? Default= T
#' @examples
#' mat <- matrix(0:11, ncol = 3)
#' col <- matrix(6:-5, ncol = 3)
#' vl_balloons_plot(mat)
#' vl_balloons_plot(x= mat, color_var = col)
#' vl_balloons_plot(x= mat,
#' color_var = col,
#' balloon_size_legend = "test",
#' balloon_col_legend = "test2")
#' vl_balloons_plot(x= mat,
#' color_var = col,
#' balloon_size_legend = "test",
#' balloon_col_legend = "test2", 
#' cex.balloons = 8,
#' main= "Check!")
#' 
#' @return Balloon plot
#' @export
#' 
vl_balloons_plot <- function(x, ...)
  UseMethod("vl_balloons_plot")


#' @describeIn vl_balloons_plot matrix_method
#' @export

vl_balloons_plot.matrix <- function(x,
                                    color_var,
                                    x_breaks,
                                    color_breaks,
                                    col= c("cornflowerblue", "lightgrey", "tomato"),
                                    main= NA, 
                                    balloon_size_legend= NA,
                                    balloon_col_legend= NA,
                                    cex.balloons= 4,
                                    auto_margins= T)
{
  # Functions
  cex_scale <- function(x)
    (x-min(x_breaks, na.rm = T))/(max(x_breaks, na.rm = T)-min(x_breaks, na.rm = T))*cex.balloons+1
  compute_balloon_radius <- function(cex= cex.balloons)
    grconvertX(0.75/2, "chars", "in")*cex/2
  
  # Checks
  if(missing(x_breaks))
    x_breaks <- range(x, na.rm= T)
  if(missing(color_var))
  {
    color_var <- matrix(rep(1, nrow(x)*ncol(x)), 
                        nrow = nrow(x),
                        ncol = ncol(x))
    color_breaks <- c(0,1,2)
  }
  color_var <- color_var[nrow(color_var):1,]
  if(!identical(dim(x), dim(color_var)))
    stop("x and color_var matrices should have identical dimensions")
  if(missing(color_breaks))
    color_breaks <- seq(min(color_var, na.rm= T), 
                        max(color_var, na.rm= T), 
                        length.out= length(col))
  Cc <- circlize::colorRamp2(color_breaks, col)
  
  
  # Scale data
  scaled <- cex_scale(x)
  if(is.null(colnames(scaled)))
    colnames(scaled) <- seq(ncol(scaled))
  if(is.null(rownames(scaled)))
    rownames(scaled) <- seq(nrow(scaled))
  scaled <- scaled[nrow(scaled):1,]
  
  # Margins
  if(auto_margins)
  {
    bot <- 0.5+max(strwidth(colnames(x), "inches"))
    left <- 0.5+max(strwidth(rownames(x), "inches"))
    top <- 0.5+strheight(main, units = "inches")
    leg.width <- strwidth(c(balloon_size_legend, balloon_col_legend), "in")
    right <- 0.5+max(c(leg.width, 0.5))
    par(mai= c(bot, left, top, right),
        xaxs= "i",
        yaxs= "i")
  }
  # Init plot
  plot.new()
  # Compute xlim
  ball.radius <- compute_balloon_radius(cex= cex.balloons+1)
  plot.size.x <- grconvertX(1, "npc", "in")-grconvertX(0, "npc", "in")
  adj.x <- (ball.radius*1.3)/plot.size.x
  ext.x <- (ncol(scaled)-1)*adj.x
  xl <- c(1-ext.x, ncol(scaled)+ext.x)
  # Compute ylim
  plot.size.y <- grconvertY(1, "npc", "in")-grconvertY(0, "npc", "in")
  adj.y <- (ball.radius*1.3)/plot.size.y
  ext.y <- (nrow(scaled)-1)*adj.y
  yl <- c(1-ext.y, nrow(scaled)+ext.y)
  # Plot window
  plot.window(xlim = xl,
              ylim = yl)
  # Title
  title(main)
  # Grid
  segments(1, 
           seq(nrow(scaled)), 
           ncol(scaled),
           seq(nrow(scaled)))
  segments(seq(ncol(scaled)),
           1,
           seq(ncol(scaled)),
           nrow(scaled))
  # Points
  col_vec <- c(color_var)
  col_vec[!is.na(col_vec)] <- Cc(col_vec[!is.na(col_vec)])
  points(col(scaled),
         row(scaled),
         cex= c(scaled),
         pch= 21,
         bg= col_vec,
         col= "black")
  # Axes
  axis(1,
       at= seq(ncol(scaled)),
       labels = colnames(x),
       las= 2,
       lwd= 0,
       line= -0.5)
  axis(2,
       at= seq(nrow(scaled)),
       labels = rownames(scaled),
       las= 2,
       lwd= 0,
       line= -0.5)

  # Size legend
  x.leg <- grconvertX(1, "npc", "in")+grconvertX(0.5, "line", "in")+ball.radius
  x.leg.title <- grconvertX(x.leg-ball.radius, "in", "user")
  x.leg.text <- grconvertX(x.leg+ball.radius, "in", "user")
  x.leg <- grconvertX(x.leg, "in", "user")
  b.top <- grconvertY(nrow(scaled), "user", "in")-grconvertY(0.5, "line", "in")-ball.radius
  scale.values <- rev(axisTicks(x_breaks, log= F, nint = 4))
  scale.cex <- cex_scale(scale.values)
  scale.radius <- compute_balloon_radius(cex= scale.cex)
  scale.y <- sapply(seq(scale.cex[-1]), function(i) scale.radius[i]+scale.radius[i+1]+0.05)
  scale.y <- grconvertY(b.top-cumsum(c(0, scale.y)), "in", "user")
  points(rep(x.leg, length(scale.y)),
         scale.y, 
         cex= scale.cex, 
         xpd= T)
  text(x.leg.text, 
       scale.y, 
       scale.values,
       xpd= T, 
       pos= 4, 
       cex= 0.6,
       offset = 0.25)
  text(x.leg.title, 
       nrow(scaled), 
       balloon_size_legend,
       xpd= T, 
       pos= 4,
       offset = 0.5)
  
  # Color legend
  x.left <- grconvertX(x.leg, "user", "in")-ball.radius/2
  x.right <- x.left+grconvertX(1, "line", "in")
  y.leg <- grconvertY(min(scale.y), "user", "in")-grconvertY(2, "line", "in")
  y.top <- grconvertY(min(scale.y), "user", "in")-grconvertY(3, "line", "in")
  y.bot <- y.top-grconvertY(5, "line", "in")
  x.left <- grconvertX(x.left, "in", "user")
  y.bot <- grconvertY(y.bot, "in", "user")
  y.leg <- grconvertY(y.leg, "in", "user")
  x.right <- grconvertX(x.right, "in", "user")
  y.top <- grconvertY(y.top, "in", "user")
  col.im <- matrix(Cc(seq(max(color_breaks), min(color_breaks), length.out = 100)), ncol= 1)
  rasterImage(col.im,
              x.left, 
              y.bot, 
              x.right,
              y.top, 
              xpd= T)
  tick.lab <- axisTicks(range(color_breaks), log= F, nint = 4)
  tick.y <- y.bot+(y.top-y.bot)*(tick.lab-min(color_breaks))/(max(color_breaks)-min(color_breaks))
  text(rep(x.right, length(tick.y)), 
       tick.y,
       labels = tick.lab,
       cex= 0.6,
       pos= 4, 
       xpd= T, 
       offset = 0.25)
  text(x.leg.title,
       y.leg,
       labels = balloon_col_legend,
       pos= 4, 
       xpd= T, 
       offset = 0.5)
}

