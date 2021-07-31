#' bw Average tracks
#'
#' Plots average tracks for a set bw files around (potentially) several sets of peaks
#'
#' @param dat Data.table containing the variables to plot
#' @param col Colors to use for bars. Wrapped with colorRampPalette
#' @param connectionColDT Data.table containing V1/V2/Cc columns with V1/V2 containing all possible transitions. See examples below
#' @examples 
#' if connectionColDT is set to null, the data.table is built internatlly using:
#'   computeConnectionCol <- function(dat, col)
#' {
#'   .c <- CJ(unlist(dat), unlist(dat), unique = T)
#'   .c[, Cc:= col]
#'   return(.c)
#' }
#' @return An object that can be used with the vl_average_bw_track_plot_only() function.
#' @export



vl_alluvial_plot <- function(dat,
                             col= c("tomato", "lightgrey", "cornflowerblue"),
                             connectionColDT= NULL)
{
  # Build connections color table
  computeConnectionCol <- function(dat, col)
  {
    .c <- CJ(unlist(dat), unlist(dat), unique = T)
    .c[, Cc:= col]
    return(.c)
  }
  if(is.null(connectionColDT))
  {
    N_cmb <- length(unique(unlist(dat)))^2
    Cc <- colorRampPalette(col)(N_cmb)
    Cc <- adjustcolor(Cc, 0.6)
    connectionColDT <- computeConnectionCol(dat, 
                                            col = Cc) 
  }
  if(!is.data.table(connectionColDT) | !identical(c("V1", "V2", "Cc"), names(connectionColDT)))
    stop("connectionColDT should be a data.table. See ?vl_alluvial_plot")
  
  # Compute constant
  xpos <- seq(0, 1, length.out= length(dat)*2)
  Cc <- colorRampPalette(col)(length(unique(unlist(dat))))
  cCc <- connectionColDT
  
  # PLOT
  plot.new()
  for(i in seq(dat))
  {
    classes <- table(dat[[i]])
    
    ypos <- c(0, cumsum(classes)/nrow(dat))
    
    x1 <- xpos[(i-1)*2+1]
    x2 <- xpos[i*2]
    x3 <- xpos[i*2+1]
    
    rect(xleft = x1, 
         xright = x2, 
         ybottom = ypos[-length(ypos)],
         ytop = ypos[-1],
         col= Cc)
    text(x= mean(c(x1, x2)),
         y= ypos[-length(ypos)]+diff(ypos)/2,
         labels= paste0(names(classes), "\n(", classes, ")"),
         xpd= T)
    text(x= mean(c(x1, x2)),
         y= 1,
         labels= names(dat)[i], 
         xpd= T, 
         pos=3)
    
    cols <- na.omit(names(dat)[i:(i+1)])
    if(length(cols)==2)
    {
      pol <- dat[, .N, keyby= c(cols)]
      pol[, ybot1:= cumsum(N)/sum(N)]
      pol[, ytop1:= c(0, ybot1[-.N])]
      setkeyv(pol, names(pol)[2])
      pol[, ybot2:= cumsum(N)/sum(N)]
      pol[, ytop2:= c(0, ybot2[-.N])]
      pol <- merge(pol, cCc, 
                   by.x= names(pol)[1:2],
                   by.y= c("V1", "V2"))
      pol[, {
        polygon(c(x2, x2, x3, x3),
                c(ytop1, ybot1, ybot2, ytop2),
                border= NA,
                col= Cc[1])
      }, (pol)]
    }
  }
}