#' Get bounding boxes for repelled text labels 
#' 
#' @param x,y (numeric) vectors representing x and y coordinates of points
#' @param labels (character) vector of labels
#' @param cex (numeric) size of text labels
#' @param xlim numeric vector representing the limits on the x axis like
#'   \code{c(xmin, xmax)}
#' @param ylim numeric vector representing the limits on the y axis like
#'   \code{c(ymin, ymax)}
#'  ‘force_push’ ‘force_pull’ ‘transform_coord’
#' @param point_size (numeric) approximate size of points in npc relative cooordinates. Can be small,
#'   is used in conjunction with "point_padding" options
#' @param box_padding_x (numeric) extra margin around text bounding boxes
#' @param box_padding_y (numeric) extra margin around text bounding boxes
#' @param point_padding_x (numeric) extra margin around data points
#' @param point_padding_y (numeric) extra margin around data points
#' @param force_push (numeric) magnitude of the push force (defaults to \code{1e-5})
#' @param force_pull (numeric) magnitude of the pull force (defaults to \code{1e-5})
#' @param transform_coord (character) the coordinate system used for calculating
#'   bounding boxes. The default, "npc", is to transfrom all cordinates to
#'   'Normalised Parent Coordinates'. When calling externally, use "native"
#'   instead, see examples.
#' @param ... other arguments passed to the function
#' 
#' @examples
#' data("mtcars")
#' mtcars$car <- rownames(mtcars)
#' 
#' # with base graphics
#' # ---------------------------
#' plot(mtcars$mpg ~ mtcars$wt, pch = 19)
#' 
#' coords <- get_repelled_labels(
#'   x = mtcars$wt, y = mtcars$mpg, 
#'   labels = mtcars$car,
#'   cex = 0.7, point_size = 0.01,
#'   xlim = range(mtcars$wt),
#'   ylim = range(mtcars$mpg),
#'   box_padding_x = 0.05,
#'   box_padding_y = 0.05,
#'   point_padding_x = 0.01,
#'   point_padding_y = 0.2,
#'   force_push = 1e-05,
#'   force_pull = 1e-05,
#'   transform_coord = "native"
#' )
#' 
#' with(coords, rect(x1_box, y1_box, x2_box, y2_box, col = "white"))
#' with(coords, text(x, y, labels = label, cex = 0.7))
#' 
#' # with ggplot2
#' # ---------------------------
#' 
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(y = mpg, x = wt)) +
#'   geom_point() +
#'   geom_label(aes(x = coords$x, y = coords$y, label = coords$label))
#' }
#' 
#' @export
get_repelled_labels <- function(x,
                                y,
                                labels,
                                cex = 1,
                                xlim = c(0, 1),
                                ylim = c(0, 1),
                                point_size = 0.001,
                                box_padding_x = 0.005,
                                box_padding_y = 0.005,
                                point_padding_x = 0.01,
                                point_padding_y = 0.01,
                                force_push = 1e-05,
                                force_pull = 1e-05,
                                transform_coord = "npc")
{
  # I implemented this step to allow different sizes for different labels
  if(length(cex)==1)
    cex <- rep(cex, length(x))
  
  # set unit for string dimensions
  str_unit <- ifelse(transform_coord == "npc", "inches", "user")
  
  # functions to convert units
  convert_x <- function(x, from = "inches", to = transform_coord) {
    if (from == "user") from = "native"
    grid::convertX(grid::unit(x, from), unitTo = to, valueOnly = TRUE)
  }
  convert_y <- function(y, from = "inches", to = transform_coord) {
    if (from == "user") from = "native"
    grid::convertY(grid::unit(y, from), unitTo = to, valueOnly = TRUE)
  }
  
  # convert all input units to desired coord system; default: "npc" (relative coordinates)
  x = convert_x(x, from = "native")
  y = convert_y(y, from = "native")
  
  # calculate extra margin around text for the bounding boxes
  box_margin_x <- sapply(cex, function(x) {
    convert_x(strwidth("M", str_unit, x), from = str_unit)
  })
  box_margin_y <- sapply(cex, function(x) {
    convert_x(strheight("M", str_unit, x), from = str_unit)
  })
  
  # calculate bounding box of text labels
  box_x = sapply(seq(labels), function(i) {
    convert_x(strwidth(labels[i], str_unit, cex[i]), from = str_unit) + box_margin_x[i]
  })
  box_y = sapply(seq(labels), function(i) {
    convert_y(strheight(labels[i], str_unit, cex[i]), from = str_unit) + box_margin_y[i]
  })
  
  # browser()
  # x and y positions as dataframe
  posdf <- data.frame(x = x, y = y)
  
  # bounding box coordinates for each label
  boxdf <- data.frame(
    x1 = x - box_x / 2 - box_padding_x,
    y1 = y - box_y / 2 - box_padding_y,
    x2 = x + box_x / 2 + box_padding_x,
    y2 = y + box_y / 2 + box_padding_y
  )
    
  # Compute boxes
  moved <- latticetools::repel_boxes(data_points = as.matrix(posdf),
                                     point_size = rep(point_size, nrow(posdf)),
                                     point_padding_x = point_padding_x,
                                     point_padding_y = point_padding_y,
                                     boxes = as.matrix(boxdf),
                                     xlim = xlim, ylim = ylim,
                                     hjust = rep(0, nrow(posdf)),
                                     vjust = rep(0, nrow(posdf)),
                                     force_push = force_push,
                                     force_pull = force_pull,
                                     max_time = 3, 
                                     max_overlaps = 30, 
                                     max_iter = 100000,
                                     direction = "both", 
                                     verbose = FALSE)
    
  coords <- cbind(posdf, moved)
  names(coords) <- c("x_orig", "y_orig", "x", "y", "overlapping")
  coords$x1_box <- coords$x - box_x / 2
  coords$y1_box <- coords$y - box_y / 2
  coords$x2_box <- coords$x + box_x / 2
  coords$y2_box <- coords$y + box_y / 2
  
  # convert units back to native
  coords[c("x_orig", "x", "x1_box", "x2_box")] <- apply(
    coords[c("x_orig", "x", "x1_box", "x2_box")], 2, function(x) {
      convert_x(x, from = transform_coord, to = "native")
    })
  coords[c("y_orig", "y", "y1_box", "y2_box")] <- apply(
    coords[c("y_orig", "y", "y1_box", "y2_box")], 2, function(y) {
      convert_y(y, from = transform_coord, to = "native")
    })
  
  # add labels
  coords$label <- labels
  
  # return coordinates
  return(coords)
}

#' Title
#'
#' @param x see ?plot()
#' @param y see ?plot()
#' @param labels labels to be plotted
#' @param label.sel Numeric or boolean vector indicating which labels should be plotted
#' @param label.cex Expansion factor for labels 
#' @param label.col Labels color
#' @param xlab X axis label
#' @param ylab Y axis label
#' @param xlim Scatterplot xlim. default to range(x)
#' @param ylim Scatterplot xlim. default to range(y)
#' @param rect.draw Draw rectangle?
#' @param rect.col Rectangles color
#' @param rect.col Segments colors
#' @param point_size 	(numeric) approximate size of points, defined as a fraction of char size. Can be small, is used in conjunction with "point_padding". default= 0.001
#' @param box_padding_x (numeric) extra margin around text bounding boxes, defined as a fraction of char size. Default= 0.005
#' @param box_padding_y (numeric) extra margin around text bounding boxes, defined as a fraction of char size. Default= 0.005
#' @param point_padding_x (numeric) extra margin around data points, defined as a fraction of char size. Default= 0.01
#' @param point_padding_y (numeric) extra margin around data points, defined as a fraction of char size. Default= 0.01
#' @param force_push (numeric) magnitude of the push force (defaults to 1e-5)
#' @param force_pull (numeric) magnitude of the pull force (defaults to 1e-5)
#' @param add Should the label be added to existing plot? Default= F
#' @param frame Should the frame be plotted? Deault= F
#' @param ... Extra arguments to be passed to plot()
#'
#' @return Plots a ggrepel-like scatterplot. Returns plotting object as data.table
#' @export
#'
#' @examples
#' data("mtcars")
#' mtcars$car <- rownames(mtcars)
#' vl_repelScatterplot(x = mtcars$wt,
#' y = mtcars$mpg,
#' labels = mtcars$car,
#' cex = 0.7)

vl_repelScatterplot <- function(x,
                                y= NULL,
                                col= "black",
                                labels,
                                label.sel= rep(T, length(x)),
                                label.cex= 0.7,
                                label.col= "black",
                                xlab= ifelse(is.null(y), "Index", deparse1(substitute(x))),
                                ylab= ifelse(is.null(y), deparse1(substitute(x)), deparse1(substitute(y))),
                                xlim= if(add) par("usr")[1:2] else range(x),
                                ylim= if(add) par("usr")[3:4] else range(y),
                                rect.draw= T,
                                rect.col= col,
                                seg.col= col,
                                noOverlap= F,
                                point_size = 0.001,
                                box_padding_x = 0.005,
                                box_padding_y = 0.005,
                                point_padding_x = 0.01,
                                point_padding_y = 0.01,
                                force_push = 1e-05,
                                force_pull = 1e-05,
                                add= F,
                                frame= F,
                                ...)
{
  if(is.null(y))
  {
    y <- x
    x <- seq(y)
  }
  if(!add)
    plot.default(x= x, 
                 y= y,
                 xlab= xlab,
                 ylab= ylab,
                 frame= frame,
                 col= col,
                 ...)
  
  # Use latticetools to compute positions
  dat <- data.table(x,
                    y,
                    labels,
                    label.cex,
                    rect.col,
                    label.col,
                    seg.col)
  dat <- dat[(label.sel)]
  coords <- vlfunctions::get_repelled_labels(x = dat$x,
                                             y = dat$y,
                                             labels = dat$labels,
                                             cex = dat$label.cex,
                                             point_size = point_size,
                                             xlim = xlim,
                                             ylim = ylim,
                                             box_padding_x = box_padding_x,
                                             box_padding_y = box_padding_y,
                                             point_padding_x = point_padding_x,
                                             point_padding_y = point_padding_y,
                                             force_push = force_push,
                                             force_pull = force_pull,
                                             transform_coord = "native")
  coords <- as.data.table(coords)
  coords <- cbind(coords,
                  dat[, rect.col:seg.col])
  if(noOverlap && any(!coords$overlapping))
  {
    warning(paste0(sum(!coords$overlapping, na.rm= T), " labels were overlapping and were removed from the plot!"))
    coords <- coords[!(overlapping)]
  }
  # Compute clipping masks for segments
  coords[, clip_x1:= par("usr")[1]]
  coords[, clip_x2:= par("usr")[2]]
  coords[, clip_y1:= par("usr")[3]]
  coords[, clip_y2:= par("usr")[4]]
  coords[, max.y.diff:= ((y2_box-y1_box)/2)*(abs(x_orig-x)/((x2_box-x1_box)/2))]
  coords[(y_orig-y)>=max.y.diff, clip_y1:= y2_box]
  coords[(y_orig-y)<=(-max.y.diff), clip_y2:= y1_box]
  coords[abs(y_orig-y)<=max.y.diff & x_orig>=x2_box, clip_x1:= x2_box]
  coords[abs(y_orig-y)<=max.y.diff & x_orig<=x1_box, clip_x2:= x1_box]
  # Plot segments
  coords[, seg_length.x:= min(abs(c(x1_box, x2_box)-x_orig)), .(x1_box, x2_box, x_orig)] # Check if segments are long enough
  coords[, seg_length.y:= min(abs(c(y1_box, y2_box)-y_orig)), .(y1_box, y2_box, y_orig)]
  coords[, seg_length:= seg_length.x>strwidth("M", cex= label.cex/3) | seg_length.y>strheight("M", cex= label.cex/3)]
  if(any(coords$seg_length))
  {
    coords[(seg_length), 
           {
             clip(clip_x1[1], clip_x2[1], clip_y1[1], clip_y2[1])
             segments(x[1], y[1], x_orig[1], y_orig[1], col= seg.col[1])
           }, .(x, y, x_orig, y_orig, clip_x1, clip_x2, clip_y1, clip_y2, seg.col)]
  }
  clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4])
  # Plot rectangles
  if(rect.draw)
    with(coords, rect(x1_box, y1_box, x2_box, y2_box, col = rect.col, border= NA))
  # Plot labels
  with(coords, text(x, y, labels = label, col = label.col, cex = label.cex))
  # Return object
  invisible(coords)
}
