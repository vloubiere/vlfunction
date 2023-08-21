#' Title
#'
#' @param x see ?plot()
#' @param y see ?plot()
#' @param labels labels to be plotted
#' @param label.cex Expansion factor for labels 
#' @param label.col Labels color
#' @param xlim Scatterplot xlim. default to range(x)
#' @param ylim Scatterplot xlim. default to range(y)
#' @param rect.draw Draw rectangle?
#' @param rect.col Rectangles color
#' @param points_size 
#' @param box_padding_x 
#' @param box_padding_y 
#' @param point_padding_x 
#' @param point_padding_y 
#' @param ... Extra arguments to be passed to plot()
#'
#' @return Plots a ggrepel-like scatterplot. Returns plotting object as data.table
#' @export
#'
#' @examples
#' data("mtcars")
#' mtcars$car <- rownames(mtcars)
#' vl_repelLabels(x = mtcars$wt,
#' y = mtcars$mpg,
#' labels = mtcars$car,
#' cex = 0.7)

vl_repelLabels <- function(x,
                           y,
                           labels,
                           label.cex= 0.7,
                           label.col= "black",
                           xlim= range(x),
                           ylim= range(y),
                           rect.draw= T,
                           rect.col= adjustcolor("lightgrey", 0.5),
                           point_size= 0.01,
                           box_padding_x = 0,
                           box_padding_y = 0,
                           point_padding_x = 0,
                           point_padding_y = 0,
                           ...)
{
  plot(x, y, ...)
  # USe latticetools to compute positions
  coords <- latticetools::get_repelled_labels(x = x,
                                              y = y,
                                              labels = labels,
                                              cex = label.cex,
                                              point_size = point_size,
                                              xlim = xlim,
                                              ylim = ylim,
                                              box_padding_x = box_padding_x,
                                              box_padding_y = box_padding_y,
                                              point_padding_x = point_padding_x,
                                              point_padding_y = point_padding_y,
                                              force_push = 1e-05,
                                              force_pull = 1e-05,
                                              transform_coord = "native")
  coords <- as.data.table(coords)
  coords[, rect.col:= rect.col]
  coords[, label.col:= label.col]
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
  coords[, {
    clip(clip_x1[1], clip_x2[1], clip_y1[1], clip_y2[1])
    segments(x[1], y[1], x_orig[1], y_orig[1])
  }, .(x, y, x_orig, y_orig, clip_x1, clip_x2, clip_y1, clip_y2)]
  clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4])
  # Plot rectangles
  if(rect.draw)
    with(coords, rect(x1_box, y1_box, x2_box, y2_box, col = rect.col, border= NA))
  # Plot labels
  with(coords, text(x, y, labels = label, col = label.col, cex = 0.7))
  # Return object
  invisible(coords)
}
